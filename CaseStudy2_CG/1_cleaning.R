library(readxl)
library(dplyr)
library(data.table)
library(getopt)
library(plyr)
library(lme4)
library(pbapply)
library(parallel)
library(sva)
library(RUVIIIC)
library(tibble)
library(diann)

rm(list = ls())
setwd("./CaseStudy2_CG/data/")

PQM <- "MaxLFQ"
PQM_Dir <- "lfq"


## Clinical data -----------------------
## lb
raw <- read.csv("rawfiles/clinical/CGZ_lb_vs chipscreen/lb.txt")
des <- read_excel("rawfiles/clinical/dataset_spec.xlsx", sheet = 2)
output_path <- "clindata.csv"

cols <- c("SUBJID", "VISIT", "VISITDY", "CTLBYN", "ablfl", "ablfl02", "PARAMCD", "AVAL1C")

# ## vs
# raw <- read.csv("rawfiles/clinical/CGZ_lb_vs chipscreen/vs.txt")
# des <- read_excel("rawfiles/clinical/dataset_spec.xlsx", sheet = 1)
# output_path <- "clindata.csv"
# 
# cols <- c("SUBJID", "VISIT", "VISITDY", "ablfl", "ablfl02", "PARAMCD", "AVAL1N")

## long to wide
long <- raw %>%
  select(all_of(cols)) %>%
  reshape2::melt(., id = cols, na.rm = TRUE)

vars <- cols[!cols %in% c("PARAMCD", "AVAL1N", "AVAL1C")]
form <- formula(paste(paste0(vars, collapse = " + "), "PARAMCD", sep = " ~ "))
wide <- reshape2::dcast(long, form, value.var = "AVAL1C")

## impute NAs
wide <- wide %>%
  filter(!is.na(VISITDY)) %>%
  mutate(ID = paste(SUBJID, VISIT, VISITDY, sep = "_"))

wide_new_list <- pblapply(unique(wide$ID), function(id) {
  wide_new_i <- wide %>%
    filter(ID %in% id) %>%
    select(!any_of(c("ID", cols))) %>%
    apply(., 2, function(x) {
      if (is.numeric(x)) {
        x <- median(x, na.rm = TRUE)
      } else if (is.character(x)) {
        x <- unique(na.omit(x))
      }
      return(x)
    }) %>%
    t %>%
    as.data.frame %>%
    mutate(ID = id)

  return(wide_new_i)
})

wide_new <- wide_new_list %>%
  rbindlist(.) %>%
  tidyr::separate(., ID, c("Subject.ID", "Week", "Day"), sep = "_") %>%
  select(Subject.ID, Week, Day, everything())

## annotate with weeks and IDs
wide_df <- wide_new %>%
  mutate_at("Week", ~ stringr::str_extract(., "(?<=\\().+?(?=周\\))")) %>%
  mutate_at("Week", ~ ifelse(is.na(.), "Baseline", paste("Week", ., sep = "")))

## save
fwrite(wide_df, output_path)


## Metadata ---------------------------
## QC info
meta_sample <- read_excel("rawfiles/metadata/样本信息汇总.xlsx")
meta_order <- read.csv("rawfiles/metadata/数据产出顺序+清洗日期.csv")
meta_qc <- read.csv("rawfiles/metadata/QC样本揭盲.csv")

meta_exp <- meta_sample %>%
  dplyr::rename(Plasma = Sample, Tube = `Ep管编号`) %>%
  left_join(., meta_order, by = "Tube") %>%
  left_join(., meta_qc, by = "Plasma") %>%
  dplyr::rename(Order = LoadingOrder) %>%
  mutate(Cleaning = sapply(Date, function(x) if (x < 20210223) return("I") else if (x > 20210306) return("III") else return("II"))) %>%
  mutate(Hemolyzed = sapply(`溶血标注`, function(x) ifelse(is.na(x), "", "Yes"))) %>%
  mutate(Type = sapply(Label, function(x) ifelse(is.na(x), "Study sample", paste("QC sample (", x, ")", sep = ""))))

## Study info
meta_subject <- read.csv("rawfiles/metadata/20200606_Chiglitazar_multiomics_metadata_zs.csv")
meta_std <- meta_subject %>%
  filter(Reported == "Yes") %>%
  select(Subject.ID, Project.ID, Site.ID, Randomized.ID, DOB, Age, Sex, DOS, Arm, Extended.Treatment, Training.Validation, Weight, BloodSample, Plasma.Week1, Plasma.Week12, Plasma.Week24) %>%
  reshape2::melt(., id = c("Subject.ID", "Project.ID", "Site.ID", "Randomized.ID", "DOB", "Age", "Sex", "DOS", "Arm", "Extended.Treatment", "Training.Validation", "Weight", "BloodSample"), na.rm = TRUE) %>%
  dplyr::rename(Plasma = value, Blood = BloodSample)

meta_all <- meta_exp %>%
  left_join(., meta_std, by = "Plasma") %>%
  arrange(Order) %>%
  mutate(ID = sapply(Tube, function(x) paste("Exp", x, sep = ""))) %>%
  mutate(Week = sapply(variable, function(x) if (x %in% "Plasma.Week1") return("Baseline") else if (x %in% "Plasma.Week24") return("Week24") else return(""))) %>%
  # mutate(Group = rep(c("I-1", "I-2", "I-3", "II-1", "II-2", "III-1", "III-2", "III-3"), c(192, 207, 104, 268, 198, 196, 166, 164))) %>%
  mutate(Group = rep(c("I-1", "I-2", "I-3", "II-1", "II-2", "II-3", "III-1", "III-2"), c(144, 255, 104, 120, 196, 150, 196, 330))) %>%
  mutate_all(as.character) %>%
  rename_with(~ gsub(".ID", "", .x, fixed = TRUE)) %>%
  select(ID, Plasma, Project, Site, Randomized, Box, Position, Tube, Cleaning, Date, Order, Group, Type, Subject, DOB, Age, Sex, DOS, Weight, Arm, Extended.Treatment, Week, Blood, Hemolyzed, Training.Validation)

meta_all[is.na(meta_all)] <- ""

fwrite(meta_all, "metadata.csv", row.names = FALSE)

# ## Mapping info
# mapped <- data.table::fread("expfiles/pep2pro/v202304/Protein_Data_latest/expdata_log2_f81257_s1476.csv")
# mapped <- mapped %>%
#   tidyr::separate_rows(Peptide, sep = "; |;") %>%
#   tidyr::separate_rows(Protein.Ids, sep = "; |;") %>%
#   dplyr::rename(Accession = Protein.Ids) %>%
#   distinct(Peptide, Accession)
# 
# fwrite(mapped, "expfiles/mapped.csv")


## functions for correction and quantification. --------------------
ratio_by_pm <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "library"
  
  df_quant_wide <- df_quant_wide %>%
    select(all_of(c("Feature", df_meta$library)))
  
  batches <- unique(df_meta$Cleaning)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    df_meta_i <- df_meta %>%
      filter(Cleaning %in% batch)
    
    df_meta_i_pm <- df_meta %>%
      filter(Cleaning %in% batch, grepl("PM", Type))
    
    match_batch <- which(samples %in% df_meta_i$library)
    match_pm <- which(samples %in% df_meta_i_pm$library)
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(PMmean = apply(.[match_pm], 1, mean)) %>%
      mutate_at(match_batch, ~ (. - PMmean)) %>%
      select(!PMmean)
    
  }
  
  return(df_corrected)
}

median_centering <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "library"
  
  df_quant_wide <- df_quant_wide %>%
    select(all_of(c("Feature", df_meta$library)))
  
  batches <- unique(df_meta$Cleaning)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    df_meta_i <- df_meta %>%
      filter(Cleaning %in% batch)
    
    match_batch <- which(samples %in% df_meta_i$library)
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(all_median = apply(.[match_batch], 1, median, na.rm = TRUE)) %>%
      mutate_at(match_batch, ~ (. - all_median)) %>%
      select(!all_median)
    
  }
  
  return(df_corrected)
}

combat_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  data_mode <- model.matrix(~ as.factor(Type), data = df_meta)
  
  df_corrected <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0, log(.))) %>%
    as.matrix() %>%
    ComBat(., df_meta$Cleaning, data_mode) %>%
    as.data.frame() %>%
    mutate_all(~ na_if(., 0)) %>%
    rownames_to_column(feature_header)
  
  return(df_corrected)
}

ruv3c_solve <- function(df_quant_wide, df_meta, ctrl_features, k) {
  
  feature_header <- colnames(df_quant_wide)[1]
  data_mode <- model.matrix(~ Type - 1, data = df_meta)
  
  colnames(df_quant_wide)[1] <- "feature"
  vary_features <- setdiff(df_quant_wide$feature, ctrl_features)
  
  df_corrected <- df_quant_wide %>%
    column_to_rownames("feature") %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, NA, log(.))) %>%
    t %>%
    RUVIII_C(k, ., data_mode,
             toCorrect = vary_features,
             controls = ctrl_features) %>%
    t %>%
    as.data.frame %>%
    mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    rownames_to_column(feature_header)
  
  return(df_corrected)
}

correct_by_loess <- function(y, x) {
  
  data_input <- data.frame(value = y, order = x)
  
  values <- data_input$value
  orders <- data_input$order
  
  if (length(values) > 2) {
    
    drift_loess <- loess(
      formula = value ~ order,
      data = data_input,
      control = loess.control(surface = "direct"),
      degree = 2,
      normalize = FALSE)
    
    drift <- predict(drift_loess, orders)
    values_final <- values - drift + median(values)
    
  } else {
    
    values_final <- values
    
  }
  
  return(values_final)
}

loess_solve <- function(df_quant_wide, df_meta) {
  
  colnames(df_quant_wide)[1] <- "Feature"
  feature_n <- nrow(df_quant_wide)
  
  batches <- unique(df_meta$Cleaning)
  df_quant_wide <- df_quant_wide %>%
    mutate_if(is.numeric, ~ ifelse(. == 0, NA, log(.)))
  
  loess_tables <- mclapply(1:feature_n, function(i) {
    
    loess_tables_i <- mclapply(batches, function(batch_id) {
      
      data_test_j <- df_quant_wide[i, ] %>%
        reshape2::melt(., id.vars = 1, verbose = FALSE) %>%
        as.data.frame %>%
        na.omit %>%
        dplyr::rename(ID = variable) %>%
        left_join(., df_meta, by = "ID") %>%
        filter(Cleaning %in% batch_id) %>%
        mutate(value_loess = correct_by_loess(value, Order))
      
      return(data_test_j)
    }, mc.cores = 2)
    
    data_test_i <- rbindlist(loess_tables_i)
    
    return(data_test_i)
  }, mc.cores = 2)
  
  data_loess <- data.table::rbindlist(loess_tables)
  
  expr_loess <- reshape2::dcast(data_loess, Feature ~ ID, value.var = "value_loess")
  
  return(expr_loess)
}

correct_by_x <- function(df_quant_wide, df_meta,
                         method = c("ratio_by_pm",
                                    "median_centering",
                                    "combat",
                                    "RUV-III-C",
                                    "loess"),
                         ctrl_features = NULL,
                         k=3) {
  
  if (method == "ratio_by_pm") {
    df_corrected <- ratio_by_pm(df_quant_wide, df_meta)
  } else if (method == "median_centering") {
    df_corrected <- median_centering(df_quant_wide, df_meta)
  } else if (method == "combat") {
    df_corrected <- combat_solve(df_quant_wide, df_meta)
  } else if (method == "RUV-III-C") {
    df_corrected <- ruv3c_solve(df_quant_wide, df_meta, ctrl_features, k)
  } else if (method == "loess") {
    df_corrected <- loess_solve(df_quant_wide, df_meta)
  }
  
  return(df_corrected)
}

maxlfq_solve <- function(quantities,
                         peptides,
                         samples,
                         margin = -10.0001) {
  
  .Call('_diann_maxlfq_solve',
        PACKAGE = 'diann',
        quantities, peptides, samples, margin)
}

maxlfq_filter <- function(df_quant,
                          protein.header = "protein_group_ids",
                          peptide.header = "peptide_id",
                          sample.id.header = "library",
                          intensity.header = "intensity") {
  
  colnames(df_quant)[match(intensity.header, colnames(df_quant))] <- "value"
  
  all_samples <- df_quant %>%
    pull(sym(sample.id.header)) %>%
    unique
  
  if (length(all_samples) > 1) {
    
    all_pairs <- all_samples %>%
      as.character %>%
      combn(., 2) %>%
      t %>%
      as.data.frame
    
    df_quant_fil <- all_pairs %>%
      merge(., df_quant,
            by.x = "V1",
            by.y = sample.id.header) %>%
      merge(., df_quant,
            by.x = c("V2", peptide.header),
            by.y = c(sample.id.header, peptide.header)) %>%
      mutate(paired = ifelse(value.x == 0 | value.y == 0, NA, 1))
    
  } else {
    
    df_quant_fil <- data.frame(paired = 1)
    
  }
  
  if (nrow(df_quant_fil) > 1) {
    
    df_quant_wide <- df_quant_fil %>%
      reshape2::dcast(., formula(paste("V1 + V2 ~ ", peptide.header, sep = "")),
                      value.var = "paired") %>%
      select(where(function(x) length(which(!is.na(x))) >= 2)) %>%
      filter(apply(., 1, function(x) length(which(!is.na(x))) >= 4))
    
    if (nrow(df_quant_wide) > 0) {
      
      keeped_peptides <- colnames(df_quant_wide)[-(1:2)]
      keeped_samples <- df_quant_wide %>%
        select(V1, V2) %>%
        as.matrix %>%
        as.vector %>%
        unique
      
      return(list(peptide_ids = keeped_peptides,
                  sample_ids = keeped_samples))
    }
  }
}

maxlfq_core <- function(df_quant, keeped_peptides, keeped_samples,
                        protein.header = "protein_group_ids",
                        peptide.header = "peptide_id",
                        sample.id.header = "library",
                        intensity.header = "intensity") {
  
  quant_validated_j <- df_quant %>%
    filter(eval(parse(text = paste(sample.id.header, "%in% keeped_samples"))),
           eval(parse(text = paste(peptide.header, "%in% keeped_peptides")))) %>%
    reshape2::acast(., formula(paste(peptide.header, "~", sample.id.header)),
                    value.var = intensity.header) %>%
    log
  
  quant_validated_j[is.infinite(quant_validated_j)] <- NA
  
  lfq_values <- maxlfq_solve(as.vector(quant_validated_j),
                             nrow(quant_validated_j),
                             ncol(quant_validated_j))
  
  names(lfq_values) <- sort(keeped_samples)
  
  return(lfq_values)
  
}

pep2pro_by_x <- function(df_quant, df_mapping,
                         method = c("MaxLFQ", "iBAQ", "TopPep3")) {
  
  df_mapping <- df_mapping %>%
    mutate(peptide_id = as.numeric(as.factor(sequence))) %>%
    mutate(protein_group_id = as.numeric(as.factor(protein_names)))
  
  df_quant <- df_quant %>%
    reshape2::melt(., id = 1,
                   variable.name = "library",
                   value.name = "intensity") %>%
    left_join(df_mapping, ., by = "sequence",
              relationship = "many-to-many") %>%
    na.omit
  
  all_samples <- sort(unique(df_quant$library))
  all_proteins <- unique(df_quant$protein_group_id)
  
  if (method == "MaxLFQ") {
    
    lfq_tables <- mclapply(all_proteins, function(i) {
      
      df_quant_i <- df_quant %>%
        filter(protein_group_id == i)
      
      # keeped_list <- maxlfq_filter(df_quant = df_quant_i)
      # 
      # final_peptides <- keeped_list[[1]]
      # final_samples <- keeped_list[[2]]
      
      final_peptides <- unique(df_quant_i$peptide_id)
      final_samples <- unique(df_quant_i$library)
      
      lfq_i <- rep(NA, length(all_samples))
      names(lfq_i) <- all_samples
      
      if (length(final_peptides) > 1 & length(final_samples) > 1) {
        
        lfq_values <- maxlfq_core(df_quant = df_quant_i,
                                  keeped_peptides = final_peptides,
                                  keeped_samples = final_samples)
        
        lfq_i[match(names(lfq_values), names(lfq_i))] <- lfq_values
        
      }
      
      quant_lfq_i <- data.frame(protein_group_id = i,
                                library = names(lfq_i),
                                lfq = lfq_i)
      
      return(quant_lfq_i)
      
    }, mc.cores = 4)
    
    df_lfq <- data.table::rbindlist(lfq_tables)
    
    df_lfq_wide <- reshape2::dcast(df_lfq,
                                   protein_group_id ~ library,
                                   value.var = "lfq")
    
    df_final <- df_mapping %>%
      distinct(protein_group_id, protein_names) %>%
      left_join(., df_lfq_wide, by = "protein_group_id") %>%
      select(!protein_group_id)
    
  }
  
  if (method == "iBAQ") {
    
    ibaq_tables <- mclapply(all_proteins, function(i) {
      
      df_quant_i <- df_quant %>%
        filter(protein_group_id == i) %>%
        group_by(library, protein_names) %>%
        dplyr::summarise(ibaq = sum(intensity), .groups = "keep") %>%
        mutate_at("ibaq", ~ ifelse(. == 0, NA, log(.)))
      
      return(df_quant_i)
      
    }, mc.cores = 4)
    
    df_ibaq <- data.table::rbindlist(ibaq_tables)
    
    df_ibaq_wide <- reshape2::dcast(df_ibaq,
                                    protein_names ~ library,
                                    value.var = "ibaq")
    
    df_final <- df_ibaq_wide
  }
  
  if (method == "TopPep3") {
    
    top3_tables <- mclapply(all_proteins, function(i) {
      
      df_quant_i <- df_quant %>%
        filter(protein_group_id == i) %>%
        filter(intensity != 0) %>%
        group_by(peptide_id) %>%
        dplyr::summarise(total_intensity = sum(intensity))
      
      if (nrow(df_quant_i) >= 3) {
        
        kept_peps <- df_quant_i %>%
          filter(total_intensity >= sort(total_intensity, decreasing = TRUE)[3]) %>%
          pull(peptide_id)
        
        df_top3_i <- df_quant %>%
          filter(protein_group_id == i) %>%
          filter(peptide_id %in% kept_peps) %>%
          group_by(library, protein_names) %>%
          dplyr::summarise(top3 = mean(intensity, na.rm = TRUE), .groups = "keep") %>%
          mutate_at("top3", ~ ifelse(. == 0, NA, log(.)))
        
        return(df_top3_i)
      }
    }, mc.cores = 4)
    
    df_top3 <- data.table::rbindlist(top3_tables)
    
    df_top3_wide <- reshape2::dcast(df_top3,
                                    protein_names ~ library,
                                    value.var = "top3")
    
    df_final <- df_top3_wide
  }
  
  return(df_final)
}


## Peptide-level data: 约20分钟-----------------------------------
meta_raw <- fread("metadata.csv")
expr_raw0 <- fread("rawfiles/peptide/Sequence_merge.csv")

## check column names
dict_colnames <- sapply(colnames(expr_raw0), function(x) {
    if (x == "Area._A49.x") {
      x <- "Area._A95"
    } else if (x == "Area._A49.y") {
      x <- "Area._A49"
    }
    if (grepl("Area", x)) {
      x_split <- strsplit(x, "_", fixed = T) %>% unlist
      x_final <- paste("Exp", x_split[2], sep = "")
    } else if (grepl("Exp", x)) {
      x_final <- x
    } else {
      x_final <- "Feature"
    }
  })

setdiff(meta_raw$ID, dict_colnames)

## log-e transformation: remove peptides with 0 quantity
expr_log <- expr_raw0 %>%
  plyr::rename(dict_colnames) %>%
  tibble::column_to_rownames("Feature") %>%
  mutate_all(~ na_if(., 0)) %>%
  log %>%
  tibble::rownames_to_column("Feature")

setdiff(meta_raw$ID, colnames(expr_log))

fwrite(expr_log, "expfiles/peptide/expdata_log.csv", row.names = FALSE)

## quantile normalization
expr_quantile <- expr_log %>%
  column_to_rownames("Feature") %>%
  limma::normalizeQuantiles(.,ties = T)

expr_quantile <- expr_quantile %>%
  rownames_to_column("Feature")

fwrite(expr_quantile, "expfiles/peptide/expr_quantile.csv", row.names = FALSE)

## loess curve fitting: 约3小时
df_quantile_exp <- expr_quantile %>%
  mutate_if(is.numeric, exp) %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))

df_loess <- correct_by_x(df_quant_wide = df_quantile_exp,
                         df_meta = meta_raw,
                         method = "loess")

fwrite(df_loess, file = paste("expfiles/peptide/expdata_loess.csv", sep = ""), row.names = FALSE)


## Protein-level data: 约2小时 ---------------------
mapping <- fread("mapped.csv")
meta_raw <- fread("metadata.csv")
expr_log <- fread("expfiles/peptide/expdata_log.csv")

## protein quantification
df_sum_final <- expr_log %>%
  mutate_if(is.numeric, exp)

colnames(df_sum_final)[1] <- "sequence"
colnames(mapping) <- c("sequence", "protein_names")

expr_pro <- pep2pro_by_x(df_quant = df_sum_final,
                         df_mapping = mapping,
                         method = PQM)

expr_pro_path <- paste("expfiles/protein/raw_mapped/", PQM_Dir, "/expdata_log.csv", sep = "")
fwrite(expr_pro, expr_pro_path, row.names = FALSE)

expr_pro_fot <- expr_pro %>%
  mutate_if(is.numeric, ~ ifelse(. <= -1e+06, NA, .)) %>%
  mutate_if(is.numeric, exp) %>%
  mutate_if(is.numeric, ~ (.) * 10 ^ (6) / sum(., na.rm = TRUE)) %>%
  mutate_if(is.numeric, log)

expr_pro_fot_path <- paste("expfiles/protein/raw_mapped/", PQM_Dir, "/expdata_log_fot.csv", sep = "")
fwrite(expr_pro_fot, expr_pro_fot_path, row.names = FALSE)

## quantile normalization
expr_pro <- fread(expr_pro_path)
colnames(expr_pro)[1] <- "Feature"

expr_quantile <- expr_pro %>%
  column_to_rownames("Feature") %>%
  limma::normalizeQuantiles(., ties = T)

expr_quantile <- expr_quantile %>%
  rownames_to_column("Feature")

expr_quantile_path <- paste("expfiles/protein/raw_mapped/", PQM_Dir, "/expdata_quantile.csv", sep = "")
fwrite(expr_quantile, expr_quantile_path, row.names = FALSE)

## loess curve fitting: 约1小时
df_quantile_exp <- expr_quantile %>%
  mutate_if(is.numeric, exp) %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))

df_loess <- correct_by_x(df_quant_wide = df_quantile_exp,
                         df_meta = meta_raw,
                         method = "loess")

expr_loess_path <- paste("expfiles/protein/raw_mapped/", PQM_Dir, "/expdata_loess.csv", sep = "")
fwrite(df_loess, file = expr_loess_path, row.names = FALSE)


## one-way ANOVA. 约1分钟 ------------------------
data_list <- list(expr_log, expr_pro)
types <- c("peptide", "protein")
anova_pep_tables <- pblapply(1:2, function(i) {
  
  df_quant <- data_list[[i]]
  colnames(df_quant)[1] <- "feature"
  
  df_sum_i <- df_quant %>%
    mutate_if(is.numeric, ~ ifelse(. <= -1e+06, NA, .)) %>%
    na.omit %>%
    select(all_of(c("feature", meta_raw$ID[grepl("PM", meta_raw$Type)])))
  
  all_features <- df_sum_i$feature
  
  anova_tables_i <- mclapply(all_features, function(feature_id) {
    
    df_anova_j <- df_sum_i %>%
      filter(feature %in% feature_id) %>%
      reshape2::melt(., id = 1, variable.name = "ID", na.rm = TRUE) %>%
      left_join(., meta_raw, by = "ID") %>%
      rstatix::anova_test(value ~ Cleaning)
    
    df_anova_j$feature <- feature_id
    
    return(df_anova_j)
    
  }, mc.cores = 3)
  
  df_anova_i <- anova_tables_i %>%
    rbindlist(.) %>%
    mutate(type = types[i])
  
  return(df_anova_i)
  
})

df_anova <- anova_pep_tables %>%
  rbindlist(.) %>%
  select(type, feature, everything())

df_anova_path <- paste("expfiles/protein/raw_mapped/", PQM_Dir, "/anova.csv", sep = "")
fwrite(df_anova, df_anova_path, row.names = FALSE)


## Peptide-level correction. 约3小时 ------------------------
mapping <- fread("mapped.csv")
meta_raw <- fread("metadata.csv")
expr_loess <- fread("expfiles/peptide/expdata_loess.csv")
anova_dt <- fread(paste("expfiles/protein/raw_mapped/", PQM_Dir, "/anova.csv", sep = ""))

colnames(mapping) <- c("sequence", "protein_names")

df_anova <- anova_dt %>%
  filter(type %in% "peptide", p < .05)

correct_methods <- c("combat", "median_centering", "ratio_by_pm", "RUV-III-C")
label_names <- c("combat", "median", "ratio", "ruv")

## 约3小时
correct_pep_tables <- pblapply(1:4, function(i) {
  
  if (correct_methods[i] %in% "RUV-III-C") {
    
    ctrl_peptides <- df_anova %>%
      pull(feature) %>%
      unique
    
  } else {
    
    ctrl_peptides <- NULL
    
  }
  
  df_sum_i <- expr_loess %>%
    mutate_if(is.numeric, exp) %>%
    mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))
  
  meta_i <- meta_raw %>%
    filter(ID %in% colnames(df_sum_i))
  
  df_i <- correct_by_x(df_quant_wide = df_sum_i,
                       df_meta = meta_i,
                       method = correct_methods[i],
                       ctrl_features = ctrl_peptides, k = 3)
  
  fwrite(df_i, file = paste("expfiles/peptide/expdata_",
                            label_names[i], ".csv", sep = ""), row.names = FALSE)
  
  return(df_i)
})

## 约1小时
correct_pro_tables <- pblapply(1:4, function(j) {
  
  df_j <- fread(paste("expfiles/peptide/expdata_",
                      label_names[j], ".csv", sep = ""), showProgress = FALSE)
  
  df_j_tmp <- df_j %>%
    mutate_if(is.numeric, exp)
  
  colnames(df_j_tmp)[1] <- "sequence"
  
  df_loess_tmp <- expr_loess %>%
    filter(Feature %in% df_j_tmp$sequence)
  
  df_j_tmp[is.na(df_loess_tmp)] <- NA
  
  df_final <- pep2pro_by_x(df_quant = df_j_tmp,
                           df_mapping = mapping,
                           method = PQM)
  
  fwrite(df_final, file = paste("expfiles/protein/final_mapped/", PQM_Dir, "/expdata_",
                          label_names[j], ".csv", sep = ""), row.names = FALSE)
  
})


## Protein-level correction. 约40分钟 ------------------------
meta_raw <- fread("metadata.csv")
expr_loess <- fread(paste("expfiles/protein/raw_mapped/", PQM_Dir, "/expdata_loess.csv", sep = ""))
anova_dt <- fread(paste("expfiles/protein/raw_mapped/", PQM_Dir, "/anova.csv", sep = ""))

df_anova <- anova_dt %>%
  filter(type %in% "protein", p < .05)

df_uncorrected <- expr_loess %>%
  mutate_if(is.numeric, exp) %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))

correct_methods <- c("combat", "median_centering", "ratio_by_pm", "RUV-III-C")
label_names <- c("combat", "median", "ratio", "ruv")

correct_pro_tables <- pblapply(1:4, function(j) {
  
  if (correct_methods[j] %in% "RUV-III-C") {
    
    ctrl_proteins <- df_anova %>%
      pull(feature) %>%
      unique
    
  } else {
    
    ctrl_proteins <- NULL
  }
  
  meta_i <- meta_raw %>%
    filter(ID %in% colnames(df_uncorrected))
  
  df_i <- correct_by_x(df_quant_wide = df_uncorrected,
                       df_meta = meta_i,
                       method = correct_methods[j],
                       ctrl_features = ctrl_proteins,
                       k = 3)
  
  expr_loess_tmp <- expr_loess %>%
    filter(Feature %in% df_i$Feature)
  
  df_i[is.na(expr_loess_tmp)] <- NA
  
  fwrite(df_i, file = paste("expfiles/protein/raw_mapped/", PQM_Dir,
                            "/expdata_", label_names[j], ".csv", sep = ""))
  
  return(df_i)
})
