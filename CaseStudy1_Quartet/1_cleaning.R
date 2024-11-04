# ' Title: Data cleaning; Batch correction.
# ' Author: Qiaochu Chen

library(parallel)
library(pbapply)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(diann)
library(sva)
library(harmony)
library(RUVIIIC)


rm(list = ls())
setwd('~/Desktop/毕业课题/manuscript-batch/CaseStudy1_Quartet/data/')


## metadata ----------------------------------
rm_order1 <- rep(c("FDU_Quartet_Peptide_D5_20171106", "FDU_Quartet_Peptide_D6_20171106", "FDU_Quartet_Peptide_F7_20171106", "FDU_Quartet_Peptide_M8_20171106"), 3)
rm_order2 <- rep(c("FDU_Quartet_Peptide_D5_20171106", "FDU_Quartet_Peptide_D6_20171106", "FDU_Quartet_Peptide_F7_20171106", "FDU_Quartet_Peptide_M8_20171106"), each = 3)
nvg_samples <- apply(cbind(rep(c("D5", "D6", "F7", "M8"), each = 3), rep(1:3, 4)), 1, function(x) paste0(x, collapse = "_"))

meta <- data.frame(
  entry1 = c(
    sapply(64:75, function(x) paste("LFQ intensity", x)),
    sapply(1:12, function(x) paste("LFQ intensity", x)),
    sapply(nvg_samples, function(x) paste("LFQ intensity", x)),
    sapply(76:87, function(x) paste("LFQ intensity", x)),
    apply(cbind(1:12, c("05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16")), 1, function(x) paste("LFQ intensity ", x[2],"_XM_SLM-20200706-DIA_", x[1], sep = "")),
    sapply(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), function(x) paste("LFQ intensity", x))),
  entry2 = c(
    sapply(64:75, function(x) paste("Intensity", x)),
    sapply(1:12, function(x) paste("Intensity", x)),
    sapply(nvg_samples, function(x) paste("Intensity", x)),
    sapply(76:87, function(x) paste("Intensity", x)),
    apply(cbind(1:12, c("05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16")), 1, function(x) paste("Intensity ", x[2],"_XM_SLM-20200706-DIA_", x[1], sep = "")),
    sapply(c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"), function(x) paste("Intensity", x))),
  rm = c(rep(rm_order1, 5), rm_order2),
  omics = "Shotgun-Proteomics",
  platform = "ThermoFisher",
  pattern = rep(c("DDA", "DDA", "DDA", "DIA", "DIA", "DIA"), each = 12),
  lab = rep(c("APT", "FDU", "NVG", "APT", "FDU", "BGI"), each = 12),
  instrument = rep(c("QE", "Lumos", "QEHFX", "QEHFX", "Lumos", "QEHF"), each = 12),
  sample = c(
    rep(rep(c("D5", "D6", "F7", "M8"), 3), 2),
    rep(c("D5", "D6", "F7", "M8"), each = 3),
    rep(rep(c("D5", "D6", "F7", "M8"), 3), 3)),
  replicate = c(
    rep(rep(1:3, each = 4), 2), 
    rep(1:3, 4),
    rep(rep(1:3, each = 4), 3)),
  database = "Uniprot",
  tool = "MaxQuant") %>%
  mutate(name = apply(., 1, function(x) paste0(x[9:10], collapse = "_"))) %>%
  mutate(batch = sapply(rep(1:6, each = 12), function(x) paste("DAT", x, sep = ""))) %>%
  mutate(library = apply(., 1, function(x) paste0(x[c(14, 9:10)], collapse = "_"))) %>%
  select(library, batch, name, everything())

write.csv(meta, "expfiles/meta.csv", row.names = F)


## get peptide-level data and mapping info. ---------------
## functions
input_mq_files <- function(dir) {
  
  meta <- fread("./expfiles/meta.csv")
  
  label_i <- unlist(strsplit(dir, split = "\\/"))[4]
  lab_i <- unlist(strsplit(label_i, split = "_"))[1]
  pattern_i <- unlist(strsplit(label_i, split = "_"))[2]
  meta_i <- meta %>%
    filter(lab %in% lab_i) %>%
    filter(pattern %in% pattern_i) %>%
    mutate(experiment = sapply(entry1, function(x) gsub("LFQ intensity ", "", x)))
  
  if (grepl("BGI_DIA", dir)) {
    meta_i <- meta_i %>%
      mutate(experiment = as.numeric(factor(experiment)))
  } else if (grepl("APT|FDU_DDA", dir)) {
    meta_i <- meta_i %>%
      mutate(experiment = as.numeric(experiment))
  }
  
  evi_file <- paste(dir, "/evidence.txt", sep = "")
  txt_evidence <- fread(evi_file, showProgress = FALSE)
  
  output_list <- list(meta_i, txt_evidence)
  
  return(output_list)
}
preprocess_evidence <- function(df_evidence, df_meta) {
  
  df_meta <- df_meta %>%
    select(library, batch, experiment)
  
  df_evidence <- df_evidence %>%
    rename_with(., ~ tolower(gsub(" ", "_", .x, fixed = TRUE)))
  
  df_mapping <- df_evidence %>%
    select(protein_group_ids:peptide_id, sequence, proteins:protein_names) %>%
    mutate(batch = df_meta$batch[1]) %>%
    distinct()
  
  df_sum <- df_evidence %>%
    group_by(sequence, experiment) %>%
    dplyr::summarise(intensity = sum(intensity, na.rm = TRUE), .groups = "keep") %>%
    ungroup() %>%
    left_join(., df_meta, by = "experiment") %>%
    select(sequence, library, intensity)
  
  df_final <- list(df_mapping, df_sum)
  
  return(df_final)
}

## inputs
all_files <- list.files("./rawfiles/MaxQuant", full.names = T)
all_dirs <-  all_files[!grepl("\\/raw$", all_files)]

meta <- fread("./expfiles/meta.csv")

## evidence.txt
cleaned_lists <- pblapply(all_dirs, function(dir) {
  
  data_list <- input_mq_files(dir)
  
  meta_i <- data_list[[1]]
  txt_evidence <- data_list[[2]]
  
  cleaned_data_list <- preprocess_evidence(txt_evidence, meta_i)
  
  return(cleaned_data_list)
  
})

df_mapping_final <- cleaned_lists %>%
  sapply(., function(x) x[1]) %>%
  rbindlist(.)

df_sum_final <- cleaned_lists %>%
  sapply(., function(x) x[2]) %>%
  rbindlist(.) %>%
  reshape2::dcast(., sequence ~ library, value.var = "intensity") %>%
  mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))

## save
fwrite(df_sum_final, "./expfiles/peptide/expdata_intensity.csv")
fwrite(df_mapping_final, "./expfiles/peptide/pep2prot.csv")

df_log_final <- df_sum_final %>%
  mutate_if(is.numeric, ~ ifelse(. == 0, NA, log(.)))

fwrite(df_log_final, "./expfiles/peptide/expdata_log.csv")


## functions for correction and quantification. --------------------
rm(list = ls())
ratio_by_d6 <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "library"
  batches <- unique(df_meta$batch)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    match_batch <- which(grepl(batch, samples))
    match_d6 <- match_batch[which(grepl("D6", samples[match_batch]))]
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(D6mean = apply(.[match_d6], 1, mean)) %>%
      mutate_at(match_batch, ~ (. - D6mean)) %>%
      select(!D6mean)
    
  }
  
  return(df_corrected)
}

median_centering <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "library"
  batches <- unique(df_meta$batch)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    match_batch <- which(grepl(batch, samples))
    
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

mean_centering <- function(df_quant_wide, df_meta) {
  
  colnames(df_meta)[1] <- "library"
  batches <- unique(df_meta$batch)
  samples <- colnames(df_quant_wide)
  
  df_corrected <- df_quant_wide
  
  for (batch in batches) {
    
    match_batch <- which(grepl(batch, samples))
    
    df_corrected <- df_corrected %>%
      as.data.frame() %>%
      mutate_at(match_batch, ~ (.) * 10 ^ (6)/ sum(.)) %>%
      mutate_at(match_batch, ~ ifelse(. == 0, NA, log(.))) %>%
      mutate(all_mean = apply(.[match_batch], 1, mean, na.rm = TRUE)) %>%
      mutate_at(match_batch, ~ (. - all_mean)) %>%
      select(!all_mean)
    
  }
  
  return(df_corrected)
}

combat_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  data_mode <- model.matrix(~ as.factor(sample), data = df_meta)
  
  df_corrected <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0, log(.))) %>%
    as.matrix() %>%
    ComBat(., df_meta$batch, data_mode) %>%
    as.data.frame() %>%
    mutate_all(~ na_if(., 0)) %>%
    rownames_to_column(feature_header)
  
  return(df_corrected)
}

sva_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  
  df_log_wide <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0.01, log(.))) %>%
    as.matrix()
  
  df_meta <- df_meta %>%
    column_to_rownames("library") %>%
    select(sample, batch)
  
  data_mode <- model.matrix(~ as.factor(sample), data = df_meta)
  data_mode0 <- model.matrix(~ 1, data = df_meta)
  
  sv_num <- num.sv(df_log_wide, data_mode, method = "be")
  sv_obj <- sva(df_log_wide, data_mode, data_mode0, n.sv = sv_num)
  
  sv_mode <- cbind(data_mode, sv_obj$sv)
  hat <- solve(t(sv_mode) %*% sv_mode) %*% t(sv_mode)
  beta <- (hat %*% t(df_log_wide))
  
  num_type <- ncol(data_mode)
  sva_be <- t((sv_mode[, -c(1:num_type)]) %*% beta[-c(1:num_type), ])
  sva_corrected <- (df_log_wide - sva_be) %>%
    as.data.frame() %>%
    rownames_to_column(feature_header)
  
  sva_corrected[df_quant_wide == 0] <- NA
  
  return(sva_corrected)
  
}

harmony_solve <- function(df_quant_wide, df_meta) {
  
  feature_header <- colnames(df_quant_wide)[1]
  
  df_log_wide <- df_quant_wide %>%
    column_to_rownames(feature_header) %>%
    mutate_all(~ (.) * 10 ^ (6)/ sum(.)) %>%
    mutate_all(~ ifelse(. == 0, 0, log(.))) %>%
    as.matrix()
  
  harmony_corrected <- df_log_wide %>%
    HarmonyMatrix(., df_meta, vars_use = "batch", do_pca = FALSE) %>%
    as.data.frame() %>%
    rownames_to_column(feature_header)
  
  harmony_corrected[df_quant_wide == 0] <- NA
  
  return(harmony_corrected)
  
}

ruv3c_solve <- function(df_quant_wide, df_meta, ctrl_features, k) {
  
  feature_header <- colnames(df_quant_wide)[1]
  data_mode <- model.matrix(~ sample - 1, data = df_meta)
  
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

correct_by_x <- function(df_quant_wide, df_meta,
                         method = c("ratio_by_d6",
                                    "median_centering",
                                    "mean_centering",
                                    "combat",
                                    "sva",
                                    "harmony",
                                    "RUV-III-C"),
                         ctrl_features = NULL,
                         k=3) {
  
  if (method == "ratio_by_d6") {
    df_corrected <- ratio_by_d6(df_quant_wide, df_meta)
  } else if (method == "median_centering") {
    df_corrected <- median_centering(df_quant_wide, df_meta)
  } else if (method == "mean_centering") {
    df_corrected <- mean_centering(df_quant_wide, df_meta)
  } else if (method == "combat") {
    df_corrected <- combat_solve(df_quant_wide, df_meta)
  } else if(method == "sva") {
    df_corrected <- sva_solve(df_quant_wide, df_meta)
  } else if(method == "harmony") {
    df_corrected <- harmony_solve(df_quant_wide, df_meta)
  } else if(method == "RUV-III-C") {
    df_corrected <- ruv3c_solve(df_quant_wide, df_meta, ctrl_features, k)
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
    filter(protein_names != "") %>%
    separate_rows(protein_names, sep = ";") %>%
    distinct(sequence, protein_names) %>%
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
        dplyr::summarise(ibaq = sum(intensity)) %>%
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


## no batch effect correction. 约15分钟/10分钟/13分钟 --------------------------
## inputs
meta <- fread("./expfiles/meta.csv")
df_sum_final <- fread("./expfiles/peptide/expdata_intensity.csv")
df_mapping_final <- fread("./expfiles/peptide/pep2prot.csv")

# ## peptides: balanced design
# design_scenario <- "balanced"
# 
# df_pep <- df_sum_final

## peptides: counfounded design
design_scenario <- "confounded"

df_pep <- df_sum_final %>%
  select(sequence, contains("D6"),
         starts_with("DAT1_D5"), starts_with("DAT4_D5"),
         starts_with("DAT2_F7"), starts_with("DAT5_F7"),
         starts_with("DAT3_M8"), starts_with("DAT6_M8"))

# ## peptides: random design
# design_scenario <- "random"
# 
# set.seed(2023)
# index6 <- which(grepl("D6", colnames(df_sum_final)))
# index578 <- setdiff(2:73, index6)
# index <- c(sample(index578, 36), index6)
# 
# df_pep <- df_sum_final %>%
#   select(sequence, index)

## proteins
pro_tables <- pblapply(c("iBAQ", "TopPep3", "MaxLFQ"), function(method_id) {
  
  df_mapping_i <- df_mapping_final
  
  df_final <- pep2pro_by_x(df_quant = df_pep,
                           df_mapping = df_mapping_i,
                           method = method_id)
  
  return(df_final)
})

folder_names <- c("ibaq", "top3", "lfq")
save_pro <- pblapply(1:3, function(i) {
  
  df_uncorrected <- pro_tables[[i]]
  
  fwrite(df_uncorrected, file = paste("./expfiles/protein/",
                                      design_scenario,
                                      "/raw_mapped/",
                                      folder_names[i],
                                      "/expdata_log.csv",
                                      sep = ""))
  
  df_fot <- df_uncorrected %>%
    mutate_if(is.numeric, ~ ifelse(. <= -1e+06, NA, .)) %>%
    mutate_if(is.numeric, exp) %>%
    mutate_if(is.numeric, ~ (.) * 10 ^ (6) / sum(., na.rm = TRUE)) %>%
    mutate_if(is.numeric, log)
  
  fwrite(df_fot, file = paste("./expfiles/protein/",
                              design_scenario,
                              "/raw_mapped/",
                              folder_names[i],
                              "/expdata_log_fot.csv",
                              sep = ""))
  
})


## two-way ANOVA. 约2分钟 -----------
data_list <- c(list(df_sum_final), pro_tables)

pep2pro_methods <- c("", "iBAQ", "TopPep3", "MaxLFQ")
feature_types <- c("peptide", "protein", "protein", "protein")

anova_tables <- pblapply(1:4, function(i) {
  
  df_quant <- data_list[[i]]
  colnames(df_quant)[1] <- "feature"
  
  df_sum_i <- df_quant %>%
    mutate_if(is.numeric, ~ ifelse(. <= -1e+06 | . == 0, NA, .)) %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 1))
  
  all_features <- df_sum_i$feature
  
  anova_tables_i <- mclapply(all_features, function(feature_id) {
    
    df_anova_j <- df_sum_i %>%
      filter(feature %in% feature_id) %>%
      reshape2::melt(., id = 1, variable.name = "library",
                     value.name = "intensity") %>%
      tidyr::separate(., library, c("batch", "sample", "replicate")) %>%
      rstatix::anova_test(intensity ~ sample * batch)
    
    df_anova_j$feature <- feature_id
    
    return(df_anova_j)
    
  }, mc.cores = 3)
  
  df_anova_i <- anova_tables_i %>%
    rbindlist(.) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    mutate(type = feature_types[i])
  
  return(df_anova_i)
  
})

## combine & save
df_anova <- anova_tables %>%
  rbindlist(.) %>%
  select(type, quantitation, feature, everything())

fwrite(df_anova, file = paste("../results/tables/1_",
                              design_scenario, "_anova.csv",
                              sep = ""))


## peptide-level batch effect correction. 约2小时 --------------------
df_anova <- fread("./../results/tables/1_balanced_anova.csv")

sub_anova <- df_anova %>%
  filter(type %in% "peptide") %>%
  filter(p < .05 & Effect %in% "batch")

df_anova_filter <- df_anova %>%
  filter(type %in% "peptide") %>%
  filter(feature %in% sub_anova$feature) %>%
  group_by(feature) %>%
  dplyr::summarise(count = length(which(p > .05 & grepl("sample", Effect)))) %>%
  filter(count == 2)

pep2pro_methods <- c("iBAQ", "TopPep3", "MaxLFQ")
folder_names <- c("ibaq", "top3", "lfq")

correct_methods <- c("ratio_by_d6", "median_centering", "combat", "RUV-III-C")
label_names <- c("ratio", "median", "combat", "ruv")

## 约10分钟/5分钟/7分钟
correct_pep_tables <- pblapply(1:4, function(i) {
  
  if (correct_methods[i] %in% "RUV-III-C") {
  
    ctrl_peptides <- df_anova_filter %>%
      pull(feature) %>%
      unique
       
  } else {
    
    ctrl_peptides <- NULL
    
  }
  
  meta_i <- meta %>% filter(library %in% colnames(df_pep))
  df_i <- correct_by_x(df_pep, meta_i, correct_methods[i], ctrl_peptides, k=6)
  
  fwrite(df_i, file = paste("./expfiles/peptide/",
                            design_scenario, "/expdata_",
                            label_names[i], ".csv", sep = ""))
  
  return(df_i)
})

## 约40分钟/20分钟/30分钟
correct_pro_tables <- pblapply(1:3, function(i) {
  
  df_mapping_i <- df_mapping_final
  
  tables_i <- mclapply(1:4, function(j) {

    df_j <- correct_pep_tables[[j]] %>%
      mutate_if(is.numeric, exp)
    
    df_j_final <- pep2pro_by_x(df_quant = df_j,
                               df_mapping = df_mapping_i,
                               method = pep2pro_methods[i])
    
    fwrite(df_j_final, file = paste("./expfiles/protein/",
                                    design_scenario, "/final_mapped/",
                                    folder_names[i], "/expdata_",
                                    label_names[j], ".csv", sep = ""))
    
  }, mc.cores = 2)
})


## protein-level batch effect correction. 约30秒钟 --------------------
df_anova <- fread("./../results/tables/1_balanced_anova.csv")

sub_anova <- df_anova %>%
  filter(type %in% "protein") %>%
  filter(p < .05 & Effect %in% "batch")

df_anova_filter <- df_anova %>%
  filter(type %in% "protein") %>%
  filter(feature %in% sub_anova$feature) %>%
  group_by(quantitation, feature) %>%
  dplyr::summarise(count = length(which(p > .05 & grepl("sample", Effect)))) %>%
  filter(count == 2)

pep2pro_methods <- c("iBAQ", "TopPep3", "MaxLFQ")
folder_names <- c("ibaq", "top3", "lfq")

correct_methods <- c("ratio_by_d6", "median_centering", "combat", "RUV-III-C")
label_names <- c("ratio", "median", "combat", "ruv")

correct_pro_tables <- pblapply(1:3, function(i) {
  
  df_uncorrected <- fread(paste("./expfiles/protein/",
                                design_scenario, "/raw_mapped/",
                                folder_names[i], "/expdata_log.csv",
                                sep = ""))
  
  df_uncorrected <- df_uncorrected %>%
    dplyr::rename(feature = protein_names) %>%
    mutate_if(is.numeric, exp) %>%
    mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .))
  
  tables_i <- mclapply(1:4, function(j) {
    
    if (correct_methods[j] %in% "RUV-III-C") {
      
      noNA_proteins <- df_uncorrected %>%
        mutate_if(is.numeric, ~ ifelse(. == 0, NA, .)) %>%
        na.omit %>%
        pull(feature) %>%
        unique
      
      ctrl_proteins <- df_anova_filter %>%
        filter(quantitation %in% pep2pro_methods[i]) %>%
        pull(feature) %>%
        unique
      
      ctrl_proteins <- intersect(ctrl_proteins, noNA_proteins)
      
    } else {
      
      ctrl_proteins <- NULL
    }
    
    meta_i <- meta %>% filter(library %in% colnames(df_uncorrected))
    df_i <- correct_by_x(df_uncorrected, meta_i, correct_methods[j], ctrl_proteins, k=6)
    
    fwrite(df_i, file = paste("./expfiles/protein/",
                              design_scenario, "/raw_mapped/",
                              folder_names[i], "/expdata_",
                              label_names[j], ".csv", sep = ""))
    
    return(df_i)
  }, mc.cores = 4)
})
