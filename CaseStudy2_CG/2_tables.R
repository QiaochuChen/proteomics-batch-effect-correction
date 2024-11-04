## Created on Jun 18th, 2023
## Updated on Jul 5th, 2024
## Author: Qiaochu Chen

library(parallel)
library(pbapply)
library(dplyr)
library(tibble)
library(reshape2)
library(data.table)
library(readxl)
library(pROC)

rm(list = ls())


## input --------------------
setwd("~/Desktop/毕业课题/manuscript-batch/CaseStudy2_CG/")
meta <- read.csv("data/metadata.csv")

## peptide level
pep_files <- paste("data/expfiles/peptide",
                   "/expdata_",
                   c("log_fot", "combat", "median", "ratio", "ruv"),
                   ".csv",
                   sep = "")

## protein level
uniprot <- read.csv("../../database/uniprot/uniprot_20220823.csv")
pro_files <- paste("data/expfiles/protein",
                   rep(rep(c("/raw_mapped", "/final_mapped", "/raw_mapped"), c(1, 4, 4)), 3),
                   rep(c("/ibaq", "/lfq", "/top3"), each = 9),
                   "/expdata_",
                   rep(c("log_fot", rep(c("combat", "median", "ratio", "ruv"), 2)), 3),
                   ".csv",
                   sep = "")

data_list <- mclapply(pro_files, data.table::fread)


## global ----------------------
source("~/Desktop/毕业课题/utils/PCA.R")
source("~/Desktop/毕业课题/utils/DEP.R")

setwd("~/Desktop/毕业课题/manuscript-batch/CaseStudy2_CG/results/tables/")

# ## 这19个IDs：基于log2水平80%PCs(前28个PCs)后Cook距离大于四倍平均值
# ## 2023年4月初计算了先去除这19个样本后批次校正的数据。
# ## 2023年4月底计算全部1495个样本批次校正的数据。
# rm_ids <- c("ExpB19", "ExpC86", "ExpD43", "ExpE6", "ExpF65", "ExpE30",
#             "ExpE31", "ExpG52", "ExpI14", "ExpI49", "ExpI75", "ExpI85",
#             "ExpI93", "ExpJ1", "ExpK3", "ExpH34", "ExpL85", "ExpP21", "ExpQ42")

## PCA检出异常
rm_ids <- c("ExpF25", "ExpB19")

data_levels <- rep(c("Uncorrected", rep("Peptide-corrected", 4), rep("Protein-corrected", 4)), 3)

pep2pro_methods <- rep(c("iBAQ", "MaxLFQ", "TopPep3"), each = 9)

correct_methods <- rep(c("Log transformed", rep(c("ComBat", "Median centering", "Ratio", "RUV-III-C"), 2)), 3)


## outliers: Cook's distance --------------------
expr <- fread("../../data/expfiles/peptide/expdata_log.csv")
colnames(expr)[1] <- "Feature"

expr[is.na(expr)] <- 0
na_cutoff <- 1 * (ncol(expr) - 1)

exprdata_t <- expr %>%
  column_to_rownames("Feature") %>%
  filter(apply(., 1, function(x) length(which(x != 0)) >= na_cutoff)) %>%
  t

dim(exprdata_t)

metadata <- meta %>%
    column_to_rownames("ID") %>%
    mutate(Plate = sapply(Tube, function(x) unlist(strsplit(x, split = "[0-9]"))[1])) %>%
    dplyr::select(Type, Project, Site, Plate, Cleaning, Date, Type, Sex, Week)

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        group = "Cleaning", center = TRUE, scale = TRUE,
                        biplot = FALSE, plot = FALSE)

colnames(pca_results$pcs_props)[which(pca_results$pcs_props[3, ] > .8)][1]

sub_pca <- pca_results$pcs_values %>%
  dplyr::select(library, any_of(sapply(1:28, function(x) paste("PC", x, sep = "")))) %>%
  left_join(., meta, by = c("library" = "ID")) %>%
  dplyr::select(Order, starts_with("PC"))

expr_mode <- lm(Order ~ ., data = sub_pca)
vect_cook <- cooks.distance(expr_mode)
y_breaks <- 4 * mean(vect_cook)

sub_cook <- data.frame(cookd = vect_cook, Order = names(vect_cook)) %>%
  mutate(Order = as.integer(Order)) %>%
  left_join(., meta, by = "Order") %>%
  mutate(Outlier = sapply(cookd, function(x) ifelse(x > y_breaks, "Yes", "No"))) %>%
  mutate(Outlier = factor(Outlier, levels = c("Yes", "No")))

length(sub_cook$Outlier[sub_cook$Outlier == "Yes"])

data.table::fwrite(sub_cook, "1_cook.csv")


## CV: qc samples -----------------
calculate_cv <- function(df_quant_wide, pattern = "PM", df_meta = meta) {
  
  match_pat <- which(grepl(pattern, df_meta$Type))
  
  cv_values <- apply(df_quant_wide, 1, function(x) {
    
    target_x <- as.numeric(x[match_pat])
    target_x <- na.omit(target_x)
    
    if (length(target_x) >= 3) {
      cv <- sd(target_x) / mean(target_x) * 100
    } else{
      cv <- NA
    }
    
    return(cv)
  })
  
  return(cv_values)
}

cv_tables <- mclapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  expr <- expr %>%
    select(feature, all_of(meta$ID))
  
  sub_cv <- expr %>%
    mutate_if(is.numeric, exp) %>%
    mutate(cv_pm = calculate_cv(., "PM")) %>%
    mutate(cv_p10 = calculate_cv(., "P10")) %>%
    mutate(cv_p11 = calculate_cv(., "P11")) %>%
    mutate(level = data_levels[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    select(level, correction, quantitation, feature, starts_with("cv")) %>%
    filter(apply(., 1, function(x) length(which(is.na(x))) <= 1))
  
  return(sub_cv)
  
}, mc.cores = 3)

sub_cv_final <- cv_tables %>%
  rbindlist(.) %>%
  mutate_all(~ ifelse(is.infinite(.), NA, .))

fwrite(sub_cv_final, "./1_cv.csv")


## Rank --------------------------
rank_tables <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  sub_rank <- expr %>%
    tibble::column_to_rownames("feature") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 1495)) %>%
    mutate(mean = apply(., 1, mean, na.rm = TRUE)) %>%
    arrange(desc(mean)) %>%
    mutate(rank = 1:nrow(.)) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    tibble::rownames_to_column("feature") %>%
    select(feature, level, correction, quantitation, mean, rank)
  
  return(sub_rank)
})

sub_rank_final <- rbindlist(rank_tables)

fwrite(sub_rank_final, "./1_rank.csv")


## Count --------------------
pep2prot_dt <- fread("../../data/mapped.csv")
cv_dt <- fread("./1_cv.csv")

colnames(pep2prot_dt) <- c("sequence", "protein_names")

## 约1分钟
count_tables <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "Feature"
  
  colnm <- colnames(expr)[2:ncol(expr)]
  count <- apply(expr[, 2:ncol(expr)], 2, function(x) length(x[!is.na(x)]))
  
  sub_count <- expr %>%
    reshape2::melt(., id = 1, variable.name = "ID", na.rm = TRUE) %>%
    left_join(., meta, by = "ID") %>%
    group_by(Type) %>%
    summarise(count = length(unique(Feature))) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(sub_count)
})

sub_count_final <- count_tables %>%
  rbindlist(.)

sub_count_filter <- cv_dt %>%
  group_by(level, correction, quantitation) %>%
  summarise_at(vars(cv_pm:cv_p11), ~ length(which(. < 30))) %>%
  reshape2::melt(., id = 1:3, variable.name = "Type", value.name = "filter_count") %>%
  mutate_at("Type", ~ paste("QC sample (", toupper(gsub("cv_", "", .)), ")", sep = ""))

sub_count_total <- pep2prot_dt %>%
  filter(protein_names != "") %>%
  tidyr::separate_rows(protein_names, sep = ";") %>%
  summarise(total_count = length(unique(protein_names)))

sub_count_summary <- sub_count_final %>%
  right_join(., sub_count_filter, by = c("Type", "level", "correction", "quantitation")) %>%
  mutate(total_count = sub_count_total$total_count) %>%
  mutate(`Missed` = (total_count - count) / total_count * 100,
         `Unqualified` = (count - filter_count) / total_count * 100,
         `Qualified` = filter_count / total_count * 100) %>%
  mutate_if(is.numeric, ~ round(., digits = 2)) %>%
  select(Type, level, correction, quantitation, everything())

fwrite(sub_count_summary, "./1_count.csv")


## UMAP: all samples ------------------
umap_objs <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  metadata <- meta %>%
    filter(!Type %in% "QC sample (PM)") %>%
    filter(!ID %in% rm_ids) %>%
    column_to_rownames("ID")
  
  exprdata_t <- expr %>%
    column_to_rownames("feature") %>%
    select(any_of(rownames(metadata))) %>%
    na.omit %>%
    t
  
  ## UMAP分析样本
  umap_results <- main_umap(exprdata_t = exprdata_t, metadata = metadata,
                            group = "Type", center = TRUE, scale = TRUE,
                            pct_threshold = .6, plot = FALSE)
  
  return(umap_results)
  
})

umap_tables <- pblapply(1:27, function(i) {
  
  umap_results <- umap_objs[[i]]
  
  sub_umap_values <- umap_results$table %>%
    dplyr::rename(UMAP1 = V1, UMAP2 = V2) %>%
    select(any_of(c("sample", "library", colnames(meta), "UMAP1", "UMAP2"))) %>%
    mutate(level = data_levels[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    mutate(correction = correct_methods[i])
  
  return(sub_umap_values)
  
})

sub_umap_final <- rbindlist(umap_tables)

fwrite(sub_umap_final, "2_umap.csv")


## PCA: all samples ------------------
cv_dt <- fread("1_cv.csv")
sub_cv <- cv_dt %>%
  mutate(n = apply(., 1, function(x) length(which(as.numeric(x[5:7]) <= 30)))) %>%
  mutate(qualified = ifelse(n == 3, "Yes", "No")) %>%
  select(feature, level, correction, quantitation, qualified)

rank_dt <- fread("1_rank.csv")

## 按Qualified/rank分类
pca_objs <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  
  sub_cv_tmp <- sub_cv %>%
    filter(level %in% data_levels[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    left_join(., rank_dt, by = c("feature",  "level", "correction", "quantitation"))
  
  pca_objs_by_x <- list()
  
  for (all_id in c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ ")) {
    
    if (all_id %in% "Qualified") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(qualified %in% "Yes")
      
    } else if (all_id %in% "All") {
      
      sub_cv_x <- sub_cv_tmp
      
    } else if (all_id %in% "1 ~ 200") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(rank <= 200)
      
    } else if (all_id %in% "200 ~ 500") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(rank <= 500 & rank >200)
      
    } else if (all_id %in% "500 ~ ") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(rank < 500)
      
    }
    
    metadata <- meta %>%
      filter(!ID %in% rm_ids) %>%
      column_to_rownames("ID")
    
    exprdata_t <- expr %>%
      filter(feature %in% sub_cv_x$feature) %>%
      column_to_rownames("feature") %>%
      select(any_of(rownames(metadata))) %>%
      na.omit %>%
      t
    
    if (ncol(exprdata_t) > 1) {
      
      ## PCA分析QC样本
      pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "Cleaning",
                              biplot = FALSE, snr = FALSE, plot = FALSE)
      
      pca_results <- c(pca_results, list(i = i), list(all = all_id))
      
    } else {
      
      pca_results <- NULL
    }
    
    pca_objs_by_x <- c(pca_objs_by_x, list(pca_results))
  }
  
  return(pca_objs_by_x)
})

sub_pca_tables <- pblapply(1:27, function(i) {
  
  pca_objs_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results <- pca_objs[[i]]
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_pca_values_by_x <- pca_results_x$pcs_values %>%
        select(any_of(c("sample", "library", colnames(meta), "PC1", "PC2"))) %>%
        mutate(all = all_id) %>%
        mutate(level = data_levels[i]) %>%
        mutate(correction = correct_methods[i]) %>%
        mutate(quantitation = pep2pro_methods[i])
      
      return(sub_pca_values_by_x)
      
    }
  })
  
  sub_pca_values <- rbindlist(pca_objs_by_x)
  
  return(sub_pca_values)
  
})

sub_props_tables <- pblapply(1:27, function(i) {
  
  sub_props_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results <- pca_objs[[i]]
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_pca_props_by_x <- pca_results_x$pcs_props %>%
        t %>% as.data.frame %>%
        mutate(all = all_id) %>%
        mutate(level = data_levels[i]) %>%
        mutate(correction = correct_methods[i]) %>%
        mutate(quantitation = pep2pro_methods[i]) %>%
        tibble::rownames_to_column("PC")
      
      return(sub_pca_props_by_x)
    }
  })
  
  sub_props <- rbindlist(sub_props_by_x)
  
  return(sub_props)
  
})

sub_pca_final <- rbindlist(sub_pca_tables)
sub_pca_props_final <- rbindlist(sub_props_tables)

fwrite(sub_pca_final, "2_pca.csv")
fwrite(sub_pca_props_final, "2_pca_props.csv")


## PCA: QC samples ------------------
cv_dt <- fread("1_cv.csv")
sub_cv <- cv_dt %>%
  mutate(n = apply(., 1, function(x) length(which(as.numeric(x[5:7]) <= 30)))) %>%
  mutate(qualified = ifelse(n == 3, "Yes", "No")) %>%
  select(feature, level, correction, quantitation, qualified)

rank_dt <- fread("1_rank.csv")

# 约30分钟: 按Qualified/rank分类
pca_objs <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  
  sub_cv_tmp <- sub_cv %>%
    filter(level %in% data_levels[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    left_join(., rank_dt, by = c("feature",  "level", "correction", "quantitation"))
  
  pca_objs_by_x <- list()
  
  for (all_id in c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ ")) {
    
    if (all_id %in% "Qualified") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(qualified %in% "Yes")
      
    } else if (all_id %in% "All") {
      
      sub_cv_x <- sub_cv_tmp
      
    } else if (all_id %in% "1 ~ 200") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(rank <= 200)
      
    } else if (all_id %in% "200 ~ 500") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(rank <= 500 & rank >200)
      
    } else if (all_id %in% "500 ~ ") {
      
      sub_cv_x <- sub_cv_tmp %>%
        filter(rank < 500)
      
    }
    
    metadata <- meta %>%
      filter(Type %in% c("QC sample (PM)", "QC sample (P10)", "QC sample (P11)")) %>%
      filter(!ID %in% rm_ids) %>%
      column_to_rownames("ID")
    
    exprdata_t <- expr %>%
      filter(feature %in% sub_cv_x$feature) %>%
      column_to_rownames("feature") %>%
      select(any_of(rownames(metadata))) %>%
      na.omit %>%
      t
    
    if (ncol(exprdata_t) > 1) {
      
      ## PCA分析QC样本
      pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "Type",
                              biplot = FALSE, snr = TRUE, plot = FALSE)
      
      pca_results <- c(pca_results, list(i = i), list(all = all_id))
      
    } else {
      
      pca_results <- NULL
    }
    
    pca_objs_by_x <- c(pca_objs_by_x, list(pca_results))
  }
  
  return(pca_objs_by_x)
})

sub_pca_tables <- pblapply(1:27, function(i) {
  
  pca_objs_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results <- pca_objs[[i]]
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_pca_values_by_x <- pca_results_x$pcs_values %>%
        select(any_of(c("sample", "library", colnames(meta), "PC1", "PC2"))) %>%
        mutate(all = all_id) %>%
        mutate(level = data_levels[i]) %>%
        mutate(correction = correct_methods[i]) %>%
        mutate(quantitation = pep2pro_methods[i])
      
      return(sub_pca_values_by_x)
      
    }
  })
  
  sub_pca_values <- rbindlist(pca_objs_by_x)
  
  return(sub_pca_values)
  
})

sub_props_tables <- pblapply(1:27, function(i) {
  
  sub_props_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results <- pca_objs[[i]]
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_pca_props_by_x <- pca_results_x$pcs_props %>%
        t %>% as.data.frame %>%
        mutate(all = all_id) %>%
        mutate(level = data_levels[i]) %>%
        mutate(correction = correct_methods[i]) %>%
        mutate(quantitation = pep2pro_methods[i]) %>%
        tibble::rownames_to_column("PC")
      
      return(sub_pca_props_by_x)
    }
  })
  
  sub_props <- rbindlist(sub_props_by_x)
  
  return(sub_props)
  
})

sub_snr_tables <- pblapply(1:27, function(i) {
  
  pca_objs_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results <- pca_objs[[i]]
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_snr_by_x <- pca_results_x$snr_results %>%
        mutate(all = all_id) %>%
        mutate(level = data_levels[i]) %>%
        mutate(correction = correct_methods[i]) %>%
        mutate(quantitation = pep2pro_methods[i])
      
      return(sub_snr_by_x)
      
    }
  })
  
  sub_snr <- rbindlist(pca_objs_by_x)
  
  return(sub_snr)
  
})

sub_pca_qc_final <- rbindlist(sub_pca_tables)
sub_pca_qc_props_final <- rbindlist(sub_props_tables)
sub_snr_final <- rbindlist(sub_snr_tables)

fwrite(sub_snr_final, "3_pca_snr.csv")
fwrite(sub_pca_qc_final, "3_pca_qc.csv")
fwrite(sub_pca_qc_props_final, "3_pca_qc_props.csv")


## Batch related info: weighted PCs Spearman correlated with the batch ------
pca_dt <- fread("../tables/2_pca.csv")
pca_props_dt <- fread("../tables/2_pca_props.csv")

brinfo_list <- pblapply(1:27, function(i) {
  
  pcs_props <- pca_props_dt %>%
    filter(PC %in% c("PC1", "PC2")) %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i])
  
  all_ids <- unique(pcs_props$all)
  
  brinfo_list_by_x <- mclapply(all_ids, function(all_id) {
    
    pcs_props_x <- pcs_props %>%
      filter(all %in% all_id) %>%
      pull(`Proportion of Variance`)
    
    sp_cors <- pca_dt %>%
      filter(all %in% all_id) %>%
      filter(level %in% data_levels[i]) %>%
      filter(quantitation %in% pep2pro_methods[i]) %>%
      filter(correction %in% correct_methods[i]) %>%
      select(Order, PC1, PC2) %>%
      cor(., method = "spearman") %>%
      .[1, 2:3]
    
    brinfo <- sum(pcs_props_x * abs(sp_cors))
    
    sub_brinfo_x <- data.frame(all = all_id,
                               level = data_levels[i],
                               quantitation = pep2pro_methods[i],
                               correction = correct_methods[i],
                               brinfo = brinfo)
    
    return(sub_brinfo_x)
    
  })
  
  sub_brinfo_tmp <- rbindlist(brinfo_list_by_x)
  
  return(sub_brinfo_tmp)
})

sum_brinfo <- rbindlist(brinfo_list)

fwrite(sum_brinfo, "../tables/2_brinfo.csv")


## PVCA ----------------
## 约40分钟
sub_pvca_tables <- pblapply(1:27, function(i) {

  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"

  expr[is.na(expr)] <- 0
  na_cutoff <- 1 * (ncol(expr) - 1)

  exprdata_t <- expr %>%
    column_to_rownames("feature") %>%
    filter(apply(., 1, function(x) length(which(x != 0)) >= na_cutoff)) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    t

  metadata <- meta %>%
    filter(!ID %in% rm_ids) %>%
    column_to_rownames("ID") %>%
    mutate(Plate = sapply(Tube, function(x) unlist(strsplit(x, split = "[0-9]"))[1])) %>%
    dplyr::select(Type, Project, Site, Plate, Cleaning, Date, Sex, Week)

  exprdata_t <- exprdata_t[rownames(metadata), ]

  if (ncol(exprdata_t) > 0) {
    
    pvca_results <- main_pvca(exprdata_t = exprdata_t, metadata = metadata, plot = FALSE)
    
    sub_pvca <- pvca_results$table %>%
      mutate(level = data_levels[i]) %>%
      mutate(quantitation = pep2pro_methods[i]) %>%
      mutate(correction = correct_methods[i])
    
    return(sub_pvca)
  }
})

sub_pvca_final <- rbindlist(sub_pvca_tables)

fwrite(sub_pvca_final, "3_pvca.csv")


## DEPs of study samples at baseline ----------------
## 约6分钟: Training set & Validation set
dep_tables <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  dep_tables_i <- mclapply(1:2, function(j) {
    
    metadata <- meta %>%
      filter(Week %in% "Baseline") %>%
      filter(ID %in% colnames(expr)) %>%
      filter(Training.Validation == j) %>%
      mutate_at("Sex", ~ ifelse(. %in% "F", "Female", "Male"))
    
    exprdata <- expr %>%
      column_to_rownames("feature") %>%
      select(all_of(metadata$ID)) %>%
      mutate_all(~ ifelse(is.na(.), 0, 2 ^ (.))) %>%
      filter(apply(., 1, function(x) length(which(x != 0)) > 0))
    
    group <- factor(metadata$Sex, ordered = TRUE)
    
    results_tmp <- main_dep(exprdata, group,
                            test_method = "t",
                            plot_volcano = FALSE,
                            silent = TRUE)
    
    results1 <- results_tmp$data %>%
      dplyr::mutate(Training.Validation = j)
    
    return(results1)
  })
  
  results_all <- dep_tables_i %>%
    data.table::rbindlist(.) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(results_all)
})

sub_deps_baseline <- data.table::rbindlist(dep_tables)

write.csv(sub_deps_baseline, "4_deps_sex_baseline_perTV.csv", row.names = FALSE)

## 约6分钟
dep_tables <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  metadata <- meta %>%
    filter(Week %in% "Baseline") %>%
    filter(ID %in% colnames(expr)) %>%
    mutate_at("Sex", ~ ifelse(. %in% "F", "Female", "Male"))
  
  exprdata <- expr %>%
    column_to_rownames("feature") %>%
    select(all_of(metadata$ID)) %>%
    mutate_all(~ ifelse(is.na(.), 0, 2 ^ (.))) %>%
    filter(apply(., 1, function(x) length(which(x != 0)) > 0))
  
  group <- factor(metadata$Sex, ordered = TRUE)
  
  results_tmp <- main_dep(exprdata, group,
                          test_method = "t",
                          plot_volcano = FALSE,
                          silent = TRUE)
  
  results_all <- results_tmp$data %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(results_all)
})

sub_deps_baseline <- data.table::rbindlist(dep_tables)

write.csv(sub_deps_baseline, "4_deps_sex_baseline.csv", row.names = FALSE)


## Reference dataset of Sex DEPs ----------------
dep_tables <- pblapply(c(1, 10, 19), function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  dep_tables_i <- mclapply(c("I", "II", "III"), function(j) {
    
    metadata <- meta %>%
      filter(Week %in% "Baseline") %>%
      filter(Cleaning %in% j) %>%
      filter(ID %in% colnames(expr)) %>%
      mutate_at("Sex", ~ ifelse(. %in% "F", "Female", "Male"))
    
    exprdata <- expr %>%
      column_to_rownames("feature") %>%
      select(all_of(metadata$ID)) %>%
      mutate_all(~ ifelse(is.na(.), 0, 2 ^ (.))) %>%
      filter(apply(., 1, function(x) length(which(x != 0)) > 0))
    
    group <- factor(metadata$Sex, ordered = TRUE)
    
    results_tmp <- main_dep(exprdata, group,
                            test_method = "t",
                            plot_volcano = FALSE,
                            silent = TRUE)
    
    results1 <- results_tmp$data %>%
      dplyr::mutate(Cleaning = j)
    
    return(results1)
  })
  
  results_all <- dep_tables_i %>%
    data.table::rbindlist(.) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(results_all)
})

sub_deps_baseline_perbatch <- data.table::rbindlist(dep_tables)

write.csv(sub_deps_baseline_perbatch, "4_deps_sex_baseline_perbatch.csv", row.names = FALSE)

ref_deps <- sub_deps_baseline_perbatch %>%
  filter(p.adj < .05) %>%
  reshape2::dcast(., feature + Cleaning ~ quantitation, value.var = "estimate") %>%
  mutate(direction = apply(., 1, function(x) {
    if (length(which(as.numeric(x[3:5]) > 0)) == 3) {
      a <- "up"
    } else if (length(which(as.numeric(x[3:5]) < 0)) == 3) {
      a <- "down"
    } else {
      a <- "unknown"
    }
    return(a)
  })) %>%
  filter(direction %in% c("up", "down")) %>%
  reshape2::dcast(., feature ~ Cleaning, value.var = "direction") %>%
  mutate(vote = apply(., 1, function(x) 3 - sum(is.na(x)))) %>%
  filter(vote >= 1) %>%
  pull(feature) %>%
  unique

sub_ref <- sub_deps_baseline_perbatch %>%
  filter(feature %in% ref_deps) %>%
  group_by(feature) %>%
  summarise(log2FC.ref = mean(estimate)) %>%
  mutate(omic = "Protein", pair = "Female/Male") %>%
  mutate(change.ref = ifelse(log2FC.ref > 0, "up", "down")) %>%
  select(omic, feature, pair, log2FC.ref, change.ref)

fwrite(sub_ref, "4_deps_sex_baseline_ref.csv")


## RC & RMSE: Sex ------------------
ref_dt <- fread("4_deps_sex_baseline_ref.csv")
dep_dt <- fread("4_deps_sex_baseline.csv")

sub_ref <- ref_dt %>%
  select(!omic)

sub_test <- dep_dt %>%
  filter(feature %in% sub_ref$feature) %>%
  filter(p.adj < .05) %>%
  mutate(pair = paste(group1, group2, sep = "/"))

cv_dt <- fread("1_cv.csv")
sub_cv <- cv_dt %>%
  mutate(n = apply(., 1, function(x) length(which(as.numeric(x[5:7]) <= 30)))) %>%
  mutate(qualified = ifelse(n == 3, "Yes", "No")) %>%
  select(feature, level, correction, quantitation, qualified)

rank_dt <- fread("1_rank.csv")

sub_test <- sub_test %>%
  left_join(., sub_cv, by = c("feature", "level", "correction", "quantitation")) %>%
  left_join(., rank_dt, by = c("feature", "level", "correction", "quantitation"))

metrics_tables1 <- pblapply(1:27, function(i) {
  
  sub_test_i <- sub_test %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    select(feature, pair, qualified, estimate)
  
  sub_comb_i <- sub_test_i %>%
    right_join(., sub_ref, by = c("feature", "pair")) %>%
    mutate(class = apply(., 1, function(x) {
      if (is.na(x[4])) label <- "FN"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) > 0)) label <- "TP"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) <= 0)) label <- "FP"
      return(label)
    }))
  
  metrics_by_pair <- mclapply(c("All", "Qualified"), function(all_id) {
    
    if (all_id %in% "Qualified") {
      sub_comb_x <- sub_comb_i %>%
        filter(qualified %in% "Yes")
    } else {
      sub_comb_x <- sub_comb_i
    }
    
    sub_stat_i <- sub_comb_x %>%
      group_by(pair, class) %>%
      summarise(count = length(feature), .groups = "keep")
    
    metrics_by_pair_x <- mclapply(c("Female/Male"), function(pair_id) {
      
      sub_comb_i_j <- sub_comb_x %>%
        filter(pair %in% pair_id) %>%
        na.omit
      
      if (nrow(sub_comb_i_j) >= 3) {
        cor_list <- cor.test(formula = ~ log2FC.ref + estimate,
                             data = sub_comb_i_j,
                             method = "pearson")
        rc_value <- ifelse(cor_list$p.value < 0.05, cor_list$estimate, NA)
        
        rmse_value <- sum((sub_comb_i_j$log2FC.ref - sub_comb_i_j$estimate) ^ 2) / nrow(sub_comb_i_j)
        rmse_value <- sqrt(rmse_value)
        
      } else {
        rc_value <- NA
        rmse_value <- NA
        
      }
      
      metrics_dt_j <- data.frame(all = ifelse(all_id %in% "Qualified", "Qualified", "All"),
                                 pair = pair_id, rc = rc_value, rmse = rmse_value)
      
      return(metrics_dt_j)
      
    }, mc.cores = 3)
    
    metrics_dt_x <- metrics_by_pair_x %>%
      rbindlist(.) %>%
      left_join(., sub_stat_i, by = "pair") 
    
    return(metrics_dt_x)
  })
  
  metrics_dt_i <- metrics_by_pair %>%
    rbindlist(.) %>%
    mutate(scenario = "Random design") %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    select(scenario, level, correction, quantitation, all, everything())
  
  return(metrics_dt_i)
})

metrics_tables2 <- pblapply(1:27, function(i) {
  
  sub_test_i <- sub_test %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    select(feature, pair, rank, estimate)
  
  sub_comb_i <- sub_test_i %>%
    right_join(., sub_ref, by = c("feature", "pair")) %>%
    mutate(class = apply(., 1, function(x) {
      if (is.na(x[4])) label <- "FN"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) > 0)) label <- "TP"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) <= 0)) label <- "FP"
      return(label)
    }))
  
  metrics_by_pair <- mclapply(c("1 ~ 200", "200 ~ 500", "500 ~"), function(all_id) {
    
    if (all_id %in% "1 ~ 200") {
      sub_comb_x <- sub_comb_i %>%
        filter(rank <= 200)
      
    } else if (all_id %in% "200 ~ 500") {
      sub_comb_x <- sub_comb_i %>%
        filter(rank > 200 & rank <= 500)
      
    } else {
      sub_comb_x <- sub_comb_i %>%
        filter(rank > 500)
      
    }
    
    sub_stat_i <- sub_comb_x %>%
      group_by(pair, class) %>%
      summarise(count = length(feature), .groups = "keep")
    
    metrics_by_pair_x <- mclapply(c("Female/Male"), function(pair_id) {
      
      sub_comb_i_j <- sub_comb_x %>%
        filter(pair %in% pair_id) %>%
        na.omit
      
      if (nrow(sub_comb_i_j) >= 3) {
        cor_list <- cor.test(formula = ~ log2FC.ref + estimate,
                             data = sub_comb_i_j,
                             method = "pearson")
        rc_value <- ifelse(cor_list$p.value < 0.05, cor_list$estimate, NA)
        
        rmse_value <- sum((sub_comb_i_j$log2FC.ref - sub_comb_i_j$estimate) ^ 2) / nrow(sub_comb_i_j)
        rmse_value <- sqrt(rmse_value)
        
      } else {
        rc_value <- NA
        rmse_value <- NA
        
      }
      
      metrics_dt_j <- data.frame(pair = pair_id, rc = rc_value, rmse = rmse_value)
      
      return(metrics_dt_j)
      
    }, mc.cores = 3)
    
    metrics_dt_x <- metrics_by_pair_x %>%
      rbindlist(.) %>%
      mutate(all = all_id) %>%
      left_join(., sub_stat_i, by = "pair") 
    
    return(metrics_dt_x)
  })
  
  metrics_dt_i <- metrics_by_pair %>%
    rbindlist(.) %>%
    mutate(scenario = "Random design") %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    select(scenario, level, correction, quantitation, all, everything())
  
  return(metrics_dt_i)
})

metrics_dt <- metrics_tables1 %>%
  c(., metrics_tables2) %>%
  rbindlist(.) %>%
  reshape2::dcast(., scenario + level + correction + quantitation + pair + all + rc + rmse ~ class, value.var = "count") %>%
  select(!`NA`) %>%
  mutate_at(9:11, ~ ifelse(is.na(.), 0, .)) %>%
  mutate(precision = TP / (TP + FP)) %>%
  mutate(recall = TP / (TP + FN)) %>%
  mutate(f1 = 2 * TP / (2 * TP + FP + FN)) %>%
  mutate_at(12:14, ~ ifelse(is.nan(.), NA, .))

data.table::fwrite(metrics_dt, "4_with_reference_dataset.csv")


## DEPs of QC samples  ----------------
## 约3分钟
dep_tables <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  metadata <- meta %>%
    filter(Type %in% c("QC sample (P10)", "QC sample (P11)")) %>%
    filter(ID %in% colnames(expr)) %>%
    mutate_at("Type", ~ str_extract(., "(?<=\\().+(?=\\))"))
  
  exprdata <- expr %>%
    column_to_rownames("feature") %>%
    select(all_of(metadata$ID)) %>%
    mutate_all(~ ifelse(is.na(.), 0, 2 ^ (.))) %>%
    filter(apply(., 1, function(x) length(which(x != 0)) > 0))
  
  group <- factor(metadata$Type, ordered = TRUE)
  
  results_tmp <- main_dep(exprdata, group,
                          test_method = "t",
                          plot_volcano = FALSE,
                          silent = TRUE)
  
  results_all <- results_tmp$data %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(results_all)
})

sub_deps_qc <- data.table::rbindlist(dep_tables)

write.csv(sub_deps_qc, "4_deps_qc.csv", row.names = FALSE)


## Reference dataset of P10-P11 DEPs ----------------
dep_tables <- pblapply(1:27, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  dep_tables_i <- mclapply(c("I", "II", "III"), function(j) {
    
    metadata <- meta %>%
      filter(Type %in% c("QC sample (P10)", "QC sample (P11)")) %>%
      filter(Cleaning %in% j) %>%
      filter(ID %in% colnames(expr)) %>%
      mutate_at("Type", ~ str_extract(., "(?<=\\().+(?=\\))"))
    
    exprdata <- expr %>%
      column_to_rownames("feature") %>%
      select(all_of(metadata$ID)) %>%
      mutate_all(~ ifelse(is.na(.), 0, 2 ^ (.))) %>%
      filter(apply(., 1, function(x) length(which(x != 0)) > 0))
    
    group <- factor(metadata$Type, ordered = TRUE)
    
    results_tmp <- main_dep(exprdata, group,
                            test_method = "t",
                            plot_volcano = FALSE,
                            silent = TRUE)
    
    results1 <- results_tmp$data %>%
      dplyr::mutate(Cleaning = j)
    
    return(results1)
  })
  
  results_all <- dep_tables_i %>%
    data.table::rbindlist(.) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(results_all)
})

sub_deps_qc_perbatch <- data.table::rbindlist(dep_tables)

write.csv(sub_deps_qc_perbatch, "4_deps_qc_perbatch.csv", row.names = FALSE)

sub_deps_qc_perbatch <- fread("4_deps_qc_perbatch.csv")

n1 <- sub_deps_qc_perbatch %>%
  filter(level %in% "Uncorrected") %>%
  filter(p < .05) %>%
  group_by(Cleaning) %>%
  summarise(length(unique(feature)))

n2 <- sub_deps_qc_perbatch %>%
  filter(level %in% "Uncorrected") %>%
  filter(p < .05, abs(estimate) > log2(1.2)) %>%
  group_by(Cleaning) %>%
  summarise(length(unique(feature)))

ref_deps <- sub_deps_qc_perbatch %>%
  filter(level %in% "Uncorrected") %>%
  filter(p < .05, abs(estimate) > log2(1.2)) %>%
  reshape2::dcast(., feature + level + correction + quantitation ~ Cleaning, value.var = "estimate") %>%
  mutate(direction = apply(., 1, function(x) {
    if (length(which(as.numeric(x[5:7]) > 0)) >= 3) {
      a <- "up"
    } else if (length(which(as.numeric(x[5:7]) < 0)) >= 3) {
      a <- "down"
    } else {
      a <- "unknown"
    }
    return(a)
  })) %>%
  filter(direction %in% c("up", "down")) %>%
  # reshape2::dcast(., feature ~ level + correction + quantitation, value.var = "direction") %>%
  # mutate(vote = apply(., 1, function(x) 3 - sum(is.na(x)))) %>%
  # filter(vote >= 3) %>%
  pull(feature) %>%
  unique

sub_ref <- sub_deps_qc_perbatch %>%
  filter(feature %in% ref_deps) %>%
  group_by(feature) %>%
  summarise(log2FC.ref = mean(estimate)) %>%
  mutate(omic = "Protein", pair = "P10/P11") %>%
  mutate(change.ref = ifelse(log2FC.ref > 0, "up", "down")) %>%
  select(omic, feature, pair, log2FC.ref, change.ref)

fwrite(sub_ref, "4_deps_qc_ref.csv")


## RC & RMSE: P10-P11 ------------------
ref_dt <- fread("4_deps_qc_ref.csv")
dep_dt <- fread("4_deps_qc.csv")

sub_ref <- ref_dt %>%
  select(!omic)

sub_test <- dep_dt %>%
  filter(feature %in% sub_ref$feature) %>%
  filter(p < .05, abs(estimate) > log2(1.2)) %>%
  mutate(pair = paste(group1, group2, sep = "/"))

cv_dt <- fread("1_cv.csv")
sub_cv <- cv_dt %>%
  mutate(n = apply(., 1, function(x) length(which(as.numeric(x[5:7]) <= 30)))) %>%
  mutate(qualified = ifelse(n == 3, "Yes", "No")) %>%
  select(feature, level, correction, quantitation, qualified)

rank_dt <- fread("1_rank.csv")

sub_test <- sub_test %>%
  left_join(., sub_cv, by = c("feature", "level", "correction", "quantitation")) %>%
  left_join(., rank_dt, by = c("feature", "level", "correction", "quantitation"))

metrics_tables1 <- pblapply(1:27, function(i) {
  
  sub_test_i <- sub_test %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    select(feature, pair, qualified, estimate)
  
  sub_comb_i <- sub_test_i %>%
    right_join(., sub_ref, by = c("feature", "pair")) %>%
    mutate(class = apply(., 1, function(x) {
      if (is.na(x[4])) label <- "FN"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) > 0)) label <- "TP"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) <= 0)) label <- "FP"
      return(label)
    }))
  
  metrics_by_pair <- mclapply(c("All", "Qualified"), function(all_id) {
    
    if (all_id %in% "Qualified") {
      sub_comb_x <- sub_comb_i %>%
        filter(qualified %in% "Yes")
    } else {
      sub_comb_x <- sub_comb_i
    }
    
    sub_stat_i <- sub_comb_x %>%
      group_by(pair, class) %>%
      summarise(count = length(feature), .groups = "keep")
    
    metrics_by_pair_x <- mclapply(c("P10/P11"), function(pair_id) {
      
      sub_comb_i_j <- sub_comb_x %>%
        filter(pair %in% pair_id) %>%
        na.omit
      
      if (nrow(sub_comb_i_j) >= 3) {
        cor_list <- cor.test(formula = ~ log2FC.ref + estimate,
                             data = sub_comb_i_j,
                             method = "pearson")
        rc_value <- ifelse(cor_list$p.value < 0.05, cor_list$estimate, NA)
        
        rmse_value <- sum((sub_comb_i_j$log2FC.ref - sub_comb_i_j$estimate) ^ 2) / nrow(sub_comb_i_j)
        rmse_value <- sqrt(rmse_value)
        
      } else {
        rc_value <- NA
        rmse_value <- NA
        
      }
      
      metrics_dt_j <- data.frame(all = ifelse(all_id %in% "Qualified", "Qualified", "All"),
                                 pair = pair_id, rc = rc_value, rmse = rmse_value)
      
      return(metrics_dt_j)
      
    }, mc.cores = 3)
    
    metrics_dt_x <- metrics_by_pair_x %>%
      rbindlist(.) %>%
      left_join(., sub_stat_i, by = "pair") 
    
    return(metrics_dt_x)
  })
  
  metrics_dt_i <- metrics_by_pair %>%
    rbindlist(.) %>%
    mutate(scenario = "Random design") %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    select(scenario, level, correction, quantitation, all, everything())
  
  return(metrics_dt_i)
})

metrics_tables2 <- pblapply(1:27, function(i) {
  
  sub_test_i <- sub_test %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    select(feature, pair, rank, estimate)
  
  sub_comb_i <- sub_test_i %>%
    right_join(., sub_ref, by = c("feature", "pair")) %>%
    mutate(class = apply(., 1, function(x) {
      if (is.na(x[4])) label <- "FN"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) > 0)) label <- "TP"
      if (!is.na(x[4]) & (as.numeric(x[4]) * as.numeric(x[5]) <= 0)) label <- "FP"
      return(label)
    }))
  
  metrics_by_pair <- mclapply(c("1 ~ 200", "200 ~ 500", "500 ~"), function(all_id) {
    
    if (all_id %in% "1 ~ 200") {
      sub_comb_x <- sub_comb_i %>%
        filter(rank <= 200)
      
    } else if (all_id %in% "200 ~ 500") {
      sub_comb_x <- sub_comb_i %>%
        filter(rank > 200 & rank <= 500)
      
    } else {
      sub_comb_x <- sub_comb_i %>%
        filter(rank > 500)
      
    }
    
    sub_stat_i <- sub_comb_x %>%
      group_by(pair, class) %>%
      summarise(count = length(feature), .groups = "keep")
    
    metrics_by_pair_x <- mclapply(c("P10/P11"), function(pair_id) {
      
      sub_comb_i_j <- sub_comb_x %>%
        filter(pair %in% pair_id) %>%
        na.omit
      
      if (nrow(sub_comb_i_j) >= 3) {
        cor_list <- cor.test(formula = ~ log2FC.ref + estimate,
                             data = sub_comb_i_j,
                             method = "pearson")
        rc_value <- ifelse(cor_list$p.value < 0.05, cor_list$estimate, NA)
        
        rmse_value <- sum((sub_comb_i_j$log2FC.ref - sub_comb_i_j$estimate) ^ 2) / nrow(sub_comb_i_j)
        rmse_value <- sqrt(rmse_value)
        
      } else {
        rc_value <- NA
        rmse_value <- NA
        
      }
      
      metrics_dt_j <- data.frame(pair = pair_id, rc = rc_value, rmse = rmse_value)
      
      return(metrics_dt_j)
      
    }, mc.cores = 3)
    
    metrics_dt_x <- metrics_by_pair_x %>%
      rbindlist(.) %>%
      mutate(all = all_id) %>%
      left_join(., sub_stat_i, by = "pair") 
    
    return(metrics_dt_x)
  })
  
  metrics_dt_i <- metrics_by_pair %>%
    rbindlist(.) %>%
    mutate(scenario = "Random design") %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    select(scenario, level, correction, quantitation, all, everything())
  
  return(metrics_dt_i)
})

metrics_dt <- metrics_tables1 %>%
  c(., metrics_tables2) %>%
  rbindlist(.) %>%
  reshape2::dcast(., scenario + level + correction + quantitation + pair + all + rc + rmse ~ class, value.var = "count") %>%
  select(!`NA`) %>%
  mutate_at(9:11, ~ ifelse(is.na(.), 0, .)) %>%
  mutate(precision = TP / (TP + FP)) %>%
  mutate(recall = TP / (TP + FN)) %>%
  mutate(f1 = 2 * TP / (2 * TP + FP + FN)) %>%
  mutate_at(12:14, ~ ifelse(is.nan(.), NA, .))

a <- metrics_dt %>%
  filter(all %in% "All")

data.table::fwrite(metrics_dt, "4_with_reference_dataset.csv")


## Prediction of Study samples: Sex -------------------
# ## get DEPs
# dep_table <- data.table::fread("4_deps_sex_baseline_perTV.csv")
# sex_deps <- dep_table %>%
#   filter(Training.Validation == 1) %>%
#   filter(p.adj < 0.05) %>%
#   pull(feature) %>%
#   unique

## get sex-related DEPs
dep_table <- read_excel("../../data/41591_2019_673_MOESM3_ESM.xlsx", sheet = 6, skip = 1)
sex_deps <- dep_table %>%
  mutate_at(2:12, as.numeric) %>%
  filter(q.Sex < 0.05) %>%
  mutate(Gene.Symbol = stringr::str_extract_all(ID, "^(.+)(?=\\.[0-9]{4}.+)")) %>%
  tidyr::separate_rows(Gene.Symbol, sep = "\\.") %>%
  left_join(., uniprot, by = "Gene.Symbol") %>%
  pull(Protein.ID) %>%
  unique

data_list_new <- c(data_list[c(1, 10, 19)], data_list)
data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

## set negative controls and filter by DEPs; missing values < 20%
subset_tables <- pblapply(1:30, function(i) {
  
  if (data_levels_new[i] %in% "Negative control") {
    
    set.seed(2023)
    meta_new <- meta %>%
      filter(Week %in% "Baseline") %>%
      mutate(Sex_random = sample(Sex, replace = TRUE)) %>%
      select(ID, Training.Validation, Sex_random) %>%
      dplyr::rename(Sex = Sex_random)
    
  } else {
    
    meta_new <- meta %>%
      filter(Week %in% "Baseline") %>%
      select(ID, Training.Validation, Sex)
    
  }
  
  expr <- data_list_new[[i]]
  colnames(expr)[1] <- "Accession"
  
  expr_new <- expr %>%
    select(Accession, all_of(meta_new$ID)) %>%
    # na.omit %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 1496 * .2)) %>%
    filter(Accession %in% sex_deps)
  
  expr_sub <- expr_new %>%
    reshape2::melt(., id = 1, variable.name = "ID", na.rm = TRUE) %>%
    right_join(., meta_new, by = "ID") %>%
    select(Training.Validation, ID, Sex, Accession, value) %>%
    mutate(level = data_levels_new[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(expr_sub)
})

## combine all
sub_wide <- subset_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Training.Validation + ID + Sex + level + correction + quantitation ~ Accession)

fwrite(sub_wide, "published_sex_deps_0.2_missing_values/5_subset_expr_sex_baseline.csv")

## after using JMP Pro: MCC
model_dt <- read_excel("./5_models_sex.xlsx")
sub_model <- model_dt %>%
  select(Training.Validation, level, correction, quantitation, Sex, ends_with("Predicted"))

data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

mcc_tables <- pblapply(1:30, function(i) {
  
  sub_tmp_i <- sub_model %>%
    filter(level %in% data_levels_new[i]) %>%
    filter(correction %in% correct_methods_new[i]) %>%
    filter(quantitation %in% pep2pro_methods_new[i])
  
  all_models <- c("Boosted Tree")
  
  class_tables <- mclapply(all_models, function(model_id) {
    
    sub_tmp_j <- sub_tmp_i %>%
      select(Training.Validation, level, correction, quantitation, Sex, sym(paste(model_id, "Predicted"))) %>%
      mutate(class = apply(., 1, function(x) {
        if(as.character(x[5]) %in% "F" & as.character(x[6]) %in% "F") {
          a <- "TP"
        } else if(as.character(x[5]) %in% "M" & as.character(x[6]) %in% "F") {
          a <- "FP"
        } else if(as.character(x[5]) %in% "F" & as.character(x[6]) %in% "M") {
          a <- "FN"
        } else if(as.character(x[5]) %in% "M" & as.character(x[6]) %in% "M") {
          a <- "TN"
        } else {
          a <- NA
        }
        return(a)
      })) %>%
      mutate(model = model_id) %>%
      na.omit
    
    return(sub_tmp_j)
  })
  
  sub_mcc_i <- class_tables %>%
    rbindlist(.) %>%
    reshape2::dcast(., Training.Validation + level + correction + quantitation + model ~ class,
                    value.var = "class", fun.aggregate = length) %>%
    mutate(mcc = (TP * TN - FP * FN) / sqrt((TP+FP) * (TP + FN) * (TN + FP) * (TN + FN)))
  
  
  return(sub_mcc_i)
})

sub_mcc <- rbindlist(mcc_tables)

fwrite(sub_mcc, "6_mcc.csv")

## after using JMP Pro: ROC
model_dt <- read_excel("./5_models_sex.xlsx")

sub_model <- model_dt %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(correction, sep = "_") %>%
  reshape2::dcast(., Training.Validation + Sex + ID + quantitation + correction ~ level, value.var = "Boosted Tree Prob F")

tr_va <- rep(1:2, 12)
pep2pro_methods_new <- rep(c("iBAQ", "MaxLFQ", "TopPep3"), each = 8)
correct_methods_new <- rep(rep(c("ComBat", "Median centering", "Ratio", "RUV-III-C"), each = 2), 3)

roc_results <- pblapply(1:24, function(i) {
  
  sub_tmp_i <- sub_model %>%
    filter(Training.Validation == tr_va[i]) %>%
    filter(correction %in% correct_methods_new[i]) %>%
    filter(quantitation %in% pep2pro_methods_new[i])
  
  roc_results_i <- roc(Sex ~ `Negative control` + `Uncorrected` + `Peptide-corrected` + `Protein-corrected`,
                       data = sub_tmp_i, aur = TRUE, ci = TRUE, percent = TRUE,
                       smooth = FALSE, levels = c("F", "M"), direction = ">")
  
  return(roc_results_i)
})

delong_tables <- pblapply(1:24, function(i) {
  
  roc_results_i <- roc_results[[i]]
  
  compared_lists <- combn(1:4, 2)
  
  roc_tests_i <- mclapply(1:6, function(j) {
    
    roc_test_j <- roc.test(roc_results_i[[compared_lists[1, j]]],
                           roc_results_i[[compared_lists[2, j]]])
    
    sub_test_j <- data.frame(level1 = names(roc_results_i)[compared_lists[1,j]],
                             level2 = names(roc_results_i)[compared_lists[2,j]],
                             p = roc_test_j$p.value)
    
    return(sub_test_j)
  })
  
  sub_roc_i <- roc_tests_i %>%
    rbindlist(.) %>%
    mutate(Training.Validation = tr_va[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(sub_roc_i)
})

roc_tables <- pblapply(1:24, function(i) {
  
  roc_results_i <- roc_results[[i]]
  
  roc_tables_i <- mclapply(1:4, function(j) {
    
    sub_roc_j <- data.frame(level = names(roc_results_i)[j],
                            threshold = roc_results_i[[j]]$thresholds,
                            sensitivity = roc_results_i[[j]]$sensitivities,
                            specificity = roc_results_i[[j]]$specificities)
    
    return(sub_roc_j)
  })
  
  sub_roc_i <- roc_tables_i %>%
    rbindlist(.) %>%
    mutate(Training.Validation = tr_va[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(sub_roc_i)
})

auc_tables <- pblapply(1:24, function(i) {
  
  roc_results_i <- roc_results[[i]]
  
  auc_tables_i <- mclapply(1:4, function(j) {
    
    auc <- roc_results_i[[j]]$auc
    sub_auc_j <- data.frame(level = names(roc_results_i)[j], auc = auc)
    
    return(sub_auc_j)
  })
  
  sub_auc_i <- auc_tables_i %>%
    rbindlist(.) %>%
    mutate(Training.Validation = tr_va[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(sub_auc_i)
})

sub_roc <- rbindlist(roc_tables)

sub_delong <- rbindlist(delong_tables)

sub_auc <- rbindlist(auc_tables)

saveRDS(list(AUC = sub_auc, DelongTest = sub_delong, ROC = sub_roc),
        "./6_roc_sex.rds")


## Prediction of Study samples: Age -------------------
data_list_new <- c(data_list[c(1, 10, 19)], data_list)
data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

## set negative controls and filter by missing values
subset_tables <- pblapply(1:30, function(i) {
  
  if (data_levels_new[i] %in% "Negative control") {
    
    set.seed(2023)
    meta_new <- meta %>%
      filter(Week %in% "Baseline") %>%
      mutate(Age_random = sample(Age, replace = FALSE)) %>%
      select(ID, Training.Validation, Age_random) %>%
      dplyr::rename(Age = Age_random)
    
  } else {
    
    meta_new <- meta %>%
      filter(Week %in% "Baseline") %>%
      select(ID, Training.Validation, Age)
    
  }
  
  expr <- data_list_new[[i]]
  colnames(expr)[1] <- "Accession"
  
  expr_new <- expr %>%
    select(Accession, all_of(meta_new$ID)) %>%
    # filter(Accession %in% sex_deps) %>%
    na.omit
  
  expr_sub <- expr_new %>%
    reshape2::melt(., id = 1, variable.name = "ID", na.rm = TRUE) %>%
    right_join(., meta_new, by = "ID") %>%
    select(Training.Validation, ID, Age, Accession, value) %>%
    mutate(level = data_levels_new[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(expr_sub)
})

## combine all
sub_wide <- subset_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Training.Validation + ID + Age + level + correction + quantitation ~ Accession)

fwrite(sub_wide, "5_subset_expr_age_baseline.csv")

## after using JMP Pro: R square
model_dt <- read_excel("./5_models_age.xlsx")
sub_model <- model_dt %>%
  select(Training.Validation, level, correction, quantitation,
         Age, ends_with("Predicted"), ends_with("Residual"))

data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

r2_tables <- pblapply(1:30, function(i) {
  
  sub_tmp_i <- sub_model %>%
    filter(level %in% data_levels_new[i]) %>%
    filter(correction %in% correct_methods_new[i]) %>%
    filter(quantitation %in% pep2pro_methods_new[i]) %>%
    mutate(model = "Bootstrap Forest")
  
  sub_r2_i <- sub_tmp_i %>%
    group_by(Training.Validation, level, correction, quantitation, model) %>%
    summarise(r2 = 1 - sum((Age - `Bootstrap Forest Predicted`) ^ 2) / sum((Age - mean(Age)) ^ 2), .groups = "drop")
  
  return(sub_r2_i)
})

sub_r2 <- rbindlist(r2_tables)

fwrite(sub_r2, "6_r2_age.csv")


## Prediction of Study samples: HbA1C ------------------
clin_dt <- fread("../../data/clindata.csv")
clin <- clin_dt %>%
  filter(HBA1C != "") %>%
  tidyr::separate_rows(., HBA1C, sep = "\\|") %>%
  mutate_at("HBA1C", as.numeric) %>%
  group_by(Subject.ID, Week) %>%
  summarise(HBA1C = round(mean(HBA1C), 2), .groups = "drop") %>%
  select(Subject.ID, Week, HBA1C) %>%
  na.omit

data_list_new <- c(data_list[c(1, 10, 19)], data_list)
data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

## check
set.seed(2023)
HBA1C_random <- sample(clin$HBA1C, replace = FALSE)
HBA1C_random2 <- sample(HBA1C_random, replace = FALSE)
HBA1C_random3 <- sample(HBA1C_random2, replace = FALSE)
HBA1C_random4 <- sample(HBA1C_random3, replace = FALSE)
HBA1C_true <- clin$HBA1C
df <- data.frame(HBA1C_true = HBA1C_true,
                 HBA1C_random = HBA1C_random,
                 HBA1C_random2 = HBA1C_random2,
                 HBA1C_random3 = HBA1C_random3,
                 HBA1C_random4 = HBA1C_random4)
p <- ggplot(df, aes(x = HBA1C_random4, y = HBA1C_true)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_bw() + 
  coord_fixed();p
ggsave("~/Desktop/HBA1C_trueVSrandom.pdf", p, width = 6, height = 6)


## set negative controls and filter by missing values
subset_tables <- pblapply(1:30, function(i) {
  
  if (data_levels_new[i] %in% "Negative control") {
    
    set.seed(2023)
    clin_new <- clin %>%
      mutate(HBA1C_random = sample(HBA1C, replace = FALSE)) %>%
      select(Subject.ID, Week, HBA1C_random) %>%
      dplyr::rename(HBA1C = HBA1C_random)
  } else {
    
    clin_new <- clin
  }
  
  meta_new <- meta %>%
    filter(!is.na(Subject)) %>%
    select(ID, Subject, Week, Training.Validation)
  
  expr <- data_list_new[[i]]
  colnames(expr)[1] <- "Accession"
  
  expr_new <- expr %>%
    select(Accession, all_of(meta_new$ID)) %>%
    na.omit
  
  expr_sub <- expr_new %>%
    reshape2::melt(., id = 1, variable.name = "ID", na.rm = TRUE) %>%
    left_join(., meta_new, by = "ID") %>%
    left_join(., clin_new, by = c("Subject" = "Subject.ID", "Week" = "Week")) %>%
    select(Training.Validation, ID, HBA1C, Accession, value) %>%
    mutate(level = data_levels_new[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(expr_sub)
})

## combine all
sub_wide <- subset_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Training.Validation + ID + HBA1C + level + correction + quantitation ~ Accession,
                  value.var = "value")

fwrite(sub_wide, "5_subset_expr_hba1c.csv")

## after using JMP Pro: RMSE
model_dt <- read_excel("./5_models_hba1c.xlsx")
sub_model <- model_dt %>%
  select(Training.Validation, level, correction, quantitation, HBA1C, ends_with("Predicted"))

data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

rmse_tables <- pblapply(1:30, function(i) {
  
  sub_tmp_i <- sub_model %>%
    filter(level %in% data_levels_new[i]) %>%
    filter(correction %in% correct_methods_new[i]) %>%
    filter(quantitation %in% pep2pro_methods_new[i]) %>%
    reshape2::melt(., id = 1:5, value.name = "HBA1C Predicted", variable.name = "model") %>%
    mutate_at("model", ~ gsub(" Predicted", "", .))
  
  sub_rmse_i <- sub_tmp_i %>%
    na.omit %>%
    group_by(Training.Validation, level, correction, quantitation, model) %>%
    summarise(r2 = 1 - sum((HBA1C - `HBA1C Predicted`) ^ 2) / sum((HBA1C - mean(HBA1C)) ^ 2), .groups = "drop")
    # summarise(rmse = sqrt(sum((HBA1C - `HBA1C Predicted`) ^ 2, na.rm = TRUE) / nrow(.)), .groups = "drop")
  
  return(sub_rmse_i)
})

sub_rmse <- rbindlist(rmse_tables)

fwrite(sub_rmse, "6_r2_hba1c.csv")


## Prediction of Study samples: delta HbA1C ------------------
clin_dt <- fread("../../data/clindata.csv")
clin <- clin_dt %>%
  filter(HBA1C != "") %>%
  tidyr::separate_rows(., HBA1C, sep = "\\|") %>%
  mutate_at("HBA1C", as.numeric) %>%
  reshape2::dcast(., Subject.ID ~ Week, value.var = "HBA1C", fun.aggregate = mean) %>%
  mutate(delta_HBA1C = (Week24 - Baseline) / Baseline) %>%
  select(Subject.ID, delta_HBA1C) %>%
  na.omit

data_list_new <- c(data_list[c(1, 10, 19)], data_list)
data_levels_new <- c(rep("Negative control", 3), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19)], pep2pro_methods)

## set negative controls and filter by missing values
subset_tables <- pblapply(1:30, function(i) {
  
  if (data_levels_new[i] %in% "Negative control") {
    
    set.seed(2023)
    clin_new <- clin %>%
      mutate(delta_HBA1C_random = sample(delta_HBA1C, replace = FALSE)) %>%
      select(Subject.ID, delta_HBA1C_random) %>%
      dplyr::rename(delta_HBA1C = delta_HBA1C_random)
  } else {
    
    clin_new <- clin
  }
    
  meta_new <- meta %>%
    filter(!is.na(Subject)) %>%
    select(ID, Subject, Training.Validation, Week)
  
  expr <- data_list_new[[i]]
  colnames(expr)[1] <- "Accession"
  
  expr_new <- expr %>%
    select(Accession, all_of(meta_new$ID)) %>%
    na.omit
  
  expr_delta <- expr_new %>%
    reshape2::melt(., id = 1, variable.name = "ID", na.rm = TRUE) %>%
    right_join(., meta_new, by = "ID") %>%
    reshape2::dcast(., Accession + Subject + Training.Validation ~ Week, value.var = "value",
                    fun.aggregate = mean, na.rm = TRUE) %>%
    mutate(delta_Protein = (Week24 - Baseline) / Baseline) %>%
    select(Training.Validation, Subject, Accession, delta_Protein) %>%
    right_join(., clin_new, by = c("Subject" = "Subject.ID")) %>%
    mutate(level = data_levels_new[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(expr_delta)
})

## combine all
sub_wide <- subset_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Training.Validation + Subject + delta_HBA1C + level + correction + quantitation ~ Accession,
                  value.var = "delta_Protein")

fwrite(sub_wide, "5_subset_expr_hba1c_delta.csv")

