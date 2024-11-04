## Created on May 13th, 2023
## Updated on Apr 25th, 2024
## Author: Qiaochu Chen

library(parallel)
library(pbapply)
library(readxl)
library(dplyr)
library(tibble)
library(reshape2)
library(data.table)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)

rm(list = ls())


## input ---------------------
setwd("~/Desktop/毕业课题/manuscript-batch/CaseStudy1_Quartet/")
meta <- read.csv("data/expfiles/meta.csv")

## protein level
all_files <- paste("data/expfiles/protein",
                   rep(rep(c("/balanced", "/confounded"), each = 27)),
                   rep(rep(c("/raw_mapped", "/final_mapped", "/raw_mapped"), c(1, 4, 4)), 6),
                   rep(rep(c("/ibaq", "/lfq", "/top3"), each = 9), 2),
                   "/expdata_",
                   rep(c("log_fot", rep(c("combat", "median", "ratio", "ruv"), 2)), 6),
                   ".csv",
                   sep = "")

data_list <- mclapply(all_files, function(x) data.table::fread(x))

## global ----------------------
setwd("~/Desktop/毕业课题/manuscript-batch/CaseStudy1_Quartet/results/tables/")

source("~/Desktop/毕业课题/utils/PCA.R")
source("~/Desktop/毕业课题/utils/DEP.R")

batches <- unique(meta$batch)
names(batches) <- unique(apply(meta, 1, function(x) paste(x[10], x[9], sep = "_")))

data_levels <- rep(rep(c("Uncorrected", "Peptide-corrected", "Protein-corrected"), c(1, 4, 4)), 6)

scenarios <- rep(c("Balanced design", "Confounded design"), each = 27)

correct_methods <- rep(c("Log transformed", rep(c("ComBat", "Median centering", "Ratio", "RUV-III-C"), 2)), 6)

pep2pro_methods <- rep(rep(c("iBAQ", "MaxLFQ", "TopPep3"), each = 9), 2)


## Peptide-level diagnosis ---------------------
expr_pep <- fread("../../data/expfiles/peptide/expdata_log.csv")

## count
sub_count <- expr_pep %>%
  reshape2::melt(., id = 1, variable.name = "library", na.rm = TRUE) %>%
  tidyr::separate(., library, c("batch", "sample", "replicate")) %>%
  group_by(sample, batch) %>%
  summarise(count = length(unique(sequence)))

fwrite(sub_count, "1_count_peptide_level.csv")

## CV
sub_cv <- expr_pep %>%
  reshape2::melt(., id = 1, variable.name = "library", na.rm = TRUE) %>%
  tidyr::separate(., library, c("batch", "sample", "replicate")) %>%
  mutate_at("value", exp) %>%
  group_by(sequence, sample) %>%
  mutate(cv_global = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE) * 100) %>%
  ungroup %>%
  group_by(sequence, sample, batch) %>%
  mutate(cv_intrabatch = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE) * 100) %>%
  distinct(sequence, sample, batch, cv_intrabatch, cv_global)

fwrite(sub_cv, "1_cv_peptide_level.csv")

## PCA
exprdata_t <- expr_pep %>%
  column_to_rownames("sequence") %>%
  filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
  mutate_all(~ ifelse(is.infinite(.), NA, .)) %>%
  na.omit %>%
  t

metadata <- meta %>%
  filter(library %in% rownames(exprdata_t)) %>%
  column_to_rownames("library")

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        center = TRUE, scale = TRUE, group = "sample",
                        biplot = FALSE, dictGroups = metadata$sample,
                        snr = TRUE, plot = FALSE)

saveRDS(pca_results, "./1_pca_peptide_level.rds")


## extract retention time
mq_files <- list.files("../../data/rawfiles/MaxQuant", full.names = T, recursive = T)
mq_evidence <-  mq_files[grepl("evidence", mq_files)]
mq_evidence <- mq_evidence[c(1, 4, 6, 2, 4, 3)]

evidence_tables <- pblapply(1:6, function(i) {
  
  evi_dt <- fread(mq_evidence[i], showProgress = FALSE)
  
  evi_dt_i <- evi_dt %>%
    mutate(batch = batches[i]) %>%
    tidyr::separate_rows(`Protein names`, sep = ";") %>%
    select(batch, Sequence, `Protein names`,
           `Retention time`, `Retention length`,
           `m/z`,  Mass, Charge) %>%
    distinct(batch, Sequence, `Protein names`, .keep_all = TRUE)
  
  return(evi_dt_i)
  
})

sub_rt <- evidence_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Sequence + `Protein names` ~ batch, value.var = "Retention time") %>%
  dplyr::rename_if(is.numeric, ~ paste("RT", ., sep = "_"))

sub_rl <- evidence_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Sequence + `Protein names` ~ batch, value.var = "Retention length") %>%
  mutate_all(~ ifelse(. %in% "NaN", NA, .)) %>%
  dplyr::rename_if(is.numeric, ~ paste("RL", ., sep = "_"))

sub_mz <- evidence_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., Sequence + `Protein names` ~ batch, value.var = "m/z") %>%
  dplyr::rename_if(is.numeric, ~ paste("MZ", ., sep = "_"))

sub_evi <- sub_rt %>%
  full_join(., sub_rl, by = c("Sequence", "Protein names"), relationship = "many-to-many") %>%
  full_join(., sub_mz, by = c("Sequence", "Protein names"), relationship = "many-to-many")

data.table::fwrite(sub_evi, "1_evidence_peptide_level.csv")


## Calculate Pearson correlation coefficient --------------------
pcc_tables <- pblapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  cor_mat <- expr %>%
    column_to_rownames("feature") %>%
    select(!contains("D6")) %>%
    cor(., use = "pairwise.complete.obs")
  
  cor_mat[upper.tri(cor_mat, diag = TRUE)] <- NA
  
  sub_pcc <- cor_mat %>%
    as.data.frame %>%
    rownames_to_column("library1") %>%
    reshape2::melt(., id = 1, variable.name = "library2", value.name = "PCC", na.rm = TRUE) %>%
    mutate(group1 = str_extract(library1,"(?<=DAT[1-6]_).+?(?=_)")) %>%
    mutate(group2 = str_extract(library2,"(?<=DAT[1-6]_).+?(?=_)")) %>%
    mutate(class = ifelse(group1 == group2, "Intra", "Inter")) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(sub_pcc)
  
})

sub_pcc_final <- rbindlist(pcc_tables)

fwrite(sub_pcc_final, "./2_pcc.csv")


## Rank --------------------------
rank_tables <- pblapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  sub_rank <- expr %>%
    tibble::column_to_rownames("feature") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 72)) %>%
    mutate(mean = apply(., 1, mean, na.rm = TRUE)) %>%
    arrange(desc(mean)) %>%
    mutate(rank = 1:nrow(.)) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    tibble::rownames_to_column("feature") %>%
    select(feature, scenario, level, correction, quantitation, mean, rank)
  
  return(sub_rank)
})

sub_rank_final <- rbindlist(rank_tables)

fwrite(sub_rank_final, "./2_rank.csv")


## Calculate CV --------------------
calculate_cv <- function(df_quant_wide, pattern = "D5") {
  
  match_pat <- which(grepl(pattern, colnames(expr)))
  
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

cv_tables <- mclapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  sub_cv <- expr %>%
    tibble::column_to_rownames("feature") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 72)) %>%
    tibble::rownames_to_column("feature") %>%
    melt(., id = 1, variable.name = "library", na.rm = TRUE) %>%
    tidyr::separate(., library, c("batch", "sample", "replicate")) %>%
    mutate_at("value", exp) %>%
    group_by(feature, sample) %>%
    summarise(cv = sd(value, na.rm = TRUE) / mean(value, na.rm = TRUE) * 100) %>%
    reshape2::melt(., id = 1:2, variable.name = "type", na.rm = TRUE) %>%
    reshape2::dcast(., feature ~ type + sample, value.var = "value") %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    select(feature, scenario, level, correction, quantitation, everything())
  
  return(sub_cv)
  
}, mc.cores = 4)

sub_cv_final <- cv_tables %>%
  rbindlist(.) %>%
  mutate_all(~ ifelse(is.infinite(.), NA, .)) %>%
  mutate_if(is.numeric, ~ round(., digits = 3))

fwrite(sub_cv_final, "./2_cv.csv")


## Count --------------------
pep2prot_dt <- fread("../../data/expfiles/peptide/pep2prot.csv")
cv_dt <- fread("./2_cv.csv")

## 约5分钟
count_tables <- pblapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  sub_count <- expr %>%
    tibble::column_to_rownames("feature") %>%
    filter(apply(., 1, function(x) sum(is.na(x)) < 72)) %>%
    tibble::rownames_to_column("feature") %>%
    reshape2::melt(., id = 1, variable.name = "library", na.rm = TRUE) %>%
    tidyr::separate(., library, c("batch", "sample", "replicate")) %>%
    group_by(sample) %>%
    summarise(count = length(unique(feature))) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(sub_count)
})

sub_count_final <- rbindlist(count_tables)

sub_count_filter <- cv_dt %>%
  group_by(scenario, level, correction, quantitation) %>%
  summarise_at(vars(cv_D5:cv_M8), ~ length(which(. < 30))) %>%
  reshape2::melt(., id = 1:4, variable.name = "sample", value.name = "filter_count") %>%
  mutate_at("sample", ~ gsub("cv_", "", .))

sub_count_total <- pep2prot_dt %>%
  filter(protein_names != "") %>%
  tidyr::separate_rows(protein_names, sep = ";") %>%
  summarise(total_count = length(unique(protein_names)))

sub_count_summary <- sub_count_final %>%
  left_join(., sub_count_filter, by = c("sample", "scenario",
                                        "level", "correction", "quantitation")) %>%
  mutate(total_count = sub_count_total$total_count) %>%
  mutate(`Missed` = (total_count - count) / total_count * 100,
         `Unqualified` = (count - filter_count) / total_count * 100,
         `Qualified` = filter_count / total_count * 100) %>%
  mutate_if(is.numeric, ~ round(., digits = 2)) %>%
  select(sample, scenario, level, correction, quantitation, everything())

fwrite(sub_count_summary, "./2_count.csv")


## Calculate PVCA -----------------------
## 约1分钟
sub_pvca_tables <- pblapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  exprdata_t <- expr %>%
    column_to_rownames("feature") %>%
    select(!contains("D6")) %>%
    filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
    mutate_all(~ ifelse(is.infinite(.), NA, .)) %>%
    # mutate_all(~ ifelse(is.na(.), log(10 ^ (-6)), .)) %>%
    na.omit %>%
    t
  
  if (ncol(exprdata_t) != 0) {
    
    metadata <- meta %>%
    filter(library %in% rownames(exprdata_t)) %>%
    column_to_rownames("library") %>%
    dplyr::select(lab, instrument, sample)
  
    exprdata_t <- exprdata_t[rownames(exprdata_t) %in% rownames(metadata), ]
  
    # ## 判断多重共线性
    # exprdata_t_cor <- cor(exprdata_t)
    # exprdata_t_kappa <- kappa(exprdata_t_cor, exact = TRUE)
    # 
    # if (exprdata_t_kappa > 1000) {
    # 
    #   exprdata_t_eigen <- eigen(exprdata_t_cor)
    #   exprdata_t_values <- exprdata_t_eigen$values
    # 
    #   exprdata_t <- exprdata_t[, -which(exprdata_t_values > 1)]
    # 
    # }
  
    pvca_results <- main_pvca(exprdata_t = exprdata_t, metadata = metadata, plot = FALSE)
    
    sub_pvca <- pvca_results$table %>%
      mutate(scenario = scenarios[i]) %>%
      mutate(level = data_levels[i]) %>%
      mutate(correction = correct_methods[i]) %>%
      mutate(quantitation = pep2pro_methods[i])
  
    return(sub_pvca)

  }
  
})

sub_pvca_final <- data.table::rbindlist(sub_pvca_tables)

data.table::fwrite(sub_pvca_final, "3_pvca.csv")


## Calculate PCA and SNR --------------------
cv_dt <- fread("2_cv.csv")
sub_cv <- cv_dt %>%
  mutate(n = apply(., 1, function(x) length(which(as.numeric(x[6:9]) <= 30)))) %>%
  mutate(qualified = ifelse(n == 4, "Yes", "No")) %>%
  select(feature, level, scenario, correction, quantitation, qualified)

rank_dt <- fread("2_rank.csv")

# 约1小时: 按Qualified/rank分类
pca_objs <- pblapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  sub_cv_tmp <- sub_cv %>%
    filter(level %in% data_levels[i]) %>%
    filter(scenario %in% scenarios[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    left_join(., rank_dt, by = c("feature", "scenario", "level", "correction", "quantitation"))
  
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
    
    exprdata_t <- expr %>%
      filter(feature %in% sub_cv_x$feature) %>%
      column_to_rownames("feature") %>%
      # select(!contains("D6")) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0)) %>%
      mutate_all(~ ifelse(is.infinite(.), NA, .)) %>%
      # mutate_all(~ ifelse(is.na(.), log(10 ^ (-6)), .)) %>%
      na.omit %>%
      t
    
    if (ncol(exprdata_t) > 1) {
      
      metadata <- meta %>%
        filter(library %in% rownames(exprdata_t)) %>%
        column_to_rownames("library")
      
      # ## 判断多重共线性
      # exprdata_t_cor <- cor(exprdata_t)
      # exprdata_t_kappa <- kappa(exprdata_t_cor, exact = TRUE)
      # 
      # if (exprdata_t_kappa > 1000) {
      #   
      #   exprdata_t_eigen <- eigen(exprdata_t_cor)
      #   exprdata_t_values <- exprdata_t_eigen$values
      #   
      #   exprdata_t <- exprdata_t[, -which(exprdata_t_values > 1)]
      #   
      # }
      
      ## PCA分析样本
      pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                              center = TRUE, scale = TRUE, group = "sample",
                              biplot = FALSE, dictGroups = metadata$sample,
                              snr = TRUE, plot = FALSE)
      
      pca_results <- c(pca_results, list(i = i), list(all = all_id))
      
    } else {
      
      pca_results <- NULL
    }
    
    pca_objs_by_x <- c(pca_objs_by_x, list(pca_results))
  }

  return(pca_objs_by_x)
})

sub_pca_tables <- pblapply(1:54, function(i) {
  
  pca_results <- pca_objs[[i]]
  
  pca_objs_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_pca_values_by_x <- pca_results_x$pcs_values %>%
        select(any_of(c("sample", "library", colnames(meta), "PC1", "PC2"))) %>%
        mutate(all = all_id)
      
      return(sub_pca_values_by_x)
      
    }
  })
  
  sub_pca_values <- pca_objs_by_x %>%
    rbindlist(.) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(sub_pca_values)
  
})

sub_props_tables <- pblapply(1:54, function(i) {
  
  pca_results <- pca_objs[[i]]
  
  sub_props_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_pca_props_by_x <- pca_results_x$pcs_props %>%
        t %>% as.data.frame %>%
        mutate(all = all_id) %>%
        tibble::rownames_to_column("PC")
      
      return(sub_pca_props_by_x)
    }
  })
  
  sub_props <- sub_props_by_x %>%
    rbindlist(.) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(sub_props)
  
})

sub_snr_tables <- pblapply(1:54, function(i) {
  
  pca_results <- pca_objs[[i]]
  
  pca_objs_by_x <- mclapply(c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "), function(all_id) {
    
    pca_results_x <- pca_results[[match(all_id, c("All", "Qualified", "1 ~ 200", "200 ~ 500", "500 ~ "))]]
    
    if (!is.null(pca_results_x)) {
      
      sub_snr_by_x <- pca_results_x$snr_results %>%
        mutate(all = all_id)
      
      return(sub_snr_by_x)
      
    }
  })
  
  sub_snr <- pca_objs_by_x %>%
    rbindlist(.) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(sub_snr)
  
})

sub_pca_qc_final <- data.table::rbindlist(sub_pca_tables)
sub_pca_qc_props_final <- data.table::rbindlist(sub_props_tables)
sub_snr_final <- data.table::rbindlist(sub_snr_tables)

data.table::fwrite(sub_snr_final, "3_pca_snr.csv")
data.table::fwrite(sub_pca_qc_final, "3_pca_qc.csv")
data.table::fwrite(sub_pca_qc_props_final, "3_pca_qc_props.csv")


## Calculate DEPs: D5, F7, M8 --------------------
## As a whole: 约1小时30分钟
dep_tables <- pblapply(1:54, function(i) {
  
  expr <- data_list[[i]]
  colnames(expr)[1] <- "feature"
  
  sample_pairs <- list(c("D5", "F7"), c("D5", "M8"), c("F7", "M8"))
  
  dep_list <- mclapply(sample_pairs, function(sample_pair) {
    
    metadata <- meta %>%
      filter(sample %in% sample_pair) %>%
      filter(library %in% colnames(expr))
    
    exprdata <- expr %>%
      column_to_rownames("feature") %>%
      select(all_of(metadata$library)) %>%
      mutate_all(~ ifelse(is.na(.)|. == Inf, 0, exp(.))) %>%
      filter(apply(., 1, function(x) sd(x, na.rm = TRUE) != 0))
    
    group <- factor(metadata$sample, ordered = TRUE)
    
    dep_list <- main_dep(exprdata, group, na_threshold = 0,
                         test_method = "t", p_adjust = "BH",
                         plot_volcano = FALSE,
                         silent = TRUE)
    
    dep_table <- dep_list$data
    
    return(dep_table)
    
  }, mc.cores = 3)
  
  dep_dt <- dep_list %>%
    rbindlist(.) %>%
    # filter(adj.P.Val < 0.05) %>% ## 如果使用limma计算差异表达的话
    filter(p < 0.05) %>%
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i])
  
  return(dep_dt)
  
})

sub_dep <- rbindlist(dep_tables)

data.table::fwrite(sub_dep, "4_dep_by_sample.csv")


## Calculate RC & RMSE: reference dataset from Yaqing ------------------
mapping_dt <- fread("../../data/expfiles/peptide/pep2prot.csv")

df_mapping <- mapping_dt %>%
  filter(protein_names != "") %>%
  tidyr::separate_rows(protein_names, sep = ";") %>%
  distinct(gene_names, protein_names)

ref_dt <- read_excel("../../data/rawfiles/Quartet_reference_datasets/Supplementary Table 3. Reference datasets of differentially expressed features.xlsx")

sub_ref <- ref_dt %>%
  filter(omic %in% "Protein") %>%
  left_join(., df_mapping, by = c("feature" = "gene_names"), relationship = "many-to-many") %>%
  na.omit %>%
  distinct(protein_names, pair, log2FC.ref, change.ref)

dep_dt <- fread("4_dep_by_sample.csv")

sub_test <- dep_dt %>%
  filter(feature %in% sub_ref$protein_names) %>%
  filter(p < .05) %>%
  mutate(pair = paste(group1, group2, sep = "/"))

cv_dt <- fread("2_cv.csv")
sub_cv <- cv_dt %>%
  mutate(n = apply(., 1, function(x) length(which(as.numeric(x[6:9]) <= 30)))) %>%
  mutate(qualified = ifelse(n == 4, "Yes", "No")) %>%
  select(feature, level, scenario, correction, quantitation, qualified)

rank_dt <- fread("2_rank.csv")

sub_test <- sub_test %>%
  left_join(., sub_cv, by = c("feature", "scenario", "level", "correction", "quantitation")) %>%
  left_join(., rank_dt, by = c("feature", "scenario", "level", "correction", "quantitation"))

metrics_tables1 <- pblapply(1:54, function(i) {
  
  sub_test_i <- sub_test %>%
    filter(scenario %in% scenarios[i]) %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    select(feature, pair, qualified, estimate)
  
  sub_comb_i <- sub_test_i %>%
    right_join(., sub_ref, by = c("feature" = "protein_names", "pair" = "pair")) %>%
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
    
    metrics_by_pair_x <- mclapply(c("D5/F7", "D5/M8", "F7/M8"), function(pair_id) {
      
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
    mutate(scenario = scenarios[i]) %>%
    mutate(level = data_levels[i]) %>%
    mutate(correction = correct_methods[i]) %>%
    mutate(quantitation = pep2pro_methods[i]) %>%
    select(scenario, level, correction, quantitation, all, everything())
  
  return(metrics_dt_i)
})

metrics_tables2 <- pblapply(1:54, function(i) {
  
  sub_test_i <- sub_test %>%
    filter(scenario %in% scenarios[i]) %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) %>%
    select(feature, pair, rank, estimate)
  
  sub_comb_i <- sub_test_i %>%
    right_join(., sub_ref, by = c("feature" = "protein_names", "pair" = "pair")) %>%
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
    
    metrics_by_pair_x <- mclapply(c("D5/F7", "D5/M8", "F7/M8"), function(pair_id) {
      
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
    mutate(scenario = scenarios[i]) %>%
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


## Prediction of Study samples: Label -------------------
## get DEPs
dep_table <- data.table::fread("4_dep_by_sample.csv")
sub_deps <- dep_table %>%
  filter(p.adj < 0.05, abs(estimate) > 1)

data_list_new <- c(data_list[c(1, 10, 19, 28, 37, 46)], data_list)
scenarios_new <- c(scenarios[c(1, 10, 19, 28, 37, 46)], scenarios)
data_levels_new <- c(rep("Negative control", 6), data_levels)
correct_methods_new <- c(correct_methods[c(1, 10, 19, 28, 37, 46)], correct_methods)
pep2pro_methods_new <- c(pep2pro_methods[c(1, 10, 19, 28, 37, 46)], pep2pro_methods)

## set negative controls and filter by missing values
subset_tables <- pblapply(1:60, function(i) {
  
  if (data_levels_new[i] %in% "Negative control") {
    
    set.seed(2023)
    meta_new <- meta %>%
      mutate(sample_random = sample(sample, replace = TRUE)) %>%
      select(library, sample_random) %>%
      dplyr::rename(sample = sample_random)
    
    target_deps <- sub_deps %>%
      filter(scenario %in% scenarios_new[i]) %>%
      filter(level %in% "Uncorrected") %>%
      filter(correction %in% correct_methods_new[i]) %>%
      filter(quantitation %in% pep2pro_methods_new[i]) %>%
      pull(feature)
    
  } else {
    
    meta_new <- meta %>%
      select(library, sample)
    
    target_deps <- sub_deps %>%
      filter(scenario %in% scenarios_new[i]) %>%
      filter(level %in% data_levels_new[i]) %>%
      filter(correction %in% correct_methods_new[i]) %>%
      filter(quantitation %in% pep2pro_methods_new[i]) %>%
      pull(feature)
  }
  
  expr <- data_list_new[[i]]
  colnames(expr)[1] <- "feature"
  
  expr_new <- expr %>%
    select(feature, any_of(meta_new$library)) %>%
    filter(feature %in% target_deps) %>%
    na.omit
  
  expr_sub <- expr_new %>%
    reshape2::melt(., id = 1, variable.name = "library", na.rm = TRUE) %>%
    left_join(., meta_new, by = "library") %>%
    select(library, sample, feature, value) %>%
    mutate(scenario = scenarios_new[i]) %>%
    mutate(level = data_levels_new[i]) %>%
    mutate(correction = correct_methods_new[i]) %>%
    mutate(quantitation = pep2pro_methods_new[i])
  
  return(expr_sub)
})

## combine all
sub_wide <- subset_tables %>%
  rbindlist(.) %>%
  reshape2::dcast(., library + sample + scenario + level + correction + quantitation ~ feature)

fwrite(sub_wide, "5_subset_expr_sample.csv")


## GO + KEGG -----------------
uniprot <- read.csv("../../../../database/uniprot/uniprot_20220823.csv")
dep_dt <- read.csv("./4_dep_by_sample.csv")

sub_db <- uniprot %>%
  filter(Gene.Symbol != "") %>%
  distinct(Protein.Fullname, Gene.Symbol)

sub_dep <- dep_dt %>%
  filter(p.adj < .05 & abs(estimate) > 1) %>%
  inner_join(sub_db, ., by = c("Protein.Fullname" = "feature"))

## 约1小时
go_tables <- pblapply(1:54, function(i) {

  sample_pairs <- list(c("D5", "F7"), c("D5", "M8"), c("F7", "M8"))
  
  go_list <- mclapply(sample_pairs, function(sample_pair) {
    
    deps_i <- sub_dep %>%
      filter(scenario %in% scenarios[i]) %>%
      filter(level %in% data_levels[i]) %>%
      filter(quantitation %in% pep2pro_methods[i]) %>%
      filter(correction %in% correct_methods[i]) %>%
      filter(group1 %in% sample_pair[1], group2 %in% sample_pair[2]) %>%
      pull(Gene.Symbol) %>%
      unique
    
    if (length(deps_i) > 1) {
  
      genes_i <- suppressMessages({
          suppressWarnings({
          deps_i %>%
          bitr(.,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db) %>%
          pull(ENTREZID)
        })
      })
  
      ego_ALL <- enrichGO(gene = genes_i, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                          ont = "ALL", pAdjustMethod = "BH", minGSSize = 1,
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                          readable = TRUE)
  
      sub_go_j <- ego_ALL %>%
        as.data.frame %>%
        mutate(scenario = scenarios[i]) %>%
        mutate(level = data_levels[i]) %>%
        mutate(correction = correct_methods[i]) %>%
        mutate(quantitation = pep2pro_methods[i]) %>%
        mutate(group1 = sample_pair[1], group2 = sample_pair[2]) %>%
        select(pattern, scenario, level, correction, quantitation, group1, group2, everything())
      
      return(sub_go_j)
    }

  })
  
  sub_go <- rbindlist(go_list)
  
  return(sub_go)
    
})

sub_go_final <- rbindlist(go_tables)

fwrite(sub_go_final, "5_go_all.csv")

## 约2分钟
kegg_tables <- pblapply(1:54, function(i) {
  
  sample_pairs <- list(c("D5", "F7"), c("D5", "M8"), c("F7", "M8"))
  
  deps_tmp <- sub_dep %>%
    filter(scenario %in% scenarios[i]) %>%
    filter(level %in% data_levels[i]) %>%
    filter(quantitation %in% pep2pro_methods[i]) %>%
    filter(correction %in% correct_methods[i]) 
  
  if (nrow(deps_tmp) !=0 ) {
    
    kegg_list <- mclapply(sample_pairs, function(sample_pair) {
      
      deps_i <- deps_tmp %>%
        filter(group1 %in% sample_pair[1], group2 %in% sample_pair[2]) %>%
        pull(Gene.Symbol) %>%
        unique
      
      if (length(deps_i) > 1) {
        
        genes_i <- suppressMessages({
          suppressWarnings({
            deps_i %>%
              bitr(.,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db) %>%
              pull(ENTREZID)
          })
        })
        
        ekegg <- suppressMessages({
          suppressWarnings({enrichKEGG(gene = genes_i, keyType = "kegg", #organism = "human",
                                       qvalueCutoff = 0.2, pvalueCutoff = 0.05)
          })
        })
        
        if (!is.null(ekegg)) {
          
          sub_kegg_j <- ekegg %>%
            as.data.frame %>%
            mutate(group1 = sample_pair[1], group2 = sample_pair[2]) %>%
            select(group1, group2, everything())
          
          return(sub_kegg_j)
          
        }
      }
    })
    
    sub_kegg <- kegg_list %>%
      rbindlist(.) %>%
      mutate(scenario = scenarios[i]) %>%
      mutate(level = data_levels[i]) %>%
      mutate(correction = correct_methods[i]) %>%
      mutate(quantitation = pep2pro_methods[i])
    
    if (nrow(sub_kegg) > 1) return(sub_kegg)
    
  }
  
})

sub_kegg_final <- rbindlist(kegg_tables)

fwrite(sub_kegg_final, "5_kegg.csv")


