## Created on Mar 14th, 2024
## Updated on Jan 27th, 2025
## Author: Qiaochu Chen
## Title: To combine all figures

library(parallel)
library(pbapply)
library(readxl)
library(data.table)
library(dplyr)
library(tibble)
library(reshape2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(gg.gap)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggradar)
library(VennDiagram)
library(UpSetR)

rm(list = ls())

## global ----------------------
source("~/Desktop/utils/PCA.R")
source("~/Desktop/utils/DEP.R")

setwd("~/Desktop/PGx/Project_BEC_Proteomics/")

meta1 <- read.csv("./CaseStudy1_Quartet/data/expfiles/meta.csv")
meta2 <- read.csv("./CaseStudy2_CG/data/metadata.csv")

data_levels <- rep(rep(c("Uncorrected", "Peptide-corrected", "Protein-corrected"), c(1, 4, 4)), 9)

scenarios <- rep(c("Balanced design", "Confounded design", "Random design"), each = 27)

correct_methods <- rep(c("Log transformed", rep(c("Ratio", "Median centering", "RUV-III-C", "ComBat"), 2)), 9)

pep2pro_methods <- rep(rep(c("MaxLFQ", "iBAQ", "TopPep3"), each = 9), 3)

dictLabelsScenario <- rep(c("Multi-lab (Balanced)", "Multi-lab (Confounded)", "Multi-batch"), each = 27)
names(dictLabelsScenario) <- scenarios

dictLabelsBatch <- c("L1-1", "L2-1", "L3", "L1-2", "L2-2", "L4")
names(dictLabelsBatch) <- paste("DAT", 1:6, sep = "")

dictLabelsName <- 1:81
names(dictLabelsName) <- paste(scenarios, data_levels, correct_methods, pep2pro_methods, sep = "_")

dictColorsSample1 <- c("#4CC3D9", "#7BC8A4", "#FFC65D", "#F16745")
names(dictColorsSample1) <- c("D5", "D6", "F7", "M8")

dictColorsSample2 <- c("#BEAED4", "#FDC086", "#7FC97F", "#6BAED6")
names(dictColorsSample2) <- c("QC sample (P10)", "QC sample (P11)", "QC sample (PM)", "Study sample")

dictColorsBatch1 <- c("#A65628", "#FF7F00", "#FDBF6F", "#FFFF33", "#F781BF", "#999999")
names(dictColorsBatch1) <- unique(meta1$batch)

dictShapesBatch1 <- 12:17
names(dictShapesBatch1) <- unique(meta1$batch)

dictColorsBatch2 <- c("#FB9A99", "#FDBF6F", "#1F78B4")
names(dictColorsBatch2) <- c("I", "II", "III")

dictColorsLevel <- c("#A50F15", "#984EA3", "#377EB8")
names(dictColorsLevel) <- unique(data_levels)

dictColorsMethod <- c("#FFFF99", "#7FC97F", "#BEAED4", "#386CB0", "#FDC086")
names(dictColorsMethod) <- unique(correct_methods)

dictColorsPep2pro <- c("#66C2A5", "#8DA0CB","#FC8D62")
names(dictColorsPep2pro) <- unique(pep2pro_methods)

dictColorsScenario <- c("#A6D854", "#FFD92F", "pink")
names(dictColorsScenario) <- unique(scenarios)

dictColorsSex <- c("F" = "#e98989", "M" = "#78b7ee")
dictColorsWeek <- c("Baseline" = "#1B9E77", "Week24" = "#D95F02")

dictColorsTV <- c("#7db7e6", "#e6e384")
names(dictColorsTV) <- 1:2

dictColorsClass <- c("#b80d0d", "#d97c11", "#70c402", "#0f9115")
names(dictColorsClass) <- c("Bad", "Fair", "Good", "Great")

dictClassNames <- c("Training", "Validation")
names(dictClassNames) <- 1:2

dictLabelsML <- c("Decision Tree", "Bootstrap", "Boosted", "KNN", "SVM", "Stepwise", "Lasso")
names(dictLabelsML) <- c("Decision Tree", "Bootstrap Forest", "Boosted Tree",
                             "K Nearest Neighbors", "Support Vector Machines",
                             "Fit Stepwise", "Generalized Regression Lasso")

## 2023年12月重新按照清洗批次分组
# x_breaks <- c(144, 399, 503, 623, 819, 969, 1165)
x_breaks <- c(503, 969)

# ## 这19个IDs：基于log2水平80%PCs(前28个PCs)后Cook距离大于四倍平均值
# ## 2023年4月初计算了先去除这19个样本后批次校正的数据。
# ## 2023年4月底计算全部1495个样本批次校正的数据。
# rm_ids <- c("ExpB19", "ExpC86", "ExpD43", "ExpE6", "ExpF65", "ExpE30",
#             "ExpE31", "ExpG52", "ExpI14", "ExpI49", "ExpI75", "ExpI85",
#             "ExpI93", "ExpJ1", "ExpK3", "ExpH34", "ExpL85", "ExpP21", "ExpQ42")
# ## 这14个IDs：基于UMAP发现的离群值。
# rm_ids <- c("ExpB19", "ExpB48", "ExpB55", "ExpB81", "ExpD45", "ExpD38", "ExpD76", "ExpE81",
#             "ExpF25", "ExpG69", "ExpK18", "ExpK21", "ExpK35", "ExpK40")

rm_ids <- c("ExpF25", "ExpB19", "ExpP72")


## Basic information (Count): Quartet ---------------------
## peptide-level
count_dt <- fread("./CaseStudy1_Quartet/results/tables/1_count_peptide_level.csv")

mean(count_dt$count)
sd(count_dt$count)
range(count_dt$count)

p_count_pep1 <- ggplot(count_dt, aes(x = sample, y = count)) +
  stat_summary(aes(fill = sample), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .15) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, margin = unit(c(0, .3, 0, 0), "cm")),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 20)) +
  labs(y = "Number of peptides") +
  scale_fill_manual(values = dictColorsSample1) +
  facet_grid(cols = vars(batch), labeller = as_labeller(dictLabelsBatch)) +
  # scale_x_discrete(labels = dictLabelsBatch) +
  scale_y_continuous(n.breaks = 10);p_count_pep1

## Basic information (CV): Quartet ---------------------
## CV
cv_dt <- fread("./CaseStudy1_Quartet/results/tables/1_cv_peptide_level.csv")

sub_cv <- cv_dt %>%
  select(sequence, sample, starts_with("cv")) %>%
  reshape2::melt(., id = 1:2, variable.name = "type", na.rm = TRUE) %>%
  mutate_at("type", ~ gsub("cv_", "", .)) %>%
  mutate_at("type", ~ Hmisc::capitalize(gsub("intra", "intra-", .))) %>%
  mutate(name = paste(sample, type)) %>%
  mutate_at("name", ~ factor(., levels = unique(name)))

mean(sub_cv$value[sub_cv$type %in% "Intra-batch"])
sd(sub_cv$value[sub_cv$type %in% "Intra-batch"])
mean(sub_cv$value[sub_cv$type %in% "Global"])
sd(sub_cv$value[sub_cv$type %in% "Global"])

p_cv_pep1 <- ggplot(sub_cv, aes(x = sample, y = value)) +
  geom_violin(aes(fill = sample, color = sample), alpha = .8) +
  geom_boxplot(fill = "white", color = "black", width = .1, outlier.size = .5) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, margin = unit(c(0, .3, 0, 0), "cm")),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16)) +
  labs(y = "Coefficient of variation (%)") +
  facet_grid(cols = vars(type)) +
  scale_color_manual(values = dictColorsSample1) +
  scale_fill_manual(values = dictColorsSample1) +
  scale_y_continuous(n.breaks = 10)


## Batch effect diagnosis (Intensity): Quartet ---------------------
## peptide-level
expr_dt <- fread("./CaseStudy1_Quartet/data/expfiles/peptide/expdata_log.csv")

sub_intensity <- expr_dt %>%
  reshape2::melt(., id = 1, na.rm = TRUE, variable.name = "library") %>%
  left_join(., meta1, by = "library")

p_intensity_pep1 <- ggplot(sub_intensity, aes(x = sample, y = value)) +
  geom_boxplot(aes(fill = sample), outlier.size = .5) +
  theme_bw() +
  theme(legend.position = "none",
        title = element_text(size = 16),
        text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20, margin = unit(c(0, .3, 0, 0), "cm")),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(size = 20)) +
  labs(title = "Peptide-level",
       subtitle = "",
       y = "Peptide intensity\n(Log transformed)") +
  scale_fill_manual(values = dictColorsSample1) +
  facet_grid(cols = vars(batch), labeller = as_labeller(dictLabelsBatch)) +
  # scale_x_discrete(labels = dictLabelsBatch) +
  scale_y_continuous(n.breaks = 10)

## protein-level
all_files <- list.files("./CaseStudy1_Quartet/data/expfiles/protein", recursive = TRUE, full.names = TRUE)
all_expdata_log <- all_files[grepl("log_fot", all_files)]
all_expdata_log <- all_expdata_log[grepl("balanced|confounded", all_expdata_log)]
plot_titles <- rep(unique(dictLabelsScenario)[1:2], each = 3)
plot_subtitles <- rep(c("iBAQ", "MaxLFQ", "TopPep3"), 3)

p_intensity_list1 <- pblapply(1:6, function(i) {
  
  expr_tmp <- fread(all_expdata_log[i])
  
  sub_intensity <- expr_tmp %>%
    reshape2::melt(., id = 1, na.rm = TRUE, variable.name = "library") %>%
    left_join(., meta1, by = "library")
  
  p_intensity <- ggplot(sub_intensity, aes(x = sample, y = value)) +
    geom_boxplot(aes(fill = sample), outlier.size = .5) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          text = element_text(size = 16),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20, margin = unit(c(0, .3, 0, 0), "cm")),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black"),
          strip.text = element_text(size = 20)) +
    labs(title = plot_titles[i],
         subtitle = plot_subtitles[i],
         y = "Protein quantity\n(Log transformed)") +
    scale_fill_manual(values = dictColorsSample1) +
    facet_grid(cols = vars(batch), labeller = as_labeller(dictLabelsBatch)) +
    # scale_x_discrete(labels = dictLabelsBatch) +
    scale_y_continuous(n.breaks = 10)
  
  return(p_intensity)
})


## Batch effect diagnosis (Intensity): CGZ --------------
## peptide-level
expr_dt <- fread("./CaseStudy2_CG/data/expfiles/peptide/expdata_log.csv")

sub_intensity <- expr_dt %>%
  reshape2::melt(., na.rm = TRUE, id = 1, variable.name = "ID") %>%
  group_by(ID) %>%
  summarise(intensity = mean(value, na.rm = TRUE)) %>%
  left_join(., meta2, by = "ID")

p_intensity_pep21 <- ggplot(sub_intensity, aes(x = Order, y = intensity)) +
  geom_point(aes(color = Cleaning), size = 2, alpha = .7) +
  geom_smooth(aes(group = Cleaning), color = "darkgrey") +
  geom_vline(xintercept = x_breaks, lty = 2, color = "black") +
  scale_color_manual(values = dictColorsBatch2, name = "Batch") +
  scale_x_continuous(breaks = c(seq(0, 1495, 200)[-c(4, 6)], x_breaks)) +
  theme_classic() +
  theme(title = element_text(size = 16),
        text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none") +
  labs(title = "Peptide-level",
       subtitle = "",
       y = "Peptide intensity\n(Log transformed)",
       x = "Injection order")

sub_intensity_qc <- sub_intensity %>%
  filter(grepl("QC", Type))

p_intensity_pep22 <- ggplot(sub_intensity_qc, aes(x = Order, y = intensity)) +
  geom_point(aes(color = Type), size = 5) +
  geom_smooth(aes(group = Cleaning), color = "darkgrey") +
  geom_vline(xintercept = x_breaks, lty = 2, color = "black") +
  scale_color_manual(values = dictColorsSample2, name = "Sample") +
  scale_x_continuous(breaks = c(seq(0, 1495, 200)[-c(4, 6)], x_breaks)) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.position = "none") +
  labs(y = "Peptide intensity\n(Log transformed)", x = "Injection order")

## protein-level
all_files <- list.files("./CaseStudy2_CG/data/expfiles/protein", recursive = TRUE, full.names = TRUE)
all_expdata_log <- all_files[grepl("log_fot", all_files)]
plot_titles <- rep(unique(dictLabelsScenario)[3], each = 3)
plot_subtitles <- c("iBAQ", "MaxLFQ", "TopPep3")

p_intensity_list2 <- pblapply(1:3, function(i) {
  
  expr_tmp <- fread(all_expdata_log[i])
  
  sub_intensity <- expr_tmp %>%
    reshape2::melt(., id = 1, na.rm = TRUE, variable.name = "ID") %>%
    group_by(ID) %>%
    summarise(intensity = mean(value, na.rm = TRUE)) %>%
    left_join(., meta2, by = "ID")
  
  p_intensity <- ggplot(sub_intensity, aes(x = Order, y = intensity)) +
    geom_point(aes(color = Cleaning), size = 5, alpha = .7) +
    geom_smooth(aes(group = Cleaning), color = "darkgrey") +
    geom_vline(xintercept = x_breaks, lty = 2, color = "black") +
    scale_color_manual(values = dictColorsBatch2, name = "Batch") +
    scale_x_continuous(breaks = c(seq(0, 1495, 200)[-c(4, 6)], x_breaks)) +
    theme_classic() +
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 16),
          title = element_text(size = 20),
          legend.position = "none")  +
    labs(title = plot_titles[i],
         subtitle = plot_subtitles[i],
         y = "Protein quantity\n(Log transformed)",
         x = "Injection order")
  
  return(p_intensity)
})


## Batch effect diagnosis (Peptide-level PCA): Quartet ---------------------
## peptide-level
pca_obj <- readRDS("./CaseStudy1_Quartet/results/tables/1_pca_peptide_level.rds")
pca_dt <- pca_obj$pcs_values
snr_dt <- pca_obj$snr_results
prop1 <- pca_obj$pcs_props$PC1[2]
prop2 <- pca_obj$pcs_props$PC2[2]

p_pca_pep1 <- ggplot(pca_dt,aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample, shape = batch), size = 2) +
  theme_bw() +
  theme(legend.position = "right",
        title = element_text(size = 16),
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.margin = unit(c(1,.5, .5, .5), "cm")) +
  scale_shape_manual(values = dictShapesBatch1, labels = dictLabelsBatch) +
  scale_color_manual(values = dictColorsSample1) +
  scale_y_continuous(limits = c(-20, 20)) +
  scale_x_continuous(limits = c(-70, 70)) +
  labs(color = "Sample", shape = "Batch",
       title = "Peptide-level", 
       subtitle = sprintf("SNR = %.2f", snr_dt$snr),
       x = sprintf("PC1 (%.2f%%)", prop1 * 100),
       y = sprintf("PC2 (%.2f%%)", prop2 * 100));p_pca_pep1


## Batch effect diagnosis (Protein-level PCA): Quartet --------------
## protein-level
all_pca_values <- fread("CaseStudy1_Quartet/results/tables/3_pca_qc.csv")
all_pca_props <- fread("CaseStudy1_Quartet/results/tables/3_pca_qc_props.csv")
all_pca_snr <- fread("CaseStudy1_Quartet/results/tables/3_pca_snr.csv")

sub_pca_values <- all_pca_values %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  filter(scenario %in% unique(scenarios)[1:2]) %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

sub_pca_props <- all_pca_props %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  filter(scenario %in% unique(scenarios)[1:2]) %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

sub_pca_snr <- all_pca_snr %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  filter(scenario %in% unique(scenarios)[1:2]) %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

all_dataset_groups <- unique(sub_pca_values$dataset_group)

p_pca_pro_list1 <- pblapply(all_dataset_groups, function(dataset_group_id) {
  
  sub_values_tmp <- sub_pca_values %>%
    filter(dataset_group %in% dataset_group_id)
  
  sub_props_tmp <- sub_pca_props %>%
    filter(dataset_group %in% dataset_group_id)
  
  sub_snr_tmp <- sub_pca_snr %>%
    filter(dataset_group %in% dataset_group_id)
  
  p_pca_tmp <- ggplot(sub_values_tmp,aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample, shape = batch), size = 2) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 16),
          text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1,.5, .5, .5), "cm")) +
    scale_shape_manual(values = dictShapesBatch1) +
    scale_color_manual(values = dictColorsSample1) +
    labs(color = "Sample", shape = "Batch",
         title = unlist(strsplit(dataset_group_id, split = "="))[2], 
         subtitle = sprintf("SNR = %.2f", sub_snr_tmp$snr),
         x = sprintf("PC1 (%.2f%%)", sub_props_tmp$`Proportion of Variance`[1] * 100),
         y = sprintf("PC2 (%.2f%%)", sub_props_tmp$`Proportion of Variance`[2] * 100))
  
  p_pca_tmp
  
  return(p_pca_tmp)
})

names(p_pca_pro_list1) <- all_dataset_groups


## Batch effect diagnosis (Peptide-level PCA): CGZ --------------
## peptide-level PCA: all samples
expr <- fread("./CaseStudy2_CG/data/expfiles/peptide/expdata_log.csv")
colnames(expr)[1] <- "Feature"

expr[is.na(expr)] <- 0
na_cutoff <- 1 * (ncol(expr) - 1)

exprdata_t <- expr %>%
  column_to_rownames("Feature") %>%
  filter(apply(., 1, function(x) length(which(x != 0)) >= na_cutoff)) %>%
  t

dim(exprdata_t)

metadata <- meta2 %>%
  column_to_rownames("ID") %>%
  mutate(Plate = sapply(Tube, function(x) unlist(strsplit(x, split = "[0-9]"))[1])) %>%
  dplyr::select(Type, Project, Site, Plate, Cleaning, Date, Type, Sex, Week)

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        group = "Cleaning", center = TRUE, scale = TRUE,
                        biplot = FALSE, plot = FALSE)

sub_pca_pep <- pca_results$pcs_values
sub_props_pep <- pca_results$pcs_props

p_pca_pep21 <- ggplot(sub_pca_pep,aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample), size = 2, alpha = .7) +
  theme_bw() +
  theme(legend.position = "right",
        title = element_text(size = 16),
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.margin = unit(c(1,.5, .5, .5), "cm")) +
  scale_color_manual(values = dictColorsBatch2, name = "Batch") +
  # scale_y_continuous(limits = c(-20, 20)) +
  # scale_x_continuous(limits = c(-70, 70)) +
  labs(title = "Peptide-level", 
       subtitle = "",
       x = sprintf("PC1 (%.2f%%)", sub_props_pep$PC1[2] * 100),
       y = sprintf("PC2 (%.2f%%)", sub_props_pep$PC2[2] * 100))

## peptide-level PCA: QC samples
exprdata_t <- expr %>%
  column_to_rownames("Feature") %>%
  select(all_of(meta2$ID[grepl("QC", meta2$Type)])) %>%
  filter(apply(., 1, function(x) length(which(x != 0)) >= 64)) %>%
  t

dim(exprdata_t)

metadata <- meta2 %>%
  filter(grepl("QC", Type)) %>%
  column_to_rownames("ID") %>%
  mutate(Plate = sapply(Tube, function(x) unlist(strsplit(x, split = "[0-9]"))[1])) %>%
  dplyr::select(Type, Project, Site, Plate, Cleaning, Date, Type, Sex, Week)

pca_results <- main_pca(exprdata_t = exprdata_t, metadata = metadata,
                        group = "Cleaning", center = TRUE, scale = TRUE,
                        biplot = FALSE, plot = FALSE, snr = TRUE)

sub_pca_pep <- pca_results$pcs_values
sub_props_pep <- pca_results$pcs_props
sub_snr_pep <- pca_results$snr_results

p_pca_pep22 <- ggplot(sub_pca_pep,aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Type, shape = sample), size = 5) +
  theme_bw() +
  theme(legend.position = "right",
        title = element_text(size = 16),
        text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.margin = unit(c(1,.5, .5, .5), "cm")) +
  scale_shape_manual(values = 15:17, name = "Batch") +
  scale_color_manual(values = dictColorsSample2, name = "Sample") +
  # scale_y_continuous(limits = c(-20, 20)) +
  # scale_x_continuous(limits = c(-70, 70)) +
  labs(title = "Peptide-level", 
       subtitle = sprintf("SNR = %.2f", sub_snr_pep$snr),
       x = sprintf("PC1 (%.2f%%)", sub_props_pep$PC1[2] * 100),
       y = sprintf("PC2 (%.2f%%)", sub_props_pep$PC2[2] * 100))


## Batch effect diagnosis (Protein-level PCA): CGZ QC samples --------------
## protein-level
all_pca_values <- fread("CaseStudy2_CG/results/tables/3_pca_qc.csv")
all_pca_props <- fread("CaseStudy2_CG/results/tables/3_pca_qc_props.csv")
all_pca_snr <- fread("CaseStudy2_CG/results/tables/3_pca_snr.csv")

sub_pca_values <- all_pca_values %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

sub_pca_props <- all_pca_props %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

sub_pca_snr <- all_pca_snr %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

all_dataset_groups <- unique(sub_pca_values$dataset_group)

p_pca_pro_list2 <- pblapply(all_dataset_groups, function(dataset_group_id) {
  
  sub_values_tmp <- sub_pca_values %>%
    filter(dataset_group %in% dataset_group_id)
  
  sub_props_tmp <- sub_pca_props %>%
    filter(dataset_group %in% dataset_group_id)
  
  sub_snr_tmp <- sub_pca_snr %>%
    filter(dataset_group %in% dataset_group_id)
  
  p_pca_tmp <- ggplot(sub_values_tmp,aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample, shape = Cleaning), size = 5) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 16),
          text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1,.5, .5, .5), "cm")) +
    scale_shape_manual(values = 15:17) +
    scale_color_manual(values = dictColorsSample2) +
    labs(color = "Sample", shape = "Batch",
         title = unlist(strsplit(dataset_group_id, split = "="))[2], 
         subtitle = sprintf("SNR = %.2f", sub_snr_tmp$snr),
         x = sprintf("PC1 (%.2f%%)", sub_props_tmp$`Proportion of Variance`[1] * 100),
         y = sprintf("PC2 (%.2f%%)", sub_props_tmp$`Proportion of Variance`[2] * 100))
  
  p_pca_tmp
  
  return(p_pca_tmp)
})

names(p_pca_pro_list2) <- all_dataset_groups


## Batch effect diagnosis (Protein-level PCA): CGZ all samples --------------
## protein-level
all_pca_values <- fread("CaseStudy2_CG/results/tables/2_pca.csv")
all_pca_props <- fread("CaseStudy2_CG/results/tables/2_pca_props.csv")

sub_pca_values <- all_pca_values %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

sub_pca_props <- all_pca_props %>%
  filter(all %in% "All") %>%
  filter(correction %in% "Log transformed") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("scenario", ~ dictLabelsScenario[.]) %>%
  mutate(dataset_group = paste(scenario, quantitation, sep = "="))

all_dataset_groups <- unique(sub_pca_values$dataset_group)

p_pca_all_pro_list2 <- pblapply(all_dataset_groups, function(dataset_group_id) {
  
  sub_values_tmp <- sub_pca_values %>%
    filter(dataset_group %in% dataset_group_id)
  
  sub_props_tmp <- sub_pca_props %>%
    filter(dataset_group %in% dataset_group_id)
  
  p_pca_tmp <- ggplot(sub_values_tmp,aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 2) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 16),
          text = element_text(size = 16),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          plot.margin = unit(c(1,.5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsBatch2) +
    labs(color = "Sample", shape = "Batch",
         title = unlist(strsplit(dataset_group_id, split = "="))[2], 
         subtitle = "",
         x = sprintf("PC1 (%.2f%%)", sub_props_tmp$`Proportion of Variance`[1] * 100),
         y = sprintf("PC2 (%.2f%%)", sub_props_tmp$`Proportion of Variance`[2] * 100))
  
  p_pca_tmp
  
  return(p_pca_tmp)
})

names(p_pca_all_pro_list2) <- all_dataset_groups


## Count ------------------------
count_dt1 <- fread("./CaseStudy1_Quartet/results/tables/2_count.csv")
count_dt2 <- fread("./CaseStudy2_CG/results/tables/1_count.csv")

sub_count2 <- count_dt2 %>%
  rename(sample = Type) %>%
  mutate(scenario = "Random design")

sub_count <- count_dt1 %>%
  rbind(., sub_count2) %>%
  mutate_at("level", ~ factor(., levels = unique(data_levels)))

sub_count_delta <- sub_count %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Median centering_Ratio_RUV-III-C", .)) %>%
  tidyr::separate_rows("correction", sep = "_") %>%
  reshape2::dcast(., sample + scenario + correction + quantitation ~ level, value.var = "filter_count") %>%
  mutate_at("Protein-corrected", ~. - `Uncorrected`) %>%
  mutate_at("Peptide-corrected", ~ . - `Uncorrected`) %>%
  select(!Uncorrected) %>%
  reshape2::melt(., id = 1:4) %>%
  mutate(class = ifelse(value > 0, "Increased", "Decreased")) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

nrow(sub_count_delta[sub_count_delta$scenario %in% "Random design" &
                       sub_count_delta$variable %in% "Protein-corrected",])
nrow(sub_count_delta[sub_count_delta$scenario %in% "Random design" &
                       sub_count_delta$class %in% "Increased" &
                       sub_count_delta$variable %in% "Protein-corrected",])

nrow(sub_count_delta[sub_count_delta$scenario %in% "Random design" &
                       sub_count_delta$variable %in% "Peptide-corrected",])
nrow(sub_count_delta[sub_count_delta$scenario %in% "Random design" &
                       sub_count_delta$class %in% "Increased" &
                       sub_count_delta$variable %in% "Peptide-corrected",])

p_count_bar1 <- ggplot(sub_count_delta, aes(x = variable)) +
  geom_bar(aes(fill = class), width = .7, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = c("Decreased" = "#66C2A5", "Increased" = "#E78AC3"),
                    name = "Number of qualified proteins") +
  scale_y_continuous(labels = scales::percent,
                     name = "Proportion of combinations") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

sub_count <- sub_count %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_count_bar12 <- ggplot(sub_count, aes(x = correction, y = Qualified)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Median centering", "Log transformed"),
                                 c("Ratio", "Log transformed"),
                                 c("Log transformed", "ComBat"),
                                 c("Log transformed", "RUV-III-C")),
              map_signif_level = TRUE,
              y_position = 40,
              step_increase = .05,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod) +
  scale_y_continuous(n.breaks = 10, name = "Number of qualified proteins") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_count_bar11 <- ggplot(sub_count, aes(x = quantitation, y = Qualified)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = 42,
              step_increase = .05,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro) +
  scale_y_continuous(n.breaks = 10, name = "Number of qualified proteins") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_count_bar_levelgroups1 <- ggplot(sub_count, aes(x = level, y = Qualified)) +
  geom_boxplot(aes(fill = level), width = .7) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected")),
              map_signif_level = TRUE,
              # y_position = 42,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 10, name = "Number of qualified proteins") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

sub_count_delta <- count_dt1 %>%
  rbind(., sub_count2) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Median centering_Ratio_RUV-III-C", .)) %>%
  tidyr::separate_rows("correction", sep = "_") %>%
  mutate(miss_count = total_count - count) %>%
  reshape2::dcast(., sample + scenario + correction + quantitation ~ level, value.var = "miss_count") %>%
  mutate_at("Protein-corrected", ~. - `Uncorrected`) %>%
  mutate_at("Peptide-corrected", ~ . - `Uncorrected`) %>%
  select(!Uncorrected) %>%
  reshape2::melt(., id = 1:4) %>%
  mutate(class = ifelse(value > 0, "Increased", "Decreased"))

p_count_bar2 <- ggplot(sub_count_delta, aes(x = variable)) +
  geom_bar(aes(fill = class), width = .7, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = c("Decreased" = "#66C2A5", "Increased" = "#E78AC3"),
                    name = "Number of missed proteins") +
  scale_y_continuous(labels = scales::percent,
                     name = "Proportion of combinations") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_count_bar22 <- ggplot(sub_count, aes(x = correction, y = Missed)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Median centering", "Log transformed"),
                                 c("Ratio", "Log transformed"),
                                 c("Log transformed", "ComBat"),
                                 c("Log transformed", "RUV-III-C")),
              map_signif_level = TRUE,
              y_position = 95,
              step_increase = .07,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 10),
                     name = "Number of missed proteins") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_count_bar21 <- ggplot(sub_count, aes(x = quantitation, y = Missed)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = 95,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 100, 10),
                     name = "Number of missed proteins") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_count_bar_levelgroups2 <- ggplot(sub_count, aes(x = level, y = Missed)) +
  geom_boxplot(aes(fill = level), width = .7) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected")),
              map_signif_level = TRUE,
              # y_position = 42,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(n.breaks = 10, name = "Number of missed proteins") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))


## CV ------------------------
cv_dt1 <- fread("./CaseStudy1_Quartet/results/tables/2_cv.csv")
rank_dt1 <- fread("./CaseStudy1_Quartet/results/tables/2_rank.csv")
cv_dt2 <- fread("./CaseStudy2_CG/results/tables/1_cv.csv")
rank_dt2 <- fread("./CaseStudy2_CG/results/tables/1_rank.csv")

sub_rank1 <- rank_dt1 %>%
  filter(level %in% "Uncorrected") %>%
  select(feature, scenario, quantitation, rank)

sub_cv_rank1 <- cv_dt1 %>%
  mutate(cv = apply(.[, 6:9], 1, mean, na.rm = TRUE)) %>%
  left_join(., rank_dt1, by = c("feature", "scenario", "level", "quantitation", "correction")) %>%
  select(feature, level, scenario, quantitation, correction, cv, rank)

sub_rank2 <- rank_dt2 %>%
  filter(level %in% "Uncorrected") %>%
  select(feature, quantitation, rank)

sub_cv_rank2 <- cv_dt2 %>%
  mutate(cv = apply(.[, 5:7], 1, mean, na.rm = TRUE)) %>%
  left_join(., sub_rank2, by = c("feature", "quantitation")) %>%
  mutate(scenario = "Random design") %>%
  select(feature, level, scenario, quantitation, correction, cv, rank)

sub_cv_rank <- rbind(sub_cv_rank1, sub_cv_rank2)

lg_level <- ggplot(sub_cv_rank2, aes(x = correction, y = cv)) +
  geom_col(aes(fill = level)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  theme_bw()

p_cv_list <- pblapply(unique(scenarios), function(scenario_id) {
  
  sub_tmp <- sub_cv_rank %>%
    filter(scenario %in% scenario_id) %>%
    filter(rank <= 5000) %>%
    filter(cv < 100) %>%
    select(feature, scenario, level, quantitation, correction, cv, rank) %>%
    mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
    tidyr::separate_rows(correction, sep = "_") %>%
    mutate(level = factor(level, levels = unique(level))) %>%
    mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
    mutate(correction = factor(correction, levels = unique(correct_methods)))
  
  p_title <- ggplot(sub_tmp) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 24),
          plot.margin = unit(c(.5, 1.2, 0, 2.2), "cm"))+
    facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))
  
  axis_y <- "CV (%)"
  axis_x <- "Rank"
  
  p_tmp <- ggplot(sub_tmp, aes(x = rank, y = cv)) +
    geom_point(aes(color = level), alpha = .03, shape = 16) +
    geom_smooth(aes(color = level), method = "gam", formula = y ~ s(x, bs = "cs")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 20),
          axis.text.x = element_text(size = 14),
          axis.title.y = element_text(size = 20),
          axis.text.y = element_text(size = 14),
          strip.text = element_text(size = 20),
          strip.background = element_blank(),
          plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsLevel) +
    scale_y_continuous(limits = c(0, 100), n.breaks = 5, name = axis_y) +
    scale_x_continuous(limits = c(0, 5000), name = axis_x) +
    facet_grid(cols = vars(correction), rows = vars(quantitation))
  
  p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
  
})


## PCA: Quartet -------------------------------
sub_pca_qc_final <- fread("CaseStudy1_Quartet/results/tables/3_pca_qc.csv")
sub_pca_qc_props_final <- fread("CaseStudy1_Quartet/results/tables/3_pca_qc_props.csv")
sub_pca_snr <- fread("CaseStudy1_Quartet/results/tables/3_pca_snr.csv")

lg_pca_batch <- ggplot(sub_pca_qc_final) +
  geom_point(aes(x = PC1, y = PC2, shape = batch), size = 4) +
  scale_shape_manual(values = dictShapesBatch1) +
  labs(shape = "Batch") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

lg_pca_sample <- ggplot(sub_pca_qc_final) +
  geom_point(aes(x = PC1, y = PC2, color = sample), size = 4) +
  scale_color_manual(values = dictColorsSample1) +
  labs(color = "Sample") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

lg_method <- ggplot(sub_pca_qc_final,aes(x = PC1, y = PC2)) +
  geom_col(aes(fill = correction)) +
  scale_fill_manual(values = dictColorsMethod, name = "BECA") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

lg_level <- ggplot(sub_pca_qc_final,aes(x = PC1, y = PC2)) +
  geom_col(aes(fill = level)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  theme_bw() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

p_pca_list1 <- pblapply(unique(scenarios)[1:2], function(scenario_id) {
  
  sub_pca_tmp <- sub_pca_qc_final %>%
    filter(all %in% "All") %>%
    filter(scenario %in% scenario_id)
  
  sub_snr_tmp <- sub_pca_snr %>%
    filter(all %in% "All") %>%
    filter(scenario %in% scenario_id)
  
  p_title <- ggplot(sub_pca_tmp) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 30),
          plot.margin = unit(c(.5, 1.3, 0, 1.4), "cm"))+
    facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))
  
  p_pca_list <- mclapply(unique(data_levels), function(level_id) {
    
    sub_tmp <- sub_pca_tmp %>%
      mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
      mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods))) %>%
      filter(level %in% level_id)
    
    sub_snr <- sub_snr_tmp %>%
      mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
      mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods))) %>%
      filter(level %in% level_id)
    
    p_title <- ggplot(sub_tmp) +
      theme_classic() +
      theme(line = element_blank(), 
            strip.text = element_text(size = 24),
            plot.margin = unit(c(.5, 1.4, 0, 1.35), "cm"))+
      facet_grid(cols = vars(level))
    
    p_tmp <- ggplot(sub_tmp,aes(x = PC1, y = PC2)) +
      geom_point(aes(color = sample, shape = batch), size = 5) +
      geom_text(x = 0, y = max(sub_pca_tmp$PC2) * 0.7, aes(label = sprintf("SNR = %.2f", snr)), data = sub_snr, size = 6) +
      theme_bw() +
      theme(legend.position = "none",
            title = element_text(size = 20),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.text = element_text(size = 20),
            strip.background = element_blank(),
            plot.margin = unit(c(.5,.5, .5, .5), "cm")) +
      scale_shape_manual(values = dictShapesBatch1) +
      scale_color_manual(values = dictColorsSample1) +
      scale_y_continuous(limits = c(-60, 60)) +
      labs(color = "Sample", shape = "Batch") +
      facet_grid(cols = vars(correction), rows = vars(quantitation), drop = TRUE)
    
    p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
    
    return(p_final)
  })
  
  p_pca_all <- plot_grid(plotlist = p_pca_list, ncol = 3, rel_widths = c(1, 3, 3))
  
  p_pca_final <- plot_grid(p_title, p_pca_all, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_pca_final)
  
})


## PCA: CGZ -----------------------
## all samples
pca_dt2 <- fread("./CaseStudy2_CG/results/tables/2_pca.csv")
pca_props_dt2 <- fread("./CaseStudy2_CG/results/tables/2_pca_props.csv")

sub_pca <- pca_dt2 %>%
  mutate(scenario = "Random design") %>%
  filter(all %in% "All") %>%
  filter(!library %in% rm_ids)

lg_pca_batch21 <- ggplot(sub_pca,aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample), size = 5) +
  theme_bw() +
  scale_color_manual(values = dictColorsBatch2, name = "Batch")

p_title <- ggplot(sub_pca) +
  theme_classic() +
  theme(line = element_blank(), 
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, 1.2, 0, 1.9), "cm"))+
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_pca_list <- pblapply(unique(data_levels), function(level_id) {
  
  sub_tmp <- sub_pca %>%
    mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
    filter(level %in% level_id)
  
  p_title <- ggplot(sub_tmp) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 24),
          plot.margin = unit(c(.5, 1.2, 0, 1.9), "cm"))+
    facet_grid(cols = vars(level))
  
  p_tmp <- ggplot(sub_tmp,aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample), size = 5) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          strip.text = element_text(size = 20),
          strip.background = element_blank(),
          plot.margin = unit(c(.5,.5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsBatch2, name = "Batch") +
    # scale_y_continuous(limits = c(-10, 10)) +
    facet_grid(cols = vars(correction), rows = vars(quantitation), drop = TRUE)
  
  p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
})

p_pca21 <- plot_grid(plotlist = p_pca_list, ncol = 3, rel_widths = c(1, 3, 3))

p_pca21_final <- plot_grid(p_title, p_pca21, nrow = 2, rel_heights = c(.1, 1))

## QC samples
pca_dt2 <- fread("./CaseStudy2_CG/results/tables/3_pca_qc.csv")
pca_props_dt2 <- fread("./CaseStudy2_CG/results/tables/3_pca_qc_props.csv")

sub_pca <- pca_dt2 %>%
  mutate(scenario = "Random design") %>%
  filter(all %in% "All") %>%
  filter(!library %in% rm_ids)

lg_pca_batch22 <- ggplot(sub_pca,aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = Cleaning), size = 5) +
  theme_bw() +
  scale_shape_manual(values = 15:17, name = "Batch")

lg_pca_sample22 <- ggplot(sub_pca,aes(x = PC1, y = PC2)) +
  geom_point(aes(color = sample), size = 5) +
  theme_bw() +
  scale_color_manual(values = dictColorsSample2, name = "Sample")

p_title <- ggplot(sub_pca) +
  theme_classic() +
  theme(line = element_blank(), 
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, 1.2, 0, 2), "cm"))+
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_pca_list2 <- pblapply(unique(data_levels), function(level_id) {
  
  sub_tmp <- sub_pca %>%
    mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
    filter(level %in% level_id)
  
  p_title <- ggplot(sub_tmp) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 20),
          plot.margin = unit(c(.5, 1.2, 0, 2), "cm"))+
    facet_grid(cols = vars(level))
  
  p_tmp <- ggplot(sub_tmp,aes(x = PC1, y = PC2)) +
    geom_point(aes(color = sample, shape = Cleaning), size = 5) +
    theme_bw() +
    theme(legend.position = "none",
          title = element_text(size = 20),
          strip.text = element_text(size = 20),
          strip.background = element_blank(),
          plot.margin = unit(c(.5,.5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsSample2) +
    scale_shape_manual(values = 15:17) + 
    # scale_y_continuous(limits = c(-10, 10)) +
    facet_grid(cols = vars(correction), rows = vars(quantitation), drop = TRUE)
  
  p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
})

p_pca22 <- plot_grid(plotlist = p_pca_list2, ncol = 3, rel_widths = c(1, 3, 3))

p_pca22_final <- plot_grid(p_title, p_pca22, nrow = 2, rel_heights = c(.1, 1))


## SNR ----------------------------
sub_pca_snr1 <- fread("./CaseStudy1_Quartet/results/tables/3_pca_snr.csv")
sub_pca_snr2 <- fread("./CaseStudy2_CG/results/tables/3_pca_snr.csv")

sub_pca_snr2 <- sub_pca_snr2 %>%
  mutate(scenario = "Random design") %>%
  select(Inter, Intra, snr, scenario, everything())

sub_pca_snr <- sub_pca_snr1 %>%
  rbind(., sub_pca_snr2)

p_snr_list <- pblapply(unique(scenarios), function(scenario_id) {
  
  sub_tmp0 <- sub_pca_snr %>%
    filter(all %in% "All") %>%
    filter(scenario %in% scenario_id) %>%
    mutate(level = factor(level, levels = unique(data_levels))) %>%
    mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
    mutate(correction = factor(correction, levels = unique(correct_methods))) %>%
    mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods)))
  
  p_title <- ggplot(sub_tmp0) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 30),
          plot.margin = unit(c(.5, .5, 0, 2.7), "cm")) +
    facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))
  
  sub_tmp <- sub_tmp0 %>%
    mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
    mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
    arrange(correction, level) %>%
    mutate(name = paste(level, correction, quantitation, sep = "_"))
  
  dictLabelsName_tmp <- rep(c(1, 3, 4, 6, 7, 9, 10, 12, 13), each = 3)
  names(dictLabelsName_tmp) <- unique(sub_tmp$name)
  
  sub_tmp <- sub_tmp %>%
    mutate_at("name", ~ dictLabelsName_tmp[.]) %>%
    mutate_at("name", as.factor)
  
  axis_y <- "SNR"
  
  p_tmp <- ggplot(sub_tmp, aes(x = correction, y = snr)) +
    geom_text(aes(label = sprintf("%.1f", snr),
                  hjust = sapply(level, function(x) {
                    if (as.character(x) %in% "Peptide-corrected") {
                      hjust_x <- 1.25
                    } else if (as.character(x) %in% "Protein-corrected") {
                      hjust_x <- -.25
                    } else {
                      hjust_x <- .5
                    }
                  })),
              vjust = -.6, size = 6) +
    geom_col(aes(fill = level), width = .85, position = position_dodge(width = .9)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 24),
          strip.text = element_text(size = 24),
          plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
    scale_fill_manual(values = dictColorsLevel) +
    scale_y_continuous(limits = c(-1, 30), breaks = seq(0, 30, 5), name = axis_y) +
    facet_grid(cols = vars(quantitation))
  
  p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
  
})

sub_pca_snr$snr <- round(sub_pca_snr$snr, 2)

sub_pca_snr <- sub_pca_snr %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate_at("snr", ~ round(., 2)) %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

sub_pca_snr11 <- sub_pca_snr %>%
  filter(all %in% c("All", "Qualified"))

sub_pca_snr_kw <- sub_pca_snr11 %>%
  group_by(all, correction) %>%
  kruskal_test(snr ~ level) %>%
  adjust_pvalue(method = "fdr") %>%
  mutate(significance = sapply(p.adj, function(x) {
    if(as.numeric(x) < 0.01) {
      a <- "**"
    } else if (as.numeric(x) < 0.05) {
      a <- "*"
    } else {
      a <- "NS"
    }
    return(a)
    })
  )

sub_pca_snr_anova <- sub_pca_snr11 %>%
  group_by(all, correction) %>%
  anova_test(snr ~ level) %>%
  as.data.frame %>%
  mutate(significance = ifelse(p < .05, "*", "NS"))

p_snr_box11 <- ggplot(sub_pca_snr11, aes(x = all, y = snr)) +
  geom_boxplot(aes(fill = level), width = .7, position = position_dodge(width = .8)) +
  geom_text(y = 25, aes(x = all, label = ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.2f", p))), vjust = 0,
            data = sub_pca_snr_anova, size = 8) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(-2.5, 30), n.breaks = 10, name = "SNR") +
  facet_grid(cols = vars(correction))

sub_pca_snr12 <- sub_pca_snr %>%
  filter(!all %in% c("All", "Qualified"))

sub_pca_snr_kw <- sub_pca_snr12 %>%
  group_by(all, correction) %>%
  kruskal_test(snr ~ level) %>%
  adjust_pvalue(method = "fdr") %>%
  mutate(significance = sapply(p.adj, function(x) {
    if(as.numeric(x) < 0.01) {
      a <- "**"
    } else if (as.numeric(x) < 0.05) {
      a <- "*"
    } else {
      a <- "NS"
    }
    return(a)
  })
  )

sub_pca_snr_anova <- sub_pca_snr12 %>%
  group_by(all, correction) %>%
  anova_test(snr ~ level) %>%
  as.data.frame %>%
  mutate(significance = ifelse(p < .05, "*", "NS"))

p_snr_box12 <- ggplot(sub_pca_snr12, aes(x = all, y = snr)) +
  geom_boxplot(aes(fill = level), width = .7, position = position_dodge(width = .8)) +
  geom_text(y = 25, aes(x = all, label = sprintf("p = %.4f", p)), vjust = 0,
            data = sub_pca_snr_anova, size = 8) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(-2.5, 30), n.breaks = 10, name = "SNR") +
  facet_grid(cols = vars(correction))

sub_pca_snr13 <- sub_pca_snr1 %>%
  rbind(., sub_pca_snr2) %>%
  filter(all %in% "All") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

p_snr_bar1 <- ggplot(sub_pca_snr13, aes(x = level, y = snr)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Uncorrected"),
                                 c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected")),
              map_signif_level = TRUE,
              y_position = 15,
              step_increase = .06,
              tip_length = .01,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = FALSE)) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(limits = c(-1, 20), n.breaks = 10, name = "SNR") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_snr_bar2 <- ggplot(sub_pca_snr13, aes(x = correction, y = snr)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Ratio", "Log transformed"),
                                 c("Ratio", "Median centering"),
                                 c("Ratio", "ComBat"),
                                 c("Ratio", "RUV-III-C")),
              map_signif_level = TRUE,
              y_position = 18,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod) +
  scale_y_continuous(limits = c(-1, 30), breaks = seq(0, 30, 5), name = "SNR") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_snr_bar3 <- ggplot(sub_pca_snr13, aes(x = quantitation, y = snr)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = 18,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro) +
  scale_y_continuous(limits = c(-1, 25), breaks = seq(0, 25, 5), name = "SNR") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))


## BRinfo: CGZ ------
brinfo_dt <- fread("./CaseStudy2_CG/results/tables/2_brinfo.csv")

sub_brinfo <- brinfo_dt %>%
  filter(all %in% "All") %>%
  mutate(scenario = "Random design")

lg_level <- ggplot(sub_brinfo, aes(x = correction, y = brinfo)) +
  geom_col(aes(fill = level), color = "black") +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  theme_bw()

p_title <- ggplot(sub_brinfo) +
  theme_classic() +
  theme(line = element_blank(), 
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 2.7), "cm"))+
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

sub_tmp <- sub_brinfo %>%
  mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods)))  %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))%>%
  arrange(correction, level) %>%
  mutate(name = paste(level, correction, quantitation, sep = "_"))

dictLabelsName_tmp <- rep(c(1, 3, 4, 6, 7, 9, 10, 12, 13), each = 3)
names(dictLabelsName_tmp) <- unique(sub_tmp$name)

sub_tmp <- sub_tmp %>%
  mutate_at("name", ~ dictLabelsName_tmp[.]) %>%
  mutate_at("name", as.factor)

axis_y <- "BRinfo (%)"

p_tmp <- ggplot(sub_tmp, aes(x = correction, y = brinfo * 100)) +
  geom_text(aes(label = sprintf("%.1f", brinfo * 100),
                hjust = sapply(level, function(x) {
                  if (as.character(x) %in% "Peptide-corrected") {
                    hjust_x <- 1.25
                  } else if (as.character(x) %in% "Protein-corrected") {
                    hjust_x <- -.25
                  } else {
                    hjust_x <- .5
                  }
                })),
            vjust = -.6, size = 6) +
  geom_col(aes(fill = level), width = .85, position = position_dodge(width = .9)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 3), name = axis_y) +
  facet_grid(cols = vars(quantitation))

p_brinfo <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))

sub_brinfo <- brinfo_dt %>%
  mutate(scenario = "Random design") %>%
  mutate_at("brinfo", ~ round(., 2)) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

sub_brinfo_11 <- sub_brinfo %>%
  filter(all %in% c("All", "Qualified"))

sub_brinfo_anova <- sub_brinfo_11 %>%
  group_by(all, correction) %>%
  anova_test(brinfo ~ level) %>%
  as.data.frame %>%
  mutate(significance = ifelse(p < .05, "*", "NS"))

p_brinfo_box11 <- ggplot(sub_brinfo_11, aes(x = all, y = brinfo)) +
  geom_boxplot(aes(fill = level), width = .7, position = position_dodge(width = .8)) +
  geom_text(y = .2, aes(x = all, label = ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.2f", p))), vjust = 0,
            data = sub_brinfo_anova, size = 8) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, .25), n.breaks = 10, name = "BRinfo") +
  facet_grid(cols = vars(correction))

sub_brinfo_12 <- sub_brinfo %>%
  filter(!all %in% c("All", "Qualified"))

sub_brinfo_anova <- sub_brinfo_12 %>%
  group_by(all, correction) %>%
  anova_test(brinfo ~ level) %>%
  as.data.frame %>%
  mutate(significance = ifelse(p < .05, "*", "NS"))

p_brinfo_box12 <- ggplot(sub_brinfo_12, aes(x = all, y = brinfo)) +
  geom_boxplot(aes(fill = level), width = .7, position = position_dodge(width = .8)) +
  geom_text(y = .2, aes(x = all, label = sprintf("p = %.4f", p)), vjust = 0,
            data = sub_brinfo_anova, size = 8) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, .25), n.breaks = 10, name = "BRinfo") +
  facet_grid(cols = vars(correction))

sub_brinfo_13 <- brinfo_dt %>%
  filter(all %in% c("All")) %>%
  mutate(scenario = "Random design") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

p_brinfo_bar1 <- ggplot(sub_brinfo_13, aes(x = level, y = brinfo * 100)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Uncorrected"),
                                 c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected")),
              map_signif_level = TRUE,
              y_position = 12,
              step_increase = .1,
              tip_length = .02,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(limits = c(0, 16), n.breaks = 10, name = "BRinfo (%)") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_brinfo_bar2 <- ggplot(sub_brinfo_13, aes(x = correction, y = brinfo * 100)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Ratio", "Log transformed"),
                                 c("Median centering", "Log transformed"),
                                 c("ComBat", "Log transformed"),
                                 c("RUV-III-C", "Log transformed")),
              map_signif_level = TRUE,
              y_position = 12,
              step_increase = .1,
              tip_length = .02,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod) +
  scale_y_continuous(limits = c(0, 17), n.breaks = 10, name = "BRinfo (%)") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_brinfo_bar3 <- ggplot(sub_brinfo_13, aes(x = quantitation, y = brinfo * 100)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = 4,
              step_increase = .1,
              tip_length = .02,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro) +
  scale_y_continuous(limits = c(0, 8), n.breaks = 10, name = "BRinfo (%)") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))


## DEPs: ----------------------------
ref_dt1 <- fread("./CaseStudy1_Quartet/results/tables/4_with_reference_dataset.csv")
ref_dt2 <- fread("./CaseStudy2_CG/results/tables/4_with_reference_dataset.csv")

ref_dt <- rbind(ref_dt1, ref_dt2)

ref_dt <- ref_dt %>%
  filter(all %in% "All")

sub_tp_delta <- ref_dt %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Median centering_Ratio_RUV-III-C", .)) %>%
  tidyr::separate_rows("correction", sep = "_") %>%
  reshape2::dcast(., pair + scenario + correction + quantitation ~ level, value.var = "TP") %>%
  mutate_at("Protein-corrected", ~. - `Uncorrected`) %>%
  mutate_at("Peptide-corrected", ~ . - `Uncorrected`) %>%
  select(!Uncorrected) %>%
  reshape2::melt(., id = 1:4) %>%
  mutate(class = sapply(value, function(x) {
    if (as.numeric(x) > 0) {
      l <- "Increased"
    } else if (as.numeric(x) == 0) {
      l <- "Not changed"
    } else if (as.numeric(x) < 0) {
      l <- "Decreased"
    }
  }))

nrow(sub_tp_delta[sub_tp_delta$scenario %in% "Confounded design" &
                    sub_tp_delta$variable %in% "Protein-corrected",])
nrow(sub_tp_delta[sub_tp_delta$scenario %in% "Confounded design" &
                    !(sub_tp_delta$class %in% "Decreased") &
                    sub_tp_delta$variable %in% "Protein-corrected",])

p_tp_bar <- ggplot(sub_tp_delta, aes(x = variable)) +
  geom_bar(aes(fill = class), width = .7, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = c("Decreased" = "#66C2A5",
                               "Increased" = "#E78AC3",
                               "Not changed" = "grey"),
                    name = "TP DEPs") +
  scale_y_continuous(labels = scales::percent,
                     name = "Proportion of combinations") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

sub_fp_delta <- ref_dt %>%
  mutate_at("FP", ~ ifelse(is.na(.), 0, .)) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Median centering_Ratio_RUV-III-C", .)) %>%
  tidyr::separate_rows("correction", sep = "_") %>%
  reshape2::dcast(., pair + scenario + correction + quantitation ~ level, value.var = "FP") %>%
  mutate_at("Protein-corrected", ~. - `Uncorrected`) %>%
  mutate_at("Peptide-corrected", ~ . - `Uncorrected`) %>%
  select(!Uncorrected) %>%
  reshape2::melt(., id = 1:4) %>%
  mutate(class = sapply(value, function(x) {
    if (as.numeric(x) > 0) {
      l <- "Increased"
    } else if (as.numeric(x) == 0) {
      l <- "Not changed"
    } else if (as.numeric(x) < 0) {
      l <- "Decreased"
    }
  }))

nrow(sub_fp_delta[sub_fp_delta$scenario %in% "Confounded design" &
                    sub_fp_delta$variable %in% "Protein-corrected",])
nrow(sub_fp_delta[sub_fp_delta$scenario %in% "Confounded design" &
                    sub_fp_delta$class %in% "Increased" &
                    sub_fp_delta$variable %in% "Protein-corrected",])

p_fp_bar <- ggplot(sub_fp_delta, aes(x = variable)) +
  geom_bar(aes(fill = class), width = .7, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = c("Decreased" = "#66C2A5",
                               "Increased" = "#E78AC3",
                               "Not changed" = "grey"),
                    name = "FP DEPs") +
  scale_y_continuous(labels = scales::percent,
                     name = "Proportion of combinations") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

sub_fn_delta <- ref_dt %>%
  mutate_at(8:10, ~ ifelse(is.na(.), 0, .)) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Median centering_Ratio_RUV-III-C", .)) %>%
  tidyr::separate_rows("correction", sep = "_") %>%
  reshape2::dcast(., pair + scenario + correction + quantitation ~ level, value.var = "FN") %>%
  mutate_at("Protein-corrected", ~. - `Uncorrected`) %>%
  mutate_at("Peptide-corrected", ~ . - `Uncorrected`) %>%
  select(!Uncorrected) %>%
  reshape2::melt(., id = 1:4) %>%
  mutate(class = sapply(value, function(x) {
    if (as.numeric(x) > 0) {
      l <- "Increased"
    } else if (as.numeric(x) == 0) {
      l <- "Not changed"
    } else if (as.numeric(x) < 0) {
      l <- "Decreased"
    }
  }))

p_fn_bar <- ggplot(sub_fn_delta, aes(x = variable)) +
  geom_bar(aes(fill = class), width = .7, position = "fill") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = c("Decreased" = "#66C2A5",
                               "Increased" = "#E78AC3",
                               "Not changed" = "grey"),
                    name = "FN DEPs") +
  scale_y_continuous(labels = scales::percent,
                     name = "Proportion of combinations") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

sub_ref <- ref_dt %>%
  select(!rc & !rmse & !precision & !recall & !f1) %>%
  reshape2::melt(., id = 1:6) %>%
  mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_dep_bar0_list <- pblapply(c("TP", "FN", "FP"), function(class_id) {
  
  sub_ref_i <- sub_ref %>%
    filter(variable %in% class_id)
  
  p_dep_bar0 <- ggplot(sub_ref_i, aes(x = level, y = value)) +
    geom_boxplot(aes(fill = level), width = .7) +
    geom_signif(comparisons = list(c("Uncorrected", "Peptide-corrected"),
                                   c("Peptide-corrected", "Protein-corrected"),
                                   c("Uncorrected", "Protein-corrected")),
                map_signif_level = TRUE,
                step_increase = .1,
                # tip_length = .05,
                textsize = 4,
                test = "t.test",
                test.args = list(pairwise = TRUE)) +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_text(size =20),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 24),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = dictColorsLevel) +
    scale_y_continuous(n.breaks = 5, name = "Number of DEPs") +
    facet_grid(cols = vars(scenario), rows = vars(variable), scales = "free_y",
               labeller = labeller(scenario = dictLabelsScenario))
  
  return(p_dep_bar0)
})

p_dep_bar0 <- plot_grid(p_dep_bar0_list[[1]] +
                         theme(axis.text.x = element_blank()),
                       p_dep_bar0_list[[2]] +
                         theme(axis.text.x = element_blank(),
                               strip.text.x = element_blank()),
                       p_dep_bar0_list[[3]] +
                         theme(strip.text.x = element_blank(),
                               plot.margin = unit(c(.2, .2, 0, .5), units = "cm")),
                       nrow = 3, rel_heights = c(1.5, 1.3, 1.8))

p_dep_bar1_list <- pblapply(c("TP", "FN", "FP"), function(class_id) {
  
  sub_ref_i <- sub_ref %>%
    filter(variable %in% class_id)
  
  p_dep_bar1 <- ggplot(sub_ref_i, aes(x = correction, y = value)) +
    stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
    geom_signif(comparisons = list(c("Ratio", "Log transformed"),
                                   c("Ratio", "Median centering"),
                                   c("Ratio", "ComBat"),
                                   c("Ratio", "RUV-III-C")),
                map_signif_level = TRUE,
                step_increase = .1,
                textsize = 4,
                test = "t.test",
                test.args = list(pairwise = TRUE)) +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_text(size =20),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = dictColorsMethod) +
    scale_y_continuous(n.breaks = 5, name = "Number of DEPs") +
    facet_grid(cols = vars(scenario), rows = vars(variable), scales = "free_y",
               labeller = labeller(scenario = dictLabelsScenario))
  
  return(p_dep_bar1)
  
})

p_dep_bar1 <- plot_grid(p_dep_bar1_list[[1]] +
                          theme(axis.text.x = element_blank()),
                        p_dep_bar1_list[[2]] +
                          theme(axis.text.x = element_blank(),
                                strip.text.x = element_blank()),
                        p_dep_bar1_list[[3]] +
                          theme(strip.text.x = element_blank(),
                                plot.margin = unit(c(.2, .2, 0, .5), units = "cm")),
                        nrow = 3, rel_heights = c(1.5, 1.3, 1.8))

p_dep_bar2_list <- pblapply(c("TP", "FN", "FP"), function(class_id) {
  
  sub_ref_i <- sub_ref %>%
    filter(variable %in% class_id)
  
  p_dep_bar2 <- ggplot(sub_ref_i, aes(x = quantitation, y = value)) +
    stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
    geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                   c("MaxLFQ", "TopPep3"),
                                   c("iBAQ", "TopPep3")),
                map_signif_level = TRUE,
                step_increase = .1,
                textsize = 4,
                test = "t.test",
                test.args = list(pairwise = TRUE)) +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_text(size =20),
          legend.text = element_text(size = 16),
          strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
    scale_fill_manual(values = dictColorsPep2pro) +
    scale_y_continuous(n.breaks = 5, name = "Number of DEPs") +
    facet_grid(cols = vars(scenario), rows = vars(variable), scales = "free_y",
               labeller = labeller(scenario = dictLabelsScenario))
  
  return(p_dep_bar2)
})

p_dep_bar2 <- plot_grid(p_dep_bar2_list[[1]] +
                          theme(axis.text.x = element_blank()),
                        p_dep_bar2_list[[2]] +
                          theme(axis.text.x = element_blank(),
                                strip.text.x = element_blank()),
                        p_dep_bar2_list[[3]] +
                          theme(strip.text.x = element_blank(),
                                plot.margin = unit(c(.2, .2, 0, .5), units = "cm")),
                        nrow = 3, rel_heights = c(1.5, 1.3, 1.8))


## F1 score ------------------
ref_dt1 <- fread("./CaseStudy1_Quartet/results/tables/4_with_reference_dataset.csv")
ref_dt2 <- fread("./CaseStudy2_CG/results/tables/4_with_reference_dataset.csv")

ref_dt <- rbind(ref_dt1, ref_dt2)

sub_pr <- ref_dt %>%
  mutate(name = paste(scenario, level, correction, quantitation, sep = "_")) %>%
  mutate(name = dictLabelsName[name])

p_f1_list <- pblapply(unique(scenarios), function(scenario_id) {
  
  sub_tmp0 <- sub_pr %>%
    filter(scenario %in% scenario_id) %>%
    mutate(level = factor(level, levels = unique(data_levels))) %>%
    mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
    mutate(correction = factor(correction, levels = unique(correct_methods)))
  
  p_title <- ggplot(sub_tmp0) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 30),
          plot.margin = unit(c(.5, .5, 0, 3.1), "cm"))+
    facet_grid(cols = vars(scenario),
               labeller = labeller(scenario = dictLabelsScenario))
  
  sub_tmp <- sub_tmp0 %>%
    mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
    mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
    arrange(correction, level) %>%
    mutate(name = paste(level, correction, quantitation, sep = "_"))
  
  dictLabelsName_tmp <- rep(c(1, 3, 4, 6, 7, 9, 10, 12, 13), each = 3)
  names(dictLabelsName_tmp) <- unique(sub_tmp$name)
  
  sub_tmp <- sub_tmp %>%
    mutate_at("name", ~ dictLabelsName_tmp[.]) %>%
    mutate_at("name", as.factor)
  
  axis_y <- "F1 score (%)"
  
  p_tmp <- ggplot(sub_tmp, aes(x = correction, y = f1 * 100)) +
    stat_summary(aes(fill = level), fun = mean,
                 geom = "col", width = .7, position = position_dodge(width = .9)) +
    stat_summary(aes(group = level), fun.data = "mean_se",
                 geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 24),
          strip.text = element_text(size = 24),
          plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
    scale_fill_manual(values = dictColorsLevel) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10), name = axis_y) +
    facet_grid(cols = vars(quantitation))
  
  p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
  
})

p_f1_list2 <- pblapply(unique(scenarios), function(scenario_id) {
  
  sub_tmp0 <- sub_pr %>%
    filter(scenario %in% scenario_id) %>%
    mutate(level = factor(level, levels = unique(data_levels))) %>%
    mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods)))
  
  p_title <- ggplot(sub_tmp0) +
    theme_classic() +
    theme(line = element_blank(), 
          strip.text = element_text(size = 30),
          plot.margin = unit(c(.5, .5, 0, 3.1), "cm"))+
    facet_grid(cols = vars(scenario),
               labeller = labeller(scenario = dictLabelsScenario))
  
  sub_tmp <- sub_tmp0 %>%
    mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
    mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Median centering_Ratio_RUV-III-C", .)) %>%
    tidyr::separate_rows("correction", sep = "_") %>%
    arrange(correction, level) %>%
    mutate(name = paste(level, correction, quantitation, sep = "_"))
  
  dictLabelsName_tmp <- rep(c(1, 2, 3, 4), each = 9)
  names(dictLabelsName_tmp) <- unique(sub_tmp$name)
  
  sub_tmp <- sub_tmp %>%
    mutate_at("name", ~ dictLabelsName_tmp[.])
  
  axis_y <- "F1 score (%)"
  
  p_tmp <- ggplot(sub_tmp, aes(x = correction, y = f1 * 100)) +
    stat_summary(aes(group = level), fun.data = "mean_se",
                 geom = "errorbar", size = 1, width = .5, position = position_dodge(width = 0)) +
    stat_summary(aes(group = level, color = level, lty = level), fun = mean,
                 geom = "line", size = 1.5) +
    stat_summary(aes(color = level, shape = level), fun = mean,
                 geom = "point", size = 5, position = position_dodge(width = 0)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_text(size = 24),
          strip.text = element_text(size = 24),
          plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
    scale_color_manual(values = dictColorsLevel) +
    scale_linetype_manual(values = c("Uncorrected" = 3,
                                     "Peptide-corrected" = 4,
                                     "Protein-corrected" = 1)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10), name = axis_y) +
    facet_grid(cols = vars(quantitation))
  
  p_final <- plot_grid(p_title, p_tmp, nrow = 2, rel_heights = c(.1, 1))
  
  return(p_final)
  
})

sub_pr <- sub_pr %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

sub_pr11 <- sub_pr %>%
  filter(all %in% c("All", "Qualified"))

sub_f1_anova <- sub_pr11 %>%
  group_by(all, correction) %>%
  anova_test(f1 ~ level) %>%
  as.data.frame %>%
  mutate(significance = ifelse(p < .05, "*", "NS"))

p_f1_box11 <- ggplot(sub_pr11, aes(x = all, y = f1)) +
  geom_boxplot(aes(fill = level), width = .7, position = position_dodge(width = .8)) +
  geom_text(y = 1.05, aes(x = all, label = ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.2f", p))), vjust = 0,
            data = sub_f1_anova, size = 8) +
  # ggpubr::stat_anova_test(aes(group = level), size = 5) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, .1), name = "F1 score") +
  facet_grid(cols = vars(correction))

sub_pr12 <- sub_pr %>%
  filter(all %in% c("All")) %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

p_f1_bar12 <- ggplot(sub_pr12, aes(x = level, y = f1)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7, position = "dodge") +
  stat_summary(aes(group = level), fun.data = "mean_se", geom = "errorbar", width = 0.2, position = position_dodge(width = .7)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected")),
              map_signif_level = TRUE,
              y_position = 1,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, .1), name = "F1 score") +
  facet_grid(cols = vars(correction))

sub_pr2 <- ref_dt %>%
  filter(all %in% "All") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

p_f1_bar1 <- ggplot(sub_pr2, aes(x = level, y = f1)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7, position = "dodge") +
  stat_summary(aes(group = level), fun.data = "mean_se", geom = "errorbar", width = 0.2, position = position_dodge(width = .7)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected")),
              map_signif_level = TRUE,
              y_position = 1,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, .1), name = "F1 score") +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

p_f1_bar2 <- ggplot(sub_pr2, aes(x = correction, y = f1)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Ratio", "Log transformed"),
                                 c("Log transformed", "Median centering"),
                                 c("Log transformed", "ComBat"),
                                 c("Log transformed", "RUV-III-C")),
              map_signif_level = TRUE,
              y_position = 1,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod, name = "BECA") +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1, .1), name = "F1 score") +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

p_f1_bar3 <- ggplot(sub_pr2, aes(x = quantitation, y = f1)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = .95,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro, name = "PQM") +
  scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, .1), name = "F1 score") +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))


## RC -------------------
withref_dt1 <- fread("./CaseStudy1_Quartet/results/tables/4_with_reference_dataset.csv")
withref_dt2 <- fread("./CaseStudy2_CG/results/tables/4_with_reference_dataset.csv")

withref_dt <- rbind(withref_dt1, withref_dt2)

## barplot
sub_rc <- withref_dt %>%
  mutate(name = paste(scenario, level, correction, quantitation, sep = "_")) %>%
  mutate(name = dictLabelsName[name])

sub_rc11 <- sub_rc %>%
  filter(all %in% c("1 ~ 200", "200 ~ 500", "500 ~")) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

sub_rc_anova <- sub_rc11 %>%
  group_by(all, correction) %>%
  anova_test(rc ~ level) %>%
  as.data.frame %>%
  mutate(significance = ifelse(p < .05, "*", "NS"))

p_rc_box11 <- ggplot(sub_rc11, aes(x = all, y = rc)) +
  geom_boxplot(aes(fill = level), width = .7, position = position_dodge(width = .8)) +
  geom_text(y = 1.05, aes(x = all, label = ifelse(p < 0.0001, "p < 0.0001", sprintf("p = %.2f", p))), vjust = 0,
            data = sub_rc_anova, size = 8) +
  # ggpubr::stat_anova_test(aes(group = level), size = 5) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(.5, 1.15), breaks = seq(0, 1, .1), name = "Correlation with Reference") +
  facet_grid(cols = vars(correction))

sub_rc12 <- sub_rc %>%
  filter(all %in% c("All")) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

p_rc_bar12 <- ggplot(sub_rc12, aes(x = level, y = rc)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7,  position = position_dodge(width = .75)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1), name = "Correlation with Reference") +
  facet_grid(cols = vars(correction))

sub_rc2 <- sub_rc %>%
  mutate(level = factor(level, levels = unique(data_levels))) %>%
  mutate(scenario = factor(scenario, levels = unique(scenarios))) %>%
  mutate(quantitation = factor(quantitation, levels = unique(pep2pro_methods))) %>%
  mutate(correction = factor(correction, levels = unique(correct_methods)))

p_rc_bar1 <- ggplot(sub_rc2, aes(x = level, y = rc)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7, position = "dodge") +
  stat_summary(aes(group = level), fun.data = "mean_se", geom = "errorbar", width = 0.2, position = position_dodge(width = .7)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected")),
              map_signif_level = TRUE,
              y_position = .96,
              step_increase = .05,
              tip_length = .007,
              textsize = 4,
              test = "t.test") +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel, name = "Level") +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, .1), name = "Correlation with Reference") +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

p_rc_bar2 <- ggplot(sub_rc2, aes(x = correction, y = rc)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Ratio", "Log transformed"),
                                 c("Log transformed", "Median centering"),
                                 c("Log transformed", "ComBat"),
                                 c("Log transformed", "RUV-III-C")),
              map_signif_level = TRUE,
              y_position = .96,
              step_increase = .04,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod, name = "BECA") +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, .1), name = "Correlation with Reference") +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

p_rc_bar3 <- ggplot(sub_rc2, aes(x = quantitation, y = rc)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = .96,
              step_increase = .04,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 20, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro, name = "PQM") +
  scale_y_continuous(limits = c(0, 1.3), breaks = seq(0, 1, .1), name = "Correlation with Reference") +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

## heatmap
sub_ref <- withref_dt %>%
  reshape2::acast(., quantitation ~ scenario + level + correction,
                  value.var = "rc", fun.aggregate = mean, na.rm = TRUE)

sub_ref2 <- sub_ref %>%
  as.data.frame %>%
  # mutate_if(is.numeric, ~ ifelse(is.na(.), 0, .)) %>%
  as.matrix

sub_ref3 <- sub_ref2
sub_ref3[is.na(sub_ref2)] <- ""
sub_ref3[sub_ref2 < .8 & sub_ref2 > .5] <- ""
sub_ref3[sub_ref2 < .9 & sub_ref2 > .8] <- "*"
sub_ref3[sub_ref2 > .9] <- "**"

sub_ref2_annot <- colnames(sub_ref2) %>%
  sapply(., function(x) {unlist(strsplit(x, split = "_"))}) %>%
  t %>%
  as.data.frame

colnames(sub_ref2_annot) <- c("Scenario", "Level", "BECA")

ha_col <- HeatmapAnnotation(df = sub_ref2_annot,
                            col = list(Scenario = dictColorsScenario,
                                       Level = dictColorsLevel,
                                       BECA = dictColorsMethod),
                            which = "col",
                            show_legend = FALSE,
                            gap = unit(5, "points"),
                            show_annotation_name = TRUE,
                            annotation_name_side = "left",
                            gp = gpar(col = "transparent"),
                            annotation_name_gp = gpar(fontsize = 20, col = "black"),
                            simple_anno_size = unit(0.8, "cm"))

fa = rep(c("Multi-lab (Balanced)", "Multi-lab (Confounded)", "Large-scale"), each = 9)
dend1 = cluster_within_group(sub_ref2, fa)

ha <- Heatmap(matrix = sub_ref2,
              name = "RC",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(label = sub_ref3[i, j],
                          x = x, y = y, 
                          gp = gpar(col = "black", fontsize = 24))
              },
              rect_gp = gpar(col = "white", lwd = 2, fontsize = 2),
              heatmap_legend_param = list(direction = "horizontal", font = 40),
              col = colorRamp2(c(.5, .8, 1), c("#B8E186", "white", "#DE77AE")),
              column_split = 3,
              column_dend_height = unit(5, "mm"),
              column_title_gp = gpar(fontsize = 24),
              show_column_dend = TRUE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              cluster_columns = dend1,
              cluster_rows = FALSE,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 20),
              top_annotation = ha_col)
p_ha1 <- grid.grabExpr(draw(ha, heatmap_legend_side = "bottom"))


## 5-fold cross validation (R square) of sex: CGZ ---------------------------
cross_dt <- read_excel("./CaseStudy2_CG/results/tables/no_missing_values/5_models_sex_training_5fold.xlsx")

sub_cross <- cross_dt %>%
  filter(Method %in% c("Decision Tree", "Bootstrap Forest", "Boosted Tree", "K Nearest Neighbors", "Support Vector Machines")) %>%
  tidyr::separate(`Original Data`, c("level", "correction", "quantitation"), sep = ", ") %>%
  mutate_at(13:15, ~ str_extract(., "(?<==).+")) %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_title <- ggplot(sub_cross) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 3.4), "cm"))+
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

axis_y <- "R Square"

p_cross_sex <- ggplot(sub_cross, aes(x = reorder(Method, -`Entropy RSquare`), y = `Entropy RSquare`)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Uncorrected", "Negative control")),
              map_signif_level = TRUE,
              y_position = .4,
              step_increase = .1,
              tip_length = .01,
              textsize = 5,
              test = "wilcox.test") +
  geom_boxplot(aes(fill = Method)) +
  # stat_summary(aes(fill = level), fun = mean,
  #              geom = "col", width = .7, position = position_dodge(width = .9)) +
  # stat_summary(aes(group = level), fun.data = "mean_se",
  #              geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  facet_grid(cols = vars(level)) +
  scale_fill_brewer(palette = "Greens") +
  # scale_fill_manual(values = c(dictColorsLevel, "Negative control" = "#999999")) +
  scale_x_discrete(labels = dictLabelsML) +
  scale_y_continuous(n.breaks = 10, name = axis_y)

# p_cross_sex <- plot_grid(p_title, p_cross_sex, nrow = 2, rel_heights = c(.1, 1))


## MCC of sex: CGZ ---------------------------
mcc_dt <- fread("./CaseStudy2_CG/results/tables/no_missing_values/6_mcc_sex.csv")

sub_mcc <- mcc_dt %>%
  filter(Training.Validation == 2) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods)))

lg_level_mcc <- ggplot(sub_mcc, aes(x = correction, y = mcc)) +
  geom_col(aes(fill = level)) +
  scale_fill_manual(values = c("Negative control" = "#999999", dictColorsLevel))

p_title <- ggplot(sub_mcc) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 3), "cm"))+
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

axis_y <- "MCC"

p_mcc1 <- ggplot(sub_mcc, aes(x = level, y = mcc)) +
  geom_hline(yintercept = .3, color = "black", lty = 2) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Uncorrected", "Negative control")),
              map_signif_level = TRUE,
              step_increase = .2,
              tip_length = .05,
              textsize = 6,
              test = "t.test") +
  stat_summary(aes(fill = level), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = level), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  facet_grid(cols = vars(correction)) +
  scale_fill_manual(values = c(dictColorsLevel, "Negative control" = "#999999")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, .1), name = axis_y)

# p_mcc_final <- plot_grid(p_title, p_mcc, nrow = 2, rel_heights = c(.1, 1))

sub_mcc2 <- mcc_dt %>%
  filter(Training.Validation == 2, !level %in% "Negative control") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_mcc2 <- ggplot(sub_mcc2, aes(x = correction, y = mcc)) +
  stat_summary(aes(fill = correction), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = correction), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  geom_signif(comparisons = list(c("Log transformed", "Ratio"),
                                 c("Log transformed", "Median centering"),
                                 c("Log transformed", "RUV-III-C"),
                                 c("Log transformed", "ComBat")),
              map_signif_level = TRUE,
              step_increase = .2,
              tip_length = .05,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c(dictColorsMethod)) +
  scale_y_continuous(limits = c(-.1, 1), n.breaks = 10, name = axis_y)

p_mcc3 <- ggplot(sub_mcc2, aes(x = quantitation, y = mcc)) +
  stat_summary(aes(fill = quantitation), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = quantitation), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  geom_signif(comparisons = list(c("MaxLFQ", "TopPep3"),
                                 c("MaxLFQ", "iBAQ"),
                                 c("TopPep3", "iBAQ")),
              map_signif_level = TRUE,
              step_increase = .2,
              tip_length = .05,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c(dictColorsPep2pro)) +
  scale_y_continuous(limits = c(-.1, 1), n.breaks = 10, name = axis_y)


## ROC of sex: CGZ ---------------------------
roc_obj <- readRDS("./CaseStudy2_CG/results/tables/no_missing_values/6_roc_sex.rds")

sub_roc <- roc_obj$ROC %>%
  filter(Training.Validation == 2) %>%
  mutate(x = 1 - specificity / 100, y = sensitivity / 100)

sub_delong <- roc_obj$DelongTest %>%
  filter(level1 %in% c("Uncorrected") & level2 %in% c("Peptide-corrected", "Protein-corrected")) %>%
  select(!level1) %>%
  dplyr::rename(level = level2)

sub_auc <- roc_obj$AUC %>%
  left_join(., sub_delong, by = c("Training.Validation", "level", "correction", "quantitation")) %>%
  filter(Training.Validation == 2) %>%
  mutate(y_pos = rep(c(.35, .25, .15, .05), 12))

lg_level_roc <- ggplot(sub_roc, aes(x = x, y = y)) +
  geom_path(aes(color = level), size = 1.5) +
  scale_color_manual(values = c("Negative control" = "#999999", dictColorsLevel))

p_roc <- ggplot(sub_roc, aes(x = x, y = y)) +
  geom_path(aes(color = level), size = 1.5) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype=4)+
  geom_text(x = .7, data = sub_auc, size = 6,
            mapping = aes(y = y_pos,label = ifelse(is.na(p), sprintf("%s = %.2f%%", level, auc),
                                                   sprintf("%s = %.2f%% (%.4f)", level, auc, p)))) +
  geom_text(x = .30, y = 1, data = sub_delong, size = 5,
            mapping = aes(label = sprintf("Delong's Test P value = %.4f", p))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 24),
        strip.text = element_text(size = 30),
        strip.background = element_blank(),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values = c("Negative control" = "#999999", dictColorsLevel)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  facet_grid(cols = vars(correction), rows = vars(quantitation))


## 5-fold cross validation (R square) of age: CGZ ---------------------------
cross_dt <- read_excel("./CaseStudy2_CG/results/tables/no_missing_values/5_models_age_training_5fold.xlsx")

sub_cross <- cross_dt %>%
  filter(!Method %in% c("Fit Least Squares", "Neural Boosted")) %>%
  tidyr::separate(`Original Data`, c("level", "correction", "quantitation"), sep = ", ") %>%
  mutate_at(8:10, ~ str_extract(., "(?<==).+")) %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_title <- ggplot(sub_cross) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 3), "cm"))+
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

axis_y <- "R Square"

p_cross_age <- ggplot(sub_cross, aes(x = reorder(Method, -`RSquare`), y = `RSquare`)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Uncorrected", "Negative control")),
              map_signif_level = TRUE,
              # y_position = 1.2,
              step_increase = .15,
              tip_length = .05,
              textsize = 5,
              test = "t.test") +
  geom_boxplot(aes(fill = Method)) +
  # stat_summary(aes(fill = level), fun = mean,
  #              geom = "col", width = .7, position = position_dodge(width = .9)) +
  # stat_summary(aes(group = level), fun.data = "mean_se",
  #              geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  facet_grid(cols = vars(level)) +
  scale_fill_brewer(palette = "Oranges") +
  scale_x_discrete(labels = dictLabelsML) +
  scale_y_continuous(n.breaks = 10, name = axis_y)

# p_cross_age <- plot_grid(p_title, p_cross_age, nrow = 2, rel_heights = c(.1, 1))


## R square of age: CGZ ---------------------------
model_dt <- fread("./CaseStudy2_CG/results/tables/no_missing_values/6_r2_age.csv")

sub_r2 <- model_dt %>%
  filter(Training.Validation == 2) %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "Ratio_Median centering_RUV-III-C_ComBat", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  arrange(correction, level) %>%
  mutate(name = paste(level, correction, quantitation, sep = "_"))

lg_level <- ggplot(sub_r2, aes(x = correction, y = r2)) +
  geom_col(aes(fill = level)) +
  scale_fill_manual(values = c("Negative control" = "#999999", dictColorsLevel))

p_title <- ggplot(sub_r2) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 3.4), "cm")) +
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

axis_y <- "R Square"

p_r2_age1 <- ggplot(sub_r2, aes(x = level, y = r2)) +
  stat_summary(aes(fill = level), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = level), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Uncorrected", "Negative control")),
              map_signif_level = TRUE,
              step_increase = .2,
              tip_length = .05,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c(dictColorsLevel, "Negative control" = "#999999")) +
  scale_y_continuous(limits = c(-.1, .5), n.breaks = 10, name = axis_y) +
  facet_grid(cols = vars(correction))

# p_r2_final <- plot_grid(p_title, p_r2_age, nrow = 2, rel_heights = c(.1, 1))

sub_r3 <- model_dt %>%
  filter(Training.Validation == 2, !level %in% "Negative control") %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))
  
p_r2_age2 <- ggplot(sub_r3, aes(x = correction, y = r2)) +
  stat_summary(aes(fill = correction), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = correction), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  geom_signif(comparisons = list(c("Log transformed", "Ratio"),
                                 c("Log transformed", "Median centering"),
                                 c("Log transformed", "RUV-III-C"),
                                 c("Log transformed", "ComBat")),
              map_signif_level = TRUE,
              step_increase = .2,
              tip_length = .05,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c(dictColorsMethod)) +
  scale_y_continuous(limits = c(-.1, .5), n.breaks = 10, name = axis_y)

p_r2_age3 <- ggplot(sub_r3, aes(x = quantitation, y = r2)) +
  stat_summary(aes(fill = quantitation), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = quantitation), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  geom_signif(comparisons = list(c("MaxLFQ", "TopPep3"),
                                 c("MaxLFQ", "iBAQ"),
                                 c("TopPep3", "iBAQ")),
              map_signif_level = TRUE,
              step_increase = .2,
              tip_length = .05,
              textsize = 6,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c(dictColorsPep2pro)) +
  scale_y_continuous(limits = c(-.1, .5), n.breaks = 10, name = axis_y)


## 5-fold cross validation (R square) of HBA1C: CGZ ---------------------------
cross_dt <- read_excel("./CaseStudy2_CG/results/tables/no_missing_values/5_models_hba1c_training_5fold.xlsx")

sub_cross <- cross_dt %>%
  # filter(Method %in% c("Decision Tree", "Bootstrap Forest", "Boosted Tree", "K Nearest Neighbors")) %>%
  tidyr::separate(`Original Data`, c("level", "correction", "quantitation"), sep = ", ") %>%
  mutate_at(8:10, ~ str_extract(., "(?<==).+")) %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_title <- ggplot(sub_cross) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 3), "cm"))+
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

axis_y <- "R Square"

p_cross <- ggplot(sub_cross, aes(x = reorder(Method, -`RSquare`), y = `RSquare`)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Uncorrected", "Negative control")),
              map_signif_level = TRUE,
              # y_position = 1.2,
              step_increase = .15,
              tip_length = .05,
              textsize = 5,
              test = "t.test") +
  geom_boxplot(aes(fill = Method)) +
  # stat_summary(aes(fill = level), fun = mean,
  #              geom = "col", width = .7, position = position_dodge(width = .9)) +
  # stat_summary(aes(group = level), fun.data = "mean_se",
  #              geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  facet_grid(cols = vars(level)) +
  scale_fill_brewer(palette = "PuRd")+
  scale_y_continuous(limits = c(-.1, .5), n.breaks = 10, name = axis_y);p_cross

p_cross_final <- plot_grid(p_title, p_cross, nrow = 2, rel_heights = c(.1, 1))


## R square of HBA1C: CGZ ---------------------------
model_dt <- fread("./CaseStudy2_CG/results/tables/no_missing_values/6_r2_hba1c.csv")

sub_r2 <- model_dt %>%
  filter(Training.Validation == 2) %>%
  mutate(scenario = "Random design") %>%
  mutate_at("level", ~ factor(., levels = c("Negative control", unique(data_levels)))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  arrange(correction, level) %>%
  mutate(name = paste(level, correction, quantitation, sep = "_"))

lg_level <- ggplot(sub_r2, aes(x = correction, y = r2)) +
  geom_col(aes(fill = level)) +
  scale_fill_manual(values = c("Negative control" = "#999999", dictColorsLevel))

p_title <- ggplot(sub_r2) +
  theme_classic() +
  theme(line = element_blank(),
        strip.text = element_text(size = 30),
        plot.margin = unit(c(.5, .5, 0, 3), "cm"))+
  facet_grid(cols = vars(scenario),
             labeller = labeller(scenario = dictLabelsScenario))

dictLabelsName_tmp <- rep(c(1, 3, 5, 6, 8, 9, 11, 12, 14, 15), each = 3)
names(dictLabelsName_tmp) <- unique(sub_r2$name)

sub_tmp <- sub_r2 %>%
  mutate_at("name", ~ dictLabelsName_tmp[.]) %>%
  mutate_at("name", as.factor)

axis_y <- "R Square"

p_r2 <- ggplot(sub_tmp, aes(x = level, y = r2)) +
  stat_summary(aes(fill = level), fun = mean,
               geom = "col", width = .7, position = position_dodge(width = .9)) +
  stat_summary(aes(group = level), fun.data = "mean_se",
               geom = "errorbar", width = .4, position = position_dodge(width = .9)) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected"),
                                 c("Uncorrected", "Negative control")),
              map_signif_level = TRUE,
              y_position = .8,
              step_increase = 1,
              tip_length = .2,
              textsize = 6,
              test = "wilcox.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 30, angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 24),
        plot.margin = unit(c(.2, .5, .5, .5), "cm")) +
  scale_fill_manual(values = c(dictColorsLevel, "Negative control" = "#999999")) +
  scale_y_continuous(limits = c(-.1, .2), n.breaks = 10, name = axis_y);p_r2

p_r2_final <- plot_grid(p_title, p_r2, nrow = 2, rel_heights = c(.1, 1))


## All: Radar/Corplot/PVCA --------------------------
sum_count <- count_dt1 %>%
  rbind(., sub_count2) %>%
  group_by(scenario, level, correction, quantitation) %>%
  summarise(Qualified = mean(Qualified), Missed = mean(Missed))

sum_f1 <- withref_dt %>%
  filter(all %in% "All") %>%
  group_by(scenario, level, correction, quantitation) %>%
  summarise(`F1 Score` = mean(f1, na.rm = TRUE))

sum_rc <- sub_ref %>%
  reshape2::melt(., value.name = "RC") %>%
  tidyr::separate(Var2, c("scenario", "level", "correction"), sep = "_") %>%
  rename(quantitation = Var1) %>%
  select(scenario, level, correction, quantitation, RC)

sum_snr <- sub_pca_snr1 %>%
  rbind(., sub_pca_snr2) %>%
  filter(all %in% "All") %>%
  rename(SNR = snr) %>%
  select(scenario, level, correction, quantitation, SNR)

sum_brinfo <- brinfo_dt %>%
  filter(all %in% "All") %>%
  rename(BRinfo = brinfo) %>%
  mutate(scenario = "Random design") %>%
  select(scenario, level, correction, quantitation, BRinfo)

sum_mcc <- mcc_dt %>%
  filter(Training.Validation == 2) %>%
  filter(!level %in% "Negative control") %>%
  rename(MCC = mcc) %>%
  mutate(scenario = "Random design") %>%
  select(scenario, level, correction, quantitation, MCC)

sum_r2 <- model_dt %>%
  filter(Training.Validation == 2) %>%
  filter(!level %in% "Negative control") %>%
  rename(`R Square` = r2) %>%
  mutate(scenario = "Random design") %>%
  select(scenario, level, correction, quantitation, `R Square`)

sum_all_raw <- sum_count %>%
  full_join(., sum_snr, by = c("scenario", "level", "correction", "quantitation")) %>%
  full_join(., sum_brinfo, by = c("scenario", "level", "correction", "quantitation")) %>%
  full_join(., sum_f1, by = c("scenario", "level", "correction", "quantitation")) %>%
  full_join(., sum_rc, by = c("scenario", "level", "correction", "quantitation")) %>%
  full_join(., sum_mcc, by = c("scenario", "level", "correction", "quantitation")) %>%
  full_join(., sum_r2, by = c("scenario", "level", "correction", "quantitation")) %>%
  ungroup()

sum_quantiles <- sum_all_raw %>%
  summarise_if(is.numeric, quantile, na.rm = TRUE)

qc_range_norm <- function(x, decreasing = FALSE) {
  
  quantiles <- quantile(x)
  
  x_new <- x
  
  x_new[x <= quantiles[2] & x >= quantiles[1]] <- 2.5
  
  x_new[x <= quantiles[3] & x > quantiles[2]] <- 5
  
  x_new[x <= quantiles[4] & x > quantiles[3]] <- 7.5
  
  x_new[x <= quantiles[5] & x > quantiles[4]] <- 10
  
  if (decreasing) x_new <- 10 - x_new
  
  return(x_new)
}

sum_all_scaled <- sum_all_raw %>%
  group_by(scenario) %>%
  # mutate_at("Missed", ~ qc_range_norm(., decreasing = TRUE)) %>%
  # mutate_at(c("Qualified", "F1 Score", "SNR", "RC"), ~ qc_range_norm(., decreasing = FALSE)) %>%
  mutate_at(grep("Missed|BRinfo", colnames(sum_all_raw)),
            ~ protqc::qc_linear_norm(.,
                                     min(., na.rm = TRUE),
                                     max(., na.rm = TRUE),
                                     decreasing = TRUE)) %>%
  mutate_at(grep("Qualified|SNR|F1|RC|MCC|Square", colnames(sum_all_raw)),
            ~ protqc::qc_linear_norm(.,
                                     min(., na.rm = TRUE),
                                     max(., na.rm = TRUE),
                                     decreasing = FALSE)) %>%
  ungroup() %>%
  rename_if(is.numeric, ~ paste(., "(scaled)")) %>%
  mutate(`Total score` = apply(., 1, function(x) mean(as.numeric(x[5:12]), na.rm = TRUE))) %>%
  mutate_if(is.numeric, ~ round(., digits = 3))

sum_all_final <- left_join(sum_all_raw, sum_all_scaled, by = c("scenario", "level", "correction", "quantitation"))

fwrite(sum_all_final, "./extended_data_table.csv")

## Radar plot
sum_all_by_beca <- sum_all_scaled %>%
  filter(!level %in% "Uncorrected") %>%
  rename_if(is.numeric, ~ gsub(" (scaled)", "", ., fixed = TRUE)) %>%
  reshape2::dcast(., scenario + level ~ correction + quantitation, value.var = "Total score")

dictColumnsName <- c("CB_I", "CB_M", "CB_T",
                     "MC_I", "MC_M", "MC_T",
                     "RA_I", "RA_M", "RA_T",
                     "RU_I", "RU_M", "RU_T")

names(dictColumnsName) <- colnames(sum_all_by_beca)[3:14]

sum_all_by_beca <- sum_all_by_beca %>%
  plyr::rename(dictColumnsName)

p_radar_list <- pblapply(unique(scenarios), function(scenario_id) {
  
  sum_all_by_beca_i <- sum_all_by_beca %>%
    filter(scenario %in% scenario_id) %>%
    select(!scenario)
  
  p_radar <- ggradar(sum_all_by_beca_i, values.radar = c("0%", "50%", "100%"),
                     grid.min = 0, grid.mid = 5, grid.max = 10,
                     plot.title = dictLabelsScenario[scenario_id],
                     # plot.legend = ifelse(label_id %in% "Confounded design_MaxLFQ", TRUE, FALSE),
                     plot.legend = FALSE,
                     legend.text.size = 8,
                     group.colours = dictColorsLevel[2:3],
                     fill = TRUE,
                     fill.alpha = .3)
  
  return(p_radar)

})

sum_all_by_level <- sum_all_scaled %>%
  filter(!level %in% "Uncorrected") %>%
  rename_if(is.numeric, ~ gsub(" (scaled)", "", ., fixed = TRUE)) %>%
  reshape2::dcast(., scenario + correction + quantitation ~ level, value.var = "Total score") %>%
  mutate(delta = `Protein-corrected` - `Peptide-corrected`) %>%
  mutate(class = ifelse(delta > 0 , "Yes", "No"))

p_all_percent <- ggplot(sum_all_by_level) +
  geom_bar(aes(y = class, fill = class), position = "stack") +
  scale_x_continuous(labels = ~ sprintf("%.2f%%", ./36),
                     n.breaks = 10, name = "Proportion") +
  scale_fill_manual(values = c("Yes" = "#377EB8", "No" = "#984EA3"),
                    name = "Better performances at protein level") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 20))

## Heatmap
sum_all2 <- sum_all_scaled %>%
  rename_if(is.numeric, ~ gsub(" (scaled)", "", ., fixed = TRUE)) %>%
  reshape2::melt(., id = 1:4, variable.name = "metric") %>%
  reshape2::dcast(., level + correction + quantitation ~ scenario + metric, value.var = "value") %>%
  mutate("Total score (mean)" = apply(., 1, function(x) mean(as.numeric(x[c(12, 21, 30)])))) %>%
  arrange(desc(across(ends_with("Total score (mean)")))) %>%
  mutate(class = apply(., 1, function(x) {
    if (as.numeric(x[31]) <= 5.380833) {
      a <- "Bad"
    } else if (as.numeric(x[31]) <= 6.230000 & as.numeric(x[31]) > 5.358500) {
      a <- "Fair"
    } else if (as.numeric(x[31]) <= 7.123000 & as.numeric(x[31]) > 6.230000) {
      a <- "Good"
    } else if (as.numeric(x[31]) <= 8.220667 & as.numeric(x[31]) > 7.123000) {
      a <- "Great"
    }
    return(a)
  })) %>%
  mutate_at("class", ~ factor(., levels = c("Great", "Good", "Fair", "Bad"))) %>%
  mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods))) %>%
  select(!`Total score (mean)`) %>%
  select(class, level, correction, quantitation, everything())

sum_all_matrix <- sum_all2 %>%
  mutate(label = paste(class, level, correction, quantitation, sep = "_")) %>%
  tibble::column_to_rownames("label") %>%
  select_if(is.numeric) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  as.matrix

sum_all_annot1 <- sum_all2 %>%
  mutate(label = paste(class, level, correction, quantitation, sep = "_")) %>%
  tibble::column_to_rownames("label") %>%
  select_if(is.factor) %>%
  rename_with(~ Hmisc::capitalize(.))

ha_row <- HeatmapAnnotation(df = sum_all_annot1,
                            col = list(Class = dictColorsClass,
                                       Level = dictColorsLevel,
                                       Quantitation = dictColorsPep2pro,
                                       Correction = dictColorsMethod),
                            which = "row",
                            show_legend = TRUE,
                            gap = unit(5, "points"),
                            show_annotation_name = TRUE,
                            annotation_name_side = "top",
                            gp = gpar(col = "white", lwd = 2),
                            annotation_name_gp = gpar(fontsize = 16, col = "darkgrey"),
                            annotation_legend_param = list(ncol = 1,
                                                           grid_height = unit(8, "mm"),
                                                           grid_width = unit(6, "mm"),
                                                           title_gp = gpar(fontsize = 16),
                                                           labels_gp = gpar(fontsize = 14)),
                            simple_anno_size = unit(.8, "cm"))

sum_all_annot2 <- colnames(sum_all2)[5:31] %>%
  sapply(., function(x) {unlist(strsplit(x, split = "_"))}) %>%
  t %>%
  as.data.frame

colnames(sum_all_annot2) <- c("Scenario", "Metric")

dictMetricsName <- sum_all_annot2$Metric
names(dictMetricsName) <- rownames(sum_all_annot2)

dictScenariosName <- dictLabelsScenario[sum_all_annot2$Scenario]
names(dictScenariosName) <- rownames(sum_all_annot2)

sum_all_annot2 <- sum_all_annot2 %>%
  select(Scenario)

ha_col <- HeatmapAnnotation(df = sum_all_annot2,
                            col = list(Scenario = dictColorsScenario),
                            which = "col",
                            show_legend = FALSE,
                            gap = unit(3, "points"),
                            show_annotation_name = FALSE,
                            annotation_name_side = "right",
                            gp = gpar(col = "transparent"),
                            annotation_legend_param = list(ncol = 1,
                                                           grid_height = unit(8, "mm"),
                                                           grid_width = unit(6, "mm"),
                                                           title_gp = gpar(fontsize = 16),
                                                           labels_gp = gpar(fontsize = 14)),
                            simple_anno_size = unit(.5, "cm"))

set_color <- colorRamp2(c(1, 5, 10), c("navyblue", "white", "#FB6A4A"))

ha <- Heatmap(matrix = sum_all_matrix,
              name = "Score",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.circle(x = x, y = y, 
                            r = sum_all_matrix[i, j]/7 * min(unit.c(width, height)), 
                            gp = gpar(fill = set_color(sum_all_matrix[i, j]), col = "transparent"))
              },
              show_heatmap_legend = TRUE,
              column_labels = dictMetricsName,
              column_split = rep(c("Multi-lab (Balanced)", "Multi-lab (Confounded)", "Multi-batch"), each = 9),
              column_title_gp = gpar(fontsize = 20),
              cluster_columns = FALSE,
              col = set_color,
              rect_gp = gpar(type = "none"),
              heatmap_legend_param = list(ncol = 1,
                                          grid_height = unit(8, "mm"),
                                          grid_width = unit(6, "mm"),
                                          title_gp = gpar(fontsize = 16),
                                          labels_gp = gpar(fontsize = 14)),
              show_column_names = TRUE,
              column_names_gp = gpar(fontsize = 16),
              column_names_rot = 45,
              column_names_side = "top",
              show_row_names = FALSE,
              cluster_rows = FALSE,
              left_annotation = ha_row,
              top_annotation = ha_col)

p_ha_all <- grid.grabExpr(draw(ha, merge_legend = TRUE))


## compare total scores between groups
sum_all_scaled <- sum_all_scaled %>%
  mutate_at("level", ~ factor(., levels = unique(data_levels))) %>%
  mutate_at("correction", ~ factor(., levels = unique(correct_methods))) %>%
  mutate_at("quantitation", ~ factor(., levels = unique(pep2pro_methods)))

p_all_bar1 <- ggplot(sum_all_scaled, aes(x = level, y = `Total score`)) +
  stat_summary(aes(fill = level), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Peptide-corrected", "Protein-corrected"),
                                 c("Peptide-corrected", "Uncorrected"),
                                 c("Uncorrected", "Protein-corrected")),
              map_signif_level = TRUE,
              y_position = 8,
              step_increase = .1,
              tip_length = .02,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsLevel) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 1), name = "Total score") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_all_bar2 <- ggplot(sum_all_scaled, aes(x = correction, y = `Total score`)) +
  stat_summary(aes(fill = correction), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("Log transformed", "Ratio"),
                                 c("Log transformed", "Median centering"),
                                 c("Log transformed", "RUV-III-C"),
                                 c("Log transformed", "ComBat")),
              map_signif_level = TRUE,
              y_position = 8,
              step_increase = .1,
              tip_length = .02,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsMethod) +
  scale_y_continuous(limits = c(0, 11), breaks = seq(0, 10, 1), name = "Total score") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))

p_all_bar3 <- ggplot(sum_all_scaled, aes(x = quantitation, y = `Total score`)) +
  stat_summary(aes(fill = quantitation), fun = mean, geom = "bar", width = .7) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = .2) +
  geom_signif(comparisons = list(c("MaxLFQ", "iBAQ"),
                                 c("MaxLFQ", "TopPep3"),
                                 c("iBAQ", "TopPep3")),
              map_signif_level = TRUE,
              y_position = 8,
              step_increase = .1,
              tip_length = .01,
              textsize = 4,
              test = "t.test",
              test.args = list(pairwise = TRUE)) +
  theme_classic() +
  theme(legend.position = "right",
        legend.title = element_text(size =20),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 24, margin = unit(rep(.3, 4), "cm")),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_manual(values = dictColorsPep2pro) +
  scale_y_continuous(limits = c(0, 10), n.breaks = 10, name = "Total score") +
  facet_grid(cols = vars(scenario), labeller = as_labeller(dictLabelsScenario))


## Figure 2 ----------------------
p_blank <- ggplot() + theme_void()

figure2ab <- plot_grid(p_intensity_pep1 + theme(plot.margin = unit(c(1,.5, 1.7, .5), "cm")),
                       p_pca_pep1 + theme(plot.margin = unit(c(1, 2.5, .5, .5), "cm")),
                       labels = c("a", "b"), label_size = 30,
                       ncol = 2, rel_widths = c(1, .8))

figure2c <- plot_grid(plotlist = c(p_pca_pro_list1[1:3], list(p_blank)),
                      labels = c("c"), label_size = 30,
                      ncol = 4, rel_widths = c(1, 1, 1, .4))

figure2de <- plot_grid(p_intensity_pep21 + theme(plot.margin = unit(c(1,.5, .3, .5), "cm")),
                          p_pca_pep21 + theme(plot.margin = unit(c(1, 3, .5, .5), "cm")),
                          labels = c("d", "e"), label_size = 30,
                          ncol = 2, rel_widths = c(1, .8))

figure2f <- plot_grid(plotlist = c(p_pca_all_pro_list2, list(p_blank)),
                      labels = c("f"), label_size = 30,
                      ncol = 4, rel_widths = c(1, 1, 1, .4))

figure2 <- plot_grid(figure2ab, figure2c, figure2de, figure2f,
                     nrow = 4)

ggsave("./figure2.pdf", figure2, height = 18, width = 16)

sp_figure1bc <- plot_grid(p_intensity_pep22 + theme(plot.margin = unit(c(1,.5, .3, .5), "cm")),
                       p_pca_pep22 + theme(plot.margin = unit(c(1, .3, .5, .5), "cm")),
                       labels = c("b", "c"), label_size = 30,
                       ncol = 2, rel_widths = c(1, .8))

sp_figure1d <- plot_grid(plotlist = c(p_pca_pro_list2, list(p_blank)),
                         labels = c("d"), label_size = 30,
                         ncol = 4, rel_widths = c(1, 1, 1, .4))

sp_figure1 <- plot_grid(p_count_pep1 + theme(plot.margin = unit(c(1, .5, .8, .5), "cm")),
                        sp_figure1bc,
                        sp_figure1d,
                        labels = c("a", ""), label_size = 30,
                        ncol = 1)

ggsave("./extended_figure1.pdf", sp_figure1, height = 14, width = 16)


## Figure 3 -------------------------
figure3ab <- plot_grid(p_count_bar_levelgroups1 +
                         theme(legend.position = "none",
                               plot.margin = unit(c(.5, .5, .5, 1.5), "cm"),
                               axis.ticks.x = element_blank()) +
                         scale_x_discrete(labels = c("Peptide-corrected" = "Peptide-\ncorrected",
                                                     "Protein-corrected" = "Protein-\ncorrected")),
                       p_count_bar_levelgroups2 +
                         theme(legend.position = "none",
                               plot.margin = unit(c(.5, .5, .5, 1.5), "cm"),
                               axis.ticks.x = element_blank()) +
                         scale_x_discrete(labels = c("Peptide-corrected" = "Peptide-\ncorrected",
                                                     "Protein-corrected" = "Protein-\ncorrected")),
                       ncol = 2, labels = c("a", "b"), label_size = 24)

# p_ha_legend <- plot_grid(get_legend(lg_level + guides(fill = guide_legend(direction = "horizontal"))),
#                          get_legend(lg_method + guides(fill = guide_legend(direction = "horizontal"))),
#                          nrow = 2)

figure3cd <- plot_grid(p_f1_bar1 + theme(legend.position = "none",
                                         plot.margin = unit(c(.5, .5, .5, 1.5), "cm")) +
                         scale_x_discrete(labels = c("Peptide-corrected" = "Peptide-\ncorrected",
                                                     "Protein-corrected" = "Protein-\ncorrected")),
                       p_rc_bar1 + theme(legend.position = "none",
                                         plot.margin = unit(c(.5, .5, .5, 1.5), "cm")) +
                         scale_x_discrete(labels = c("Peptide-corrected" = "Peptide-\ncorrected",
                                                     "Protein-corrected" = "Protein-\ncorrected")),
                       ncol = 2, labels = c("c", "d"), label_size = 24)

figure3ef <- plot_grid(p_f1_box11 + theme(legend.position = "none",
                                         plot.margin = unit(c(.5, .5, .5, 1.5), "cm")),
                       p_rc_box11 + theme(legend.position = "none",
                                         plot.margin = unit(c(.5, .5, .5, 1.5), "cm")), 
                       nrow = 2, rel_heights = c(1, 1),
                       labels = c("e", "f"), label_size = 24)

figure3 <- plot_grid(figure3ab,
                     figure3cd,
                     figure3ef,
                     # p_ha1,
                     # p_ha_legend + theme(plot.margin = unit(c(0, 0, 1, 0), "cm")) ,
                     nrow = 3, rel_heights = c(1, 1, 2))

ggsave("./figure3.pdf", figure3, height = 26, width = 24)

sp_figure2 <- plot_grid(p_count_bar11 + theme(legend.position = "none",
                                              plot.margin = unit(c(.5, .5, 1.8, .5), "cm")),
                        p_count_bar12 + theme(legend.position = "none",
                                              plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                        p_count_bar21 + theme(legend.position = "none",
                                              plot.margin = unit(c(.5, .5, 1.8, .5), "cm")),
                        p_count_bar22 + theme(legend.position = "none",
                                              plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                        ncol = 2, labels = c("a", "b", "c", "d"), label_size = 30)

ggsave("./extended_figure2.pdf", sp_figure2, height = 14, width = 24)

sp_figure3 <- plot_grid(p_cv_list[[1]],
                        p_cv_list[[2]],
                        p_cv_list[[3]],
                        get_legend(lg_level + guides(fill = guide_legend(direction = "horizontal"))),
                        nrow = 4, rel_heights = c(1, 1, 1, .1),
                        labels = c("a", "b", "c"), label_size = 24)

ggsave("./extended_figure3.pdf", sp_figure3, height = 24, width = 16)

sp_figure7 <- p_dep_bar0
ggsave("./extended_figure7.pdf", sp_figure7, height = 13, width = 15)

sp_figure8 <- p_dep_bar1
ggsave("./extended_figure8.pdf", sp_figure8, height = 13, width = 15)

sp_figure9 <- p_dep_bar2
ggsave("./extended_figure9.pdf", sp_figure9, height = 13, width = 15)

sp_figure5 <- plot_grid(p_f1_list[[1]] + theme(plot.margin = unit(c(.5, 1, .5, 1), "cm")),
                        p_f1_list[[2]] + theme(plot.margin = unit(c(.5, 1, .5, 1), "cm")),
                        p_f1_list[[3]] + theme(plot.margin = unit(c(.5, 1, .5, 1), "cm")),
                        get_legend(lg_level +
                                     guides(fill = guide_legend(direction = "horizontal")) +
                                     theme(legend.title = element_text(size = 30),
                                           legend.text = element_text(size = 24))),
                        nrow = 4, rel_heights = c(1, 1, 1, .1),
                        labels = c("a", "b", "c", ""), label_size = 24)

ggsave("./extended_figure5.pdf", sp_figure5, height = 26, width = 23)

sp_figure6 <- plot_grid(p_f1_bar2 + theme(plot.margin = unit(c(.5, 0, 0, .5), "cm")),
                        p_f1_bar3 + theme(plot.margin = unit(c(.5, 0, 0, .5), "cm")),
                        nrow = 2, rel_heights = c(1, .85),
                        labels = c("a", "b"), label_size = 24)

ggsave("./extended_figure6.pdf", sp_figure6, height = 10, width = 12)

sp_figure10 <- plot_grid(p_rc_bar2 + theme(plot.margin = unit(c(.5, 0, 0, .5), "cm")),
                        p_rc_bar3 + theme(plot.margin = unit(c(.5, 0, 0, .5), "cm")),
                        nrow = 2, rel_heights = c(1, .85),
                        labels = c("a", "b"), label_size = 24)

ggsave("./extended_figure10.pdf", sp_figure10, height = 10, width = 12)


## Figure 4 -------------------------
p_legend <- plot_grid(get_legend(lg_pca_batch + theme(legend.direction = "horizontal",
                                                      legend.title = element_text(size = 30),
                                                      legend.text = element_text(size = 24))),
                      get_legend(lg_pca_sample + theme(legend.direction = "horizontal",
                                                       legend.title = element_text(size = 30),
                                                       legend.text = element_text(size = 24))),
                      ncol = 2)

figure4a <- plot_grid(p_pca_list1[[1]] + theme(plot.margin = unit(c(.5, 0, .5, 1), "cm")),
                      p_legend,
                      nrow = 2, rel_heights = c(1, .1),
                      labels = c("a", ""), label_size = 30)

figure4bc <- plot_grid(p_snr_bar1 +
                         scale_x_discrete(labels = c("Peptide-corrected" = "Peptide-\ncorrected",
                                                     "Protein-corrected" = "Protein-\ncorrected")),
                       p_brinfo_bar1 +
                         scale_x_discrete(labels = c("Peptide-corrected" = "Peptide-\ncorrected",
                                                     "Protein-corrected" = "Protein-\ncorrected")),
                       ncol = 2, rel_widths = c(3, 1),
                       labels = c("b", "c"), label_size = 30)

figure4de <- plot_grid(p_snr_box11 + theme(legend.position = "none"),
                       p_brinfo_box11 + theme(legend.position = "none"),
                       nrow = 2, labels = c("d", "e"), label_size = 30)

figure4 <- plot_grid(figure4a + theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                     figure4bc + theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                     figure4de + theme(plot.margin = unit(c(.5, .5, 0, .5), "cm")),
                     get_legend(lg_level + labs(fill = "Level") + 
                                  theme(legend.direction = "horizontal",
                                        legend.title = element_text(size = 30),
                                        legend.text = element_text(size = 24))),
                     nrow = 4, rel_heights = c(2, 1.1, 2, .1))

ggsave("./figure4.pdf", figure4, height = 28, width = 24)

p_legend2 <- plot_grid(get_legend(lg_pca_batch22 + theme(legend.direction = "horizontal",
                                                      legend.title = element_text(size = 30),
                                                      legend.text = element_text(size = 24))),
                      get_legend(lg_pca_sample22 + theme(legend.direction = "horizontal",
                                                       legend.title = element_text(size = 30),
                                                       legend.text = element_text(size = 24))),
                      ncol = 2)

sp_figure11 <- plot_grid(p_pca_list1[[2]] + theme(plot.margin = unit(c(.5, 0, .5, 1), "cm")),
                        p_legend,
                        p_pca21_final + theme(plot.margin = unit(c(.5, 0, .5, 1), "cm")),
                        get_legend(lg_pca_batch21 + theme(legend.direction = "horizontal",
                                                          legend.title = element_text(size = 30),
                                                          legend.text = element_text(size = 24))),
                        p_pca22_final + theme(plot.margin = unit(c(.5, 0, .5, 1), "cm")),
                        p_legend2,
                        nrow = 6, rel_heights = c(1, .1, 1, .1, 1, .1),
                        labels = c("a", "", "b", "", "c", ""), label_size = 30)

ggsave("./extended_figure11.pdf", sp_figure11, height = 30, width = 24)

sp_figure12cd <- plot_grid(p_brinfo_bar2 + theme(plot.margin = unit(c(.5, 0, 0, 1.2), "cm")),
                          p_brinfo_bar3 + theme(plot.margin = unit(c(.5, 0, 2, 1.2), "cm")),
                          ncol = 2, labels = c("c", "d"), label_size = 24)

sp_figure12 <- plot_grid(p_snr_bar2 + theme(legend.position = "none",
                                            plot.margin = unit(c(.5, 0, 0, 1.2), "cm")),
                        p_snr_bar3 + theme(legend.position = "none",
                                           plot.margin = unit(c(.5, 0, 0, 1.2), "cm")),
                        sp_figure12cd + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
                        nrow = 3, rel_heights = c(1.2, 1, 1.2),
                        labels = c("a", "b", ""), label_size = 24)

ggsave("./extended_figure12.pdf", sp_figure12, height = 15, width = 12)


## Figure 5 -------------------------
figure5 <- plot_grid(p_cross_sex + theme(plot.margin = unit(c(.5, 0, 0, 1), "cm")),
                     p_cross_age + theme(plot.margin = unit(c(.5, 0, 0, 1), "cm")),
                     p_mcc1 + theme(plot.margin = unit(c(.5, 0, .5, 1.5), "cm"),
                                   axis.text.x = element_blank()),
                     p_r2_age1 + theme(plot.margin = unit(c(.5, 0, .8, 1), "cm"),
                                  axis.text.x = element_blank()),
                     get_legend(lg_level_mcc + theme(legend.direction = "horizontal",
                                                     legend.title = element_text(size = 30),
                                                     legend.text = element_text(size = 24))),
                     nrow = 5, rel_heights = c(1, 1, .8, .8, .1),
                     labels = c("a", "b", "c", "d", ""), label_size = 30)

ggsave("./figure5.pdf", figure5, height = 25, width = 20)

sp_figure13 <- plot_grid(p_mcc2 + theme(plot.margin = unit(c(.5, 0, 0, 0), "cm")),
                         p_mcc3 + theme(plot.margin = unit(c(.5, 0, 2.5, .5), "cm")),
                         p_r2_age2 + theme(plot.margin = unit(c(.5, 0, 0, 0), "cm")),
                         p_r2_age3 + theme(plot.margin = unit(c(.5, 0, 2.5, .5), "cm")),
                         ncol = 2, rel_widths = c(1.5, 1),
                         labels = c("a", "b", "c", "d"), label_size = 30)
ggsave("./extended_figure13.pdf", sp_figure13, height = 16, width = 20)


## Figure 6 -------------------------
p_bg <- ggplot(sub_auc) +  theme_void()

figure6a <- plot_grid(p_bg, p_bg, p_radar_list[[1]], p_bg, p_radar_list[[2]], p_bg, p_radar_list[[3]], p_bg,
                      ncol = 8, rel_widths = c(.1, .1, 1, .1, 1, .1, 1, .5))

figure6 <- plot_grid(p_bg, figure6a,
                     p_all_percent + theme(plot.margin = unit(c(.5, 0, .5, .5), "cm")),
                     p_bg, p_ha_all,
                     nrow = 5, rel_heights = c(.1, 1, 1, .1, 3.5),
                     labels = c("a", "", "b", "c", ""), label_size = 24)

ggsave("./figure6.pdf", figure6, height = 22, width = 15)

sp_figure14 <- plot_grid(p_all_bar1 + theme(legend.position = "none",
                                            plot.margin = unit(c(.5, 0, 0, 1), "cm")),
                         p_all_bar2 + theme(legend.position = "none",
                                            plot.margin = unit(c(.5, 0, 0, 1), "cm")),
                         p_all_bar3 + theme(legend.position = "none",
                                            plot.margin = unit(c(.5, 0, 0, 1), "cm")),
                         nrow = 3, rel_heights = c(1, 1, .85),
                         labels = c("a", "b", "c"), label_size = 24)

ggsave("./extended_figure14.pdf", sp_figure14, height = 16, width = 12)


## Abstract -----------
all_results <- read_excel("./submitted_20241108/Extended_Data_Table_Evaluation_Summary.xlsx")

sub_snr <- all_results %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Ratio_Median centering_RUV-III-C", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  reshape2::dcast(., scenario + correction + quantitation ~ level, value.var = "SNR") %>%
  mutate(delta1 = `Protein-corrected` - `Peptide-corrected`,
         delta2 = `Protein-corrected` - `Uncorrected`) %>%
  mutate(pass = ifelse(delta1 > 0, "Yes", "No"))

table(sub_snr$pass)

sub_total <- all_results %>%
  mutate_at("correction", ~ ifelse(. %in% "Log transformed", "ComBat_Ratio_Median centering_RUV-III-C", .)) %>%
  tidyr::separate_rows(., correction, sep = "_") %>%
  reshape2::dcast(., scenario + correction + quantitation ~ level, value.var = "Total score") %>%
  mutate(delta1 = `Protein-corrected` - `Peptide-corrected`,
         delta2 = `Protein-corrected` - `Uncorrected`) %>%
  mutate(pass = ifelse(delta1 > 0, "Yes", "No"))

table(sub_total$pass)

df_mapping_final <- fread("./expfiles/peptide/pep2prot.csv")

stat_mapping <- df_mapping_final %>%
  filter(!protein_names %in% "") %>%
  group_by(batch, protein_group_ids) %>%
  summarise(n_peps = length(unique(peptide_id)))

median(stat_mapping$n_peps)
mean(stat_mapping$n_peps)
quantile(stat_mapping$n_peps)

df_mapping_final <- fread("../../CaseStudy2_CG/data/mapped.csv")

stat_mapping <- df_mapping_final %>%
  filter(!Accession %in% "") %>%
  group_by(Accession) %>%
  summarise(n_peps = length(unique(Peptide)))

median(stat_mapping$n_peps)
mean(stat_mapping$n_peps)
quantile(stat_mapping$n_peps)

