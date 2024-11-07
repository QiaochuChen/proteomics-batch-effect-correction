library(rstatix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(limma)
library(edgeR)
library(DEqMS)

calculate_p <- function(exprdata, group, test_method, p_adjust, na_threshold, transform_to_log2) {

  na_cut <- na_threshold * ncol(exprdata)

  exprdf <- exprdata %>%
    as.data.frame %>%
    filter(apply(., 1, function(x) length(which(as.numeric(x) == 0)) <= na_cut))
  
  if (test_method == "limma") transform_to_log2 <- FALSE
  
  if (test_method == "deqms") transform_to_log2 <- TRUE
  
  if (transform_to_log2 == TRUE) {
    exprdata <- exprdata %>%
      mutate_if(is.numeric, ~ ifelse(. == 0, NA, log2(.)))
  }

  if (test_method == "wilcox") {
    
    names(group) <- colnames(exprdf)
    
    exprdf2long <- exprdf %>%
      tibble::rownames_to_column("feature") %>%
      reshape2::melt(., id = 1) %>%
      mutate(sample = group[variable])
    
    all_features <- exprdf2long %>%
      group_by(feature, sample) %>%
      summarise(sd = sd(value, na.rm = TRUE)) %>%
      filter(sd != 0) %>%
      pull(feature) %>%
      unique
    
    p_tables <- mclapply(all_features, function(feature_id) {
      
      exprdf2long_i <- exprdf2long %>%
        filter(feature %in% feature_id) %>%
        select(variable, sample, value)
      
      w_result <- exprdf2long_i %>%
        wilcox_test(value ~ sample, p.adjust.method = "none", detailed = TRUE) %>%
        mutate(feature = feature_id)
      
      return(w_result)
      
    }, mc.cores = 2)
    
    final_results <- p_tables %>%
      data.table::rbindlist() %>%
      adjust_pvalue(method = p_adjust) %>%
      select(feature, everything())
    
  } else if (test_method == "t") {

    names(group) <- colnames(exprdf)
    
    exprdf2long <- exprdf %>%
      tibble::rownames_to_column("feature") %>%
      reshape2::melt(., id = 1) %>%
      mutate(sample = group[variable])
    
    all_features <- exprdf2long %>%
      group_by(feature, sample) %>%
      summarise(sd = sd(value, na.rm = TRUE)) %>%
      filter(sd != 0) %>%
      pull(feature) %>%
      unique

    p_tables <- mclapply(all_features, function(feature_id) {

      exprdf2long_i <- exprdf2long %>%
        filter(feature %in% feature_id) %>%
        select(variable, sample, value)

      t_result <- exprdf2long_i %>%
        t_test(value ~ sample, p.adjust.method = "none", detailed = TRUE) %>%
        mutate(feature = feature_id)

      return(t_result)
    }, mc.cores = 2)

    final_results <- p_tables %>%
      data.table::rbindlist() %>%
      adjust_pvalue(method = p_adjust) %>%
      select(feature, everything())

  } else if (test_method == "limma") {

    dge <- DGEList(counts = exprdf)
    design <- model.matrix(~ group)

    keep <- filterByExpr(dge, design)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)

    v <- voom(dge, design, plot = F)
    fit <- lmFit(v, design)

    fit <- eBayes(fit)
    p_results <- topTable(fit, coef = ncol(design),
                          sort.by = "logFC", number = Inf,
                          adjust.method = p_adjust)

    final_results <- p_results %>%
      tibble::rownames_to_column("feature") %>%
      mutate(group1 = levels(group)[1]) %>%
      mutate(group2 = levels(group)[2])

  } else if (test_method == "deqms") {
    
    exprdf1 <- exprdf %>%
      mutate_if(is.numeric, ~ . + 1)
    
    exprmt_log2 <- exprdf1 %>%
      mutate_if(is.numeric, log2) %>%
      as.matrix
    
    exprdf_rowmins <- exprdf1 %>%
      mutate(count = rowMins(as.matrix(.))) %>%
      select(count)
    
    class <- factor(as.numeric(group))
    design <- model.matrix(~ 0 + class)
    
    fit1 <- lmFit(exprmt_log2, design = design)
    
    fit2_contrasts <- makeContrasts(class1 - class2, levels = design)
    fit2 <- contrasts.fit(fit1, contrasts = fit2_contrasts)
    
    fit3 <- eBayes(fit2)
    fit3$count <- exprdf_rowmins[rownames(fit3$coefficients), "count"]
    
    fit4 <- spectraCounteBayes(fit3)
    
    p_results <- outputResult(fit4, coef_col = 1)
    
    final_results <- p_results %>%
      dplyr::rename(feature = gene) %>%
      mutate(group1 = levels(group)[1]) %>%
      mutate(group2 = levels(group)[2])
    
  }

  return(final_results)

}

plot_deps <- function(final_results, labels, fc_threshold, p_threshold) {

  if (!is.null(labels)) {

    plot_results <- final_results %>%
      mutate(label = sapply(Feature, function(x) ifelse(x %in% labels, x, NA)))

  } else {

    plot_results <- final_results %>%
      mutate(dep_pass = ifelse(adj.P.Val < p_threshold &
                                 abs(logFC) > fc_threshold, "DEP", "Non-DEP")) %>%
      mutate(label = ifelse(dep_pass %in% "DEP", Feature, NA))

  }

  p <- ggplot(data = plot_results, aes(x = logFC, y = - 10 * log10(adj.P.Val))) +
    geom_point(aes(color = dep_pass), alpha = 0.8, size = 2) +
    geom_label_repel(aes(label = label), show.legend = FALSE) +
    geom_vline(xintercept = - fc_threshold, linetype="dashed",color = "red") +
    geom_vline(xintercept = fc_threshold, linetype="dashed",color = "red") +
    geom_hline(yintercept = -log10(p_threshold), linetype="dashed",color = "red") +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    scale_color_brewer(palette = "Set1") +
    labs(y = "-log10(p.adjust)", color = "");p

  return(p)

}

main_dep <- function(exprdata, group, na_threshold = 0, transform_to_log2 = TRUE,
                     test_method = "limma", p_adjust = "BH",
                     fc_threshold = 1, p_threshold = 0.05, labels = NULL,
                     plot_volcano = TRUE,
                     silent = FALSE) {

  if (!silent) print("Calculate p.values and FCs...")
  final_results <- calculate_p(exprdata, group,
                               test_method, p_adjust, na_threshold)

  if (plot_volcano) {
    if (!silent) print("A volcano plot...")
    p <- plot_deps(final_results, labels, fc_threshold, p_threshold)
    final_results <- list(data = final_results, plot = p)
  } else {
    final_results <- list(data = final_results)
  }

  return(final_results)

}
