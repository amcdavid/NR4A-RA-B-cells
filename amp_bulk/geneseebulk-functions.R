run_deseq_design = function (dge_, formula_, subset_) {
  dge_ = subset_dge(subset_, dge_)
  design(dge_) = formula_
  DESeq(dge_)
}

subset_dge = function (subset_, dge_) 
{
  if (!missing(subset_) && !is.null(subset_) && !is.na(subset_)) {
    cd = as.data.frame(colData(dge_))
    cd[["__idx__"]] = seq_len(nrow(cd))
    subset_ = parse(text = subset_, n = 1)[[1]]
    cd = dplyr::filter(cd, !!subset_)
    dge_ = dge_[, cd[["__idx__"]]]
    cdd = colData(dge_)
    for (j in seq_len(ncol(cdd))) {
      if (is.factor(cdd[, j])) 
        colData(dge_)[, j] = droplevels(cdd[, j])
    }
  }
  dge_
}

clamp = function(x, modulus = 5){
  x[x < -modulus] = -modulus
  x[x > modulus] = modulus
  x
}

modified_volc = function (dd, name, direction, heatmap_top_group = "group", heatmap_max_gene = Inf,
          subset_idx = TRUE, title = "", subtitle = "", ...) {
  if (!missing(name) && !any(name ==  resultsNames(dd))) 
    stop("bad coefficient ", name, ". Options are ", paste0(resultsNames(dd), 
                                                            collapse = ", "), ".")
  res = results(dd, name = name, ...) %>% as.data.frame() %>% 
    rownames_to_column("SYMBOL") %>% mutate(sign = sign(log2FoldChange)) %>% 
    group_by(sign) %>% mutate(p_rank = rank(pvalue), label = ifelse(p_rank < 
                                                                      50, SYMBOL, ""))
  res[is.na(res$padj), "padj"] = 1
  trans = vst(dd)
  sub = trans[filter(res, padj < 0.1, p_rank < heatmap_max_gene)$SYMBOL, 
              subset_idx]
  if (nrow(sub) > 0) {
    assay(sub, "zscore") = t(scale(t(assay(sub))))
    if (!("set" %in% names(rowData(sub)))) {
      rowData(sub)$set = NA
    }
    catmat2 = rowData(sub)[, "set", drop = FALSE]
    h = ComplexHeatmap::Heatmap(assay(sub, "zscore"), top_annotation = ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(colData(sub)[, 
                                                                                                                                         heatmap_top_group]), which = "column"), name = "Z-scored\nNormalized Exp", 
                                row_names_gp = gpar(fontsize = 4), clustering_distance_rows = "spearman", 
                                clustering_distance_columns = "spearman", column_names_gp = gpar(fontsize = 4))
    print(h)
  }
  else {
    message("No significant comparisons in ", title, "(", 
            subtitle, ").")
  }
  ggplot(res, aes(x = log2FoldChange, y = -clamp(log10(pvalue), 
                                                 20), color = padj < 0.1)) + geom_point() + geom_text_repel(data = filter(res, 
                                                                                                                          p_rank < 25), mapping = aes(label = label), size = 2) + 
    geom_vline(lty = 2, xintercept = 0) + theme_minimal() + 
    scale_color_discrete("FDR", labels = c(`TRUE` = "<10%", 
                                           `FALSE` = ">=10%")) + ggtitle(title, subtitle = subtitle)
}