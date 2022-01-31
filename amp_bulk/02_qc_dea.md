QC and Differential Expression for amp\_bulk project
================
Andrew McDavid
2022-01-31

Notes:

Knit with:

``` r
print(params)
```

    ## $output_root
    ## [1] "refined/02_qc_dea"
    ## 
    ## $features_csv
    ## NULL
    ## 
    ## $fdr_level
    ## [1] 0.1
    ## 
    ## $dea_method
    ## [1] "deseq2"
    ## 
    ## $input_root
    ## [1] "private/amp_lowinput_se.rds"
    ## 
    ## $se_subset
    ## [1] "se = se[, !is.na(se$nr4a_total_score) & se$Cell.type != 'B cell']\n                                  assays(se) = list(counts = floor(assay(se, 'counts')))\n                                "
    ## 
    ## $design_csv
    ## [1] "extradata/nr4a_design.csv"
    ## 
    ## $sample_id
    ## [1] "Sample.ID"
    ## 
    ## $terms_annotation_csv
    ## [1] "extradata/nr4a_dea_terms.csv"
    ## 
    ## $cache
    ## [1] TRUE
    ## 
    ## $organism
    ## [1] "human"
    ## 
    ## $plot_covariates
    ## [1] "Cell.type"        "nr4a_total_score"

``` r
se = readRDS(params$input_root)
```

# Cross-type associations with NR4As

``` r
genes3 = tibble(SYMBOL = c('NR4A1', 'NR4A2', 'NR4A3', 'CD69', 'CD83', 'LMNA', 'GPR183', "LY9", "CXCR4"))

genes_nr4a = left_join(genes3, rowData(se) %>% as.data.frame() %>% dplyr::select(SYMBOL, ENSEMBL))
se_b2 = se[, se$Cell.type == 'B cell']
nr4a_total_score = colMeans(assay(se_b2, 'lcpm')[genes_nr4a$ENSEMBL,])
nr4a2 = assay(se_b2, 'lcpm')[filter(genes_nr4a, SYMBOL == 'NR4A2')$ENSEMBL,]
nr4a_covariates = data.frame(Donor.ID = se_b2$Donor.ID, nr4a_total_score = scale(nr4a_total_score), nr4a2 = scale(nr4a2))

ggplot(nr4a_covariates, aes(x = nr4a_total_score, y = nr4a2, label = Donor.ID)) + geom_text()
```

![](02_qc_dea_files/figure-gfm/setup-cross-1.png)<!-- -->

``` r
colData(se) = left_join(as.data.frame(colData(se)), nr4a_covariates) %>% DataFrame()
```

``` r
se = se[, !is.na(se$nr4a_total_score) & se$Cell.type != 'B cell']
                                  assays(se) = list(counts = floor(assay(se, 'counts')))
```

``` r
designs = read_csv(params$design_csv, comment = "#")
```

``` r
if(ncol(colData(se))<1) stop("Expecting non-empty colData")
if(is.null(params$plot_covariates)){
  plot_covariates = rep(names(colData(se)), 2)[1:2] # ensure at least length two
} else{
  plot_covariates = params$plot_covariates
}
```

``` r
dge = DESeqDataSet(se=se, design = as.formula(designs$formula[1]))
```

``` r
#design(dge) = ~ groupf

#genes_plus_random = union(which(!is.na(rd$group)), sample(nrow(rd), 200))

logcpm = function(x){
    total = colSums(counts(x))
 log2( 1e6*t(counts(x))/total + 1)
}

colData(dge) = cbind(colData(dge), logcpm(dge)[, !is.na(rowData(dge)$group)])
```

``` r
dge_f = dge#[!duplicated(dge$genes$ENSEMBL) | is.na(dge$genes$ENSEMBL),]
#dge_f = dge[,dge$group == params$group]
dge_f = estimateSizeFactors(dge_f)
null = rowSums(counts(dge_f) > 5) < 3
dge_f = dge_f[!null,]
ds_trans = vst(dge_f)
```

# Differential expression testing

Currently we test 3 different models. Specifications are shown below:

``` r
knitr::kable(designs)
```

| name       | formula              | subset\_expr        | contrast |
|:-----------|:---------------------|:--------------------|:---------|
| T cells    | \~nr4a\_total\_score | Cell.type==‘T cell’ | NA       |
| Fibroblast | \~nr4a\_total\_score | Cell.type==‘Fibro’  | NA       |
| Monocytes  | \~nr4a\_total\_score | Cell.type==‘Mono’   | NA       |

``` r
dea_method = match.arg(params$dea_method, c('deseq2', 'nebula'))
```

``` r
dew = purrr::map2(setNames(designs$formula, designs$name), designs$subset_expr, function(x, y) run_deseq_design(dge_f,  as.formula(x), y))
```

``` r
dew_df = purrr::map_dfr(dew, tidy, .id = 'model')  %>% dplyr::rename(SYMBOL = gene, padj = p.adjusted, pvalue = p.value, log2FoldChange = estimate)
models_and_terms = dew_df %>% count(model, term)
knitr::kable(models_and_terms)
```

| model      | term               |     n |
|:-----------|:-------------------|------:|
| Fibroblast | Intercept          | 30521 |
| Fibroblast | nr4a\_total\_score | 30521 |
| Monocytes  | Intercept          | 30521 |
| Monocytes  | nr4a\_total\_score | 30521 |
| T cells    | Intercept          | 30521 |
| T cells    | nr4a\_total\_score | 30521 |

Models fit and terms returned.

``` r
models = read_csv(params$terms_annotation_csv)
```

``` r
count_sig = semi_join(dew_df, models) %>% group_by(SYMBOL) %>% 
    summarize(n_sig = sum(padj < .1, na.rm = TRUE), geomean_pvalue = sum(log10(pvalue), na.rm = TRUE), geovar_pvalue = var(log10(pvalue), na.rm = TRUE))

sig_contrasts = semi_join(dew_df, count_sig %>% filter(rank(geomean_pvalue)<20)) %>% left_join(models)

plt = ggplot(sig_contrasts, aes(y = SYMBOL, x = log2FoldChange, color = model)) + geom_text(aes(label = term_class))
plt + geom_vline(xintercept = 0, lty = 2)
```

![](02_qc_dea_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
models_plot = filter(models, plot_volcano)

for(m in seq_len(nrow(models_plot))){
  this_model = models_plot[[m,'model']]
  print(modified_volc(dew[[this_model]],  name = models_plot[[m, 'term']], heatmap_top_group = plot_covariates, heatmap_max_gene = 40,
                      title = sprintf("%s in %s background", models_plot[[m, 'term_class']], this_model)))
                               # ^^^ Modify as needed
}
```

![](02_qc_dea_files/figure-gfm/heatmaps-1.png)<!-- -->![](02_qc_dea_files/figure-gfm/heatmaps-2.png)<!-- -->![](02_qc_dea_files/figure-gfm/heatmaps-3.png)<!-- -->![](02_qc_dea_files/figure-gfm/heatmaps-4.png)<!-- -->![](02_qc_dea_files/figure-gfm/heatmaps-5.png)<!-- -->![](02_qc_dea_files/figure-gfm/heatmaps-6.png)<!-- -->

## Venn diagrams of results

``` r
terms_by_model = left_join(models, dew_df) %>% ungroup() %>% 
  filter(padj < params$fdr_level) %>% mutate(direction = sign(log2FoldChange))

terms_by_model_list = terms_by_model %>% split(.$term_class)

for(i in seq_along(terms_by_model_list)){
  gene_subset = terms_by_model_list[[i]] %>% select(SYMBOL, model)
  n_model = count(gene_subset, model)
  gene_subset = left_join(gene_subset, n_model) %>% 
    transmute(SYMBOL, model = glue::glue("{model}({n} genes)"))
  v = venneuler::venneuler(gene_subset)
  
  plot(v, main = names(terms_by_model_list)[i])
}
```

Genes significant at 0.1 by model group.

Venn diagram

# Gene set enrichment

``` r
library(clusterProfiler)

db = switch(params$organism[1],
            mouse = "org.Mm.eg.db",
            human = "org.Hs.eg.db")

remap = bitr(unique(dew_df$SYMBOL), 'ENSEMBL', 'ENTREZID', db, drop = FALSE)

# take subset of term, and create model below...
gsea_genes = left_join(dew_df, remap, by = c('SYMBOL' = 'ENSEMBL')) %>% 
  mutate(direction = factor(sign(log2FoldChange), labels = c('-1' = 'down', '1' = 'up'))) %>% 
  left_join(models) %>% 
  filter(padj < params$fdr_level, abs(log2FoldChange) > 1, plot_volcano) %>% 
  group_by(term, direction)

term_direction_drop = gsea_genes %>% 
  group_by(model, direction) %>%
  summarize(n = n()) %>% 
  filter(n < 10)
```

``` r
knitr::kable(term_direction_drop)
```

| model | direction |   n |
|:------|:----------|----:|

``` r
gsea_genes = gsea_genes %>% anti_join(term_direction_drop)
```

Dropping the following sets for having &lt; 10 enriched genes.

``` r
cc_go = compareCluster(ENTREZID ~ model + direction, 
                       data = gsea_genes, 
                       OrgDb = db, readable = TRUE, universe = unique(remap$ENTREZID), 
                       fun = 'enrichGO')
```

``` r
dotplot(cc_go, font.size = 6, showCategory = 12) + 
  theme(axis.text.x = element_text(angle = 90))
```

![](02_qc_dea_files/figure-gfm/clusterGO-1.png)<!-- -->![](02_qc_dea_files/figure-gfm/clusterGO-2.png)<!-- -->

``` r
all_go = rbind(cc_go@compareClusterResult,
      cc_bmt@compareClusterResult,
      cc_c7@compareClusterResult)

readr::write_csv(all_go, path = 'enrichment_by_subpop.csv')
readr::write_csv(dew_df, path = 'tests_by_subpop.csv')
```

``` r
to_write = left_join(dew_df, models_plot)
write_csv(to_write, paste0(params$output_root, '_results.csv'))
```
