rmarkdown::render('phase1_low_input.Rmd', output_format = 'github_document')

rmarkdown::render('02_qc_dea.Rmd',
                  params = list(input_root = 'private/amp_lowinput_se.rds',
                                se_subset = 
                                  "se = se[, !is.na(se$nr4a_total_score) & se$Cell.type != 'B cell']
                                  assays(se) = list(counts = floor(assay(se, 'counts')))
                                ",
                                design_csv = 'extradata/nr4a_design.csv',
                                sample_id = 'Sample.ID',
                                terms_annotation_csv = 'extradata/nr4a_dea_terms.csv',
                                cache = TRUE,
                                organism = 'human',
                                plot_covariates = c('Cell.type', 'nr4a_total_score')
                  ), output_format = 'github_document')
