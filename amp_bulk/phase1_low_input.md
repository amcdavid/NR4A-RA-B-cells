AMP phase 1 lowinput
================

Notes: Only need `des_leuk` section…

    ## Warning: package 'AnnotationDbi' was built under R version 4.1.2

## Load AMP low input data

    ## 'select()' returned 1:many mapping between keys and columns

    ## Rows: 192 Columns: 67

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (38): Donor.ID, Sample.ID, Sample.Date, Cell.type, Tissue.type, Alternat...
    ## dbl (28): Number.of.cells, cDNA.concentration, Year.of.Dx, RF.result, CRP.in...

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## Rows: 68 Columns: 7

    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): sample, anti_ccp_category
    ## dbl (5): duration_yrs, anti_ccp_or_ccp2, auto_das28, krenn_inflam, anti_ccp_num

    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## Joining, by = c("Cell.type", "Sample.ID")

    ## Joining, by = "Donor.ID"

    ## Joining, by = "Tissue.type"

    ## Joining, by = "Donor.ID"

## Severity/NRA4A1

    ## converting counts to integer mode

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 4770 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## Joining, by = "ENSEMBL"

    ## Joining, by = "SYMBOL"

![](phase1_low_input_files/figure-gfm/des_tissue-1.png)<!-- -->

## Biopsy selected genes

    ## Joining, by = "ENSEMBL"

    ## Joining, by = "SYMBOL"

    ## Joining, by = "term"

![](phase1_low_input_files/figure-gfm/tissue_selected_genes-1.png)<!-- -->

### Categorical mahanoblis

    ## converting counts to integer mode

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## estimating size factors

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## final dispersion estimates

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## fitting model and testing

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## -- replacing outliers and refitting for 4281 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## fitting model and testing

    ##   Note: levels of factors in the design contain characters other than
    ##   letters, numbers, '_' and '.'. It is recommended (but not required) to use
    ##   only letters, numbers, and delimiters '_' or '.', as these are safe characters
    ##   for column names in R. [This is a message, not a warning or an error]

    ## Joining, by = "ENSEMBL"

    ## Joining, by = "SYMBOL"

![](phase1_low_input_files/figure-gfm/des_leuk-1.png)<!-- -->

    ## Joining, by = "ENSEMBL"
    ## Joining, by = "SYMBOL"

    ## Joining, by = "term"

![](phase1_low_input_files/figure-gfm/des_leuk-2.png)<!-- -->

### Check DAS

|      | Donor.ID | Sample.ID | tissue | Auto.calculation.of.DAS28.CRP | Manual.DAS.28.CRP | anti\_ccp\_num | anti\_ccp\_or\_ccp2 |
|:-----|:---------|:----------|:-------|:------------------------------|:------------------|---------------:|--------------------:|
| S232 | 300-0512 | S232      | RAB    | see manual DAS28 CRP          | DAS28CRP: 7.58    |              3 |               340.0 |
| S68  | 301-0122 | S68       | RAA    | NA                            | NA                |              3 |               140.0 |
| S100 | 300-0481 | S100      | RAB    | 6.33                          | DAS28: 5.1        |              1 |                 1.6 |
| S108 | 300-0483 | S108      | RAB    | 5.87                          | DAS28: 6.77       |              1 |                 0.0 |
| S12  | 301-0156 | S12       | RAA    | 3.87                          | NA                |              3 |                64.4 |
| S20  | 301-0158 | S20       | RAA    | NA                            | NA                |              1 |                 0.0 |

## DAS/Tissue/B cells

![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-6.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-7.png)<!-- -->![](phase1_low_input_files/figure-gfm/unnamed-chunk-2-8.png)<!-- -->
