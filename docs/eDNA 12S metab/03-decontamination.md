Sub-sampling and Decontamination
================

**.Rmd script**

This script provides code for optional sub-sampling and decontamination
based on NTCs and blanks.

## Load libraries

``` r
library(tidyverse) # for data transformation
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   4.0.0     ✔ tibble    3.3.0
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.1.0     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2) # for plotting
library(readxl) ## for reading in excel files
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(hrbrthemes)
library(ggrepel)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(Rmisc)
```

    ## Loading required package: lattice
    ## Loading required package: plyr
    ## ------------------------------------------------------------------------------
    ## You have loaded plyr after dplyr - this is likely to cause problems.
    ## If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    ## library(plyr); library(dplyr)
    ## ------------------------------------------------------------------------------
    ## 
    ## Attaching package: 'plyr'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     arrange, count, desc, failwith, id, mutate, rename, summarise,
    ##     summarize
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     compact

``` r
library(writexl)
library(vegan) ## rareify 
```

    ## Loading required package: permute

``` r
library(patchwork)
```

    ## 
    ## Attaching package: 'patchwork'
    ## 
    ## The following object is masked from 'package:cowplot':
    ## 
    ##     align_plots

``` r
library(microDecon)

# removing scientific notation
## remove line or comment out if not desired 
options(scipen=999)
set.seed(123)
```

## Load data

For the sake of demonstration, I added a NTC to the results matrix.
Change the file name for own data.

``` r
### loading RData if applicable uncomment
# load("../results/Raw_data.RData")
# data <- ASV_table_taxID_collapsed

data <- read_xlsx("example_output/Results_initial/Results_rawreads_matrix_mockNTC.xlsx") %>%
  dplyr::rename()

meta <- read.csv("example_input/SBNMS_metadata.csv") %>%
  mutate(sampleID = gsub("-", "_", sampleID))
```

## Read Counts prior to decontamination

Visualize the read count distribution. There is one sample that has more
reads than the rest of the samples (\>200,000 in this case).

``` r
## Calculate total read numbers 
totals <- data.frame(
  sampleID = colnames(data)[4:ncol(data)],   # column names for sample counts
  nonchim  = colSums(data[, 4:ncol(data)])   # total reads per sample
) %>% left_join(., meta, by = "sampleID")

total_presub <- totals %>% ggplot(., aes(x=SampleType, y=nonchim)) +
  geom_boxplot(outlier.shape = NA, aes(color=SampleType)) +
  labs(
    x = "Sample Type", y = "Total reads pre sub-sample"
  ) +
  geom_jitter(aes(fill=SampleType), shape=21, size=1.5, alpha=0.5, width=0.25) + 
  theme_bw() +
    theme(
      legend.position = "none"
    ); total_presub
```

![](03-decontamination_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Sub-sample above 200,000

Change this cut-off to what you’d like it to be based on your data

``` r
## Compute total reads per sample
sample_depths <- colSums(data[, 4:ncol(data)])  # each column is a sample
summary(sample_depths)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     226   54462   62748   66998   76951  244276

``` r
## Define target depth
### CUSTOMIZE VALUE AS NEEDED 
target_depth <- quantile(sample_depths, 0.99)
```

Edit based on the target depth above

``` r
spp_info <- data[, 1:3]
counts <- data[, 4:ncol(data)]
counts_t <- t(counts)

counts_adj <- counts_t
high_depth <- rowSums(counts_t) > target_depth

n_high <- sum(high_depth)

if (n_high == 1) {
  counts_adj[high_depth, ] <- rrarefy(
    matrix(counts_t[high_depth, ], nrow = 1),
    sample = target_depth
  )
} else if (n_high > 1) {
  counts_adj[high_depth, ] <- rrarefy(counts_t[high_depth, ], sample = target_depth)
}
```

    ## Warning in rrarefy(counts_t[high_depth, ], sample = target_depth): function
    ## should be used for observed counts, but smallest count is 11

``` r
counts_adj <- t(counts_adj)
df_adj <- cbind(spp_info, counts_adj)
```

Confirming that only the desired samples were changed

``` r
new_totals <- data.frame(
  sampleID = colnames(df_adj)[4:ncol(df_adj)],   # column names for sample counts
  new_nonchim  = colSums(df_adj[, 4:ncol(df_adj)])   # total reads per sample
)

totals_adj <- totals %>% full_join(new_totals) %>%
  mutate(difference = new_nonchim-nonchim)
```

    ## Joining with `by = join_by(sampleID)`

``` r
## confirm no species were lost (output should be empty)
zero_rows <- which(rowSums(df_adj[, 4:ncol(df_adj)]) == 0)
```

## Visualizing afterwards

``` r
total_after_sub <- totals_adj %>% ggplot(., aes(x=SampleType, y=new_nonchim)) +
  geom_boxplot(outlier.shape = NA, aes(color=SampleType)) +
  labs(
    x = "Sample Type", y = "Total reads post sub-sample"
  ) +
  geom_jitter(aes(fill=SampleType), shape=21, size=1.5, alpha=0.5, width=0.25) + 
  theme_bw() +
    theme(
      legend.position = "none"
    ); total_after_sub
```

![](03-decontamination_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
combined_plot <- total_presub + total_after_sub + plot_layout(ncol = 2)
ggsave("example_output/Figures/Filtering_subsample_before_after.png", combined_plot, width = 10, height = 5, dpi = 300)
```

## Microdecon

<https://github.com/donaldtmcknight/microDecon>

At this point, remove positive controls from the dataset (not applicable
to this example).

``` r
## Add annotation as needed for the name of the NTC or blank
df_for_microDecon <- df_adj %>%
  dplyr::select(-(Common_name:Category)) %>%
  relocate(matches("NTC|neg|Blank"), .after = Species_name)

## negative and blanks, columns 2
## sample columns 3:144

decontaminated <- decon(data = df_for_microDecon, 
                        
                        ## number of blank columns [integer]
                        numb.blanks=1, 
                        
                        ## number of individuals per group [integer]
                        ## number of samples - 1 (still troubleshooting why this is - 1..)
                        numb.ind=c(141),
                        
                        ## Indicate if taxa columns are present (T/F)
                        taxa=T)

removed <- decontaminated$reads.removed
save(removed, file = "example_output/Results_decontaminated/Decontaminated_reads_removed.RData")

decontamin_df <- decontaminated$decon.table

## if multiple blanks, a new column called Mean.blank is created, if not the column name remains the same
decontamin_df <- decontamin_df %>% 
  #dplyr::select(-Mean.blank)
  dplyr::select(-NTC)

df_filtered <- df_adj %>% dplyr::select(Species_name:Category) %>% right_join(., decontamin_df)
```

    ## Joining with `by = join_by(Species_name)`

``` r
nrow(df_adj) ## Number spp before microdecon
```

    ## [1] 43

``` r
nrow(df_filtered) ## Number spp after microdecon
```

    ## [1] 40

Exporting dataframe

``` r
df_filtered %>% write_xlsx("example_output/Results_decontaminated/Results_read_count_decontaminated.xlsx")

df_relab <- df_filtered %>% 
  gather("sampleID", "reads", 4:last_col()) %>%
  dplyr::group_by(sampleID) %>%

  ### total
  dplyr::mutate(sample_total = sum(reads)) %>%
  dplyr::group_by(sampleID, Species_name) %>%

  ## relab
  dplyr::mutate(relab = reads/sample_total) %>% ungroup() %>%
  select(-reads, -sample_total) 

df_relab %>% write_xlsx("example_output/Results_decontaminated/Results_relative_abundance_decontaminated.xlsx")
```
