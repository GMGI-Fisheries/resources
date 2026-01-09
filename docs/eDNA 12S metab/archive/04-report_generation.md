Report Generation for eDNA service clients
================

**.Rmd script**

This script evaluates your sequence quality and taxonomic assignment
quality. Figures produced in this script can go into supplemental data
for a manuscript.

## Load libraries

Download if needed (one-time only, uncomment if needed)

``` r
#install.packages("BiocManager")
# BiocManager::install("phyloseq")
# BiocManager::install("microbiome")

#install.packages("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
```

``` r
## Core tidyverse and utilities
library(tidyverse)    # includes ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats
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
library(readxl)       # read Excel files
library(writexl)      # write Excel files
library(knitr)        # for knitting reports
library(cowplot)      # plot arrangement and themes
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(ggrepel)      # text repulsion in ggplot
library(viridis)      # color scales
```

    ## Loading required package: viridisLite

``` r
library(hrbrthemes)   # themes for ggplot
library(scales)       # scale utilities
```

    ## 
    ## Attaching package: 'scales'
    ## 
    ## The following object is masked from 'package:viridis':
    ## 
    ##     viridis_pal
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard
    ## 
    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
## Specialized analysis
library(phyloseq)     # ecological stats
library(vegan)        # ecological stats
```

    ## Loading required package: permute

``` r
library(microbiome)   # alpha diversity
```

    ## 
    ## microbiome R package (microbiome.github.com)
    ##     
    ## 
    ## 
    ##  Copyright (C) 2011-2022 Leo Lahti, 
    ##     Sudarshan Shetty et al. <microbiome.github.io>
    ## 
    ## 
    ## Attaching package: 'microbiome'
    ## 
    ## The following object is masked from 'package:vegan':
    ## 
    ##     diversity
    ## 
    ## The following object is masked from 'package:scales':
    ## 
    ##     alpha
    ## 
    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     alpha
    ## 
    ## The following object is masked from 'package:base':
    ## 
    ##     transform

``` r
library(pairwiseAdonis)
```

    ## Loading required package: cluster

``` r
library(lme4)
```

    ## Loading required package: Matrix
    ## 
    ## Attaching package: 'Matrix'
    ## 
    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(car)
```

    ## Loading required package: carData
    ## 
    ## Attaching package: 'car'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     some

``` r
library(stats)
library(funrar)       # relative abundance
library(ggsankey)     # sankey plots
library(naniar)       # missing data handling
library(ggh4x)        # advanced facetting
library(Rmisc)        # summary statistics
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
## Settings
set.seed(1234)
options(scipen = 999)
```

## Load data

User edits:  
1. Replace the 3 paths below: edit example_input to your project
specific path  
2. Confirm your sampleIDs match between metadata, results df, and
filtering stats output

``` r
## METADATA
meta <- read.csv("example_input/SBNMS_metadata.csv") %>%
  mutate(sampleID = gsub("-", "_", sampleID))

## SEQ SUMMARY
filtering_stats <- read_tsv("example_input/overall_summary.tsv", show_col_types = FALSE) %>% 
  dplyr::rename(sampleID = sample)

## RESULTS MATRIX AND LONG FORMATS
df_raw <- read_xlsx("example_output/Results_initial/Results_rawreads_matrix.xlsx")
df_raw_long <- read_xlsx("example_output/Results_initial/Results_rawreads_long_format.xlsx")
df_relative <- read_xlsx("example_output/Results_initial/Results_relative_abundance_matrix.xlsx")
df_relative_long <- read_xlsx("example_output/Results_initial/Results_relative_abundance_long_format.xlsx")

## OPTIONAL: DECONTAMINATION READS
df_decont_raw <- read_xlsx("example_output/Results_decontaminated/Results_read_count_decontaminated.xlsx")
df_decont_relative <- read_xlsx("example_output/Results_decontaminated/Results_relative_abundance_decontaminated.xlsx")

## ASV INFORMATION
asv_breakdown <- read_xlsx("example_output/Breakdown_ASV_level.xlsx")
```

Adding Taxonomic Level information for heatmap df

``` r
## GMGI DB TAX INFORMATION
taxlevels <- read_excel(
  "C:/BoxDrive/Box/Science/Fisheries/Projects/eDNA/Metabarcoding Lab Resources/Reference Databases/GMGI_Vert_Ref.xlsx") %>% 
  dplyr::select("Species_name", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "species") %>%
  distinct()

## COMMON NAME SHEET 
commonNames_annotated <- read_xlsx("example_output/Taxonomic_assignments/CommonNames_required_edited.xlsx")

## DF for heatmap
heatmap_df <- df_relative_long %>% left_join(., taxlevels, by = "Species_name")

# Loop through each row of the dataframe to add taxonomic level information from required edited worksheet 
for (i in commonNames_annotated$Species_name) {
  # Extract the current row (will do this for each ASV_ID in the choice df)
  current_row <- commonNames_annotated %>% subset(Species_name==i)
  
  # Apply filter based on the current row's condition
heatmap_df <- heatmap_df %>%
  mutate(across(c(Common_name, Category, Kingdom, Phylum, Class, Order, Family, Genus, species),
                ~case_when(Species_name == current_row$Species_name ~ current_row[[cur_column()]],
                           TRUE ~ .x)))
}

## tax list 
df_tax <- heatmap_df %>% dplyr::select(Species_name, Kingdom:Genus) %>% distinct()
```

### Determine order of categories desired for visualization

``` r
unique(df_raw$Category)
```

    ## [1] "Teleost Fish"  "Bird"          "Marine Mammal" "Livestock"    
    ## [5] "Human"         "Unassigned"

``` r
category_order <- c(
  "Human", "Livestock", "Unassigned", "Bird", "Marine Mammal", "Teleost Fish"
                    )
```

## Sequence data

### Data Transformation

No user edits.

``` r
df <- full_join(filtering_stats, meta, by = "sampleID") %>%
  # filtering out columns we don't need 
  dplyr::select(-cutadapt_reverse_complemented) %>%
  
  # removing percentage icon from cutadapt_passing_filters_percent
  dplyr::mutate(cutadapt_passing_filters_percent = gsub("%", "", cutadapt_passing_filters_percent)) %>%
  
  # confirming that all columns of interest are numerical 
  dplyr::mutate_at(vars(2:10), as.numeric) %>%
  
  # data transformation so all columns of interest are together 
  gather("measure", "value", 2:10)  

head(df)
```

    ## # A tibble: 6 × 7
    ##   sampleID    Site  Section    Substrate SampleType    measure             value
    ##   <chr>       <chr> <chr>      <chr>     <chr>         <chr>               <dbl>
    ## 1 C1_bottom   C1    Central SB Water     Bottom Water  cutadapt_total_pr…  93700
    ## 2 C1_surface  C1    Central SB Water     Surface Water cutadapt_total_pr…  83038
    ## 3 C10_surface C10   Central SB Water     Surface Water cutadapt_total_pr… 151479
    ## 4 C10a        C10   Central SB Sediment  Sediment      cutadapt_total_pr…  72926
    ## 5 C10b        C10   Central SB Sediment  Sediment      cutadapt_total_pr…  64564
    ## 6 C10c        C10   Central SB Sediment  Sediment      cutadapt_total_pr… 108626

### Plotting

Suggested webpage to choose colors: <https://coolors.co/>

``` r
### User edits:
### 1. Change paths of output to desired folder (data/figures is suggested data structure)
### 2. Change x axis and color, size based on metadata desired 
### 3. Change custom colors and sizes if desired and number of colors and sizes based on metadata variable chosen 

df %>% 
  
  ## USER EDITS IN LINE BELOW 
  ggplot(., aes(x=SampleType, y=value)) + 
  
  ## adding points in jitter format 
  geom_jitter(width=0.15, alpha=0.5, fill="#0077b6", color='black', size=1, shape=21) + 
  
  ## option for additional boxplots if desired (uncomment to add)
  #geom_boxplot() +
  
  ## using facet_wrap to create grid based on variables and factor() to order them in custom format
  facet_wrap(~factor(measure, levels=c('cutadapt_total_processed', 'cutadapt_passing_filters', 
                                       'cutadapt_passing_filters_percent', 'DADA2_input',
                                 'filtered', 'denoisedF', 'denoisedR', 'merged', 'nonchim')), scales = "free") +
  
  ## graph asthetics 
  theme_bw() +
  ylab("Number of reads") + 

  theme(panel.background=element_rect(fill='white', colour='black'),
        strip.background=element_rect(fill='white', colour='black'),
        strip.text = element_text(size = 10, face="bold"),
        legend.position = "right",
        axis.text.y = element_text(size=7, color="grey30"),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=11, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=11, face="bold"))
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave("example_output/Figures/FilteringStats_full.png", width = 11, height=9)
```

Condensed plot used in contract reporting.

``` r
### User edits:
### 1. Change paths of output to desired folder (data/figures is suggested data structure)
### 2. Change x axis and color, size based on metadata desired 
### 3. Change custom colors and sizes if desired and number of colors and sizes based on metadata variable chosen 

df %>% 
  subset(measure == "cutadapt_total_processed" | measure == "nonchim") %>%

  ## USER EDITS IN LINE BELOW 
  ggplot(., aes(x=measure, y=value)) + 
  #geom_hline(yintercept = cutoff_low, color = "grey80", linetype = "dashed") +
  
  ## adding points in jitter format 
  geom_boxplot(outlier.shape = NA, fill=NA, aes(color = measure)) +
  geom_jitter(width=0.15, shape=21, alpha=0.25, size=1.5, color = 'black', aes(fill = measure)) + 

  ## graph asthetics 
  theme_bw() +
  labs(
    y="Number of reads",
    x="Bioinformatic Step"
    ) + 
  
  ## USER EDITS IN MANUAL CODE BELOW 
  # scale_color_manual(values = c("red3", "lightblue", "purple2", "gold", "green4", "black")) +
  # scale_size_manual(values = c(21,17)) +
  scale_fill_manual(values = c("#264653", "#2a9d8f")) +
  scale_color_manual(values = c("#264653", "#2a9d8f")) +
  scale_x_discrete(labels = c("cutadapt_total_processed" = "Start", 
                              "nonchim" = "Final")) +
  
  theme(panel.background=element_rect(fill='white', colour='black'),
        strip.background=element_rect(fill='white', colour='black'),
        strip.text = element_text(size = 10, face="bold"),
        legend.position = "none",
        axis.text.y = element_text(size=10, color="black"),
        axis.text.x = element_text(size=10, color="black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=11, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=11, face="bold"))
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_FilteringStats_Condensed.png", width = 4, height=4)
```

## Plot taxonomy

### Data transformation

User edits: Change decimal places as needed.

``` r
results_summary <- df_raw_long %>% 
  mutate(Category = factor(Category, levels = category_order)) %>%
  dplyr::group_by(Category) %>%
  reframe(sum_reads = sum(reads))

general_stats <- results_summary %>% 
  dplyr::mutate(total = sum(sum_reads),
         percent = sum_reads/total*100) %>% dplyr::select(Category, percent) %>% distinct() %>%
  ## round to a certain # of decimal places 
  dplyr::mutate(across(c('percent'), round, 4))
```

    ## Warning: There was 1 warning in `dplyr::mutate()`.
    ## ℹ In argument: `across(c("percent"), round, 4)`.
    ## Caused by warning:
    ## ! The `...` argument of `across()` is deprecated as of dplyr 1.1.0.
    ## Supply arguments directly to `.fns` through an anonymous function instead.
    ## 
    ##   # Previously
    ##   across(a:b, mean, na.rm = TRUE)
    ## 
    ##   # Now
    ##   across(a:b, \(x) mean(x, na.rm = TRUE))

``` r
ASV_summary <- asv_breakdown %>%
  mutate(Category = factor(Category, levels = category_order)) %>%
  dplyr::group_by(Category) %>%
  reframe(count = n_distinct(ASV_ID))

species_summary <- df_raw_long %>%
  mutate(Category = factor(Category, levels = category_order))  %>%
  dplyr::group_by(Category) %>%
  reframe(count = n_distinct(Species_name))
```

### Graph color and order options

User edits: Uncomment/comment colors as needed

``` r
category_order
```

    ## [1] "Human"         "Livestock"     "Unassigned"    "Bird"         
    ## [5] "Marine Mammal" "Teleost Fish"

``` r
fill_colors <- c("Human" = "#FE9E20", 
                 "Livestock" = "#FCCA46", 
                 "Other" = "#FFECC2", 
                 "Unassigned" = "grey85",
                 "Bird" = "#C8D2B1",
                 "Sea Turtle" = "#A1C181",
                 "Elasmobranch" = "#97ADCB",
                 "Marine Mammal" ="#DCF1F9",
                 "Teleost Fish" = "#85B6CB"
                 )
```

### Relative abundance by category

``` r
Relab_category <- df_raw_long %>% 
  
  ## Calculate reads per category per sample ID
  dplyr::group_by(sampleID, Category) %>%
  reframe(group_sum = sum(reads)) %>%
  
  ## Calculate relative abundance of that value
  dplyr::group_by(sampleID) %>%
  mutate(
    sample_total = sum(group_sum),
    group_relab = group_sum/sample_total
  ) 

Relab_category %>%
  
ggplot(., aes(y=group_relab, x=Category)) + 
  geom_boxplot(outlier.shape = NA, aes(color=Category), fill=NA) +
  geom_jitter(aes(fill=Category), width=0.2, shape=21, color='black', size = 0.75, alpha=0.35) +
    scale_color_manual(values = fill_colors) +
    scale_fill_manual(values = fill_colors) +
    labs(color = "Category") +
    scale_x_discrete(labels = scales::label_wrap(10), limits = category_order) +
    theme_bw() + 
    xlab("Category") + ylab("Relative abundance") +
    theme(panel.background=element_rect(fill='white', colour='black'),
          legend.position = "none",
        axis.text.y = element_text(size=7, color="grey30"),
        axis.text.x = element_text(size=7), # angle=45, hjust=1
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=11, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=11, face="bold"))
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_categories_relative_abundance.png", width = 6.5, height = 4)
```

### Human specific reads by Sample Type

``` r
Relab_category %>%
  
  ## subset to include Human and metadata 
  subset(Category == "Human") %>%
  left_join(., meta, by = "sampleID") %>%
  
  ggplot(., aes(x=SampleType, y=group_relab)) +
  geom_boxplot(outlier.shape=NA, fill=NA, color = "#e76f51") +
  geom_jitter(fill="#e76f51", shape=21, alpha=0.5, width=0.2) +
  labs(
    x = "Sample Type",
    y = "Relative Abundance"
  ) +
  theme_bw() +
  theme(panel.background=element_rect(fill='white', colour='black'),
          legend.position = "none",
        axis.text.y = element_text(size=7, color="grey30"),
        axis.text.x = element_text(size=7),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=11, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=11, face="bold"))
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Human.png", width=4, height=3)
```

### Piechart

Reads Piechart. Used in contract report.

``` r
### User edits:
### 1. Change paths of output to desired folder (data/figures is suggested data structure)
### 2. Change scale brewer color if desired 

piechart <- general_stats %>%
  dplyr::mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos)) %>%
  dplyr::mutate_if(is.numeric, round, digits = 4)

unassigned_percent <- piechart %>% 
  filter(Category == "Unassigned") %>% 
  pull(percent)

piechart_reads <- general_stats %>% 
  ggplot(., aes(x="", y = percent, fill = Category)) +
  geom_col(color = "black", width=1.25) +
  geom_label_repel(data = piechart,
                   aes(y = pos, label = paste0(percent, "%")),
                   size = 3, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_colors) +
  theme_bw() +
  # labs(
  #   x = NULL, y = NULL, fill = "Category", 
  #      caption = paste0("Unassigned reads = ", unassigned_percent, "%")
  # ) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    legend.position = "none",
    panel.border = element_blank(),  # Remove panel border
    panel.grid = element_blank(),    # Remove grid lines
    axis.ticks = element_blank(),    # Remove axis ticks
    axis.text = element_blank(),     # Remove axis text
    axis.title = element_blank()     # Remove axis titles
      ) +
  ggtitle("Raw reads (%)") +
  xlab("") + ylab("") + labs(fill = "Category"); piechart_reads
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Category_breakdown_percent_rawreads.png", width=4, height=3)
```

ASV Piechart (NOT used in contract report)

``` r
### User edits:
### 1. Change paths of output to desired folder (data/figures is suggested data structure)
### 2. Change scale brewer color if desired 

piechart_ASV <- ASV_summary %>%  
  dplyr::mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))

ASV_summary %>% 
  ggplot(., aes(x="", y = count, fill = Category)) +
  geom_col(color = "black", width=1.25) +
  geom_label_repel(data = piechart_ASV,
                   aes(y = pos, label = paste0(count)),
                   size = 3, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_colors) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    panel.border = element_blank(),  # Remove panel border
    panel.grid = element_blank(),    # Remove grid lines
    axis.ticks = element_blank(),    # Remove axis ticks
    axis.text = element_blank(),     # Remove axis text
    #legend.position = "none"
    axis.title = element_blank()     # Remove axis titles
      ) +
  ggtitle("Number of ASVs") +
  xlab("") + ylab("") + labs(fill = "Category")
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Category_breakdown_ASVs.png", width=4, height=3)
```

Number of species pie chart. Used in contract report.

``` r
piechart_spp <- species_summary %>%  
  dplyr::mutate(csum = rev(cumsum(rev(count))), 
         pos = count/2 + lead(csum, 1),
         pos = if_else(is.na(pos), count/2, pos))

piechart_species_plot <- species_summary %>% 
  ggplot(., aes(x="", y = count, fill = Category)) +
  geom_col(color = "black", width=1.25) +
  geom_label_repel(data = piechart_spp,
                   aes(y = pos, label = paste0(count)),
                   size = 3, nudge_x = 1, show.legend = FALSE) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_colors) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"),
    panel.border = element_blank(),  # Remove panel border
    panel.grid = element_blank(),    # Remove grid lines
    axis.ticks = element_blank(),    # Remove axis ticks
    axis.text = element_blank(),     # Remove axis text
    axis.title = element_blank()     # Remove axis titles
  ) +
  ggtitle("Number of Taxonomic Assignments") +
  xlab("") + ylab("") + labs(fill = "Category"); piechart_species_plot
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Category_breakdown_species.png", width=4, height=3)
```

Plot together and export

``` r
plot_grid(piechart_reads, piechart_species_plot,
          ncol=2,
          align = "vh"
          )
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Category_breakdown_contract.png", width=14, height=4)
```

## Generating Heatmaps

Heatmap of all species, subset to Blanks if needed for report. Or
generate heatmap of all species.

reverse label order: scale y discrete limits reverse limits=rev

<https://coolors.co/> (hit tools on the top right hand side)

### General heatmap

``` r
heatmap_df %>% 
  
 # subset(SampleType == "Blank") %>% 
  
  ## replace zeros with NAs for plotting
  replace_with_na_all(condition = ~.x == 0.00000) %>%
  
  # Create a factor for Common_name ordered by Order within each Category
  group_by(Category) %>%
  mutate(Common_name = factor(Common_name, levels = unique(Common_name[order(Order, desc(Common_name))]))) %>%
  ungroup() %>%
  
  ## ggplot basic options (USER EDIT: X AND Y AXIS)
  ggplot(., aes(x = sampleID, y = Common_name)) +
  geom_tile(aes(fill = relab), color = "black") +
  
  ## x, y, and legend labels (USER EDITS IF DESIRED)
  ylab("Common name") +
  xlab("Site") +
  labs(fill = "Relative Abundance (%)") +
  
  ## color of the tile options; direction=1 will flip the low/high (USER EDITS IF DESIRED)
  scale_fill_gradient(na.value = "white", low = "lightskyblue2", high = "#0C4D66") + 
  
  ## facet grid with Category and project variables
  facet_grid2(Category ~ SampleType, 
              scales = "free", space = "free", 
              labeller = labeller(Category = label_wrap_gen(width = 10))) +
  
  ## graph theme options
  theme_classic() +
  theme(
    ## axis text 
    axis.text.x = element_text(angle = 90, size=6, color="grey25", hjust = 1),
    axis.text.y = element_text(colour = 'black', size = 8),
    
    ## legend text and title 
    legend.text = element_text(size = 8, color="black"),
    legend.title = element_text(margin = margin(t = 0, r = 0, b = 5, l = 0), size=10, color="black", face="bold"),
    legend.position = c(-0.4, -0.05), 
    legend.key.height = unit(5, 'mm'),
    legend.direction = "horizontal",
    legend.key.width = unit(5, 'mm'),
    legend.title.align = 0.5,
    legend.title.position = "top",
    
    ## axis titles 
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=14, face="bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=14, face="bold"),
        
    ## facet wrap labels
    strip.text.x = element_text(color = "black", face = "bold", size = 12),
    strip.text.y = element_text(color = "black", face = "bold", size = 12, angle=0),
    strip.background.y = element_blank(),
    strip.clip = "off"
    )
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
## USER EDITS WIDTH AND HEIGHT TO DESIRED   
ggsave("example_output/Figures/Relative_abundance.png", width = 22, height = 14) 
```

### Top Species List visual

Used for contract report.

``` r
top_list <- df_raw_long %>%
  
  ## remove non-target categories
  filter(!Category == "Other" & !Category == "Livestock" & !Category == "Unassigned" & !Category == "Human") %>%
  dplyr::group_by(Species_name, Common_name) %>%
  
  reframe(total = sum(reads),
            log = log10(total),
            total_M = total/1000000) %>% 
  arrange(desc(total)) %>%
  
  ## select the number desired (ie top 30 hits)
  head(30) 

top_list %>% ggplot(., aes(x = fct_reorder(Common_name, log), y = log)) +
  geom_segment(aes(xend = Common_name, yend = 0), color = "#97C1DE") +  # Lollipop stick
  geom_point(size = 3, shape=21, color='grey30', fill = "#97C1DE") +  # Lollipop head
  coord_flip() +  # Flip coordinates for horizontal lollipop chart
  labs(
       x = "",
       y = "Normalized Reads") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8, color='black'), #, face="italic"
        axis.text.x = element_text(size = 6),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=10, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=10, face="bold"),
        axis.text.x.top = element_text(size = 8, color='black', face="italic", angle = 45, hjust = 0)) +
  scale_y_continuous(
    labels = comma,
    limits = c(0, max(top_list$log) + (max(top_list$log)*0.1))  # Set the upper limit to max value + 10%
  )
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Heatmap_Top_species_log.png", width=3.5, height=6)
```

### Bubble plot

Used for contract report.

If you get an error that says something like the below, then change
“(seq(0, 1, length.out = 5))” length.out=X value to reflect the value in
“data has X” within the error. This is the number of size groups in the
legend.

Error in `[[<-.data.frame`(`*tmp*`, i, value = c(“\#0C4D66”, “\#386C86”,
: replacement has 5 rows, data has 2

``` r
df_raw_long %>% 
  ## calculate reads 
  dplyr::group_by(Species_name, Common_name, Category) %>%
  reframe(sum = sum(reads)/1000000,
          xaxis = "x") %>%
  left_join(., df_tax, by = "Species_name") %>%
  
  ## remove non-target categories
  filter(!Category == "Other" & !Category == "Livestock" & !Category == "Unassigned" & !Category == "Human") %>%
  
  # Create a factor for Common_name ordered by Order within each Category
  dplyr::group_by(Category) %>%
  mutate(Common_name = factor(Common_name, levels = unique(Common_name[order(Order, desc(Common_name))]))) %>%
  ungroup() %>%
  
  ggplot(., aes(x=xaxis, y=Common_name)) + 
  
  geom_point(aes(size=sum, fill=sum), color = 'black', shape=21) +
  scale_fill_gradient(na.value = "white", low = "lightskyblue2", high = "#0C4D66") +

  facet_grid2(Category ~ ., scales = "free", space = "free") +
  
  theme_bw() +
  labs(
    x="",
    y="",
    fill = "Reads (M)",
    size = "Reads (M)"
  ) +
  
  theme(
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.text.x = element_blank(),
    axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=12, face="bold"),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=12, face="bold"),
    ## facet wrap labels
    strip.text.x = element_text(color = "black", face = "bold", size = 12),
    #strip.text.y = element_text(color = "black", face = "bold", size = 12, angle=0),
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    strip.clip = "off",
    # Combine legends
    legend.position = "right",
    legend.box = "vertical"
  ) +
  guides(
    fill = "none",
    size = guide_legend(order = 2, reverse = TRUE,
                        override.aes = list(fill = scales::seq_gradient_pal("#0C4D66", "lightskyblue2")
                                            (seq(0, 1, length.out = 2))))
  )
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Heatmap_Species_bubbleplot.png", width=4.5, height=16)
```

## Biodiversity metrics

Creating phyloseq object with the df_relative matrix. Change to df_raw
if needed.

``` r
otu <- otu_table(df_relative %>% 
                   subset(!Category == "Unassigned") %>%
                   dplyr::select(-Common_name, -Category) %>%
                   column_to_rownames(var = "Species_name"), 
                  taxa_are_rows = T)

meta_phyloseq <- sample_data(
  meta %>% ## rownames are also needed in phyloseq meta table
  mutate(sampleID2=sampleID) %>% column_to_rownames(var = "sampleID2")
)

## Merge metadata and OTU table into one phyloseq "object"
phylo_obj <- merge_phyloseq(otu, meta_phyloseq)

## view phyloseq obj 
## expected output = otu_table() with taxa and sample numbers and sample_data() with the sample and column numbers
print(phylo_obj)
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 42 taxa and 142 samples ]
    ## sample_data() Sample Data:       [ 142 samples by 5 sample variables ]

``` r
# Ensure that your OTU table doesn't contain any NA or negative values (output should be FALSE)
any(is.na(otu_table(phylo_obj)))
```

    ## [1] FALSE

``` r
any(otu_table(phylo_obj) < 0)
```

    ## [1] FALSE

### Calculating alpha diversity and species richness

This is currently done with relative abundance values which is not the
best way.. But might be better than raw reads.

``` r
alpha_div <- estimate_richness(phylo_obj, measures = c("Shannon", "Simpson")) %>%
  rownames_to_column(var = "sampleID") %>% left_join(., meta, by = "sampleID") 
```

    ## Warning in estimate_richness(phylo_obj, measures = c("Shannon", "Simpson")): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
## Species Richness
biodiv_df <- df_raw_long %>% subset(!Category == "Unassigned") %>% dplyr::group_by(sampleID) %>%
  filter(reads > 0) %>%
  reframe(Richness = n_distinct(Species_name)) %>%
  
  ## combining this information with alphadiversity 
  left_join(., alpha_div, by = "sampleID") 

biodiv_df %>% 
  write_xlsx("example_output/Biodiversity.xlsx")
```

Plotting

``` r
## plotting
biodiv_df %>%
  ggplot(., aes(x=Richness, y=Shannon)) + 
  
  annotate("rect", xmin=quantile(biodiv_df$Richness, 0.75), xmax=Inf,
             ymin=quantile(biodiv_df$Shannon, 0.75, na.rm = TRUE), ymax=Inf,
             alpha=0.2, fill="#058C42") +

  geom_point(fill = "#97C1DE", color='black', shape=21, alpha=0.5, size=2) + 
  
  labs(
    x = "Sample Richness",
    y = "Sample Evenness (Shannon Index)"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=8, color="grey30"),
        axis.text.x = element_text(size=8, color="grey30"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size=10, face="bold"),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size=10, face="bold")
  )
```

![](04-report_generation_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
ggsave("example_output/Figures/SampleReport_Biodiversity.png", width = 5.5, height = 5)
```
