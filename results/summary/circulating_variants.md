Circulating SARS-CoV-2 RBD variants
================
8/5/2020

  - [Setup](#setup)
  - [Analyzing amino acid diversity in GISAID Spike
    sequences](#analyzing-amino-acid-diversity-in-gisaid-spike-sequences)

This notebook analyzes RBD variants that have been sampled in isolates
within the current SARS-CoV-2 pandemic. It outputs a table of mutants
with the number of times each has been sampled among GISAID sequences,
which can be loaded into other analyses.

## Setup

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra","bio3d","seqinr","knitr")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages])
}
#load packages
lapply(packages, library, character.only=T)
```

    ## [[1]]
    ## [1] "yaml"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
    ## [7] "methods"   "base"     
    ## 
    ## [[2]]
    ## [1] "data.table" "yaml"       "stats"      "graphics"   "grDevices" 
    ## [6] "utils"      "datasets"   "methods"    "base"      
    ## 
    ## [[3]]
    ##  [1] "forcats"    "stringr"    "dplyr"      "purrr"      "readr"     
    ##  [6] "tidyr"      "tibble"     "ggplot2"    "tidyverse"  "data.table"
    ## [11] "yaml"       "stats"      "graphics"   "grDevices"  "utils"     
    ## [16] "datasets"   "methods"    "base"      
    ## 
    ## [[4]]
    ##  [1] "gridExtra"  "forcats"    "stringr"    "dplyr"      "purrr"     
    ##  [6] "readr"      "tidyr"      "tibble"     "ggplot2"    "tidyverse" 
    ## [11] "data.table" "yaml"       "stats"      "graphics"   "grDevices" 
    ## [16] "utils"      "datasets"   "methods"    "base"      
    ## 
    ## [[5]]
    ##  [1] "bio3d"      "gridExtra"  "forcats"    "stringr"    "dplyr"     
    ##  [6] "purrr"      "readr"      "tidyr"      "tibble"     "ggplot2"   
    ## [11] "tidyverse"  "data.table" "yaml"       "stats"      "graphics"  
    ## [16] "grDevices"  "utils"      "datasets"   "methods"    "base"      
    ## 
    ## [[6]]
    ##  [1] "seqinr"     "bio3d"      "gridExtra"  "forcats"    "stringr"   
    ##  [6] "dplyr"      "purrr"      "readr"      "tidyr"      "tibble"    
    ## [11] "ggplot2"    "tidyverse"  "data.table" "yaml"       "stats"     
    ## [16] "graphics"   "grDevices"  "utils"      "datasets"   "methods"   
    ## [21] "base"      
    ## 
    ## [[7]]
    ##  [1] "knitr"      "seqinr"     "bio3d"      "gridExtra"  "forcats"   
    ##  [6] "stringr"    "dplyr"      "purrr"      "readr"      "tidyr"     
    ## [11] "tibble"     "ggplot2"    "tidyverse"  "data.table" "yaml"      
    ## [16] "stats"      "graphics"   "grDevices"  "utils"      "datasets"  
    ## [21] "methods"    "base"

``` r
knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#read in config file
config <- read_yaml("config.yaml")

#read in file giving concordance between RBD numbering and SARS-CoV-2 Spike numbering
RBD_sites <- read.csv(config$RBD_sites,stringsAsFactors = F)

#make output directory
if(!file.exists(config$circulating_variants_dir)){
  dir.create(file.path(config$circulating_variants_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /app/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_haswellp-r0.3.1.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] knitr_1.25        seqinr_3.6-1      bio3d_2.3-4      
    ##  [4] gridExtra_2.3     forcats_0.4.0     stringr_1.4.0    
    ##  [7] dplyr_0.8.3       purrr_0.3.3       readr_1.3.1      
    ## [10] tidyr_1.0.0       tibble_2.1.3      ggplot2_3.2.1    
    ## [13] tidyverse_1.2.1   data.table_1.12.6 yaml_2.2.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.10        haven_2.1.1      lattice_0.20-38 
    ##  [5] colorspace_1.4-1 vctrs_0.2.0      generics_0.0.2   htmltools_0.4.0 
    ##  [9] rlang_0.4.1      pillar_1.4.2     glue_1.3.1       withr_2.1.2     
    ## [13] modelr_0.1.5     readxl_1.3.1     lifecycle_0.1.0  munsell_0.5.0   
    ## [17] gtable_0.3.0     cellranger_1.1.0 rvest_0.3.4      evaluate_0.14   
    ## [21] parallel_3.6.1   broom_0.5.2      Rcpp_1.0.2       scales_1.0.0    
    ## [25] backports_1.1.5  jsonlite_1.6     hms_0.5.2        digest_0.6.22   
    ## [29] stringi_1.4.3    ade4_1.7-13      grid_3.6.1       cli_1.1.0       
    ## [33] tools_3.6.1      magrittr_1.5     lazyeval_0.2.2   crayon_1.3.4    
    ## [37] pkgconfig_2.0.3  zeallot_0.1.0    MASS_7.3-51.4    xml2_1.2.2      
    ## [41] lubridate_1.7.4  assertthat_0.2.1 rmarkdown_1.16   httr_1.4.1      
    ## [45] rstudioapi_0.10  R6_2.4.0         nlme_3.1-141     compiler_3.6.1

Read in table of mutation effects on expression and ACE2 binding

``` r
mutants <- data.table(read.csv(file=config$mut_bind_expr,stringsAsFactors = F))

# #rename mutants site indices to prevent shared names with RBD_sites, simplifying some downstream calculations that cross-index these tables
# setnames(mutants, "site_RBD", "RBD_site");setnames(mutants, "site_SARS2", "SARS2_site")
```

## Analyzing amino acid diversity in GISAID Spike sequences

We constructed an alignment of all Spike sequences available on GISAID
as of 4 August, 2020. On the EpiCoV page, under downloads, one of the
pre-made options is a fasta of all Spike sequences isolated thus far,
which is updated each day. I have downloaded this file, unzipped,
replaced spaces in fasta headers with underscores, and aligned
sequences. We load in this alignment using the `read.fasta` function of
the `bio3d` package, and trim the alignment to RBD residues. We remove
sequecnes from non-human isolates (e.g. bat, pangolin, “environment”,
mink, etc.) and sequences with gap `-` characters, and then iterate
through the alignment and save any observed mutations. We then filter
mutations based on rationale below, and add counts of filtered
observations for each mutation as an ‘nobs’ colum in our overall mutants
data table.

We filter out any mutations that were *only* observed on sequences with
large numbers of missing `X` characters – from my initial pass, I saw
some singleton amino acid variants which would require \>1 nt change,
and these were only found in a single sequence with many X amino acid
characters (the first half of the sequence was well determined, but the
second half was all X’s, with the annotated “differences” being within
short stretches between Xs with determiined amino acids), which made me
realize I needed to be careful not only of sequences rich in gap “-”
characters, but also ambiguous “X” characters. However, I didn’t want to
remove all sequences with undetermined characters off the bat, because
another pattern I saw is that for isolates bearing the N439K mutation,
\>10 are well determined across the entire RBD, but \~80 have many X
characters (in part of the RBD that is *not* near the N439K sequence
call). So, my preference would be to believe a mutation observed in an
X-rich sequence *if the variant in the sequence is observed in at least
one variant that does not contain an X*, but not believe mutations that
are *only* observed in X-rich sequences. (I noticed this issue with
N439K, but this is not the only mutation like this which is observed on
0X sequences at least once but other times on sequences with X
characters.) That is the filtering I therefore do below.

``` r
alignment <- bio3d::read.fasta(file=config$GISAID_alignment, rm.dup=T)
```

    ## [1] " ** Duplicated sequence id's: Spike|hCoV-19/Wuhan/WIV04/2019|2019-12-30|EPI_ISL_402124|Original|hCoV-19^^Hubei|Human|Wuhan_Jinyintan_Hospital|Wuhan_Institute_of_Virology|Wuhan_Institute_of_Virology|China **"

``` r
#remove non-human samples
keep <- grep("Human",alignment$id);  alignment$ali <- alignment$ali[keep,]; alignment$id <- alignment$id[keep]

#remove columns that are gaps in first reference sequence
alignment$ali <- alignment$ali[,alignment$ali[1,]!="-"]

alignment_RBD <- alignment; alignment_RBD$ali <- alignment$ali[,RBD_sites$site_SARS2]

#check that the first sequence entry matches our reference RBD sequence
stopifnot(sum(!(alignment_RBD$ali[1,] == RBD_sites[,"amino_acid_SARS2"]))==0)

#remove sequences have gaps, as the amino acid calls may be generally unreliable
remove <- c()
for(i in 1:nrow(alignment_RBD$ali)){
  if(sum(alignment_RBD$ali[i,]=="-") > 0){remove <- c(remove,i)}
}

alignment_RBD$ali <- alignment_RBD$ali[-remove,];alignment_RBD$id <- alignment_RBD$id[-remove]

#output all mutation differences from the WT/reference RBD sequence
#I do this by iterating over rows and columns of the alignment matrix which is STUPID but effective
variants_vec <- c()
isolates_vec <- c()
for(j in 1:ncol(alignment_RBD$ali)){
  #print(i)
  for(i in 1:nrow(alignment_RBD$ali)){
    if(alignment_RBD$ali[i,j] != alignment_RBD$ali[1,j] & !(alignment_RBD$ali[i,j] %in% c("X","-"))){
      variants_vec <- c(variants_vec, paste(alignment_RBD$ali[1,j],j,alignment_RBD$ali[i,j],sep=""))
      isolates_vec <- c(isolates_vec, alignment_RBD$id[i])
    }
  }
}

#remove any mutations that are *only* observed in X-rich sequences of dubious quality (keep counts in X-rich sequences if they were observed in at least one higher quality isolate)
#make a data frame that gives each observed mutation, the isolate it was observed in, and the number of X characters in that sequence. Also, parse the header to give the country/geographic division of the sample
variants <- data.frame(isolate=isolates_vec,mutation=variants_vec)
for(i in 1:nrow(variants)){
  variants$number_X[i] <- sum(alignment_RBD$ali[which(alignment_RBD$id == variants[i,"isolate"]),]=="X")
  variants$geography[i] <- strsplit(as.character(variants$isolate[i]),split="/")[[1]][2]
}
#filter the sequence set for mutations observed in at least one X=0 background
variants_filtered <- data.frame(mutation=unique(variants[variants$number_X==0,"mutation"])) #only keep variants observed in at least one sequence with 0 X
for(i in 1:nrow(variants_filtered)){
  variants_filtered$n_obs[i] <- sum(variants$mutation == variants_filtered$mutation[i]) #but keep counts for any sequence with observed filtered muts
  variants_filtered$n_geography[i] <- length(unique(variants[variants$mutation == variants_filtered$mutation[i],"geography"]))
  variants_filtered$list_geography[i] <- list(list(unique(variants[variants$mutation == variants_filtered$mutation[i],"geography"])))
}

#add count to mutants df
mutants[,nobs:=0]
mutants[,ngeo:=0]
mutants[,geo_list:=as.list(NA)]
for(i in 1:nrow(mutants)){
  if(mutants$mutation_RBD[i] %in% variants_filtered$mutation){
    mutants$nobs[i] <- variants_filtered[variants_filtered$mutation==mutants$mutation_RBD[i],"n_obs"]
    mutants$ngeo[i] <- variants_filtered[variants_filtered$mutation==mutants$mutation_RBD[i],"n_geography"]
    mutants$geo_list[i] <- variants_filtered[variants_filtered$mutation==mutants$mutation_RBD[i],"list_geography"]
  }
}
```

We see 1658 amino acid polymorphisims within the 74869 sequences
uploaded in GISAID, which represents 184 of our 3819 measured missense
mutants. In the table below, we can see that many of these mutations are
observed only one or a few times, so there may still be unaccounted for
sequencinig artifacts, which we tried to account for at least minimally
with some filtering above.

``` r
kable(table(mutants[mutant!=wildtype & mutant!="*",nobs]),col.names=c("mutation count","frequency"))
```

| mutation count | frequency |
| :------------- | --------: |
| 0              |      3635 |
| 1              |        89 |
| 2              |        29 |
| 3              |        12 |
| 4              |        11 |
| 5              |         4 |
| 6              |         4 |
| 7              |         2 |
| 8              |         6 |
| 9              |         5 |
| 10             |         1 |
| 11             |         3 |
| 13             |         2 |
| 14             |         1 |
| 16             |         2 |
| 18             |         1 |
| 23             |         1 |
| 24             |         1 |
| 27             |         2 |
| 32             |         1 |
| 38             |         1 |
| 40             |         1 |
| 57             |         1 |
| 76             |         1 |
| 104            |         1 |
| 182            |         1 |
| 517            |         1 |

For curiosity’s sake, here are tables giving mutations that were seen
\>20 times, and those seen any number of times with measured binding
effects \>0.1.

| Mutation | expr, lib1 | expr, lib2 | expression effect | bind, lib1 | bind, lib2 | binding effect | number of GISAID sequences | number locations |
| :------- | ---------: | ---------: | ----------------: | ---------: | ---------: | -------------: | -------------------------: | ---------------: |
| V341I    |       0.09 |       0.33 |              0.21 |       0.03 |       0.01 |           0.02 |                         27 |                3 |
| A344S    |     \-0.59 |     \-0.64 |            \-0.62 |     \-0.15 |     \-0.13 |         \-0.14 |                         40 |                8 |
| V367F    |         NA |       0.74 |              0.74 |       0.02 |       0.13 |           0.07 |                         57 |               13 |
| P384L    |     \-0.03 |     \-0.04 |            \-0.03 |     \-0.02 |       0.04 |           0.01 |                         23 |                8 |
| N439K    |     \-0.33 |     \-0.36 |            \-0.35 |       0.11 |     \-0.02 |           0.04 |                        517 |                6 |
| S477N    |       0.02 |       0.10 |              0.06 |       0.02 |       0.09 |           0.06 |                        182 |                5 |
| T478I    |     \-0.14 |     \-0.18 |            \-0.16 |     \-0.05 |     \-0.02 |         \-0.04 |                        104 |                1 |
| P479S    |     \-0.17 |     \-0.23 |            \-0.20 |     \-0.01 |     \-0.05 |         \-0.03 |                         76 |                2 |
| V483A    |       0.01 |       0.17 |              0.09 |       0.00 |     \-0.05 |         \-0.03 |                         38 |                2 |
| G485R    |     \-0.44 |     \-0.64 |            \-0.54 |     \-0.14 |     \-0.23 |         \-0.18 |                         32 |                2 |
| A520S    |     \-0.06 |     \-0.05 |            \-0.05 |     \-0.04 |     \-0.04 |         \-0.04 |                         27 |                8 |
| A522V    |     \-0.02 |     \-0.12 |            \-0.07 |       0.03 |     \-0.09 |         \-0.03 |                         24 |                8 |

| Mutation | expr, lib1 | expr, lib2 | expression effect | bind, lib1 | bind, lib2 | binding effect | number of GISAID sequences | number locations |
| :------- | ---------: | ---------: | ----------------: | ---------: | ---------: | -------------: | -------------------------: | ---------------: |
| Y453F    |     \-0.13 |     \-0.04 |            \-0.08 |       0.19 |       0.31 |           0.25 |                          4 |                2 |
| N501Y    |     \-0.06 |     \-0.22 |            \-0.14 |       0.09 |       0.38 |           0.24 |                         16 |                3 |

What are the sequences on which these stronger affinity-enhancing muts
have been sampled?

``` r
mutants[bind_avg>0.2 & nobs>0,]
```

    ##    site_RBD site_SARS2 wildtype mutant mutation mutation_RBD bind_lib1
    ## 1:      123        453        Y      F    Y453F        Y123F      0.19
    ## 2:      171        501        N      Y    N501Y        N171Y      0.09
    ##    bind_lib2 bind_avg expr_lib1 expr_lib2 expr_avg nobs ngeo geo_list
    ## 1:      0.31     0.25     -0.13     -0.04    -0.08    4    2   <list>
    ## 2:      0.38     0.24     -0.06     -0.22    -0.14   16    3   <list>

``` r
alignment_RBD$id[alignment_RBD$ali[,171]=="Y"]
```

    ##  [1] "Spike|hCoV-19/Australia/VIC2130/2020|2020-06-19|EPI_ISL_480729|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ##  [2] "Spike|hCoV-19/Australia/VIC2151/2020|2020-06-18|EPI_ISL_480765|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ##  [3] "Spike|hCoV-19/Australia/VIC2066/2020|2020-06-12|EPI_ISL_480756|Original|hCoV-19^^Victoria|Human|Microbiological_Diagnostic_Unit_-_Public_Health_Laboratory_(MDU-PHL)|MDU-PHL|Schultz|Australia"
    ##  [4] "Spike|hCoV-19/Australia/VIC2112/2020|2020-06-17|EPI_ISL_480718|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ##  [5] "Spike|hCoV-19/Australia/VIC2089/2020|2020-06-12|EPI_ISL_480701|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ##  [6] "Spike|hCoV-19/Brazil/PE-IAM19/2020|2020-04-07|EPI_ISL_500467|Original|hCoV-19^^Pernambuco|Human|LACEN/PE|WallauLab|Wallau|Brazil"                                                              
    ##  [7] "Spike|hCoV-19/Australia/VIC2126/2020|2020-06-18|EPI_ISL_480735|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ##  [8] "Spike|hCoV-19/USA/NY-NYUMC836/2020|2020-04-21|EPI_ISL_456109|Original|hCoV-19^^New_York|Human|NYU_Langone_Health|Departments_of_Pathology_and_Medicine|Maurano|USA"                            
    ##  [9] "Spike|hCoV-19/Australia/VIC2065/2020|2020-06-15|EPI_ISL_480755|Original|hCoV-19^^Victoria|Human|Microbiological_Diagnostic_Unit_-_Public_Health_Laboratory_(MDU-PHL)|MDU-PHL|Schultz|Australia"
    ## [10] "Spike|hCoV-19/Australia/VIC2020/2020|2020-06-03|EPI_ISL_480662|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ## [11] "Spike|hCoV-19/Australia/VIC2067/2020|2020-06-12|EPI_ISL_480757|Original|hCoV-19^^Victoria|Human|Microbiological_Diagnostic_Unit_-_Public_Health_Laboratory_(MDU-PHL)|MDU-PHL|Schultz|Australia"
    ## [12] "Spike|hCoV-19/Australia/VIC2090/2020|2020-06-14|EPI_ISL_480702|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ## [13] "Spike|hCoV-19/Australia/VIC2160/2020|2020-06-16|EPI_ISL_480760|Original|hCoV-19^^Victoria|Human|Microbiological_Diagnostic_Unit_-_Public_Health_Laboratory_(MDU-PHL)|MDU-PHL|Schultz|Australia"
    ## [14] "Spike|hCoV-19/Australia/VIC2150/2020|2020-06-17|EPI_ISL_480766|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ## [15] "Spike|hCoV-19/Australia/VIC2149/2020|2020-06-17|EPI_ISL_480764|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"
    ## [16] "Spike|hCoV-19/Australia/VIC2102/2020|2020-06-15|EPI_ISL_480710|Original|hCoV-19^^Victoria|Human|Victorian_Infectious_Diseases_Reference_Laboratory_(VIDRL)|VIDRL_and_MDU-PHL|Schultz|Australia"

``` r
alignment_RBD$id[alignment_RBD$ali[,123]=="F"]
```

    ## [1] "Spike|hCoV-19/South_Africa/KRISP-0303/2020|2020-06-25|EPI_ISL_487348|Original|hCoV-19^^KZN|Human|NHLS-IALCH|KRISP|Tulio_de_Oliveira|South_Africa"                                    
    ## [2] "Spike|hCoV-19/Switzerland/190022_544_H2/2020|2020-07-07|EPI_ISL_500887|Original|hCoV-19^^Zug|Human|Viollier_AG|Department_of_Biosystems_Science_and_Engineering|Nadeau|Switzerland"  
    ## [3] "Spike|hCoV-19/Switzerland/190026_546_F10/2020|2020-07-07|EPI_ISL_500890|Original|hCoV-19^^Bern|Human|Viollier_AG|Department_of_Biosystems_Science_and_Engineering|Nadeau|Switzerland"
    ## [4] "Spike|hCoV-19/South_Africa/KRISP-0267/2020|2020-06-17|EPI_ISL_482854|Original|hCoV-19^^KZN|Human|Molecular_Diagnostics_Services_(MDS)|KRISP|de_Oliveira|South_Africa"

Illustrate table showing any circulating mutations found at any our
sites of significant antigenic escape across all mAbs.

``` r
sig_sites <- read.csv(config$significant_escape_sites)

kable(mutants[site_SARS2 %in% sig_sites$site & nobs>0,c("mutation","bind_avg","expr_avg","nobs","ngeo","geo_list")],
      col.names=c("Mutation","binding effect","expression effect","number of GISAID sequences", "number locations","locations"))
```

| Mutation | binding effect | expression effect | number of GISAID sequences | number locations | locations                                                                                     |
| :------- | -------------: | ----------------: | -------------------------: | ---------------: | :-------------------------------------------------------------------------------------------- |
| C361T    |         \-0.61 |            \-1.08 |                          4 |                1 | list(“Panama”)                                                                                |
| N370S    |         \-0.02 |            \-0.17 |                         11 |                2 | list(c(“Portugal”, “England”))                                                                |
| A372S    |         \-0.19 |            \-0.45 |                          1 |                1 | list(“Wales”)                                                                                 |
| A372T    |         \-0.16 |            \-0.28 |                          2 |                1 | list(“England”)                                                                               |
| T376I    |         \-0.07 |            \-0.30 |                          5 |                4 | list(c(“England”, “USA”, “Wales”, “Brazil”))                                                  |
| K378N    |         \-0.02 |            \-0.26 |                          3 |                1 | list(“England”)                                                                               |
| K378R    |           0.08 |              0.16 |                          1 |                1 | list(“Shanghai”)                                                                              |
| V382E    |         \-0.52 |            \-1.79 |                          1 |                1 | list(“France”)                                                                                |
| V382L    |         \-0.05 |            \-0.25 |                          8 |                6 | list(c(“England”, “Bangladesh”, “Senegal”, “Brazil”, “Australia”, “USA”))                     |
| P384L    |           0.01 |            \-0.03 |                         23 |                8 | list(c(“England”, “Scotland”, “Belgium”, “Wales”, “USA”, “South\_Africa”, “Russia”, “India”)) |
| P384S    |         \-0.09 |            \-0.38 |                          8 |                7 | list(c(“USA”, “Portugal”, “Singapore”, “Scotland”, “Netherlands”, “Latvia”, “South\_Africa”)) |
| R408I    |         \-0.09 |            \-0.46 |                         13 |                3 | list(c(“Egypt”, “England”, “India”))                                                          |
| A411D    |         \-0.69 |            \-2.02 |                          1 |                1 | list(“Iran”)                                                                                  |
| A411S    |         \-0.14 |            \-0.73 |                          8 |                3 | list(c(“England”, “USA”, “Wuhan”))                                                            |
| K417N    |         \-0.45 |              0.10 |                          4 |                3 | list(c(“Northern\_Ireland”, “Germany”, “England”))                                            |
| K417R    |         \-0.17 |              0.08 |                          1 |                1 | list(“Israel”)                                                                                |
| A435S    |         \-0.07 |            \-0.36 |                          3 |                3 | list(c(“Netherlands”, “Israel”, “Finland”))                                                   |
| V445I    |         \-0.01 |            \-0.04 |                          2 |                1 | list(“England”)                                                                               |
| G446A    |         \-0.26 |            \-0.26 |                          1 |                1 | list(“England”)                                                                               |
| G446S    |         \-0.20 |            \-0.40 |                          2 |                1 | list(“England”)                                                                               |
| G446V    |         \-0.27 |            \-0.48 |                          8 |                5 | list(c(“Australia”, “Wales”, “England”, “Scotland”, “Finland”))                               |
| N450K    |           0.04 |            \-0.01 |                          1 |                1 | list(“England”)                                                                               |
| L452R    |           0.02 |              0.32 |                          2 |                2 | list(c(“Netherlands”, “Denmark”))                                                             |
| A475V    |         \-0.14 |            \-0.21 |                          9 |                3 | list(c(“USA”, “Australia”, “England”))                                                        |
| E484A    |         \-0.07 |            \-0.23 |                          2 |                2 | list(c(“Spain”, “Northern\_Ireland”))                                                         |
| E484D    |         \-0.38 |            \-0.19 |                          2 |                2 | list(c(“Thailand”, “Germany”))                                                                |
| E484K    |           0.06 |            \-0.10 |                          6 |                4 | list(c(“Spain”, “Sweden”, “England”, “USA”))                                                  |
| E484Q    |           0.03 |            \-0.08 |                         16 |                3 | list(c(“Wales”, “India”, “USA”))                                                              |
| G485R    |         \-0.18 |            \-0.54 |                         32 |                2 | list(c(“Australia”, “China”))                                                                 |
| G485S    |         \-0.20 |            \-0.11 |                          1 |                1 | list(“England”)                                                                               |
| F490L    |         \-0.11 |            \-0.35 |                          4 |                3 | list(c(“Australia”, “USA”, “Singapore”))                                                      |
| F490S    |           0.00 |            \-0.10 |                          8 |                1 | list(“England”)                                                                               |
| S494L    |         \-0.35 |            \-1.02 |                          2 |                2 | list(c(“England”, “Switzerland”))                                                             |
| S494P    |           0.00 |            \-0.02 |                         13 |                8 | list(c(“Wales”, “USA”, “Singapore”, “India”, “England”, “Scotland”, “Sweden”, “Spain”))       |
| P499H    |         \-0.47 |            \-0.51 |                          1 |                1 | list(“Egypt”)                                                                                 |
| P499S    |         \-0.23 |            \-0.24 |                          1 |                1 | list(“England”)                                                                               |

Output table giving nobs and geographic spread for circulating mutants.

``` r
write.csv(mutants[,.(site_RBD,site_SARS2,wildtype,mutant,mutation,nobs,ngeo)],file=config$circulating_variants,row.names=F)
```
