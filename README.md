
<!-- README.md is generated from README.Rmd. Please edit that file -->

# subsamplePhyloseq

<!-- badges: start -->
<!-- badges: end -->

This package allows for the **repeated subsampling** of an existing
phyloseq object as described in [Schloss
2023](https://www.biorxiv.org/content/10.1101/2023.06.23.546312v1) and
the averaging of common alpha- and beta-diversity metrics across the
subsamples (i.e., perform rarefaction). The repeated
subsampling/rarefaction provides an approach to control for uneven
sequencing effort for the alpha- and beta-diversity metrics returned. A
list of subsampled phyloseq objects is also returned allowing for
additional metrics to be computed and averaged across the subsamples as
required. Parallel computation is implemented using the doParallel and
foreach packages.

## Installation

You can install the development version of subsamplePhyloseq from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Nick243/subsamplePhyloseq")
```

## Example

This basic example highlights the core functionality of the
subsamplePhyloseq package. First, an example phyloseq object is obtained
from the phyloseq package. Then the subsample_phyloseq() function is run
to perform the repeated subsampling and averaging of diversity metrics.
Finally, an example ordination is performed on the average Bray_Curtis
dissimilarity.

Larger values for the depth of subsampling and number of repeats should
be used when analyzing data from an actual experiment.

Run parallel::detectCores() to determine the number of available cores.
For projects with a large number of samples or library sizes consider
using using parallel::detectCores() - 1.

``` r
library(subsamplePhyloseq)
library(phyloseq)
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.1.2

#Create example data
data("GlobalPatterns")
ps <- phyloseq::subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Skin"))
ps <- phyloseq::subset_taxa(ps, phyloseq::taxa_sums(ps) > 0)

x <- phyloseq::taxa_sums(ps)
keepTaxa <- (x / sum(x)) > 1e-3
ps <- phyloseq::prune_taxa(keepTaxa, ps)


#Repeatedly subsample phyloseq object
rep_ps <- subsample_phyloseq(ps_object = ps, repeats = 2, depth = 10000, n_cores = 2)
#> subsample_phyloseq() assumes taxa comprise the rows of the OTU table.
#>           If the function returns a distance matrix of taxa, and you wanted pairwise distances
#>           for each sample, then transpose the OTU table before running.
rep_ps
#> $ps_obj
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 131 taxa and 7 samples ]
#> sample_data() Sample Data:       [ 7 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 131 taxa by 7 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 131 tips and 130 internal nodes ]
#> 
#> $repeats
#> [1] 2
#> 
#> $subsample_depth
#> [1] 10000
#> 
#> $adiv
#> # A tibble: 7 × 6
#>   id      Observed Shannon Simpson InvSimpson Fisher
#>   <chr>      <dbl>   <dbl>   <dbl>      <dbl>  <dbl>
#> 1 F21Plmr       55    3.09   0.932      14.8    7.67
#> 2 M11Fcsw       78    2.87   0.897       9.75  11.6 
#> 3 M11Plmr       84    2.62   0.852       6.74  12.7 
#> 4 M31Fcsw       82    3.19   0.907      10.8   12.2 
#> 5 M31Plmr       57    2.87   0.875       8.01   7.99
#> 6 TS28          80    3.54   0.956      22.6   11.8 
#> 7 TS29          71    2.87   0.898       9.83  10.3 
#> 
#> $bray_curtis
#>         M31Fcsw M11Fcsw M31Plmr M11Plmr F21Plmr    TS28    TS29
#> M31Fcsw 0.00000 0.53055 0.99365 0.98545 0.99390 0.62490 0.77225
#> M11Fcsw 0.53055 0.00000 0.99750 0.99020 0.99755 0.78465 0.88605
#> M31Plmr 0.99365 0.99750 0.00000 0.80835 0.36615 0.99700 0.99645
#> M11Plmr 0.98545 0.99020 0.80835 0.00000 0.72415 0.98970 0.99000
#> F21Plmr 0.99390 0.99755 0.36615 0.72415 0.00000 0.99715 0.99655
#> TS28    0.62490 0.78465 0.99700 0.98970 0.99715 0.00000 0.60170
#> TS29    0.77225 0.88605 0.99645 0.99000 0.99655 0.60170 0.00000
#> 
#> $jaccard
#>           M31Fcsw   M11Fcsw   M31Plmr   M11Plmr   F21Plmr      TS28      TS29
#> M31Fcsw 0.0000000 0.2921196 0.7860031 0.5950719 0.8138306 0.2810236 0.3094628
#> M11Fcsw 0.2921196 0.0000000 0.8318966 0.6293000 0.8589275 0.3955670 0.3659498
#> M31Plmr 0.7860031 0.8318966 0.0000000 0.4278923 0.3284314 0.7918919 0.8416052
#> M11Plmr 0.5950719 0.6293000 0.4278923 0.0000000 0.4145995 0.6041667 0.6712094
#> F21Plmr 0.8138306 0.8589275 0.3284314 0.4145995 0.0000000 0.8201754 0.8644962
#> TS28    0.2810236 0.3955670 0.7918919 0.6041667 0.8201754 0.0000000 0.2600241
#> TS29    0.3094628 0.3659498 0.8416052 0.6712094 0.8644962 0.2600241 0.0000000
#> 
#> $horn
#>           M31Fcsw   M11Fcsw   M31Plmr   M11Plmr   F21Plmr      TS28      TS29
#> M31Fcsw 0.0000000 0.2623441 0.9981900 0.9971549 0.9979632 0.4238101 0.8862384
#> M11Fcsw 0.2623441 0.0000000 0.9991828 0.9981622 0.9992998 0.6073205 0.9368171
#> M31Plmr 0.9981900 0.9991828 0.0000000 0.8829716 0.2869414 0.9988917 0.9992453
#> M11Plmr 0.9971549 0.9981622 0.8829716 0.0000000 0.7791160 0.9972771 0.9982294
#> F21Plmr 0.9979632 0.9992998 0.2869414 0.7791160 0.0000000 0.9986796 0.9988284
#> TS28    0.4238101 0.6073205 0.9988917 0.9972771 0.9986796 0.0000000 0.5692798
#> TS29    0.8862384 0.9368171 0.9992453 0.9982294 0.9988284 0.5692798 0.0000000
#> 
#> $subsampled_list
#> $subsampled_list[[1]]
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 131 taxa and 7 samples ]
#> sample_data() Sample Data:       [ 7 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 131 taxa by 7 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 131 tips and 130 internal nodes ]
#> 
#> $subsampled_list[[2]]
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 131 taxa and 7 samples ]
#> sample_data() Sample Data:       [ 7 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 131 taxa by 7 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 131 tips and 130 internal nodes ]

#Generate example ordination
my_pcoa <- cmdscale(rep_ps$bray_curtis, eig = TRUE)
lab1 <- paste0("Axis.1  [", round(100 * my_pcoa$eig[1]/sum(my_pcoa$eig), 1), "%]")
lab2 <- paste0("Axis.2  [", round(100 * my_pcoa$eig[2]/sum(my_pcoa$eig), 1), "%]")

plot_ordination(ps, my_pcoa, type = "samples", color = "SampleType") +
    labs(x = lab1, y = lab2) +
    geom_point(size = 3)
```

<img src="man/figures/README-example-1.png" width="100%" />
