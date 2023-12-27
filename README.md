
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
Finally, an example ordination is performed on the average Bray-Curtis
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
#> 1 F21Plmr       60    3.10   0.933      15.0    8.48
#> 2 M11Fcsw       78    2.86   0.894       9.47  11.4 
#> 3 M11Plmr       83    2.59   0.849       6.64  12.4 
#> 4 M31Fcsw       80    3.19   0.906      10.6   12.0 
#> 5 M31Plmr       60    2.88   0.878       8.19   8.40
#> 6 TS28          74    3.55   0.957      23.1   10.9 
#> 7 TS29          72    2.85   0.897       9.67  10.6 
#> 
#> $bray_curtis
#>         M31Fcsw M11Fcsw M31Plmr M11Plmr F21Plmr    TS28    TS29
#> M31Fcsw 0.00000 0.52195 0.99320 0.98680 0.99315 0.63370 0.78095
#> M11Fcsw 0.52195 0.00000 0.99695 0.99050 0.99645 0.78425 0.88810
#> M31Plmr 0.99320 0.99695 0.00000 0.81775 0.37625 0.99675 0.99655
#> M11Plmr 0.98680 0.99050 0.81775 0.00000 0.72965 0.99060 0.99075
#> F21Plmr 0.99315 0.99645 0.37625 0.72965 0.00000 0.99645 0.99610
#> TS28    0.63370 0.78425 0.99675 0.99060 0.99645 0.00000 0.61455
#> TS29    0.78095 0.88810 0.99655 0.99075 0.99610 0.61455 0.00000
#> 
#> $jaccard
#>           M31Fcsw   M11Fcsw   M31Plmr   M11Plmr   F21Plmr      TS28      TS29
#> M31Fcsw 0.0000000 0.2826087 0.7826087 0.6079710 0.7835082 0.3334876 0.3458976
#> M11Fcsw 0.2826087 0.0000000 0.8189117 0.6340580 0.8046668 0.3739130 0.3873728
#> M31Plmr 0.7826087 0.8189117 0.0000000 0.4150908 0.3629880 0.8089839 0.8214143
#> M11Plmr 0.6079710 0.6340580 0.4150908 0.0000000 0.3751776 0.6244851 0.6537481
#> F21Plmr 0.7835082 0.8046668 0.3629880 0.3751776 0.0000000 0.8042911 0.8171490
#> TS28    0.3334876 0.3739130 0.8089839 0.6244851 0.8042911 0.0000000 0.2178131
#> TS29    0.3458976 0.3873728 0.8214143 0.6537481 0.8171490 0.2178131 0.0000000
#> 
#> $horn
#>           M31Fcsw   M11Fcsw   M31Plmr   M11Plmr   F21Plmr      TS28      TS29
#> M31Fcsw 0.0000000 0.2464569 0.9976838 0.9980811 0.9972698 0.4343585 0.8904716
#> M11Fcsw 0.2464569 0.0000000 0.9983735 0.9991068 0.9985446 0.6036073 0.9427666
#> M31Plmr 0.9976838 0.9983735 0.0000000 0.8918437 0.2942485 0.9988168 0.9993138
#> M11Plmr 0.9980811 0.9991068 0.8918437 0.0000000 0.7889263 0.9981794 0.9986362
#> F21Plmr 0.9972698 0.9985446 0.2942485 0.7889263 0.0000000 0.9978479 0.9988216
#> TS28    0.4343585 0.6036073 0.9988168 0.9981794 0.9978479 0.0000000 0.5842069
#> TS29    0.8904716 0.9427666 0.9993138 0.9986362 0.9988216 0.5842069 0.0000000
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
