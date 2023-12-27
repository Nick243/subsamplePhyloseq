#' subsample_phyloseq
#'
#' This function performs repeated subsampling of an existing phyloseq object and returns the average across subsamples for common alpha- and beta-diversity metrics. Subsampling is performed without replacement and parallel computations supported. A list of subsampled phyloseq objects is returned allowing for additional metrics to be computed and averaged across the subsamples as required.
#'
#' @param ps_object An existing phyloseq object to repeatedly subsample.
#' @param repeats Number of repeats (defaults to n=100).
#' @param depth Library size that each sample will be subsampled to. Sampling is performed without replacement and a value required to encourage reporting of the depth (defaults to 50k).
#' @param n_cores Number of cores to use for parallel computations (default n=1).
#' @param ...
#'
#' @return A list with the following elements:
#' @return ps_object: The input ps_object (i.e., unsubsampled object)
#' @return repeats: The number of repeated subsamples
#' @return subsample_depth: The library size selected for subsampling
#' @return adiv: Alpha-diversity metrics averaged across subsamples
#' @return bray_curtis: Bray-Curtis dissimilarity averaged across subsamples
#' @return jaccard: Jaccard index averaged across subsamples
#' @return horn: Horn index averaged across subsamples
#' @return subsampled_list: returns a list of the subsampled phyloseq objects
#'
#' @export
#'
#' @examples
subsample_phyloseq <- function(ps_object, repeats = 100, depth = 50000, n_cores = 1, ...){

  message("subsample_phyloseq() assumes taxa comprise the rows of the OTU table.
          If the function returns a distance matrix of taxa, and you wanted pairwise distances
          for each sample, then transpose the OTU table before running.")

  ps <- ps_object

  #Perform repeated subsampling
  my_cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
  doParallel::registerDoParallel(cl = my_cluster)

  '%dopar%' <- foreach::'%dopar%'
  my_list <- vector(mode = "list", length = repeats)
  my_list <- foreach::foreach(i = 1:repeats, .packages = "phyloseq") %dopar% {
    my_list[[i]] <- phyloseq::rarefy_even_depth(ps, sample.size = depth, replace = F, trimOTUs = T, verbose = F)
  }
  parallel::stopCluster(my_cluster)


  #Obtain mean alpha diversity estimates
  adiv_list <- lapply(my_list, phyloseq::estimate_richness, measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher"))
  adiv_df <- do.call(rbind.data.frame, adiv_list)

  adiv_mu_df <- adiv_df |>
    dplyr::mutate(id = rep(phyloseq::sample_names(ps), repeats)) |>
    dplyr::group_by(id) |>
    dplyr::summarise(Observed = round(mean(Observed), 0),
                     Shannon = mean(Shannon),
                     Simpson = mean(Simpson),
                     InvSimpson = mean(InvSimpson),
                     Fisher = mean(Fisher))


  #Obtain mean beta diversity estimates
  '%>%' <- magrittr::'%>%'
  calc_dist <- function(i, my_method, my_bin){
    my_otu <- data.frame(phyloseq::otu_table(i)) %>%
      t(.) %>%
      data.frame(.)

    my_dist <- vegan::vegdist(my_otu, method = my_method, binary = my_bin, diag = T, upper = T, na.rm = FALSE)
    my_mat <- as.matrix(my_dist)

    return(my_mat)
  }

  get_dist_means <- function(dist_list){
    dist_bind <- do.call(cbind, dist_list)

    dist_array <- array(dist_bind, dim=c(dim(dist_list[[1]]), length(dist_list)))
    mu_dist <- apply(dist_array, c(1, 2), mean, na.rm = TRUE)

    colnames(mu_dist) <- phyloseq::sample_names(ps)
    rownames(mu_dist) <- phyloseq::sample_names(ps)

    return(mu_dist)
  }

  my_bc_list <- lapply(my_list, calc_dist, my_method = "bray", my_bin = F)
  my_jd_list <- lapply(my_list, calc_dist, my_method = "jaccard", my_bin = T)
  my_hn_list <- lapply(my_list, calc_dist, my_method = "horn", my_bin = F)

  bc_mu_df <- get_dist_means(my_bc_list)
  jd_mu_df <- get_dist_means(my_jd_list)
  hn_mu_df <- get_dist_means(my_hn_list)


  #Creating output list
  out_list <- list(ps_obj = ps_object,
                   repeats = repeats,
                   subsample_depth = depth,
                   adiv = adiv_mu_df,
                   bray_curtis = bc_mu_df,
                   jaccard = jd_mu_df,
                   horn = hn_mu_df,
                   subsampled_list = my_list)

  return(out_list)
}
