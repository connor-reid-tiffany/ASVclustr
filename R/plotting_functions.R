#' Bar plot of cluster means.
#'
#' @description Creates a barplot of cluster means. Useful for checking if cluster patterns are unique and if dispersion is low. Can add ggplot2 aesthetic layers.
#' @param melted_asv_list A data.frame. The output of melt_asv_list.
#' @param time_var A string. The metadata variable for timepoint.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr n
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 position_dodge
#' @importFrom stats sd
#' @return A ggplot2 object
#' @export
#' @examples
#' cluster_means_bar(melted_asv_list, time_var = "Timepoint")
#'
cluster_means_bar <- function(melted_asv_list, time_var){

  #establish objects
  timepoint <- cluster <- Abundance <- se <- NULL

  #remove NA values, rename timepoint variable
  means <- melted_asv_list[!is.na(melted_asv_list[,"Abundance"]),]
  colnames(means)[colnames(means) == time_var] <- "timepoint"

  #calculate means by cluster and timepoint
  means <- means %>%
    group_by(timepoint, cluster) %>%
    summarise(mean = mean(Abundance), n = n(), sd = sd(Abundance), se= sd(Abundance)/sqrt(n()))

  #sort timepoint in ascending order
  sorted_timepoint <- paste(sort(as.integer(levels(as.factor(means$timepoint)))))
  means$timepoint <- factor(means$timepoint, levels = sorted_timepoint)

  #plot
  plot <- ggplot(means, aes(x=timepoint, y=mean, fill=as.character(cluster))) +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                  width=.2,                    #Width of the error bars
                  position=position_dodge(.9)) + facet_wrap(.~cluster, scales = "free")

  return(plot)

}



#' Plot a cluster phytree.
#'
#' @description Wrapper function for ggtree that produces a tree with branches and tips colored by cluster.
#' can use ggtree aesthetic layers.
#' @param asv_list An asv_list. Must have an h_clust element.
#' @param tree A tree of phylo class produced with the seqmat given to the asv_list constructor function.
#' @param cladogram A boolean. TRUE produces a cladogram where branch lengths are equal, FALSE produces a
#' phylogenetic tree. Default is FALSE.
#' @param layout A string. See ?ggtree::ggtree for layout options. Default is "circular".
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tree
#' @importFrom ggtree geom_tippoint
#' @importFrom dplyr left_join
#' @importFrom methods isClass
#' @return A ggtree object
#' @export
#' @examples
#' cluster_phytree(asv_list = asv_list, tree = phangtree$tree)
#'
#' cluster_phytree(asv_list = asv_list, tree = phangtree$tree, cladogram = TRUE, layout = "radial")
#'
cluster_phytree <- function(asv_list, tree, cladogram=FALSE, layout="circular"){

  #establish objects
  cluster <- NULL

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be a list object with class ASVclustr")
  }

  if (is.null(asv_list$h_clust)){
    stop("asv_list must have an h_clust element to use function cluster_phytree")
  }

  #create ggtree object
  if (cladogram==FALSE){

    ggtree_ob <- ggtree(tree, layout = layout)

  } else if (cladogram==TRUE){

    ggtree_ob <- ggtree(tree, layout = layout, branch.length = "none")

  }

  #join Cluster data into ggtree data element by ASV value to color by cluster
  colnames(ggtree_ob$data)[colnames(ggtree_ob$data) == "label"] <- "OTU"
  ggtree_ob$data <- left_join(ggtree_ob$data, asv_list$h_clust, "OTU")
  ggtree_ob$data$cluster <- as.character(ggtree_ob$data$cluster)

  #plot the tree
  ggtree_ob + geom_tippoint(aes(color = cluster), size = 2) + geom_tree(aes(color = cluster))

}


#' Investigate taxa abundance by cluster.
#'
#' @description Aggregates taxa by sample, cluster, and timepoint by a specified taxonomic level.
#' @param melted_asv_list A data.frame. Produced from the function melt_asv_list.
#' @param tax_level A string. The taxonomic rank you want to aggregate by.
#' @param time_var A string. The timepoint variable in meta.
#' @param sam_var A string. The sample variable in meta.
#' @param plot_boxplot Optional. A boolean, defualt is FALSE. If FALSE,
#' returns a data.frame. Otherwise returns a ggplot2 boxplot
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_boxplot
#' @return A data.frame if plot_boxplot = FALSE, or a boxplot if plot_boxplot = TRUE
#' @export
#' @examples
#'
#' taxa_by_cluster(melted_asv_list = melted_asv_list, tax_level = Class, time_var = Timepoint,
#' sam_var = Sample)
#'
#' taxa_by_cluster(melted_asv_list = melted_asv_list, tax_level = Class, time_var = Timepoint,
#' sam_var = Sample, plot_boxplot = TRUE)
#'
#'
taxa_by_cluster <- function(melted_asv_list, tax_level, time_var, sam_var, plot_boxplot = FALSE){

  #rename variables of interest
  colnames(melted_asv_list)[colnames(melted_asv_list) == tax_level] <- "tax_level"
  colnames(melted_asv_list)[colnames(melted_asv_list) == sam_var] <- "sample"
  colnames(melted_asv_list)[colnames(melted_asv_list) == time_var] <- "timepoint"

  #group data.frame by variables, then sum abundance by these variables
  melted_asv_list <- melted_asv_list %>%

    group_by(tax_level, sample, timepoint, cluster) %>%

    summarise(taxa_cluster_sum = sum(Abundance))

  #return a new dataframe
  if (plot_boxplot == FALSE){

    return(melted_asv_list)

    #or return a plot of choice if plot_geom parameter is used
  } else if (plot_boxplot == TRUE){

    ggplot(data = melted_asv_list,
           aes(x = timepoint, y = taxa_cluster_sum, fill = tax_level, color = tax_level)) +
      geom_boxplot(alpha = 1/4, size = 1, outlier.alpha = 1) + facet_wrap(.~cluster, scales = "free")

  }

}
