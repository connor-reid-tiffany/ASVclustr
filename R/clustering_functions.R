#' Hierarchical clustering of ASVs.
#'
#' @description Calculates the cosine similarity, converts to an angular distance matrix, and then either
#' returns a dendrogram if no k is set, or returns an asv_list with an h_clust element.
#'
#' @param asv_list An asv_list.
#' @param agg_method A string of either "ward", "ward.D2", "complete", "single", "average", "mcquitty", "median", or "centroid".
#' "ward.D2" is the default. See ?hclust for more details.
#' @param k Optional. An integer. The number of clusters you wish to group your ASVs into. If k is missing asv_clustr will return
#' a dendrogram plot instead
#' @param td A string of either "none", "both", or "td_only". Default is "none":
#'           "none" will use only abundance data to perform cosine similarity calculation.
#'           "both" will use both the abundance data and time rate derivative data.
#'           "td_only" will use only the time rate derivative data.
#' @importFrom wordspace dist.matrix
#' @importFrom stats hclust
#' @importFrom methods isClass
#' @importFrom stats as.dist
#' @importFrom stats cutree
#' @importFrom stats filter
#' @return An asv_list object, or a dendrogram if no k is set
#' @export
#' @examples
#' clustr_dendro  <- asv_clustr(asv_list = asv_list)
#' plot(clustr_dendro, labels = FALSE)
#' rect.hclust(tree = clustr_dendro, border ="red", k=4)
#'
#' asv_list <- asv_clustr(asv_list = asv_list,  k=4)
#'
#' asv_list <- asv_clustr(asv_list = asv_list,  k=4, td = "both")
#'
#' asv_list <- asv_clustr(asv_list = asv_list,  k=4, td = "td_only")
#'
#'
asv_clustr <- function(asv_list, agg_method = "ward.D2", td = "none", k){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be of classes list and ASVclustr")
  }

  if (td == "none"){
    seqmat <- asv_list$seqmat

  } else if (td == "td_only"){

    if (is.null(asv_list$seqmat_td)){
      stop("asv_list must have a seqmat_td element to use parameter td = 'td_only'")
      }

    seqmat <- asv_list$seqmat_td

  } else if (td == "both"){

    if (is.null(asv_list$seqmat_td)){
      stop("asv_list must have a seqmat_td element to use parameter td = 'both'")
    }

    seqmat <- rbind(asv_list$seqmat, asv_list$seqmat_td)
  }

  #compute cosine similarity and angular distance
  cosine_sim <- dist.matrix(M = seqmat, convert = TRUE, byrow = FALSE)

  #perform hierarchical clustering
  if (missing(k)){
    h_clust <- hclust(d = as.dist(cosine_sim),
                      method = agg_method)
    return(h_clust)
  } else if (!missing(k)){
    h_clust <- as.data.frame(cutree(hclust(d = as.dist(cosine_sim),
                                           method = agg_method),k = k))
    #format cluster dataframe
    colnames(h_clust)[1] <- "cluster"
    h_clust$OTU <- rownames(h_clust)

    asv_list$h_clust <- h_clust
    return(asv_list)

  }
}




#' Calculation of time rate derivatives.
#'
#' @description Computes the time rate derivatives for every ASV in every sample and adds them as an element in your asv_list object asv_list.
#'
#' @export
#' @param asv_list An asv_list.
#' @param sam_var A string. The metadata variable for samples.
#' @param time_var A string. The metadata variable for timepoints.
#' @importFrom methods isClass
#' @return An asv_list.
#' @examples
#'
#' asv_list <- calc_td(asv_list = asv_list,  sam_var = "Sample", time_var = 'Timepoint')
#'
calc_td <- function(asv_list, sam_var, time_var){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be an object of class list and ASVclustr")
  }
  #define variables
  seqmat_td <- as.data.frame(asv_list$seqmat)
  meta <- asv_list$meta

  seqmat_td_list <- order_seqmat(seqmat = seqmat_td, meta = meta, sam_var = sam_var, time_var = time_var)

  seqmat_td_list <- lapply(seqmat_td_list, function(x)
    x[!(names(x) %in% c("SampleID","sample"))])

  #This function will compute the time derivative for every sequential time interval for each
  #ASV and will be iterated across each data frame in the list using lapply
  #f'(x) ?????? (f(x+h) - f(x))/h
  #this function can take samples with unequal length time vectors
  #diff iterates the difference down the rows within a column

  compute_derivatives <- function(x, time_var){
    #coerce the dataframe to a matrix within diff to index
    #all ASV columns
    d <- diff(as.matrix(x[!(names(x) %in% c(time_var))]))/diff(x[,ncol(x)])
    return(d)

  }
  #iterate compute derivatives across each data frame in the list.
  seqmat_derivatives_list <- lapply(seqmat_td_list, compute_derivatives,
                                    time_var = time_var)
  #change the rownames in seqtab_derivatives_list
  seqmat_derivatives_list <- lapply(seqmat_derivatives_list,
                                    function(x){ rownames(x) <- paste0(rownames(x), ".","TD"); x})
  #rbind the two lists together, element by element

  seqmat_derivatives <- do.call("rbind", seqmat_derivatives_list)

  #remove timepoint_var column in each dataframe of seqtab_order_list

  seqmat_derivatives <- seqmat_derivatives[,colnames(seqmat_derivatives)!= "timepoint"]

  #convert na values to 0
  conv_na_zero <- function(x){
    x[is.na(x)] <- 0
    return(x)
  }
  seqmat_derivatives <- apply(seqmat_derivatives, 2, conv_na_zero)
  asv_list$seqmat_td <- seqmat_derivatives

  return(asv_list)
}
#' Compare population dynamics between groups.
#'
#' @description Compares mean population dynamics of either clusters to determine if an experimental variable or batch introduces noise,
#' or individual ASVs to determine if an independent variable effects population dynamics.
#'
#' @param asv_list An asv_list.
#' @param sam_var A string. The metadata variable for samples.
#' @param time_var A string. The metadata variable for timepoints.
#' @param independent_var A string. The metadata variable for the experimental variable(s).
#' @param batch Optional. A string. A metadata variable for a batch effect. The batch effect must be separate
#' from the independent variable, i.e. levels of the batch variable must be distributed among the levels of the idependent variable.
#' @param by_asv A boolean. If TRUE, difference in integral means will be computed for every ASV. If FALSE
#' difference in integral means will be computed for each cluster. Default is FALSE
#' @param rescale A boolean. If TRUE, each ASV for each sample will be rescaled to a maximum value of 1 and
#' minimum value of 0. Otherwise raw read counts will be used. Default is FALSE
#' @importFrom dplyr left_join
#' @importFrom dplyr distinct
#' @importFrom reshape2 melt
#' @importFrom coin kruskal_test
#' @importFrom coin pvalue
#' @importFrom stats p.adjust
#' @importFrom methods isClass
#' @return An asv_list.
#' @export
#' @examples
#'
#' \dontshow{
#' new_meta <- asv_list$meta
#'
#' new_meta$batch <- c("A","B","C","D")
#'
#' asv_list$meta <- new_meta
#'
#' }
#' pd_cage_effect <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint",independent_var = "Treatment", batch="batch")
#'
#' pd_cage_effect_rescale <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", batch="batch", rescale = TRUE)
#'
#' pd_cage_effect_ASV <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", batch="batch", by_asv = TRUE)
#'
#' pd_cage_effect_ASV_rescale <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint",independent_var = "Treatment",
#' batch="batch", by_asv = TRUE, rescale = TRUE)
#'
#' pd_Treatment <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment")
#'
#' pd_Treatment_rescale <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", rescale = TRUE)
#'
#' pd_Treatment_ASV <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", by_asv = TRUE)
#'
#' pd_Treatment_ASV_rescale <- compare_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", by_asv = TRUE, rescale = TRUE)
#'
#'
compare_pd <- function(asv_list, sam_var,time_var,independent_var,batch, by_asv=FALSE, rescale=FALSE){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be an object of class list and ASVclustr")
  }

  #define variables
  seqmat <- as.data.frame(asv_list$seqmat)

  meta <- asv_list$meta

  seqmat$SampleID <- rownames(seqmat)

  #create a list of dataframes where each dataframe is one sample, and order them by timepoint
  seqmat_order_list <- order_seqmat(seqmat = seqmat, meta = meta, sam_var = sam_var, time_var = time_var)

  #function to extract timepoint column
  extract_time_var <- function(x){
    x <- x[,"timepoint"]
    return(x)
  }

  time_var_order_list <- lapply(seqmat_order_list, extract_time_var)

  #remove all metadata variables from each dataframe in the list except for the timepoint variable
  seqmat_order_list <- lapply(seqmat_order_list, function(x)
    x[!(names(x) %in% c("SampleID","sample", "timepoint"))])

  #normalize ASVs to fit within the same scale (removes influence of abundance differences)
  if (rescale==TRUE){

    rescale <- function(x){
      x <- (as.matrix(x)-min(as.matrix(x)))/(max(as.matrix(x))-min(as.matrix(x)))
      return(x)
    }

  #rescale the data by iteration
    seqmat_order_list <- lapply(seqmat_order_list, function(x) sapply(x,  rescale))
    #convert any NA values to 0

    conv_na_zero <- function(x){
      x[is.na(x)] <- 0
      return(x)
    }

    #recombine individual sample dataframes back into one dataframe
    seqmat_order_list <- lapply(seqmat_order_list, conv_na_zero)
    seqmat_order_list <- Map(cbind, seqmat_order_list, time_var_order_list)
    seqmat_order_list <- lapply(seqmat_order_list, as.data.frame)

    #the timepoint column lost its name value and must be renamed
    rn_last_col <- function(x){
      colnames(x)[ncol(x)] <- "timepoint"
      return(x)
    }

    seqmat_order_list <- lapply(seqmat_order_list, rn_last_col)

  } else if (rescale==FALSE){
    seqmat_order_list <- Map(cbind, seqmat_order_list, time_var_order_list)
    seqmat_order_list <- lapply(seqmat_order_list, as.data.frame)

    rn_last_col <- function(x){
      colnames(x)[ncol(x)] <- "timepoint"
      return(x)
    }

    seqmat_order_list <- lapply(seqmat_order_list, rn_last_col)
    seqmat_order_list <- seqmat_order_list
  }
  #compute time rate derivatives and also return the timepoint column for integration
   compute_derivatives <- function(x){

    d <- cbind(diff(as.matrix(x[!(names(x) %in% c("timepoint"))]))/diff(x[,ncol(x)]), "time_interval" = diff(x[,ncol(x)]))
    return(d)
  }

  #iterate over list to compute all time rate derivatives and time intervals
  seqmat_derivatives_list <- lapply(seqmat_order_list, compute_derivatives)
  seqmat_derivatives_list <- lapply(seqmat_derivatives_list, as.data.frame)
  #function to compute the integrals using the euclidean norm
  compute_integrals <- function(x, time_interval){

    x <- as.matrix(x[,!(names(x) %in% c(time_interval))])^2 * x[,ncol(x)]
    x <- colSums(as.matrix(x))^(1/2)
    return(x)
  }
  #iterate over list of derivatives that have  a time_interval vector!
  seqmat_integrals_list <- lapply(seqmat_derivatives_list, compute_integrals,
                                  time_interval = "time_interval")
  #column bind into one dataframe, make an OTU  column that can be matched
  #to cluster info
  seqmat_integrals <- as.data.frame(do.call("cbind", seqmat_integrals_list))
  seqmat_integrals$OTU <- rownames(seqmat_integrals)
  #Join cluster data
  seqmat_integrals <- left_join(seqmat_integrals, asv_list$h_clust,by = "OTU")
  seqmat_integrals <- melt(seqmat_integrals,  c("OTU", "cluster"))

  #test for significant difference in variance of ASV curves in each cluster
  #between independent variables
  if (missing(batch)){
    #remove unecessary vectors from meta, rename variables for modelling, merge with seqmat for modeling
    meta_2 <- meta[,names(meta) %in% c(sam_var,independent_var)]
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == "variable"] <- sam_var
    seqmat_integrals<- left_join(seqmat_integrals, meta_2,sam_var)

    colnames(seqmat_integrals)[colnames(seqmat_integrals) == independent_var] <- "independent_var"
    seqmat_integrals[,"independent_var"] <- as.factor(seqmat_integrals[,"independent_var"])
    #remove duplicates
    seqmat_integrals <- distinct(seqmat_integrals)

    #split the data into a list of dataframes of either OTUs or clusters
    if (by_asv==FALSE){

      if (is.null(asv_list$h_clust)){
        stop("asv_list must have an h_clust element if parameter by_asv=FALSE")
      }

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$cluster))

    } else if (by_asv==TRUE){
      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$OTU))
    }

    #function to compute either mann whitney u test or kruskal test if testing for multiple groups
    compute_mw_test <- function(x){

      x <- kruskal_test(data = x, value ~ independent_var)
      return(x)

    }
    #iterate across each element and correct p value for multiple tests
    seqmat_mwu_list <- lapply(seqmat_integrals_list, compute_mw_test)
    seqmat_mwu_list <- lapply(seqmat_mwu_list, pvalue)
    seqmat_mwu_list <- p.adjust(p = seqmat_mwu_list, method = "BH")

  } else if (!missing(batch)){

    if (any(check_batch(meta = meta[,c(independent_var, batch)],
                        independent_var = independent_var, batch = batch) == TRUE)){
    stop("batch effect must be distributed among multiple levels of independent variable")

    }

    #remove unecessary vectors from meta, rename variables for modelling, merge with seqmat for modeling
    meta_2 <- meta[,names(meta) %in% c(sam_var,independent_var, batch)]
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == "variable"] <- sam_var
    seqmat_integrals<- left_join(seqmat_integrals, meta_2,sam_var)


    colnames(seqmat_integrals)[colnames(seqmat_integrals) == batch] <- "batch"
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == independent_var] <- "independent_var"
    seqmat_integrals[,"independent_var"] <- as.factor(seqmat_integrals[,"independent_var"])
    seqmat_integrals[,"batch"] <- as.factor(seqmat_integrals[,"batch"])
    #remove duplicate rows
    seqmat_integrals <- distinct(seqmat_integrals)

    #split into a list of dataframes by ASV or cluster
    if (by_asv==FALSE){

      if (is.null(asv_list$h_clust)){
        stop("asv_list must have an h_clust element if parameter by_asv=FALSE")
      }

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$cluster))

      } else if (by_asv==TRUE){

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$OTU))

      }

    #function to compute a stratified mann whitney u test or kruskal test for multiple groups
    compute_vanelteren <- function(x){
      x <- kruskal_test(data = x, value ~ independent_var | batch)
      return(x)
    }

    #iterate to compute stratified test and correct for multiple testing
    seqmat_mwu_list <- lapply(seqmat_integrals_list, compute_vanelteren)
    seqmat_mwu_list <- lapply(seqmat_mwu_list, pvalue)
    seqmat_mwu_list<- p.adjust(p = seqmat_mwu_list, method = "BH")

  }

  return(seqmat_mwu_list)

}
