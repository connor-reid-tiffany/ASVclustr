#' Creation of an ASVclustr list object.
#'
#' @description Constructor function that creates a list with the added class ASVclustr. The list needs to contain at least a sequence matrix, taxa matrix, and
#' metadata data.frame.
#'
#' @param seqmat A numeric matrix of sequencing reads, with ASVs as columns and unique sample identifiers as rownames.
#' The correct format is that of a dada2 sequence count matrix.
#' @param taxmat A character matrix of taxonomic information, with rownames matching seqmat column names, and taxonomic
#' rankings (Kingdom, Phylum, Class etc) as columns. The correct format is that of a dada2 taxonomy matrix.
#' @param meta A data frame with metadata for your sequencing samples. Rownames should be equal and in identical order to the seqmat.
#' It also msut contain the following: a column of non unique sample identifiers (i.e. mouse number, plant number, etc),
#' a column of time points as numeric, a column of independent variable(s) or grouped independent variables is optional.
#'
#' @return A list with seqmat, taxmat, meta. additionally adds the class ASVclustr
#'
#' @examples
#' asv_list <- asv_list(seqmat = seqmat, taxmat = taxmat, meta = meta)
#'
#' @export
asv_list <- function(seqmat, taxmat, meta){

  if (!is.matrix(seqmat) & !is.numeric(seqmat)){
    stop("seqmat must be a numeric matrix")
  }

  if (!is.matrix(taxmat) & !is.character(taxmat)){
    stop("taxmat must be a character matrix")
  }

  if (!is.data.frame(meta) & !identical(rownames(seqmat), rownames(meta))){
    stop("meta must be a data.frame with identical rownames to seqmat")
  }


  asv_list <- list(seqmat = seqmat, taxmat = taxmat, meta = meta)


  class(asv_list) <- append(class(asv_list), "ASVclustr")
  return(asv_list)

}



#' Melt an asv_list into a long form data.frame for plotting functions.
#'
#' @description Creates a data.frame from an asv list in long-form with each variable as a single vector.
#'
#' @param asv_list An asv_list.
#' @param ratio A boolean. TRUE will convert reads into a ratio for each sample. Default is FALSE.
#' @param rescale A boolean. TRUE will rescale reads for each ASV in each sample to a value between 0 (min) and 1 (max).
#' @param sam_var A string. The name of the sample variable in metadata.
#' @importFrom reshape2 melt
#' @importFrom dplyr left_join
#' @return A data.frame in longform with every variable as a single vector for use in plotting functions.
#'
#' @examples
#' melted_asv <- melt_asv_list(asv_list = asv_list, sam_var = "Sample", time_var = "Timepoint")
#'
#' @export
melt_asv_list <- function(asv_list, ratio=FALSE, rescale=FALSE, sam_var, time_var){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be a list object with class ASVclustr")
  }

  if (ratio==FALSE){
    #'access seqmat and melt into long form of 3 columns. SampleID, ASV, and
    seqmat <- as.data.frame(asv_list$seqmat)
  } else if (ratio==TRUE){
    seqmat <- t(asv_list$seqmat)
    seqmat <- apply(seqmat, FUN = function(x) x/sum(x), 2)
    seqmat <- as.data.frame(t(seqmat))
  }
  #' rescale values between 0 and 1
  if (rescale==TRUE){
    meta <- asv_list$meta

    rescale <- function(x){
      x <- (as.matrix(x)-min(as.matrix(x)))/(max(as.matrix(x))-min(as.matrix(x)))
      return(x)
    }

    seqmat_order_list <- order_seqmat(seqmat = seqmat, meta = meta, sam_var = sam_var, time_var = time_var)
    #'remove all metadata variables from each data.frame in the list except for the timepoint variable
    seqmat_order_list <- lapply(seqmat_order_list, function(x)
      x[!(names(x) %in% c("SampleID","sample", "timepoint"))])

    seqmat_order_list <- lapply(seqmat_order_list, function(x) sapply(x,  rescale))
    #' convert NA to 0
    conv_na_zero <- function(x){
      x[is.na(x)] <- 0
      return(x)
    }

    seqmat_order_list <- lapply(seqmat_order_list, conv_na_zero)
    seqmat <- do.call("rbind", seqmat_order_list)
    rownames(seqmat) <- rownames(seqmat_order)
    seqmat <- as.data.frame(seqmat)

    } else  if (rescale==FALSE){
    seqmat <- seqmat
    meta <- asv_list$meta
    meta$SampleID <- rownames(meta)
  }
  #'create SampleID column to left join with taxa and metadata and then melt to longform
  seqmat$SampleID <- rownames(seqmat)
  seqmat_melt <- melt(seqmat)

  colnames(seqmat_melt)[2] <- "OTU"
  colnames(seqmat_melt)[3] <- "Abundance"

  taxa <- as.data.frame(asv_list$taxmat)
  taxa$OTU <- as.factor(rownames(taxa))

  h_clust <- asv_list$h_clust

  #'left join to meta, then to taxa, then to cluster
  asv_melt <- left_join(seqmat_melt, meta, "SampleID")
  asv_melt <- left_join(asv_melt, taxa, "OTU")

  if (!is.null(asv_list$h_clust)){
    h_clust <- asv_list$h_clust
    asv_melt <- left_join(asv_melt, h_clust, "OTU")
  }

  return(asv_melt)
}





#' Sum by cluster
#'
#' @description Sums ASVs by cluster.
#' @param asv_list An asv list with an h_clust element.
#' @param abs_abund Optional. A data.frame with matching rownames to seqmat and meta within the asv_list. The first
#' column must be a numeric vector of total abundance values for each unique sample.
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#'
#' @return A data.frame in longform with normalized cluster abundances for each sample.
#'
#' @examples
#'
#' clustersum <- sum_by_cluster(asv_list = asv_list)
#'
#' clustersum_total <- sum_by_cluster(asv_list = asv_list, abs_abund = abs_abund)
#'
#' @export
sum_by_cluster <- function(asv_list, abs_abund){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be a list object with class ASVclustr")
  }

  if (is.null(asv_list$h_clust)){
    stop("asv_list must have an h_clust element to use function sum_by_cluster")
  }

  #'transform counts to ratio
  seqmat <- asv_list$seqmat

  seqmat <- t(apply(seqmat, 1, FUN = function(x) x/sum(x)))

  #'pull out metadata
  meta <- asv_list$meta
  meta$SampleID <- rownames(meta)

  if (!missing(abs_abund)){

  if (!is.numeric(abs_abund[,1])){
  stop("First column of abs_abund data.frame must be a numeric vector of total abundance values for each sample")
  }
  #'creat sampleID column to match normalized abundance values to
  abs_abund_vect <- abs_abund[,1]
  #'multiply values in rows by absolute abundance
  seqmat <- seqmat * abs_abund_vect
  seqmat <- as.data.frame(seqmat)

  } else if (missing(abs_abund)){

  seqmat <- as.data.frame(seqmat)
  }
  #'remake SampleID
  seqmat$SampleID <-rownames(seqmat)
  #'melt
  seqmat <- melt(seqmat)
  colnames(seqmat)[2] <- "OTU"
  colnames(seqmat)[3] <- "Abundance"
  #'join to hclust df by OTU
  seqmat <- left_join(seqmat, asv_list$h_clust, by = "OTU")
  seqmat <- left_join(seqmat, meta, "SampleID")
  #'Aggregate OTUs by cluster

  seqmat <- seqmat %>%
    group_by(cluster,SampleID) %>%
    summarise(cluster_abundance = sum(Abundance))

  #'merge with metadata
  seqmat <- left_join(seqmat, meta, "SampleID")

  return(seqmat)
}


#' Remove samples from seqmat and meta in an asv_list.
#' @description Either removes samples by a variable level, or removes samples by name.
#' @param asv_list An asv_list object.
#' @param variable Optional. A string of a variable with levels you want to remove.
#' @param variable_level A string or character vector. If the variable parameter is used, a string of the variable level you wish to remove.
#' Otherwise, a character vector of sample names you wish to remove.
#' @importFrom dplyr left_join
#'
#' @return An asv_list object.
#'
#' @examples
#'
#' asv_list_pruned <- remove_samples(asv_list=asv_list,variable_level = c("EVelazquez0717_480.1.B", "EVelazquez0717_480.1.C"))
#'
#' asv_list_dirty <- remove_samples(asv_list = asv_list, variable = "Treatment", variable_level = "Dirty-A")
#'
remove_samples <- function(asv_list, variable, variable_level){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be of classes list and ASVclustr")
  }

  #' establish variables
  seqmat <- as.data.frame(asv_list$seqmat)
  meta <- asv_list$meta

  seqmat$SampleID <- rownames(seqmat)
  meta$SampleID <- rownames(meta)

  #'remove samples within variable level from seqmat and meta
  if (!missing(variable)){
    meta2 <- meta[,c("SampleID", variable)]
    seqmat <- left_join(seqmat, meta2, "SampleID")

    seqmat_sub <- seqmat[!seqmat[,variable] %in% variable_level,]
    rownames(seqmat_sub) <- seqmat_sub$SampleID
    seqmat_sub <- seqmat_sub[,!names(seqmat_sub) %in% c("SampleID", variable)]

    meta <- meta[!meta[,variable]%in% variable_level,]

  #' remove samples by name
  } else if (missing(variable)){

    seqmat_sub <- seqmat[!seqmat[,"SampleID"] %in% variable_level,]
    meta <- meta[!meta[,"SampleID"] %in% variable_level,]
    rownames(seqmat_sub) <- seqmat_sub$SampleID
  }

  #' remove factor levels corresponding to variable levels that were removed
  seqmat_sub[] <- lapply(seqmat_sub, function(x) if(is.factor(x)) factor(x) else x)
  meta[] <- lapply(meta,function(x) if(is.factor(x)) factor(x) else x)

  rownames(meta) <- meta$SampleID
  meta <- meta[,!names(meta) %in% "SampleID"]
  seqmat_sub <- seqmat_sub[,!names(seqmat_sub) %in% 'SampleID']

  seqmat_sub <- seqmat_sub[, colSums(seqmat_sub)!= 0] #'this code removes ASVs that are 0 across all samples after subsetting
  seqmat_sub <- as.matrix(seqmat_sub)
  asv_list$seqmat <- seqmat_sub
  asv_list$meta <- meta
  return(asv_list)
}

#' Subset ASVs by read count and presence.
#'
#' @description Subsets to ASVs that reach a required threshold of reads and or presence across samples.
#' @param asv_list An asv_list object.
#' @param sam_var A string. The name of the sample variable in metadata.
#' @param sam_threshold An integer. Threshold for the number of samples in which an ASV is present.
#' @param abund_threshold Optional. An integer. Threshold for the amount of reads present for an ASV in the dataset.
#' @importFrom dplyr left_join
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_if
#' @importFrom magrittr %>%
#' @return an asv_list object
#'
#' @examples
#'
#' asv_list_subset_abund <- remove_samples(asv_list=asv_list,variable_level = c("EVelazquez0717_480.1.B", "EVelazquez0717_480.1.C"))
#'
#' asv_list_subset_<- remove_samples(asv_list = asv_list, variable = "Treatment", variable_level = "Dirty-A")
#'
#'
subset_ASVs <- function(asv_list, sam_var, sam_threshold, abund_threshold){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be of classes list and ASVclustr")
  }
  #'define variables
  seqmat <- as.data.frame(asv_list$seqmat)
  meta <- asv_list$meta

  seqmat$SampleID <- rownames(seqmat)
  meta$SampleID <- rownames(meta)

  join <- left_join(seqmat,meta, "SampleID")

  colnames(join)[colnames(join) == sam_var] <- "sample"

  #'sum by sample to determine abundance in each sample
  join <- join %>%
    group_by(sample) %>%
    summarise_if(is.numeric, base::sum)

  remove_cols <- colnames(meta)

  join <- join[,!names(join) %in% c(remove_cols, "sample","SampleID")]

  if(missing(abund_threshold)){

    #'keep ASVs present in n samples
    above_zero_counts <- data.frame(NonZeroCounts=colSums(join!=0)[colSums(join!=0)!=0])

    above_zero_counts$OTU <- rownames(above_zero_counts)

    above_zero_counts <- above_zero_counts[above_zero_counts[,1] >= sam_threshold,]

    keep_OTUs <- above_zero_counts$OTU

    seqmat <- seqmat[,names(seqmat) %in% keep_OTUs]

  } else if (!missing(abund_threshold)){

    #'remove ASVs below a certain abundance threshold within all samples
    above_thresh_counts <- as.data.frame(colSums(join >= abund_threshold))

    above_thresh_counts$OTU <- rownames(above_thresh_counts)
    #' keeo ASVs present in n samples
    above_thresh_counts <- above_thresh_counts[above_thresh_counts[,1] >= sam_threshold,]

    keep_OTUs <- above_thresh_counts$OTU

    seqmat <- seqmat[,names(seqmat) %in% keep_OTUs]
  }

  asv_list$seqmat <- as.matrix(seqmat)
  return(asv_list)
}
