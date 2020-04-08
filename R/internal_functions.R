#' Order seqmat by sample and timepoint
#'
#'
#' @description Internal function to order a seqmat object by time point and
#' convert to a list of dataframes with each sample as an element.
#' @param seqmat seqmat from an asv_list
#' @param meta metadata from an asv_list
#' @param sam_var sample variable from metadata
#' @param time_var timepoint variable from metadata
#' @importFrom dplyr left_join
#' @importFrom dplyr arrange
#' @return a list of seqmats for every sample ordered by timepoint
#' @export
#' @examples
#' seqmat <- as.data.frame(asv_list$seqmat)
#' meta <- asv_list$meta
#' seqmat_order_list <- order_seqmat(seqmat = seqmat, meta = meta, sam_var = "Sample",
#' time_var = "Timepoint")
#'
#'
order_seqmat <- function(seqmat, meta, sam_var,time_var){

  timepoint <- NULL

  #add SampleID column to join meta with seqmat, then subset to required columns
  meta$SampleID <- rownames(meta)
  meta <- meta[,c(sam_var,time_var,"SampleID")]

  #add column to seqtab to join the seqtab df with metadata
  seqmat$SampleID <- rownames(seqmat)
  seqmat <-left_join(seqmat, meta, "SampleID")

  #convert the timepoint variable to numeric, rename variables, then order the dataframe
  #by Sample and timepoint
  seqmat[,time_var] <- as.numeric(as.character(seqmat[,time_var]))
  colnames(seqmat)[colnames(seqmat) == sam_var] <- "sample"
  colnames(seqmat)[colnames(seqmat) == time_var] <- "timepoint"

  seqmat_order <- arrange(seqmat, sample,timepoint)
  rownames(seqmat_order) <- seqmat_order$SampleID

  #split the dataframe into a list of dataframes, where each dataframe is one Sample
  seqmat_order_list <- split(seqmat_order,
                             f = as.factor(seqmat_order[,"sample"]))

  return(seqmat_order_list)

}

#' Check metadata for independence of batch effect
#'
#'
#' @description Internal function to produce a logical vector with an element for each batch effect factor level.
#' TRUE means that the batch level is inseparable from the independent variable (i.e. not distributed among multiple
#' levels of the independent variable)
#' @param meta A data.frame. Metadata from an asv_list
#' @param independent_var A string. Independent variable from metadata
#' @param batch A string. Batch effect variable from metadatanew
#' @return a list of seqmats for every sample ordered by timepoint
#' @export
#' @examples
#'
#' meta <- asv_list$meta
#' true_false_vect <- check_batch(meta = meta, independent_var = "Treatment", batch = "cage")
#'

check_batch <- function(meta, independent_var, batch){

  #subset meta to variables needed, split into list of data.frames by batch level
  independent_fact <- as.factor(meta[,independent_var])
  meta <- meta[,c(independent_var,batch)]
  meta <- split(meta, f = as.factor(meta[,batch]))

  #remove missing factor levels, reconvert elements to data.frames from matrices
  meta <- lapply(meta, function(x) sapply(x[], function (x) if(is.factor(x)) factor(x) else x))
  meta <- lapply(meta, as.data.frame)

  #create function to measure the length of factor levels of independent variable in each data.frame
  level_lengths <- function(x, independent_var){

    length(levels(as.factor(x[,independent_var])))

  }

  #iterate function across list elements
  lengths <- sapply(meta, level_lengths, independent_var)


  #create logical vector from resulting vector. If the length is 1 it means there is only
  #one independent variable level for its respective batch, which means the batch level is inseparable
  ifelse(lengths == length(levels(independent_fact)), TRUE, FALSE)

}
