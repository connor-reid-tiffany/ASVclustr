#' Order seqmat by sample and timepoint
#'
#'
#' @description Internal function to order a seqmat object by time point and convert to a list of dataframes with each sample as an element.
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

  meta$SampleID <- rownames(meta)
  meta2 <- as.data.frame(meta[,c(sam_var,time_var,"SampleID")])
  #add column to seqtab to join the seqtab df with metadata

  seqmat$SampleID <- rownames(seqmat)
  seqmat <-left_join(seqmat, meta2, "SampleID")
  #convert the timepoint variable to numeric, then order the dataframe by Sample and timepoint

  seqmat[,time_var] <- as.numeric(as.character(seqmat[,time_var]))
  colnames(seqmat)[colnames(seqmat) == sam_var] <- "sample"
  colnames(seqmat)[colnames(seqmat) == time_var] <- "timepoint"
  sample_var_2 = "sample"
  timepoint_var_2 = "timepoint"
  seqmat_order <- arrange(seqmat, sample,timepoint)
  rownames(seqmat_order) <- seqmat_order$SampleID
  #split the dataframe into a list of dataframes, where each dataframe is one Sample

  seqmat_order_list <- split(seqmat_order,
                             f = as.factor(seqmat_order[,sample_var_2]))
}


