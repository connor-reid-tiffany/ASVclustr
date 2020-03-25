#' An example asv_list
#'
#' A list containing the following : seqmat (a numeric matrix), taxmat (a character matrix),
#' meta (a data.frame), seqmat_td (a numeric matrix), and h_clust (a data.frame).
#' @format A list with 5 elements.
#'
"asv_list"


#' An example melted_asv_list data.frame
#'
#'A dataframe in longform. Suitable for use with ggplot2.
#' @format 3168 observations of 17 variables
#'
"melted_asv_list"

#' An example sequence count matrix
#'
#' A numeric matrix of sequence read counts.
#' @format A numeric matrix with 44 rows and 323 columns
#'
"seqmat"

#' An example metadata data.frame
#'
#'A data.frame with metadata for the samples (rows) in seqmat
#' @format A data.frame with 44 observations and 7 variables
#'
"meta"

#' An example taxonomy matrix
#'
#' A character matrix with taxonomy data for the ASVs (columns) in seqmat
#' @format A character matrix with 378 rows and 6 columns
#'
"taxmat"

#' An example absolute abundance data.frame
#'
#'A data.frame where the first column is a numeric vector of absolute abundance values
#' for each sample in seqmat (rownames are identical to seqmat)
#' @format A data.rame with 44 observations and 1 variable
#'
"abs_abund"

#' An example phylo tree object
#'
#' A pml object with a phylo element named 'tree' that was generated from the seqmat example data by
#' making a 16S gene alignment using the package DECIPHER and making a phylogenetic tree from the
#' alignment using the package phangorn.
#' @format A pml object with 22 elements.
#'
"phangtree"

