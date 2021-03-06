% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_functions.R
\name{check_batch}
\alias{check_batch}
\title{Check metadata for independence of batch effect}
\usage{
check_batch(meta, independent_var, batch)
}
\arguments{
\item{meta}{A data.frame. Metadata from an asv_list}

\item{independent_var}{A string. Independent variable from metadata}

\item{batch}{A string. Batch effect variable from metadatanew}
}
\value{
a list of seqmats for every sample ordered by timepoint
}
\description{
Internal function to produce a logical vector with an element for each batch effect factor level.
TRUE means that the batch level is inseparable from the independent variable (i.e. not distributed among multiple
levels of the independent variable)
}
\examples{

meta <- asv_list$meta
true_false_vect <- check_batch(meta = meta, independent_var = "Treatment", batch = "cage")

}
