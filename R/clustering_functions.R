#' Hierarchical clustering of ASVs.
#'
#' @description Calculates the cosine similarity, converts to an angular distance matrix, and then either
#' returns a dendrogram if no k is set, or returns an asv_list with an h_clust element.
#'
#' @param asv_list An asv_list.
#' @param agg_method A string of either "ward", "ward.D2", "complete", "single", "average", "mcquitty", "median", or "centroid".
#' "ward.D2" is the default. See ?hclust for more details.
#' @param k Optional. An integer. The number of clusters you wish to group your ASVs into. If a k is not chosen asv_clustr will
#' return a dendrogram plot instead
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
asv_clustr <- function(asv_list, agg_method = "ward.D2", td = "none", k = NULL){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){
    stop("asv_list must be of classes list and ASVclustr")
  }

  switch(td,

         both = {

           if (is.null(asv_list$seqmat_td)){
             stop("asv_list must have a seqmat_td element to use parameter td = 'both'")
           }

           seqmat <- rbind(asv_list$seqmat, asv_list$seqmat_td)
         },

         td_only = {

           if (is.null(asv_list$seqmat_td)){
             stop("asv_list must have a seqmat_td element to use parameter td = 'td_only'")
           }

           seqmat <- asv_list$seqmat_td

         },

         none = seqmat <- asv_list$seqmat,

         stop("input must be one of both, td_only, or none")
  )

  #compute cosine similarity and angular distance
  cosine_sim <- dist.matrix(M = seqmat, convert = TRUE, byrow = FALSE)

  #perform hierarchical clustering
  if (is.null(k)){

    h_clust <- hclust(d = as.dist(cosine_sim),
                      method = agg_method)
    return(h_clust)

  } else if (!is.null(k)){

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
#' @description Computes the time rate derivatives for every ASV in every sample and adds them as an
#' element in an asv list object.
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

  #arrange by timepoint in ascending order, convert into a list of data.frames
  #where each df is one sample
  seqmat_td_list <- order_seqmat(seqmat = seqmat_td, meta = meta,
                                 sam_var = sam_var, time_var = time_var)

  #remove character vectors from each dataframe
  seqmat_td_list <- lapply(seqmat_td_list, function(x)
    x[,!(names(x) %in% c("SampleID","sample"))])

  #This function will compute the time derivative for every sequential time interval for each
  #ASV and will be iterated across each data frame in the list using lapply
  #f'(x) = (f(x+h) - f(x))/h
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
  #change the rownames in seqtab_derivatives_list to avoid repeated rownames with seqmat
  seqmat_derivatives_list <- lapply(seqmat_derivatives_list,
                                    function(x){ rownames(x) <- paste0(rownames(x), ".","TD"); x})

  #rbind the two lists together, element by element
  seqmat_derivatives <- do.call("rbind", seqmat_derivatives_list)

  #remove timepoint_var column in each dataframe of seqtab_order_list
  seqmat_derivatives <- seqmat_derivatives[,colnames(seqmat_derivatives)!= "timepoint"]

  #convert na values to 0
  seqmat_derivatives <- apply(seqmat_derivatives, 2, function (x) ifelse(is.na(x), 0, x))

  #add derivatives to asv list
  asv_list$seqmat_td <- seqmat_derivatives

  return(asv_list)

}

#' Compute population dynamics across samples
#'
#' @description Computes the integral of the timecourse for each ASV in each sample and returns this as a nested list
#' in as_list
#' @param asv_list An asv_list.
#' @param sam_var A string. The metadata variable for samples.
#' @param time_var A string. The metadata variable for timepoints.
#' @param independent_var A string. The metadata variable for the experimental variable(s).
#' @param batch Optional. A string. A metadata variable for a batch effect. The batch effect must be separate
#' from the independent variable, i.e. levels of the batch variable must be distributed among the levels of
#' the idependent variable.
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
#' pd_cage_effect <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint",independent_var = "Treatment", batch="batch")
#'
#' pd_cage_effect_rescale <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", batch="batch", rescale = TRUE)
#'
#' pd_cage_effect_ASV <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", batch="batch", by_asv = TRUE)
#'
#' pd_cage_effect_ASV_rescale <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint",independent_var = "Treatment",
#' batch="batch", by_asv = TRUE, rescale = TRUE)
#'
#' pd_Treatment <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment")
#'
#' pd_Treatment_rescale <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", rescale = TRUE)
#'
#' pd_Treatment_ASV <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", by_asv = TRUE)
#'
#' pd_Treatment_ASV_rescale <- compute_pd(asv_list = asv_list, sam_var = "Sample",
#' time_var = "Timepoint", independent_var = "Treatment", by_asv = TRUE, rescale = TRUE)
#'
#'
compute_pd <- function(asv_list, sam_var, time_var, independent_var = NULL, batch = NULL, by_asv=FALSE, rescale=FALSE){

  if (!isClass(asv_list, Class = c("list", "ASVclustr"))){

    stop("asv_list must be an object of class list and ASVclustr")

  }

  if (by_asv==FALSE && is.null(asv_list$h_clust)){

    stop("asv_list must have an h_clust element if parameter by_asv=FALSE")

  }

  if (is.null(independent_var) && !is.null(batch)){

    stop("If testing only one covariate, only use the parameter independent_var, even if this covariate is a batch effect.
    Use independent_var and batch when testing the affect of a batch on another independent covariate")

  }

  if (!is.null(independent_var)){

  if (!is.null(batch) && any(check_batch(meta = asv_list$meta[,c(independent_var, batch)],
                      independent_var = independent_var, batch = batch) == FALSE)){

    stop("batch effect must be distributed among multiple levels of independent variable")

  }

}
  #define variables
  seqmat <- as.data.frame(asv_list$seqmat_norm)
  meta <- asv_list$meta
  seqmat$SampleID <- rownames(seqmat)

  #create a list of dataframes where each dataframe is one sample, and order them by timepoint
  seqmat_order_list <- order_seqmat(seqmat = seqmat, meta = meta,
                                    sam_var = sam_var, time_var = time_var)

  #extract timepoint column
  time_var_order_list <- lapply(seqmat_order_list, function(x) x[,"timepoint"])

  #remove all metadata variables from each dataframe in the list except for the timepoint variable
  seqmat_order_list <- lapply(seqmat_order_list, function(x)
    x[,!(names(x) %in% c("SampleID","sample", "timepoint"))])

  #normalize ASVs to fit within the same scale (removes influence of abundance differences)
  if (rescale==TRUE){

  #rescale the data
    seqmat_order_list <- lapply(seqmat_order_list, function(x) sapply(x,
                function (x)(as.matrix(x)-min(as.matrix(x)))/(max(as.matrix(x))-min(as.matrix(x)))))

  #convert any NA values to 0
  seqmat_order_list <- lapply(seqmat_order_list, function (x) ifelse(is.na(x), 0, x))

  } else if (rescale==FALSE){

    seqmat_order_list <- seqmat_order_list

  }
  #recombine individual sample data.frames back into one data.frame
  seqmat_order_list <- Map(cbind, seqmat_order_list, time_var_order_list)
  seqmat_order_list <- lapply(seqmat_order_list, as.data.frame)

  #compute time rate derivatives and also return the timepoint column for integration
   compute_derivatives <- function(x){

    d <- cbind(diff(as.matrix(x[,-ncol(x)]))/diff(x[,ncol(x)]), "time_interval" = diff(x[,ncol(x)]))
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
   if (by_asv==FALSE){

     seqmat_integrals <- left_join(seqmat_integrals, asv_list$h_clust,by = "OTU")
     seqmat_integrals <- melt(seqmat_integrals,  c("OTU", "cluster"))


   } else if (by_asv==TRUE){

     seqmat_integrals <- melt(seqmat_integrals, "OTU")

  }
  #test for significant difference in variance of ASV curves in each cluster
  #between independent variables
  if (is.null(batch)){
    #remove unecessary vectors from meta, rename variables for modelling, merge with seqmat for modeling
    if (!is.null(independent_var)){
    meta <- meta[,names(meta) %in% c(sam_var,independent_var)]
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == "variable"] <- sam_var
    seqmat_integrals <- left_join(seqmat_integrals, meta,sam_var)

    colnames(seqmat_integrals)[colnames(seqmat_integrals) == independent_var] <- "independent_var"
    seqmat_integrals[,"independent_var"] <- as.factor(seqmat_integrals[,"independent_var"])

    #remove duplicate rows
    seqmat_integrals <- distinct(seqmat_integrals)

    } else if (is.null(independent_var)){

      seqmat_integrals <- seqmat_integrals

    }
    #split the data into a list of dataframes of either OTUs or clusters
    if (by_asv==FALSE){

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$cluster))

    } else if (by_asv==TRUE){

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$OTU))

    }

  } else if (!is.null(batch)){

    #remove unecessary vectors from meta, rename variables for modelling, merge with seqmat for modeling
    meta <- meta[,names(meta) %in% c(sam_var,independent_var, batch)]
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == "variable"] <- sam_var

    seqmat_integrals<- left_join(seqmat_integrals, meta,sam_var)
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == batch] <- "batch"
    colnames(seqmat_integrals)[colnames(seqmat_integrals) == independent_var] <- "independent_var"

    #convert to factors because its a required class for the method
    seqmat_integrals[,"independent_var"] <- as.factor(seqmat_integrals[,"independent_var"])
    seqmat_integrals[,"batch"] <- as.factor(seqmat_integrals[,"batch"])

    #remove duplicate rows
    seqmat_integrals <- distinct(seqmat_integrals)

    #split into a list of dataframes by ASV or cluster
    if (by_asv==FALSE){

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$cluster))

      } else if (by_asv==TRUE){

      seqmat_integrals_list <- split(seqmat_integrals, f = as.factor(seqmat_integrals$OTU))

      }

  }

  asv_list$seqmat_integrals <- seqmat_integrals_list

  return(asv_list)

}


#' compare_pd
#' @description Compare cluster/ ASV integral means using a binomial-gamma hurdle model on median of ratios normalized
#' sequencing counts
#' @param asv_list an asv_list object
#' @param batch a boolean. TRUE if stratifying out a batch effect from another covariate, otherwise FALSE
#' @param calc_CI a boolean. TRUE will calculate 95% confidence intervals using non parametric bootstrapping
#' @param boot_k an integer. default is NULL. If calc_CI is TRUE, boot_k is the number of bootstrap replicates
#' to perform. if sample sizes are insufficient for number of bootstraps, a vector of ASVs/ clusters that had
#' insufficient sample size will be returned
#' @param groups a character vector. default is NULL. a character vector of strings denoting levels of an
#' independent covariate. the first string in the vector will be used as the reference group that all
#' other groups will be compared to.
#' @importFrom stats glm
#' @importFrom stats plogis
#' @importFrom broom tidy
#' @importFrom tibble tibble
#' @importFrom boot boot
#' @importFrom boot boot.ci
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom magrittr %>%
#' @return a dataframe with coefficients for the binomial, gamma, and hurdle models, means, p values and adjusted p values
#' and confidence intervals if calc_CI is TRUE
#' @export
#' @examples
#'

compare_pd <- function(asv_list, batch=FALSE, calc_CI=FALSE, boot_k=NULL,groups=NULL){

  #set object
  seqmat_integrals_list <- asv_list$seqmat_integrals

  #function to remove ASVs that are zero across all samples in a variable level
  remove_zero_ASVs <- function(x,variable){
    #create variable for whether or not an ASV is present in a sample
    x$non_zero <- ifelse(x$value > 0, 1, 0)
    #for the gamma distribution to model the continuous response variable of ASV population
    #during a timecourse (quantified by the integral of the abundance curve), we remove samples which
    #have an integral of zero
    non_zero_only <- subset(x, non_zero==1)
    non_zero_only[] <- lapply(non_zero_only, function(x) if(is.factor(x)) factor(x) else x)
    #this control flow step sets dataframes within a list to NULL if one factor level is all 0 values
    if (length(levels(non_zero_only[,variable])) <= length(levels(x[,variable])) - 1) {

      return(NULL)

    } else if (length(levels(non_zero_only[,variable])) == length(levels(x[,variable]))){

      return(x)

    }

  }

  #remove ASVs where there is a variable level that contains all 0 values for the integral if comparison is chosen
  if (!is.null(groups)){

    if (batch == FALSE){

      seqmat_integrals_list <- lapply(seqmat_integrals_list, remove_zero_ASVs, variable = "independent_var")

    } else if (batch == TRUE){

      seqmat_integrals_list <- lapply(seqmat_integrals_list, remove_zero_ASVs, variable = "independent_var")

      seqmat_integrals_list <- lapply(seqmat_integrals_list, remove_zero_ASVs, variable = "batch")

    }

    seqmat_integrals_list <- seqmat_integrals_list[lapply(seqmat_integrals_list, is.null)==FALSE]

  } else if (is.null(groups)){

    seqmat_integrals_list <- lapply(seqmat_integrals_list, function(x) cbind(x, non_zero = ifelse(x$value > 0, 1, 0)))

  }

  #set variable names for bootstrap repetition check function
  if (calc_CI == TRUE){

    if (is.null(groups)){

      variable <- NULL

    } else if (!is.null(groups)){

      if (batch == FALSE){

        variable <- "independent_var"

      } else if (batch == TRUE){

        variable <- "batch"

      }

    }

    boot_k <- boot_k

    #returns TRUE if sample size is sufficient for number of selected bootstraps, FALSE if not
    validate_reps <- function(x,variable){

      #remove rows which are zero
      x <- x[x$value > 0,]

      if (is.null(variable)){

        #group by variable, then tally sample size by group
        x <- x %>%
          tally()

      } else if (!is.null(variable)){

        #set column name for grouping
        colnames(x)[colnames(x) == variable] <- "variable"
        #group by variable, then tally sample size by group
        x <- x %>%
          group_by(variable) %>%
          tally()

      }

      #return logical of whether sample size is too low for R bootstraps
      x <- boot_k < choose(2*x$n - 1, x$n)

      x

    }

    #iterate on seqmat_integrals
    boot_check_list <- lapply(seqmat_integrals_list, validate_reps, variable=variable)

    boot_check <- boot_check_list[lapply(boot_check_list, function(x) any(x == FALSE))==TRUE]

    if (length(boot_check) > 0){

      message("One or more response variables have an insufficient sample size for selected number of bootstraps.
              Returning a vector of response variables with insufficient sample size")

      return(asv_vector <- names(boot_check))

    } else if (length(boot_check) == 0){


    }

  }

  #flow control for creating a model formula
  if (is.null(groups)){

    bin_model <- non_zero ~ 1
    gamma_model <- value ~ 1

  } else if (!is.null(groups)){

    if (batch == FALSE){

      bin_model <- non_zero ~ independent_var
      gamma_model <- value ~ independent_var

    } else if (batch == TRUE){

      bin_model <- non_zero ~ independent_var + batch
      gamma_model <- value ~ independent_var + batch

    }

  }

  #if a comparison vector is provided subset the variable to those two levels
  if (is.null(groups)){

    seqmat_integrals_list

  } else if (!is.null(groups)){

    seqmat_integrals_list <- lapply(seqmat_integrals_list, function(x) subset(x, independent_var %in% groups))
    #reorder the factor in the order provided in the comparison vector
    seqmat_integrals_list <- lapply(seqmat_integrals_list, function (x){ x$independent_var <- factor(x$independent_var, levels = groups); x})

  }



  #create a function for the model to iterate across the list
  bin_gamma_hurdle <- function(x){

    #fit a binomial logistic regression with the binary presence or absense as the response
    # to predict the probability of an ASV being present during the timecourse
    bin_mod <- glm(formula = bin_model, data = x, family = binomial(link = logit))
    #fit a gamma GLM to the integral values
    gamma_mod <- glm(formula = gamma_model, data = subset(x, non_zero==1),
                     family = Gamma(link = "log"))
    #extract binomial model coefficients
    bin_coef <- tidy(bin_mod)
    bin_coef$model <- "binomial"
    #extract gamma coefficients
    gamma_coef <- tidy(gamma_mod)
    gamma_coef$model <- "gamma"

    mod <- rbind(gamma_coef, bin_coef)
    #the probabiliy of the NULL hypothesis where ASV integrals are the same between groups is dependent on the probability of the
    #NULL hypothesis that groups do not have an effect on the presence or absence of an ASV, so the p value from the
    #binomial model is multiplied by the p value of the gamma model
    combined_pval <- tibble(p.value = bin_coef$p.value * gamma_coef$p.value, model = "hurdle",
                            term = gamma_coef$term, estimate = NA, std.error = NA, statistic = NA)
    #combine p_values and OTU into a data.frame
    mod <- rbind(mod, combined_pval)

  }


  compute_CI <- function(x, i){

    x <- x[i, ]
    #fit a binomial logistic regression with the binary presence or absense as the response
    # to predict the probability of an ASV being present during the timecourse
    bin_mod <- glm(formula = bin_model, data = x, family = binomial(link = logit))
    #fit a gamma GLM to the integral values
    gamma_mod <- glm(formula = gamma_model, data = subset(x, non_zero==1),
                     family = Gamma(link = "log"))
    #compute a CI on the model intercept, i.e. across all samples. useful for measuring dispersion
    if (is.null(variable)){
      #tidy the models into tibbles. the only observation in this case is the intercept coefficients
      bin_coef <- tidy(bin_mod)
      gamma_coef <- tidy(gamma_mod)

      #if a comparison between groups is being performed, remove the intercept coefficients to compute CI on desired model term
    } else if (!is.null(variable)){
      #tidy the models into tibbles and remove intercept observations
      bin_coef <- tidy(bin_mod)
      #bin_coef <- bin_coef[-1,]
      gamma_coef <- tidy(gamma_mod)
      #gamma_coef <- gamma_coef[-1,]

    }
    #return model estimates to original scale. binomial model used a logit link so take the inverse
    #gamma model used a log link so take the exponent.
    #take the sum of the logarithm of both estimate coefficients to compute the product. this is the statistic being bootstrapped
    bin_est <- plogis(bin_coef$estimate)
    gamma_est <- exp(gamma_coef$estimate)
    exp(log(bin_est) + log(gamma_est))

  }
  #iterate model across the list
  mod_list <- lapply(seqmat_integrals_list, bin_gamma_hurdle)
  #determine if we are comparing OTU means or cluster means
  cluster_or_OTU <- lapply(seqmat_integrals_list, function(x) colnames(x)=="cluster")
  cluster_or_OTU <- do.call("rbind", cluster_or_OTU)
  #set a variable to join the model dataframes and CI dataframes based on whether it is clusters or OTUs
  if (any(cluster_or_OTU==TRUE)==TRUE){

    name_variable <- "cluster"

  } else if (any(cluster_or_OTU==TRUE)==FALSE){

    name_variable <- "OTU"

  }
  #create list of OTU or cluster names to eventually join the mod list to the CI list
  names_df <- data.frame(names(mod_list))
  colnames(names_df)[1] <- name_variable
  names_list <- split(names_df, f = as.factor(names_df[,name_variable]))
  #combine names to mod list dataframes and then melt into longform. Then rename columns in each dataframe
  mod_list <- lapply(mod_list, `row.names<-`, NULL)
  names_list <- lapply(names_list, `row.names<-`, NULL)
  mod_list <- Map(f = cbind, mod_list, names_list)
  mod_list <- lapply(mod_list, function(x){ x$model_term <- paste(x$model, "_", x$term) ; x})

  #compute an adjusted p value for FDR for each model p value
  #compute padj for binomial model
  mod_list_binomial <- lapply(mod_list, function(x) subset(x, model=="binomial"))
  mod_list_binomial <- do.call("rbind", mod_list_binomial)
  mod_list_binomial <- split(mod_list_binomial, f = as.factor(mod_list_binomial$term))
  mod_list_binomial <- lapply(mod_list_binomial, function(x) {x$padj <- p.adjust(p = x$p.value, method = "BH"); x})
  mod_list_binomial <- do.call("rbind", mod_list_binomial)
  mod_list_binomial$term <- rownames(mod_list_binomial)
  #compute padj for gamma model
  mod_list_gamma <- lapply(mod_list, function(x) subset(x, model=="gamma"))
  mod_list_gamma <- do.call("rbind", mod_list_gamma)
  mod_list_gamma <- split(mod_list_gamma, f = as.factor(mod_list_gamma$term))
  mod_list_gamma <- lapply(mod_list_gamma, function(x) {x$padj <- p.adjust(p = x$p.value, method = "BH"); x})
  mod_list_gamma <- do.call("rbind", mod_list_gamma)
  mod_list_gamma$term <- rownames(mod_list_gamma)
  #compute padj for hurdle
  mod_list_hurdle <- lapply(mod_list, function(x) subset(x, model=="hurdle"))
  mod_list_hurdle <- do.call("rbind", mod_list_hurdle)
  mod_list_hurdle <- split(mod_list_hurdle, f = as.factor(mod_list_hurdle$term))
  mod_list_hurdle <- lapply(mod_list_hurdle, function(x) {x$padj <- p.adjust(p = x$p.value, method = "BH"); x})
  mod_list_hurdle <- do.call("rbind", mod_list_hurdle)
  mod_list_hurdle$term <- rownames(mod_list_hurdle)
  #row bind each of the padj dataframes into a single dataframe
  mod_list_padj <- rbind(mod_list_binomial, mod_list_gamma, mod_list_hurdle)
  mod_list_padj$term <- gsub(mod_list_padj$term, pattern = "\\..*", replacement = "")
  mod_list_padj$model_term <- paste(mod_list_padj$model, "_", mod_list_padj$term)
  #remove uneccessary columns and then split into a list by either cluster or OTU
  mod_list_padj <- mod_list_padj[,c("model", name_variable, "padj", "term", "model_term")]
  mod_list_padj <- split(mod_list_padj, f = as.factor(mod_list_padj[,name_variable]))
  #remove cluster/OTU column and then match the adjusted p values to the corresponding data.frames in mod_list
  mod_list_padj <- lapply(mod_list_padj, function(x) x[!names(x) %in% c(name_variable, "term", "model")])
  mod_list <- Map(left_join, mod_list, mod_list_padj, by = "model_term")

  return(mod_list)
  #compute a confidence interval using non-parametric bootstrapping
  if (calc_CI == TRUE){
    #extract model term names to use in CI_list
    term_names <- levels(as.factor(mod_list$`1`$term))
    #iterate bootstrap function across list, calculate and extract CI
    CI_list <- lapply(seqmat_integrals_list, boot, compute_CI, R = boot_k)
    #this function gets bca confidence intervals for all model terms including the intercept
    getCI <- function(x,i) {

      ci <- boot.ci(x,index=i)
      # extract info for bca
      ci <- t(sapply(ci["bca"],function(x) tail(c(x),2)))
      # combine with metadata: CI method, index
      ci_df <- cbind(i,rownames(ci),as.data.frame(ci))
      colnames(ci_df) <- c("term","method","lower","upper")
      ci_df

    }

    #suppressWarnings is used because the model wasn't a t distribution so getCI does not need to be supplied with
    #bootstrap variances
    CI_list <- suppressWarnings(lapply(CI_list, function(x) do.call(rbind,lapply(1:length(levels(as.factor(mod_list[[1]]$term))),getCI,x = x))))
    #create list of OTU or cluster names to eventually join the CI list to the model list
    names_df <- data.frame(names(CI_list))
    colnames(names_df)[1] <- name_variable
    names_list <- split(names_df, f = as.factor(names_df[,name_variable]))
    #combine names to CI list dataframes and then melt into longform. Then rename columns in each dataframe
    CI_list <- Map(f = cbind, CI_list, names_list)
    CI_list <- do.call("rbind", CI_list)
    CI_list$term <- term_names
    CI_list$term_cluster <- paste(CI_list$term, "_", CI_list[,name_variable])
    CI_list <- CI_list[,!names(CI_list) %in% c("cluster","term")]
    #join the model list and CI list by OTU/cluster names
    mod_list <- do.call("rbind", mod_list)
    mod_list$term_cluster <- paste(mod_list$term, "_", mod_list[,name_variable])
    #join CI_list to mod_list
    mod_df <- left_join(mod_list, CI_list, by = "term_cluster")
    mod_df <- mod_df[,!names(mod_df) %in% c("term_cluster", "model_term")]

  } else if (calc_CI == FALSE){

    #row bind all of the data.frames into a single tidy dataframe
    mod_df <- do.call("rbind", mod_list)

  }

  #clean up the term column using regex to remove words/ parenthesis
  mod_df$term <- gsub(mod_df$term, pattern = "\\(", replacement = "")

  mod_df$term <- gsub(mod_df$term, pattern = "\\)", replacement = "")

  mod_df$term <- gsub(mod_df$term, pattern = "independent_var", replacement = "")

  mod_df$name_variable <- mod_df[,name_variable]

  #create term + name_variable column to left join means to mod_df
  mod_df$term_namevar <- paste0(mod_df$term, "", mod_df$name_variable)

  #calculate base mean, term means, and fold changes where appropriate
  seqmat_integrals <- do.call("rbind", seqmat_integrals_list)

  if (name_variable == "cluster"){

    intercept_mean <- seqmat_integrals %>%

      group_by(cluster) %>%

      summarise(group_means = mean(value))

  } else if (name_variable == "OTU"){

    intercept_mean <- seqmat_integrals %>%

      group_by(OTU) %>%

      summarise(group_means = mean(value))

  }

  intercept_mean$term <- "Intercept"
  colnames(intercept_mean)[colnames(intercept_mean) == name_variable] <- "name_variable"
  intercept_mean$term_namevar <- paste0(intercept_mean$term, "", intercept_mean$name_variable)
  #join intercept mean to mod_df
  mod_df$intercept_means <- intercept_mean$group_means[match(mod_df$term_namevar, intercept_mean$term_namevar)]

  if (!is.null(groups)){

    #set base mean variable by taking the first element of the groups vector
    base_group <- groups[1]
    comparison_group <- groups[-1]

    mod_df$base_group <- base_group

    if (name_variable == "cluster"){

      group_means <- seqmat_integrals %>%

        group_by(cluster, independent_var) %>%

        summarise(group_means = mean(value))

    } else if (name_variable == "OTU"){

      group_means <- seqmat_integrals %>%

        group_by(OTU, independent_var) %>%

        summarise(group_means = mean(value))

    }

    base_means <- subset(group_means, independent_var==base_group)
    colnames(base_means)[colnames(base_means) == name_variable] <- "name_variable"
    group_means <- subset(group_means, independent_var %in% comparison_group)

    colnames(group_means)[colnames(group_means) == name_variable] <- "name_variable"
    group_means$term_namevar <- paste0(group_means$independent_var, "", group_means$name_variable)


    mod_df$base_mean <- base_means$group_means[match(mod_df$name_variable, base_means$name_variable)]
    mod_df$term_means <- group_means$group_means[match(mod_df$term_namevar, group_means$term_namevar)]
    mod_df$term_means[is.na(mod_df$term_means)] <- mod_df$intercept_means[!is.na(mod_df$intercept_means)]

    mod_df <- mod_df[, !names(mod_df) %in% "intercept_means"]
  }

  mod_df <- mod_df[,!names(mod_df) %in% c("term_namevar", "name_variable")]

  return(mod_df)

}

