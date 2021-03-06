#' @title dispRity object summary
#'
#' @description Creates a summary of a \code{dispRity} object.
#'
#' @param data A \code{dispRity} object.
#' @param quantiles The quantiles to display (default is \code{quantiles = c(50, 95)}; is ignored if the \code{dispRity} object is not bootstrapped).
#' @param cent.tend A function for summarising the bootstrapped disparity values (default is \code{\link[stats]{median}}).
#' @param recall \code{logical} value specifying whether to recall the \code{dispRity} parameters input (default = \code{FALSE}).
#' @param rounding Optional, a value for rounding the values in the output table (default = 2).
# ' @param results Optional, in the case of summarising a \code{\link{sequential.test}} which results to display (default = "coefficients")
#'
#' @return
#' A \code{data.frame} with:
#' \item{subsamples}{the subsample names.}
#' \item{n}{the number of elements in each subsample.}
#' \item{observed}{the observed disparity or the the observed central tendency (<cent_tend>) of disparity (\code{obs.<cent_tend>}).}
#' \item{bootstraps...}{if \code{data} is bootstrapped, the bootstrapped disparity's central tendency (\code{bs.<cent_tend>}) and the quantiles of the bootstrapped disparities (or, if \code{data} is not bootstrapped but disparity is calculated as a distribution - see \code{\link[dispRity]{dispRity}}) - the quantiles of the observed disparity are displayed).}
#' 
#' @examples
#' ## Load the disparity data based on Beck & Lee 2014
#' data(disparity)
#'
#' ## Summarising the results
#' summary(disparity) # default
#' ## Using different options
#' summary(disparity, quantiles = 75, cent.tend = mean, rounding = 8,
#'      recall = TRUE)
#' 
#' @seealso \code{\link{dispRity}}, \code{\link{plot.dispRity}}.
#'
#' @author Thomas Guillerme

## DEBUG
# source("sanitizing.R")
# source("make.metric.R")
# source("summary.dispRity_fun.R")
# data(BeckLee_mat50)
# groups <- as.data.frame(matrix(data = c(rep(1, 12), rep(2, 13), rep(3, 12), rep(4, 13)), dimnames = list(rownames(BeckLee_mat50))), ncol = 1)
# customised_subsamples <- custom.subsamples(BeckLee_mat50, groups)
# bootstrapped_data <- boot.matrix(customised_subsamples, bootstraps = 3, rarefaction = TRUE)
# subsamples <- extract.dispRity(sum_of_variances, observed = FALSE, keep.structure = TRUE, concatenate = TRUE)
# data <- sequential.test(subsamples, family = gaussian, correction = "hommel")

# data <- dispRity(bootstrapped_data, metric = variances)

# quantiles <- c(50, 95)
# cent.tend <- median
# recall <- FALSE
# match_call <- list() ; match_call$cent.tend <- "median"

summary.dispRity <- function(data, quantiles = c(50, 95), cent.tend = median, recall = FALSE, rounding){#, results = "coefficients") {

    #----------------------
    # SANITIZING
    #----------------------

    #Get call
    match_call <- match.call()
    #return(match_call)

    #cent.tend
    #Must be a function
    check.class(cent.tend, "function")
    #The function must work
    if(make.metric(cent.tend, silent = TRUE) != "level1") {
        stop(paste(match_call$cent.tend), "can not be used for measuring the central tendency.")
    }
    ## Update match_call if argument is empty
    if(is.null(match_call$cent.tend)) match_call$cent.tend <- "median"

    #recall
    check.class(recall, "logical")

    #rounding
    if(missing(rounding)) {
        #Set to default (see below)
        rounding <- "default"
    } else {
        check.class(rounding, "numeric")
    }

    #DATA
    #must be class dispRity
    check.class(data, "dispRity")
    #Check if it is a bootstrapped dispRity object
    if(is.null(data$disparity)) {
        stop("Disparity has not been calculated yet.\nUse the dispRity() function to do so.\n", sep = "")
    }
    
    #Check quantiles
    check.class(quantiles, "numeric", " must be any value between 1 and 100.")
    if(any(quantiles < 1) | any(quantiles > 100)) {
        stop("quantiles(s) must be any value between 1 and 100.")
    }

    #----------------------
    # TRANSFORMING THE DATA INTO A TABLE
    #----------------------

    ## Check if disparity is a value or a distribution
    is_distribution <- ifelse(length(data$disparity[[1]]$elements) != 1, TRUE, FALSE)

    ## Check the bootstraps
    bootstrapped <- ifelse(!is.null(data$call$bootstrap), TRUE, FALSE)

    ## Get the elements per subsamples
    elements <- lapply(data$subsamples, lapply.get.elements, bootstrapped)
    if(is.null(elements[[1]])) {
        elements <- list(nrow(data$subsamples[[1]]$elements))
    }

    ## Get the names of the subsamples
    names <- names(data$subsamples)
    if(is.null(names)) {
        names <- seq(1:length(data$subsamples))
    }

    ## Get the disparity values
    disparity_values <- lapply(data$disparity, lapply.observed)
    names(disparity_values) <- NULL

    ## Initialise the results
    summary_results <- data.frame(row.names = NULL, "subsamples" = rep(names, unlist(lapply(elements, length))), "n" = unlist(elements))

    ## Add the observed values
    if(is_distribution) {
        summary_results <- cbind(summary_results, as.vector(unlist(mapply(mapply.observed, lapply(disparity_values, cent.tend), elements))), row.names = NULL)
        names(summary_results)[3] <- paste("obs", match_call$cent.tend, sep=".")
    } else {
        summary_results <- cbind(summary_results, as.vector(unlist(mapply(mapply.observed, disparity_values, elements))), row.names = NULL)
        names(summary_results)[3] <- "obs"
    }

    if(!is.null(data$call$bootstrap)) {
        ## Calculate the central tendencies and the quantiles
        summary_results <- cbind(summary_results, matrix(unlist(lapply(data$disparity, lapply.summary, cent.tend, quantiles)), byrow = TRUE, ncol = (1+length(quantiles)*2)))
        ## Adding the labels
        names(summary_results)[4:length(summary_results)] <- c(paste("bs", match_call$cent.tend, sep="."), names(quantile(rnorm(5), probs = CI.converter(quantiles))))
    } else {
        if(is_distribution) {
            ## Calculate the quantiles
            summary_results <- cbind(summary_results, matrix(unlist(lapply(data$disparity, lapply, get.summary, quantiles = quantiles)), byrow = TRUE, ncol = (length(quantiles)*2)))
            ## Adding the labels
            names(summary_results)[4:length(summary_results)] <- c(names(quantile(rnorm(5), probs = CI.converter(quantiles))))
        }
    }

    ## Round the results (number of decimals = maximum number of digits in the output)
    summary_results <- rounding.fun(summary_results, rounding)

    #----------------------
    # OUTPUT
    #----------------------
    if(recall) print.dispRity(data)

    return(summary_results)

    #Summary sequential.test shortcut
    # if(length(class(data)) == 2) {
    #     if(class(data)[[1]] == "dispRity" && class(data)[[2]] == "seq.test") {
    #         #Results sanitizing
    #         check.class(results, "character") 
    #         #At least one must be "coefficients"
    #         if(is.na(match("coefficients", results))) {
    #             stop("At least one of the returned results must be 'coefficients'.")
    #         }
    #         #Results must be at least coefficients
    #         if(is.na(match("coefficients", results))) {
    #             results <- c(results, "coefficients")
    #         }

    #         #Creating the table results
    #         results_out <- summary.seq.test(data, quantiles, cent.tend, recall, rounding, results, match_call)

    #         #Checking if distribution
    #         is.distribution <- ifelse(length(data$models[[1]]) == 1, FALSE, TRUE)

    #         #Rounding the results
    #         if(is.distribution == FALSE) {
    #             results_out <- lapply(results_out, rounding.fun, rounding, seq.test = TRUE)
    #         } else {
    #             results_out <- lapply(results_out, lapply, rounding.fun, rounding, seq.test = TRUE)
    #         }

    #         return(results_out)
    #     }
    #     if(class(data)[1] == "dispRity" & class(data)[2] == "randtest") {
    #         #No summary
    #         return(data)
    #     }
    # }


}
