#' @name model.test
#'
#' @title Model Test
#'
#' @description Fit models of disparity change through time
#'
#' @param data A \code{dispRity} object used to test models of evolution through time.
#' @param models The model(s) of evolution to allow for changes in disparity-through-time using a homogenous or hetergenous model. A single model can be used as input or multiple modes can be supplied in a vector. If a vector with multiple modes is supplied then the model will test for shifts in modes at the time supplied by \code{time.split}. If a list is supplied different models will be fit to the data (See \bold{Details}).
#' @param pool.variance If \code{NULL} (default) the difference in variances will be calculated using Bartlett's Test (\code{\limk[stats]{bartlett.test}}) of equal variances. If there is no significant difference among variances, then variance in samples will be pooled and the same variance will be used for all samples. A significance difference will not pool variances and the original variance will be used for model-testing. If argument \code{TRUE} or \code{FALSE} are used, Bartlett's test will be ignored and the analyses will use the user-set pooling of variances.
#' @param time.split The age of the change in mode (\code{numeric}). The age is measured as the time before the most recent sample, and multiple ages can be supplied in a vector. If no age is supplied for models then all possible time shifts are fit in the model, and the highest likelihood model is returned. Note this only applies to heterogenous models (See \bold{Details}).
#' @param fixed.optima A \code{logical} value, whether to use an estimated optimum value in OU models (\code{FALSE} - default), or whether to set the OU optimum to the ancestral value (\code{TRUE}).
#' @param control.list A \code{list} of fine-tune control inputs for the optim function.
#' @param verbose \code{logical}, whether to display the model results as computed (\code{TRUE} - default).
#' 
#' @details The models are fit using maximum likelihood optimisation using the function optim. Fine-tuning of the search algorithms can be applied using the control.list argument. Models can be fit using a homogenous model with the same process applied to the entire sequence or models with time splits that represent a change in parameters or a shift in mode. For the time split model if a time value is provided, then the shift is tested at that value only. If no time shift is supplied then all shift times that allow for there to be at least 10 samples in each time bin are tested. If the sample is fewer than 30 samples long then no time splits are searched for (unless a time split is supplied by the user). Parameters are shared across different modes. For example, c("BM", "OU") would fit a model in which the process starts with a BM model and shifts toan OU process. The ancestral value at the start of the sequence and sigma squared value are shared across the models. Any combination of the following homogenous models (with the exception of 'multi.OU') can be fit to the data:
#' 
#'  \item{"BM"}{Fits a unbiased random walk model of evolution (Felsenstein 1985; Hunt 2006). The model optimises the ancestral state and the 'step-variance' (sigma-squared)}
#'  \item{"OU"}{The Ornstein-Uhlenbeck model of evolution in which the change in variance is constrained to an optimum value (Hansen 1997). In this model there are three parameters: optima, alpha, and ancestral state. The strength of attraction based on the parameter alpha and the ancestral state is estimated from the data. The optima value is estimated from the data, and this can lead to optima being found outside the known data values. If this is the case the model is similar to a trend model. If the argument, fixed.optima is set to TRUE, the model will not estimate optima but constrain it to the first value in the sequence}
#'  \item{"Trend"}{Fits a Brownian motion model with a directional component. This model is also known as the General Random Walk (Hunt 2006). This model has three parameters: the ancestral state, the 'step-variance' (sigma-squared), and the trend component
#"Stasis" in which traits evolve with variance (omega) around a mean (theta). This model is time-independent in that the model is guided only by the variance and attraction to the mean (Hunt 2006)}
#'  \item{"Stasis"}{Fits a model in which traits evolve with variance (omega) around a mean (theta). This model is time-independent in that the model is guided only by the variance and attraction to the mean (Hunt 2006)}
#'  \item{"Early Burst(EB)"}{Trait variance accumulates early in the evolution of a trait and decreases exponentially through time (Blomberg et al. 2003; Harmon et al. 2010). This model has three parameters: ancestral state, sigma-squared, and the exponential rate of decrease}
#'  \item{"multi.OU"}{Fits a model in which the value of the optima shifts at one or more time splits. The values of the 'step-variance' (sigma squared) and attraction to the optima (alpha) are shared across all the samples. This model can not be fit with other models - the multiOU system can be be fit to the model only}
#'
#' @examples
#' ## To Add
#' @seealso \code{\link{model.test.sim}}.
#'
#' @references
#' To Add: Hunt 2006, Hunt 2008, Harmon 2010, 
#' 
#' @author Mark N Puttick and Thomas Guillerme
#' @export
#' 


# source("sanitizing.R")
# source("model.test_fun.R")
# set.seed(123)
# data(BeckLee_mat99)
# data(BeckLee_ages)
# data(BeckLee_tree)
# continuousData <- time.subsamples(BeckLee_mat99, BeckLee_tree, method="continuous", time=seq(120, 0, length.out=120), model="gradual")
# data_bootstrapped <- boot.matrix(continuousData)
# data <- dispRity(data_bootstrapped, c(sum, variances))
# models <- list("BM", "OU", "multi.OU", c("BM", "Trend"))
# pool.variance=NULL
# fixed.optima=FALSE
# verbose = TRUE
# control.list=list(fnscale = -1)
# time.split <- 66

model.test <- function(data, models, pool.variance = NULL, time.split, fixed.optima = FALSE, control.list = list(fnscale = -1), verbose = TRUE) {
    
    ## data
    check.class(data, "dispRity")
    model_test_input <- select.model.list(data)

    ## models
    check.method(models, c("BM", "OU", "Trend", "Stasis", "EB"), msg = "models")
    n_models <- length(models)

    ## pool.variance
    if(!is.null(pool.variance)) {
        check.class(pool.variance, "logical")
    }

    ## time.split
    silent <- check.class(time.split, c("numeric", "integer"))

    ## fixed optima
    check.class(fixed.optima, "logical")

    ## control.list
    check.class(control.list, "list")
        
    ## use Bartlett's test of variance to decide whether to pool variance or not (not used if pool variance is specified as TRUE or FALSE before-hand)
    if(is.null(pool.variance)) {
        p_test <- bartlett.variance(model_test_input)
        #p_test <- stats::bartlett.test(model_test_input)$p.value TG: is this equivalent?
        if(p_test < 0.05) {
            pool.variance <- FALSE
            if(verbose) cat(paste0("Evidence of equal variance (Bartlett's test of equal variances p = ", signif(p_test, 3), ").\nVariance is not pooled.\n"))
        } else {
            pool.variance <- TRUE
            if(verbose) cat(paste0("Evidence of unequal variance (Bartlett's test of equal variances p = ", signif(p_test, 3), ").\nVariance is pooled.\n"))
            model_test_input <- pooled.variance(model_test_input, TRUE)
        }
    } else {
        if(pool.variance) {
            model_test_input <- pooled.variance(model_test_input, TRUE)
        }
    } 

    ## Running the models
    models_out <- lapply(1:n_models, lapply.models, models, model_test_input, time.split, verbose)
    
    ## Calculate the models AIC and AICc
    model_parameters <- sapply(models, length) - 1 + sapply(models_out, function(x) length(x[[1]]))
    model_likelihoods <- unlist(sapply(models_out, function(x) x[2]))
    sample_size <- length(model_test_input[[1]])
    aic <- (-2 * model_likelihoods) + (2 * model_parameters)
    aicc <- (-2 * model_likelihoods) + (2 * model_parameters) * (sample_size / ( sample_size - model_parameters - 1))
    delta_aicc <- aicc - min(aicc)
    weight_aicc <- exp(-0.5 * delta_aicc) / sum(exp(-0.5 * delta_aicc))
    order_aicc <- order(weight_aicc, decreasing = TRUE)

    ## Get the model names
    model_names <- sapply(models, function(x) paste(x, collapse = ":"))
    names(weight_aicc) <- names(delta_aicc) <- names(aicc) <- names(aic) <- names(models_out) <- model_names

    ## Generate the output format
    return_out <- list("aicc.models" = cbind(aicc, delta_aicc, weight_aicc), "full.details" = models_out)    
    class(return_out) <- c("dispRity", "model.test")
    return(return_out)
}    
