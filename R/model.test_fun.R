#~~~~~~~~~~~~~~~~~~~
# Data handling functions
#~~~~~~~~~~~~~~~~~~~

## Select elements for a model list from a dispRity object
select.model.list <- function(data, observed = TRUE, cent.tend = median, rarefaction) {

    if(observed) {
        ## If observed is required
        central_tendency <- unlist(extract.dispRity(data, observed = TRUE))

        ## If disparity is a single value
        if(unique(unlist(lapply(data$disparity, lapply, lapply, length))) != 1) {
            ## Calculate the variance from the disparity data
            variance <- unlist(lapply(extract.dispRity(data, observed = FALSE), lapply, var))
        } else {
            ## Extract directly the variance from the data
            variance <- sapply(data[[3]], function(x) var(data[[1]][x[[1]]]))
        }

    } else {

        ## Getting the disparity
        if(!missing(rarefaction)) {
            disparity_tmp <- extract.dispRity(data, observed = FALSE, rarefaction = rarefaction)    
        } else {
            disparity_tmp <- extract.dispRity(data, observed = FALSE)
        }
        ## Calculating the central tendency
        central_tendency <- unlist(lapply(disparity_tmp, lapply, cent.tend))
        ## Calculating the variance
        variance <- unlist(lapply(disparity_tmp, lapply, var))
    }

    ## Getting the length of the samples
    summary_table <- summary(data)
    if(!missing(rarefaction)) {
        sample_length <- rep(rarefaction, length(central_tendency))
    } else {
        sample_length <- summary_table$n[which(!is.na(summary_table[,3]))]
    }

    ## Samples
    if(data$call$subsamples[1] == "continuous") {
        subsamples <- sort(as.numeric(names(data$subsamples)))
    } else {
        subsamples <- seq(1:length(data$subsamples))
    }

    # Making sure the samples are in the right chronological order
    if(subsamples[1] < subsamples[2]) {
        subsamples <- rev(subsamples)
    }

    ## Returns the data
    return(list("central_tendency" = central_tendency,
                "variance" = variance,
                "sample_size" = sample_length,
                "subsamples" = subsamples))
}

#~~~~~~~~~~~~~~~~~~~
# Variance functions
#~~~~~~~~~~~~~~~~~~~

## Bartlett variance test
bartlett.variance <- function(model.test_input) {
    variance.pooled <- pooled.variance(model.test_input)
    total.n <- sum(model.test_input$sample_size)
    total.group.n <- length(model.test_input$variance)
    numerator <- (total.n - total.group.n) * log(variance.pooled) - (sum((model.test_input$sample_size - 1) * log(model.test_input$variance)))
    denominator <- 1 + (1 / (3 * (total.group.n -1))) * ((sum(1 / (model.test_input$sample_size - 1))) - (1 / (total.n - total.group.n)))
    test.statistic <- numerator / denominator
    pchisq(test.statistic, df = total.group.n - 1, lower.tail = FALSE)
}

## Pooling variance function
pooled.variance <- function(data.model.test, rescale.variance = FALSE)  {
    sample.size.vector <- data.model.test$sample_size - 1
    var.vector <- data.model.test$variance
    pooled_variance <- sum(sample.size.vector * var.vector) / sum(sample.size.vector)
    if (rescale.variance) {
        data.model.test_out <- data.model.test
        data.model.test_out$variance <- rep(pooled_variance, length(data.model.test$central_tendency))
        return(data.model.test_out)
    } else {
        return(pooled_variance)
    }
}

#~~~~~~~~~~~~~~~~~~~
# Model testing wrapper
#~~~~~~~~~~~~~~~~~~~

## Lapply loop for applying the models
lapply.models <- function(one_model, models, model_test_input, time.split, control.list, fixed.optima, verbose) {
    
    ## Selecting the first model
    model_type <- models[[one_model]]
    
    ## Select the time.split
    if(length(model_type) == 1 && model_type != "multi.OU") {
        time.split <- NULL
    } else {
        time.split <- time.split
    }
    
    ## Run a multi.OU
    if(is.null(time.split) && length(model_type) == 2 || is.null(time.split) && model_type == "multi.OU") {
        
        ## Select all times to test
        all_times <- rev(model_test_input[[4]])
        
        ## Running multiple models
        if(length(all_times) > 31) { #TG: I'm still not sure about these numbers... I guess it's a minimum amount of variance thing but what about the the ten_times variable after?

            ten_times <- 9 : (length(all_times) - 11)

            if(verbose) {
                cat(paste0("Running ",  model_type ," on ", length(ten_times), " shift times...\n"))
            }
            
            ## Lapply wrapper
            lapply.all.time <- function(one_time, model_test_input, ten_times, model_type, control.list, fixed.optima, verbose) {
                if(verbose) cat(paste0("    model ", one_time - 8, " of ", length(ten_times), "\n"))
                return(model.test.lik(model_test_input, model_type, time.split = one_time, control.list, fixed.optima = fixed.optima))
            }

            ## Run all the models
            model_test_all_times <- lapply(as.list(ten_times), lapply.all.time, model_test_input, ten_times, model_type, control.list, fixed.optima, verbose)

            ## Select the best model
            best_model <- which.max(sapply(model_test_all_times, function(x) x[[2]]))
            model_out <- model_test_all_times[best_model][[1]]

            if(verbose) {
                cat(paste0("    best split time found at ", ten_times[best_model], "\n"))
                cat(paste0("Done. Log-likelihood = ", round(model_out$value, digit = 3), ".\n"))
            }

        } else {
            
            if(verbose) {
                cat(paste0("Fewer than 30 samples are available - time split models was not run!\n"))
            } else {
                warning(paste0("Fewer than 30 samples are available - time split models was not run!\n"))
            }

            run_time_split <- FALSE
        }
    
    } else {

        if(verbose) cat(paste0("Running ", paste(model_type, collapse = ":"), " model..."))
        model_out <- model.test.lik(model_test_input, model.type.in = model_type, time.split, control.list, fixed.optima = fixed.optima)
        if(verbose) cat(paste0("Done. Log-likelihood = ", round(model_out$value, digit = 3), "\n"))
    }

    return(model_out)
}
