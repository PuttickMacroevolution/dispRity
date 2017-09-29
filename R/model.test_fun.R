## Lapply loop for applying the models
lapply.models <- function(one_model, models, time.split, verbose) {
    
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
        all_times <- max(model_test_input[[4]]) - model_test_input[[4]]
        

        if(length(all_times) > 31) {

            ten_times <- 9 : (length(all_times) - 11)

            cat("Running",  model_type ,"on", length(ten_times), "shift times")
            cat("\n")
            cat("model ", 1, " of ", length(ten_times))
            
            model_test_all_times <- lapply(ten_times, function(x) {
                cat("\r", "model ", x - 8, " of ", length(ten_times))
                model.test.lik(model_test_input, model_type, time.split=x, control.list, fixed.optima=fixed.optima)
            })


            best.model <- which.max(sapply(model_test_all_times, function(x) x[[2]]))
            model.return <- model_test_all_times[best.model][[1]]    

            cat(" best split time found at", ten_times[best.model])
            cat(" log-likelihood: ", model.return$value)
            cat(". Finished.")
            cat("\n")

        } else {
            
            if(verbose) {
                cat(paste0("Fewer than 30 samples are available - time split models was not run!\n"))
            } else {
                warning(paste0("Fewer than 30 samples are available - time split models was not run!"))
            }

            run_time_split <- FALSE
        }
    
    } else {
    
            cat("Running ", model_type, "model")
            model.return <- model.test.lik(model_test_input, model_type.in = model_type, time.split, control.list, fixed.optima=fixed.optima)
            cat(" log-likelihood: ", model.return$value)
            cat(" Finished.")
            cat("\n")
    }
    model.return
}
