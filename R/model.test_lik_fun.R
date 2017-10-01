#~~~~~~~~~~~~~~~~~~~
# Likelihood estimation functions
#~~~~~~~~~~~~~~~~~~~

## Model test likelihood
model.test.lik <- function(model_test_input, model.type.in, time.split, control.list, fixed.optima) {

    ## Setting up the model (all) parameters
    pooled_variance <- pooled.variance(data.model.test)
    half_life <- (model_test_input$subsamples[length(model_test_input$subsamples)] - model_test_input$subsamples[1]) / 4
    sts_params <- stasis.parameters(model_test_input, pooled_variance)
    p <- c()
    p[1] <- model_test_input$central_tendency[1]
    p[2] <- bm.parameters(model_test_input, pooled_variance)
    p[3] <- log(2) / half_life
    p[c(4, 9, 10)] <- model_test_input$central_tendency[length(model_test_input$central_tendency)]
    p[c(5, 11, 12)] <- sts_params[2]
    p[6] <- sts_params[1]
    p[7] <- trend.parameters(model_test_input, pooled_variance)[2]
    p[8] <- eb.parameters(model_test_input, pooled_variance)[2]
    if(is.null(control.list$ndeps)) {
        control.list$ndeps <- abs(p / 1000000)
        control.list$ndeps[control.list$ndeps == 0] <- 1e-08
    }
    lower_bounds <- c(NA, 1e-8, 1e-8, NA, NA, 1e-8, -100, -100, NA, NA, NA, NA)
    upper_bounds <- c(NA, 100, 100, NA, NA, 20, 100, -1e-8, NA, NA, NA, NA)
    
    ## Running the model
    model_output <- stats::optim(par = p, fn = opt.mode,  method = "L", control = control.list, model.type.in = model.type.in, time.split = time.split, data.model.test = model_test_input, lower = lower_bounds, upper = upper_bounds, fixed.optima = fixed.optima)
    
    ## Output management
    model_output_pars <- model_output[[1]]
    names(model_output_pars) <- c("ancestral state", "sigma squared", "alpha", "optima.1", "theta.1", "omega", "trend", "eb", "optima.2", "optima.3", "theta.2", "theta.3")
    model_output$par <- get.parameters(model_output_pars, model.type.in, time.split = time.split)
    return(model_output)
}


#~~~~~~~~~~~~~~~~~~~
# Parameter functions
#~~~~~~~~~~~~~~~~~~~

## BM model parameters
bm.parameters <- function (data.model.test, pooled.variance) {
    
    sample.size <- length(data.model.test$central_tendency) - 1
    round.median.sample.size <- round(median(sample.size))
    t.step <- (data.model.test$subsamples[sample.size + 1] - data.model.test$subsamples[1]) / sample.size
    epsilon <- 2 * pooled.variance / round.median.sample.size
    data.model.test.difference <- diff(data.model.test$central_tendency)
    return((1/t.step) * ((1/sample.size) * sum(data.model.test.difference ^ 2) - epsilon))
}

## Stasis model parameters
stasis.parameters <- function (data.model.test, pooled.variance) {
    
    sample.size <- length(data.model.test$central_tendency)
    var.pooled <- pooled.variance
    theta <- mean(data.model.test$central_tendency[2:sample.size])
    omega <- var(data.model.test$central_tendency[2:sample.size]) - var.pooled / median(data.model.test$sample_size)
    return(c(omega, theta))
}

## EB model parameters
eb.parameters <- function (data.model.test, pooled.variance) {
    
    sample.size <- length(data.model.test$central_tendency) - 1
    t.step <- (data.model.test$subsamples[sample.size + 1] - data.model.test$subsamples[1]) / sample.size
    epsilon <- 2 * pooled.variance / round(median(data.model.test$sample_size))
    data.model.test.difference <- diff(data.model.test$central_tendency)
    mean.difference <- mean(data.model.test.difference)
    sigma.squared.step <- (1 / t.step) * ((1 / sample.size) * sum(mean.difference ^ 2) - epsilon)
    a <- log(1e-5) / max(data.model.test$subsamples) * (1/2)
    return(c(sigma.squared.step, a))
}

## Trend model parameters
trend.parameters <- function (data.model.test, pooled.variance) {
    
    sample.size <- length(data.model.test$central_tendency) - 1
    t.step <- (data.model.test$subsamples[sample.size + 1] - data.model.test$subsamples[1]) / sample.size
    epsilon <- 2 * pooled.variance / round(median(data.model.test$sample_size))
    data.model.test.difference <- diff(data.model.test$central_tendency)
    mean.difference <- mean(data.model.test.difference)
    trend.param <- mean.difference / t.step 
    sigma.squared.step <- (1 / t.step) * ((1 / sample.size) * sum(data.model.test.difference ^ 2) - mean.difference^2 - epsilon)
    return(c(sigma.squared.step, trend.param))
}

#~~~~~~~~~~~~~~~~~~~
# Optimisation
#~~~~~~~~~~~~~~~~~~~

## Optimise model functions
opt.mode <- function(p, model.type.in, time.split, data.model.test, fixed.optima)  {
     
    ## Select the position of the time split 
    if(!is.null(time.split)) {
        time.split <- sort(sapply(time.split, function(u) which.min(abs(u - data.model.test$subsamples))))  #TG: removed the rev() here, I'm a bit confused on which side the subsamples should be (decreasing?) - I suggest we fix that once in for all in model.test_fun.R::select.model.list().
    }
    
    ## Set up the parameters 
    total_n <- length(data.model.test$subsamples)
    sample_time <- 1:total_n
    split_here_vcv <-c(1, time.split)
    split_here_2_vcv <-c(time.split - 1, total_n)

    ## Check if the model is multi.OU
    if(any(model.type.in %in% "multi.OU")) {
        split_here_vcv <- split_here_2_vcv <- NULL
        ou_mean <- c(1, time.split, max(data.model.test$subsamples))
        split_here_vcv <- c(1, split_here_vcv)
        split_here_2_vcv <- c(split_here_2_vcv, max(data.model.test$subsamples))
    }

    ## Get the total variance covariance matrix for all time points
    total_VCV <- matrix(0, nrow = total_n, ncol = total_n)

    total_mean <- c()
    optima.level.ou <- optima.level.stasis <- 1
    model.anc <- model.alpha <- NULL
    time_int <- 1 #TG: is this the time interval? Is it one by default? What if people use models with non constant time slices/bins?
    
    for(rec_model in 1:length(model.type.in)) {

        time_x <- time_int

        if(time_x == 1) {
            
            ## Get only the concerned time slices #TG: is that what it's doing?
            data_model_test_int <- lapply(data.model.test, function(k) k[sort(sample_time[split_here_vcv[time_x] : split_here_2_vcv[time_x]])])

            ## Estimate the variance covariance
            output_vcv <- est.VCV(p, data_model_test_int, model.type = model.type.in[time_x])
            
            ## Estimate the mean
            output_mean <- est.mean(p, data_model_test_int, model.type = model.type.in[time_x], optima.level.ou = optima.level.ou, optima.level.stasis = optima.level.stasis, fixed.optima = fixed.optima, est.anc = TRUE, split.time = ou_mean)

            if(model.type.in[1] == "BM") {
                est.anc <- FALSE
                model.anc <- p[1]
                time_int <- time_x + 1
            }
            
            if(model.type.in[1] == "OU") {
                optima.level.ou <- optima.level.ou + 1
                est.anc <- FALSE
                model.anc <- tail(output_mean, 1)
                time_int <- time_x + 1
            }
            
            if(model.type.in[1] == "Stasis") {
                optima.level.stasis <- optima.level.stasis + 1
                model.anc <- p[5]
                est.anc <- FALSE
                time_int <- time_x + 1
            }
            
            if(model.type.in[1] == "Trend" || model.type.in[1] == "EB") {
                model.anc <- tail(output_mean, 1)
                est.anc <- FALSE
                time_int <- time_x + 1
            }
                
    } else {
        
        data_model_test_int <- lapply(data.model.test, function(k) k[sort(sample_time[split_here_vcv[time_x] : (split_here_2_vcv[time_x])] )])    
        
        time.out <- data_model_test_int[[4]]
        time.out.diff <- diff(time.out[1:2])
        time.out.2 <- time.out - (min(time.out) - time.out.diff)
        data_model_test_int$subsamples <- time.out.2
            
        output_vcv <- est.VCV(p, data_model_test_int, model.type=model.type.in[time_x])

        output_mean <- est.mean(p, data.model.test.in=data_model_test_int, model.type=model.type.in[time_x], optima.level.ou= optima.level.ou, optima.level.stasis= optima.level.stasis, fixed.optima=fixed.optima, est.anc=est.anc, model.anc=model.anc, split.time=NULL)
                
        if(model.type.in[time_x] == "BM") {
            est.anc <- FALSE
            model.anc <- p[1]
            time_int <- time_x + 1
        }
        
        if(model.type.in[time_x] == "OU") {
            optima.level.ou <- optima.level.ou + 1
            est.anc <- FALSE
            model.anc <- p[1]
            time_int <- time_x + 1
        }
        
        if(model.type.in[time_x] == "Stasis") {
            optima.level.stasis <- optima.level.stasis + 1
            model.anc <- p[5]
            est.anc <- FALSE
            time_int <- time_x + 1
        }
        
        if(model.type.in[time_x] == "Trend") {
            model.anc <- tail(output_mean, 1)
            est.anc <- FALSE
            time_int <- time_x + 1
        }
        
        if(model.type.in[time_x] == "EB") {
            model.anc <- tail(output_mean, 1)
            est.anc <- FALSE
            time_int <- time_x + 1
        }
            
    }

    total_VCV[split_here_vcv[time_x] : (split_here_2_vcv[time_x]), split_here_vcv[time_x] : (split_here_2_vcv[time_x]) ] <- output_vcv
    total_mean <-  c(total_mean, output_mean)
    }
    
    return(mnormt::dmnorm(t(data.model.test$central_tendency), mean =  total_mean, varcov = total_VCV, log = TRUE))
}

## Variance co-variance estimator
est.VCV <- function(p, data.model.test, model.type, fixed.optima = FALSE, est.anc = TRUE, model.anc) {
    
    if(model.type == "BM" | model.type == "Trend")  {
    
        sigma.squared <- p[2]
        VCV <- sigma.squared * outer(data.model.test$subsamples, data.model.test$subsamples, FUN = pmin)
        diag(VCV) <- diag(VCV) + data.model.test$variance / data.model.test$sample_size
        return(VCV)
    }
    
    if(model.type == "OU" || model.type == "multi.OU" )  {
        
        alpha <- p[3]
        sigma.squared <- p[2]
        VCV <- outer(data.model.test$subsamples, data.model.test$subsamples, function(x, y) abs(x - y))
        VCV <- exp(-alpha * VCV)
        VCVd <- (sigma.squared / (2 * alpha)) * (1 - exp(-2 * alpha * data.model.test$subsamples))
        VCV_two <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV_two
        diag(VCV) <- VCVd + data.model.test$variance / data.model.test$sample_size
        return(VCV)
    }    

    if(model.type == "Stasis")  {
        
        omega <- p[6]
        VCV <- diag(omega + data.model.test$variance / data.model.test$sample_size)
        return(VCV)    
    }
    
    if(model.type == "EB")  {
        
        sigma.squared <- p[2]
        r.rate <- p[8]
        time.out.eb <- data.model.test$subsamples
        VCV <- outer(sigma.squared * ((exp(r.rate * time.out.eb) - 1) / r.rate), sigma.squared * ((exp(r.rate * time.out.eb) - 1) / r.rate), FUN = pmin)    
        diag(VCV) <- diag(VCV) + data.model.test$variance / data.model.test$sample_size
        return(VCV)
    }
}

## Mean estimating functions
est.mean <- function(p, data.model.test.in, model.type, optima.level.ou, optima.level.stasis, fixed.optima = FALSE, est.anc = TRUE, model.anc, split.time) {
    
    if(model.type == "BM" || model.type == "EB") {
    
        if(est.anc) {
            anc.state <- p[1]
        } else {
            anc.state <- model.anc
        }
        sample.size <- length(data.model.test.in[[1]])
        return(rep(anc.state, sample.size))
    }
    
    if(model.type == "OU") {
        
        if(est.anc) {
            anc.state <- p[1]
        } else {
            anc.state <- model.anc
        }        
        alpha <- p[3]
                
            if (fixed.optima == FALSE) {
                if(optima.level.ou == 1) optima <- p[4]
                if(optima.level.ou == 2) optima <- p[9]
                if(optima.level.ou == 3) optima <- p[10]
        } else {
            optima <- anc.state
        }
        
        mean.ou <- optima * (1 - exp(-alpha * data.model.test.in$subsamples)) + anc.state * exp(-alpha * data.model.test.in$subsamples)
        return(mean.ou)
    }
    
    if(model.type == "multi.OU") {
        
        if(est.anc) {
            anc.state <- p[1]
        } else {
            anc.state <- model.anc
        }
        
        alpha <- p[3]
        optima <- p[c(4, 9, 10)]

        ## Selecting the splits start and end
        all.splits <- length(split.time)
        start.split <- split.time[-all.splits]
        start.split[-1]  <- start.split[-1] + 1
        end.split <- split.time[-1]
        take.away <- start.split[1] - 1
        start.split <- start.split - take.away
        end.split <- end.split - take.away

        ## Number of optima
        n.optima <- all.splits - 1
          
        ## OU mean function
        ou.mean.fun <- function (anc.state, optima, alpha, time) {
            return(optima * (1 - exp(-alpha * time)) + anc.state * exp(-alpha * time))
        }

        ou_mean_splits <- sapply(1:n.optima, function(x) ou.mean.fun(anc.state, optima[x], alpha, data.model.test.in$subsamples))

        mean.ou <- c()

        for(one_optima in 1:n.optima) {
            mean.ou <- c(mean.ou, ou_mean_splits[start.split[one_optima] : end.split[one_optima] , one_optima])
        }

        return(mean.ou)
    }    
    
    if(model.type == "Stasis") {
        
        if(optima.level.stasis == 1) theta <- p[5]
        if(optima.level.stasis  == 2) theta <- p[11]
        if(optima.level.stasis  == 3) theta <- p[12]
        sample.size <- length(data.model.test.in[[1]]) 
        return(rep(theta,  sample.size)) 
    }
    
    if(model.type == "Trend") {
        
        if(est.anc) {
            anc.state <- p[1]
        } else {
            anc.state <- model.anc
        }

        trend.param <- p[7]
        sample.size <- length(data.model.test.in[[1]])
        mean.trend <- anc.state + trend.param * data.model.test.in$subsamples
    }
}

#~~~~~~~~~~~~~~~~~~~
# Output
#~~~~~~~~~~~~~~~~~~~

## Getting parameters function
get.parameters <- function(model.output.pars, models, time.split) {

    optima.level <- stasis.level <- 1

    n.models <- length(models)
    first.model <- models[1]

    if(first.model == "BM") {
        parameters.out <- model.output.pars[c(1,2)]
        }
    if(first.model == "OU") {
        parameters.out <- model.output.pars[c(1:4)]
        optima.level <- optima.level + 1
        }
    if(first.model == "Trend") {
        parameters.out <- model.output.pars[c(1:2, 7)]
        }
    if(first.model == "EB") {
        parameters.out <- model.output.pars[c(1:2, 8)]
        }
    if(first.model == "Stasis") {
        parameters.out <- model.output.pars[c(5,6)]
        stasis.level <- stasis.level + 1
        }
    if(first.model == "multi.OU") {
        if(length(time.split) == 1) parameters.out <- model.output.pars[c(1:4, 9)]
        if(length(time.split) == 2) parameters.out <- model.output.pars[c(1:4, 9:10)]
        if(length(time.split) == 3) parameters.out <- model.output.pars[c(1:4, 9:11)]
    }


if(n.models > 1) {
    
    for(y in 2:n.models) {
        second.model <- models[y]
    
        if(second.model == "BM") {
            parameters.out <- c(parameters.out, model.output.pars[2])
            }
        if(second.model == "OU") {
            opt <- c(4, 9:10)[optima.level]
            parameters.out <- c(parameters.out, model.output.pars[c(2:3, opt)])
            optima.level <- optima.level + 1
            }
        
        if(second.model == "Trend") {
            parameters.out <- c(parameters.out, model.output.pars[c(2, 7)])
            }
        
        if(second.model == "EB") {
            parameters.out <- c(parameters.out, model.output.pars[c(2, 8)])
            }
        
        if(second.model == "Stasis") {
            stasis.opt <- c(5, 11:12)[stasis.level]
            parameters.out <- c(parameters.out, model.output.pars[c(6, stasis.opt)])
            stasis.level <- stasis.level + 1
            
            }
        }
        
        too.many <- duplicated(names(parameters.out))
        if(any(too.many)) parameters.out <- parameters.out[-which(too.many)]
    }
    
    return(parameters.out)
}