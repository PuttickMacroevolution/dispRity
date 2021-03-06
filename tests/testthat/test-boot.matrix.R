## TESTING boot.matrix

context("boot.dispRity")

## Loading the data
load("test_data.Rda")
data <- test_data$ord_data_tips


## Internal functions tests
test_that("internal: bootstrap methods", {
    ## Full bootstrap selects three value
    expect_equal(
        length(boot.full(seq(1:20), 3))
        ,3)
    ## Full bootstrap selects more values than the sequence length
    expect_equal(
        length(boot.full(seq(1:20), 40))
        ,40)

    ## Single bootstrap doesn't work with less than two values
    expect_error(boot.single(seq(1:20), 1))
    ## Single boostrap selects three value
    expect_equal(
        length(boot.single(seq(1:20), 3))
        ,3)
})

## Internal functions tests
test_that("internal: bootstrap replicates", {
    data(disparity)
    subsamples <- disparity$subsamples[[1]]

    ## Select rarefactions
    expect_equal(
        select.rarefaction(subsamples, 8)
        ,list(18, 8))
    expect_equal(
        select.rarefaction(subsamples, c(3,5))
        ,list(18, 3, 5))

    ## One bootstrap replicate
    set.seed(1)
    test_silent <- replicate.bootstraps.silent(6, 5, subsamples, boot.full)
    set.seed(1)
    expect_message(test_verbose <- replicate.bootstraps.verbose(6, 5, subsamples, boot.full))

    ## Both are the same!
    expect_is(test_silent, "matrix")
    expect_is(test_verbose, "matrix")
    expect_equal(dim(test_silent), c(6,5))
    expect_equal(dim(test_verbose), c(6,5))
    expect_true(all(as.vector(test_silent) == as.vector(test_verbose)))

    ## Bootstrap replicates wrapper
    test_boot <- bootstrap.wrapper(subsamples, bootstraps = 6, rarefaction = c(3,5), boot.type.fun = boot.full, verbose = FALSE)
    expect_is(test_boot, "list")
    expect_equal(length(test_boot), 3)
    expect_equal(unlist(lapply(test_boot, dim)), c(18, 6, 3, 6, 5, 6))
})

## Testing boot.matrix with a single matrix input
bootstraps = 5
rarefaction = FALSE
boot.type = "full"
## Sanitizing
test_that("Sanitizing works correctly", {
    expect_error(
        boot.matrix(data = "a", bootstraps, rarefaction, dimensions = FALSE, verbose = FALSE, boot.type = "full")
        )
    expect_error(
        boot.matrix(data, bootstraps = FALSE, rarefaction, dimensions = FALSE, verbose = FALSE, boot.type = "full")
        )
    expect_error(
        boot.matrix(data, bootstraps = "a", rarefaction, dimensions = FALSE, verbose = FALSE, boot.type = "full")
        )
    expect_error(
        boot.matrix(data, bootstraps, rarefaction = "a", dimensions = FALSE, verbose = FALSE, boot.type = "full")
        )
    expect_error(
        boot.matrix(data, bootstraps, rarefaction, dimensions = -1, verbose = FALSE, boot.type = "full")
        )
    expect_error(
        boot.matrix(data, bootstraps, rarefaction, dimensions = FALSE, verbose = 8, boot.type = "full")
        )
    expect_error(
        boot.matrix(data, bootstraps, rarefaction, dimensions = FALSE, verbose = FALSE, boot.type = "rangers")
        )
    expect_error(
        boot.matrix(data, bootstraps, rarefaction, dimensions = FALSE, verbose = FALSE, boot.type = 2)
        )
    expect_error(
        boot.matrix(data, bootstraps, rarefaction, dimensions = FALSE, verbose = FALSE, boot.type = "full", parallel = TRUE)
        )
})

## No bootstrap (is equal to the matrix)
test_that("No bootstraps", {
    test <- boot.matrix(data, bootstraps = 0)
    expect_is(
        test
        , "dispRity")
    expect_equal(
        length(test)
        , 3)
    expect_equal(
        as.vector(test$subsamples[[1]][[1]])
        , seq(1:nrow(test$matrix)))
    expect_equal(
        test$call$dimensions
        , ncol(test$matrix))
    expect_equal(
        length(test$subsamples[[1]])
        ,1)
    expect_equal(
        length(test$subsamples)
        ,1)
    expect_is(
        test$subsamples[[1]][[1]]
        ,"matrix")
})

## No bootstrap but remove dimensions
test_that("Remove dimensions", {
    expect_equal(
        boot.matrix(data, bootstraps = 0, dimensions = 0.5)$call$dimensions
        ,24)
    expect_equal(
        boot.matrix(data, bootstraps = 0, dimensions = 24)$call$dimensions
        ,24)
})


## Bootstraps = 5
test_that("5 bootstraps", {
    data <- test_data$ord_data_tips
    test <- boot.matrix(data, bootstraps = 5)
    expect_is(
        test
        , "dispRity")
    expect_equal(
        length(test)
        , 3)
    expect_equal(
        test$call$dimensions
        , ncol(test$matrix))
    expect_equal(
        dim(test$subsamples[[1]][[1]])
        ,c(50,1))
    expect_equal(
        dim(test$subsamples[[1]][[2]])
        ,c(50,5))
    expect_equal(
        length(test$subsamples[[1]])
        ,2)
})

## Bootstraps = 5 + Rarefaction = 5
test_that("5 bootstraps, rarefaction = 5", {
    test <- boot.matrix(data, bootstraps = 5, rarefaction = 5)
    expect_is(
        test
        , "dispRity")
    expect_equal(
        length(test)
        , 3)
    expect_equal(
        test$call$dimensions
        , ncol(test$matrix))
    expect_equal(
        dim(test$subsamples[[1]][[2]])
        ,c(50,5))
    expect_equal(
        dim(test$subsamples[[1]][[3]])
        ,c(5,5))
})

## Bootstraps = 5 + Rarefaction = TRUE
test_that("5 bootstraps, rarefaction = TRUE", {
    test <- boot.matrix(data, bootstraps = 5, rarefaction = TRUE)
    expect_equal(
        length(test$subsamples[[1]])
        , 49)
    for(rare in 2:49) {
        expect_equal(
            dim(test$subsamples[[1]][[rare]])
            ,c(50-(rare-2),5))
    }
})

## Bootstraps = 5 + Rarefaction = c(5,6) + boot.type
test_that("5 bootstraps, rarefaction = 5,6, boot type", {
    test <- boot.matrix(data, bootstraps = 5, rarefaction = c(5,6), boot.type = "single")
    expect_equal(
        test$call$bootstrap[[1]]
        , 5)
    expect_equal(
        test$call$bootstrap[[2]]
        , "single")
    expect_equal(
        test$call$bootstrap[[3]]
        , c(5,6))
})


## Bootstraps = 5 + Rarefaction = c(5,6) + subsamples
test_that("5 bootstraps, rarefaction = 5,6, subsamples", {
    ordinated_matrix <- matrix(data = rnorm(90), nrow = 10, ncol = 9, dimnames = list(letters[1:10]))
    groups <- as.data.frame(matrix(data = c(rep(1,5), rep(2,5)), nrow = 10, ncol = 1, dimnames = list(letters[1:10])))
    matrix_list <- custom.subsamples(ordinated_matrix, groups)
    test <- boot.matrix(matrix_list, bootstraps = 2, rarefaction = c(4,3))
    expect_is(
        test
        , "dispRity")
    expect_equal(
        length(test)
        , 3)
    expect_equal(
        test$call$dimensions
        , ncol(test$matrix))
    expect_equal(
        length(test$subsamples)
        ,2)
    expect_equal(
        length(test$subsamples[[1]])
        ,4)
    expect_equal(
        length(test$subsamples[[2]])
        ,4)
    expect_equal(
        dim(test$subsamples[[2]][[2]])
        ,c(nrow(test$subsamples[[2]]$elements), 2))
})

