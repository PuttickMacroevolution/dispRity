#TESTING null.test

context("null.test")

# Testing data
data(BeckLee_mat50)
single_disp <- dispRity(BeckLee_mat50, metric = ellipse.volume)
groups <- as.data.frame(matrix(data = c(rep(1, nrow(BeckLee_mat50)/2),
     rep(2, nrow(BeckLee_mat50)/2)), nrow = nrow(BeckLee_mat50), ncol = 1,
     dimnames = list(rownames(BeckLee_mat50))))
multi_disp <- dispRity(boot.matrix(custom.subsamples(BeckLee_mat50, groups), bootstraps = 100), metric = c(sum, variances))


#get.from.call
test_that("get.metric.from.call works", {
    #Errors
    expect_error(
        get.metric.from.call("a")
        )
    #Right outputs
    expect_is(
        get.metric.from.call(single_disp)
        , "function")
    expect_is(
        get.metric.from.call(multi_disp)
        , "list")
    expect_equal(
        unique(unlist(lapply(get.metric.from.call(multi_disp), class)))
        , "function")
})

#make.null.model
test_that("make.null.model works", {
    #Errors
    expect_error(
        make.null.model("a", replicates = 5, null.distrib = rnorm, null.args = NULL, null.cor = NULL, scale = FALSE)
        )
    expect_error(
        make.null.model(single_disp, replicates = -3, null.distrib = rnorm, null.args = NULL, null.cor = NULL, scale = FALSE)
        )
    expect_error(
        make.null.model(single_disp, replicates = 5, null.distrib = "rnorm", null.args = NULL, null.cor = NULL, scale = FALSE)
        )
    expect_error(
        make.null.model(single_disp, replicates = 5, null.distrib = rnorm, null.args = TRUE, null.cor = NULL, scale = FALSE)
        )

    #Right output
    expect_is(
        make.null.model(single_disp, replicates = 5, null.distrib = rnorm, null.args = NULL, null.cor = NULL, scale = FALSE)
        , "numeric")
    expect_equal(
        length(make.null.model(single_disp, replicates = 5, null.distrib = rnorm, null.args = NULL, null.cor = NULL, scale = FALSE))
        , 5)

    #Handling args properly
    my_distributions1 <- unlist(replicate(16, list(rnorm, runif, rnorm), simplify = FALSE), recursive = FALSE)
    my_distributions2 <- unlist(replicate(16, list(rnorm, runif, rgamma), simplify = FALSE), recursive = FALSE)
    set.seed(1)
    my_args <- unlist(replicate(16, list(list(mean = runif(1), sd = runif(1)), list(min = 0.1, max = 0.8), list(shape = runif(1))), simplify = FALSE), recursive = FALSE)
    my_cor.matrix <- matrix(c(unlist(replicate(47, c(1, rep(0, 48)), simplify = FALSE)),1), nrow = 48, ncol = 48, byrow = FALSE)

    set.seed(1)
    test1 <- make.null.model(single_disp, replicates = 5, null.distrib = my_distributions1, null.args = NULL, null.cor = NULL, scale = FALSE)
    set.seed(1)
    test2 <- make.null.model(single_disp, replicates = 5, null.distrib = my_distributions2, null.args = my_args, null.cor = NULL, scale = FALSE)
    set.seed(1)
    test3 <- make.null.model(single_disp, replicates = 5, null.distrib = my_distributions2, null.args = my_args, null.cor = my_cor.matrix, scale = FALSE)
    set.seed(1)
    test4 <- make.null.model(single_disp, replicates = 5, null.distrib = my_distributions1, null.args = NULL, null.cor = NULL, scale = TRUE)
    set.seed(1)
    test5 <- make.null.model(single_disp, replicates = 5, null.distrib = my_distributions2, null.args = my_args, null.cor = my_cor.matrix, scale = TRUE)

    expect_is(
        test1
        , "numeric")
    expect_is(
        test2
        , "numeric")
    expect_is(
        test3
        , "numeric")
    expect_is(
        test4
        , "numeric")
    expect_is(
        test5
        , "numeric")
    expect_equal(
        length(test1)
        , 5)
    expect_equal(
        length(test2)
        , 5)
    expect_equal(
        length(test3)
        , 5)
    expect_equal(
        length(test4)
        , 5)
    expect_equal(
        length(test5)
        , 5)

    expect_equal(
        test2
        , test3)

})

#null.test
test_that("null.test works", {
    #Errors
    expect_error(
        null.test("a", replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "two-sided", scale = FALSE)
        )
    expect_error(
        null.test(single_disp, replicates = "10", null.distrib = NULL, null.args = NULL, alter = "two-sided", scale = FALSE)
        )
    expect_error(
        null.test(single_disp, replicates = 10, null.distrib = rnorm, null.args = "NULL", alter = "two-sided", scale = FALSE)
        )
    expect_error(
        null.test(single_disp, replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "something", scale = FALSE)
        )
    #Right output
    expect_is(
        null.test(single_disp, replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "two-sided", scale = FALSE)
        , "randtest")
    expect_is(
        null.test(single_disp, replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "two-sided", scale = FALSE)
        , "dispRity")
    expect_is(
        null.test(multi_disp, replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "two-sided", scale = FALSE)
        , "randtest")
    expect_is(
        null.test(single_disp, replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "two-sided", scale = FALSE)
        , "dispRity")
    expect_equal(
        unique(unlist(lapply(null.test(multi_disp, replicates = 10, null.distrib = rnorm, null.args = NULL, alter = "two-sided", scale = FALSE), class)))
        , c("randtest", "lightrandtest"))
})