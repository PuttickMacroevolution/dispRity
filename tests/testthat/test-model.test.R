#TESTING null.test

context("model.test")

## Select model data
load("model_test_data.Rda")
data <- model_test_data

test_that("select.model.list internal", {
    model_input <- select.model.list(data)

    expect_is(model_input, "list")
    expect_equal(names(model_input), c("central_tendency", "variance", "sample_size", "subsamples"))
    expect_equal(unique(unlist(lapply(model_input, class))), "numeric")

    ## MISSING OPTIONAL ARGUMENTS TEST!
})

