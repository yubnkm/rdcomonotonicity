test_that("plot_idf_region returns a ggplot object", {
    result <- make_mock_result()

    X <- cbind(
        X1 = c(0.1, 0.3, 0.6, 0.8),
        X2 = c(0.2, 0.5, 0.4, 0.9)
    )

    treatment_rule <- function(X) {
        X[, 2L] < 0.7 - 0.4 * X[, 1L]
    }

    p <- plot_idf_region(
        X = X,
        result = result,
        threshold_function = treatment_rule
    )

    expect_s3_class(p, "ggplot")
})


test_that("plot_idf_region requires exactly two covariates", {
    result <- make_mock_result()

    X <- matrix(
        seq_len(12L),
        nrow = 4L,
        ncol = 3L
    )

    treatment_rule <- function(X) {
        rep(1, nrow(X))
    }

    expect_error(
        plot_idf_region(
            X = X,
            result = result,
            threshold_function = treatment_rule
        ),
        "exactly two"
    )
})