test_that("RDD estimator works with supported kernels", {
    dat <- make_small_dgp()

    for (kernel_name in c("gaussian", "uniform", "triangular")) {

    test_that(
        paste("RDD estimator works with", kernel_name, "kernel"),
        {
            dat <- make_small_dgp()

            fit <- RDD_extrapolate_CV_band(
                Y = dat$Y,
                X = dat$X,
                D = dat$D,
                kernel = kernel_name,
                bands = c(0.4, 0.4),
                num_folds = 2,
                order = 1,
                batch_size = 100L
            )

            expect_true(is.function(fit$q0))
            expect_true(is.function(fit$q1))
            expect_true(is.function(fit$g0_func))
            expect_true(is.function(fit$g1_func))

            expect_length(fit$y0, nrow(dat$X))
            expect_length(fit$y1, nrow(dat$X))
            expect_length(fit$S, nrow(dat$X))

            expect_true(all(fit$S %in% c(0L, 1L)))
            expect_true(is.finite(fit$band0))
            expect_true(is.finite(fit$band1))
        }
    )
    }
})


test_that("RDD estimator accepts a data frame", {
    dat <- make_small_dgp()

    fit <- RDD_extrapolate_CV_band(
        Y = dat$Y,
        X = as.data.frame(dat$X),
        D = dat$D,
        kernel = "gaussian",
        bands = c(0.4, 0.4),
        num_folds = 2,
        order = 1
    )

    expect_length(fit$S, nrow(dat$X))
})


test_that("RDD estimator removes incomplete observations", {
    dat <- make_small_dgp()

    dat$Y[1L] <- NA_real_
    dat$X[2L, 1L] <- NA_real_

    fit <- RDD_extrapolate_CV_band(
        Y = dat$Y,
        X = dat$X,
        D = dat$D,
        kernel = "gaussian",
        bands = c(0.4, 0.4),
        num_folds = 2,
        order = 1
    )

    expect_length(fit$S, nrow(dat$X) - 2L)
})


test_that("RDD estimator handles more than two covariates", {
    dat <- make_small_dgp(include_x3 = TRUE)

    fit <- RDD_extrapolate_CV_band(
        Y = dat$Y,
        X = dat$X,
        D = dat$D,
        kernel = "gaussian",
        bands = c(0.5, 0.5),
        num_folds = 2,
        order = 1,
        batch_size = 100L
    )

    expect_length(fit$S, nrow(dat$X))
    expect_true(is.function(fit$q0))
    expect_true(is.function(fit$q1))
})

test_that("bandwidth cross-validation returns selected bandwidths", {
    dat <- make_small_dgp()

    fit <- RDD_extrapolate_CV_band(
        Y = dat$Y,
        X = dat$X,
        D = dat$D,
        kernel = "gaussian",
        bands = c(0.3, 0.4, 0.5),
        num_folds = 2,
        order = 1
    )

    expect_true(fit$band0 %in% c(0.3, 0.4, 0.5))
    expect_true(fit$band1 %in% c(0.3, 0.4, 0.5))
})