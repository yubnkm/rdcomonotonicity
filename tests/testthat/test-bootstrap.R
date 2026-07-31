test_that("bootstrap returns objects with expected dimensions", {
    dat <- make_small_dgp()

    fit <- RDD_extrapolate_CV_band(
        Y = dat$Y,
        X = dat$X,
        D = dat$D,
        kernel = "gaussian",
        bands = c(0.4, 0.4),
        num_folds = 2,
        order = 1,
        batch_size = 100L
    )

    set.seed(321)

    boot <- suppressMessages(
        RDD_extrapolate_bootstrap(
            result = fit,
            Y = dat$Y,
            X = dat$X,
            D = dat$D,
            kernel = "gaussian",
            num_folds = 2,
            order = 1,
            batch_size = 100L,
            points = 7L,
            bootstrap_iter = 2L,
            parallel = FALSE
        )
    )

    expect_equal(dim(boot$q0_draws), c(7L, 2L))
    expect_equal(dim(boot$q1_draws), c(7L, 2L))

    expect_equal(
        dim(boot$y0_draws),
        c(nrow(dat$X), 2L)
    )

    expect_equal(
        dim(boot$y1_draws),
        c(nrow(dat$X), 2L)
    )

    expect_equal(
        dim(boot$multipliers),
        c(nrow(dat$X), 2L)
    )

    expect_equal(boot$bootstrap_iter, 2L)
    expect_equal(boot$points, 7L)
})