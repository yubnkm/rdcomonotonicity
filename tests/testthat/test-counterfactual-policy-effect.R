test_that("counterfactual policy effect implements the policy formula", {
    result <- make_mock_result()
    bootstrap_result <- make_mock_bootstrap_result()

    X <- cbind(
        X1 = seq_len(4L),
        X2 = seq_len(4L)
    )

    D <- c(0, 0, 1, 1)

    # Observed outcomes equal the factual conditional means.
    Y <- c(0.2, 0.4, 0.7, 0.9)

    always_treat <- function(X) {
        rep(1, nrow(X))
    }

    out <- counterfactual_policy_effect(
        result = result,
        bootstrap_result = bootstrap_result,
        policy = always_treat,
        Y = Y,
        X = X,
        D = D,
        level = 0.9
    )

    # Observations 1 and 2 each gain 0.1.
    # There are three observations with S = 1.
    expect_equal(out$estimate, 0.2 / 3)
    expect_equal(out$num_affected, 2)
    expect_equal(out$affected_share, 2 / 3)

    # Bootstrap draws are identical to the original estimates here.
    expect_equal(out$conf_rad, 0)
    expect_equal(out$conf_low, out$estimate)
    expect_equal(out$conf_high, out$estimate)
})


test_that("policy effect can be calculated without bootstrap results", {
    result <- make_mock_result()

    X <- cbind(
        X1 = seq_len(4L),
        X2 = seq_len(4L)
    )

    D <- c(0, 0, 1, 1)
    Y <- c(0.2, 0.4, 0.7, 0.9)

    always_treat <- function(X) {
        rep(1, nrow(X))
    }

    out <- counterfactual_policy_effect(
        result = result,
        bootstrap_result = NULL,
        policy = always_treat,
        Y = Y,
        X = X,
        D = D
    )

    expect_equal(out$estimate, 0.2 / 3)
    expect_null(out$bootstrap_estimates)
    expect_true(is.na(out$conf_low))
    expect_true(is.na(out$conf_high))
})


test_that("invalid policy probabilities are rejected", {
    result <- make_mock_result()

    X <- cbind(
        X1 = seq_len(4L),
        X2 = seq_len(4L)
    )

    D <- c(0, 0, 1, 1)
    Y <- c(0.2, 0.4, 0.7, 0.9)

    invalid_policy <- function(X) {
        rep(1.5, nrow(X))
    }

    expect_error(
        counterfactual_policy_effect(
            result = result,
            policy = invalid_policy,
            Y = Y,
            X = X,
            D = D
        ),
        "between 0 and 1"
    )
})