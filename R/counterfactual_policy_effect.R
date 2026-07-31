counterfactual_policy_effect <- function(
    result,
    bootstrap_result = NULL,
    policy,
    Y,
    X,
    D,
    level = 0.90
) {

    if (!is.function(policy)) {
        stop("policy must be a function")
    }
    Y <- as.numeric(Y)
    X <- as.matrix(X)
    D <- as.numeric(D)

    if (!is.null(bootstrap_result) && !is.null(bootstrap_result$complete_cases)) {
        complete_cases <- as.logical(bootstrap_result$complete_cases)
    } else {
        complete_cases <- complete.cases(cbind(Y, X, D, estimand_weights))
    }
    Y <- Y[complete_cases]
    X <- X[complete_cases, , drop = FALSE]
    D <- D[complete_cases]

    n <- length(Y)
    y0 <- as.numeric(result$y0)
    y1 <- as.numeric(result$y1)

    S <- as.numeric(result$S == 1)
    identified_n <- sum(S)

    policy_prob <- as.numeric(policy(X))

    switch_to_treated <- (1 - D) * policy_prob
    switch_to_untreated <- D * (1 - policy_prob)
    switch_prob <- switch_to_treated + switch_to_untreated

    # calculate the counterfactual policy effect
    indv_contribution <- S * ( switch_to_treated * (y1 - Y) + switch_to_untreated * (y0 - Y))
    estimate <- sum(indv_contribution) / identified_n

    # # of affected among S = 1
    num_affected <- sum(S * switch_prob)
    affected_share <- num_affected / identified_n

    ## Bootstrap
    bootstrap_estimates <- NULL
    conf_rad <- NA_real_
    conf_low <- NA_real_
    conf_high <- NA_real_

    if(!is.null(bootstrap_result)) {
        y0_draws <- as.matrix(bootstrap_result$y0_draws)
        y1_draws <- as.matrix(bootstrap_result$y1_draws)
        multipliers <- as.matrix(bootstrap_result$multipliers)
        bootstrap_iter <- as.numeric(bootstrap_result$bootstrap_iter)

        bootstrap_estimates <- vapply(
            seq(bootstrap_iter),
            FUN.VALUE = numeric(1L),
            FUN = function(b){
                multiplier_b <- multipliers[, b]

                denom_b <- sum(multiplier_b * S)
                if(!is.finite(denom_b) || denom_b <= 0) return(NA_real_)

                contribution_b <- S * (
                    switch_to_treated * (y1_draws[, b] - Y) +
                    switch_to_untreated * (y0_draws[, b] - Y)
                    )
                
                numer_b <- sum(multiplier_b * contribution_b)

                return(numer_b / denom_b)
            }
        )

        valid_bootstrap_draws <- is.finite(bootstrap_estimates)

        if(!all(valid_bootstrap_draws)) {
            warning(sum(!all(valid_bootstrap_draws)), " non finite bootstrap estimates were removed")
        }

        conf_rad <- stats::quantile(
            abs(bootstrap_estimates[valid_bootstrap_draws] - estimate),
            probs = level,
            names = FALSE,
            na.rm = TRUE
        )
        conf_low <- estimate - conf_rad
        conf_high <- estimate + conf_rad
    }

    return(
        list(
            estimate = estimate,
            conf_low = conf_low,
            conf_high = conf_high,
            conf_rad = conf_rad,
            level = level,

            bootstrap_estimates = bootstrap_estimates,

            num_affected = num_affected,
            affected_share = affected_share,

            policy_prob = policy_prob,
            switch_prob = switch_prob,
            indv_contribution = indv_contribution,

            S = S,
            identified_n = identified_n,
            complete_cases = complete_cases

        )

    )




}