#' Estimate the Effect of a Counterfactual Treatment Policy
#'
#' Estimates the mean causal effect of replacing the observed deterministic
#' treatment assignment with a user-specified counterfactual treatment policy.
#' The effect is evaluated over observations whose conditional average treatment
#' effects are identified, as indicated by `result$S == 1`.
#'
#' The counterfactual policy may be deterministic or probabilistic. A
#' deterministic policy returns only zero and one, while a probabilistic policy
#' returns treatment probabilities between zero and one.
#'
#' @param result A fitted object returned by
#'   [RDD_extrapolate_CV_band()]. It must contain `y0`, `y1`, and `S`.
#'
#' @param bootstrap_result An optional bootstrap object returned by
#'   [RDD_extrapolate_bootstrap()]. When supplied, the function calculates
#'   multiplier-bootstrap confidence intervals. The default is `NULL`.
#'
#' @param policy A function that accepts the covariate matrix `X` and returns
#'   one value per observation. Returned values must lie between zero and one
#'   and represent treatment probabilities under the counterfactual policy.
#'
#' @param Y A numeric vector of observed outcomes. 
#' 
#' @param X A numeric matrix or data frame of assignment variables or covariates.
#'   Rows correspond to observations and columns correspond to covariates.
#' 
#' @param D A numeric or integer vector indicating treatment status. Values
#'   should be `1` for treated observations and `0` for untreated observations.
#' 
#' @param level A numeric scalar between zero and one specifying the confidence
#'   level. The default is `0.90`.
#'
#' @details
#' Let `p(X)` denote the treatment probability under the counterfactual policy.
#' For observations currently untreated, the function uses the estimated treated
#' potential outcome when the counterfactual policy assigns treatment. For
#' observations currently treated, it uses the estimated untreated potential
#' outcome when the counterfactual policy removes treatment.
#'
#' The reported estimand is
#'
#' \deqn{
#' \frac{1}{\sum_i S_i}
#' \sum_i S_i \left[
#' (1-D_i)p(X_i)(\widehat{Y}_i(1)-Y_i)
#' +
#' D_i(1-p(X_i))(\widehat{Y}_i(0)-Y_i)
#' \right].
#' }
#'
#' Consequently, the effect is averaged over all observations with `S = 1`,
#' rather than only over observations whose treatment status changes.
#'
#' When `bootstrap_result` is supplied, the confidence radius is the requested
#' quantile of the absolute bootstrap deviations from the original estimate.
#'
#' @return A named list containing:
#'
#' - `estimate`: estimated counterfactual policy effect.
#' - `conf_low` and `conf_high`: confidence interval endpoints.
#' - `conf_rad`: confidence radius.
#' - `level`: requested confidence level.
#' - `bootstrap_estimates`: vector of bootstrap policy-effect estimates, or
#'   `NULL` when no bootstrap object is supplied.
#' - `num_affected`: expected number of observations with `S = 1` whose
#'   treatment status changes.
#' - `affected_share`: `num_affected` divided by the number with `S = 1`.
#' - `policy_prob`: policy treatment probabilities.
#' - `switch_prob`: probability that treatment status differs from observed
#'   treatment status.
#' - `indv_contribution`: observation-level contributions to the numerator.
#' - `S`: identification indicator used in estimation.
#' - `identified_n`: number of observations with `S = 1`.
#' - `complete_cases`: complete-case indicator for the original inputs.

#' @export
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
        complete_cases <- complete.cases(cbind(Y, X, D))
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

    if (length(policy_prob) != n) {
        stop("policy(X) must return one value per row of X.")
    }

    if (anyNA(policy_prob) || any(!is.finite(policy_prob))) {
        stop("policy(X) must return finite, non-missing values.")
    }

    if (any(policy_prob < 0 | policy_prob > 1)) {
        stop("policy(X) must return values between 0 and 1.")
    }

    if (identified_n == 0L) {
        stop("No observations have S = 1; the policy effect is not estimable.")
    }

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
            warning(sum(!valid_bootstrap_draws), " non finite bootstrap estimates were removed")
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