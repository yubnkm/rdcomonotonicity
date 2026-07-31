#' Multiplier Bootstrap for RDD Extrapolation Estimates
#'
#' Computes multiplier-bootstrap draws for the estimated potential-outcome
#' mappings and imputed conditional mean potential outcomes produced by
#' [RDD_extrapolate_CV_band()].
#'
#' Independent standard exponential multipliers are applied to the observation
#' weights. The model is then re-estimated for each bootstrap draw using the
#' bandwidth pair selected in the original fit.
#'
#' @param result A fitted object returned by
#'   [RDD_extrapolate_CV_band()]. It must contain `plotting0`, `plotting1`,
#'   `band0`, and `band1`.
#'
#' @param Y A numeric vector of observed outcomes. 
#' 
#' @param X A numeric matrix or data frame of assignment variables or covariates.
#'   Rows correspond to observations and columns correspond to covariates.
#' 
#' @param D A numeric or integer vector indicating treatment status. Values
#'   should be `1` for treated observations and `0` for untreated observations.
#'
#' @param kernel A character string specifying the kernel used in the local
#'   polynomial regressions. Available choices are `"gaussian"`, `"uniform"`,
#'   and `"triangular"`.
#'
#' @param weights An optional numeric vector of observation weights. The default
#'   assigns weight one to every observation.
#'
#' @param num_folds A positive integer giving the number of folds used for
#'   cross-validation. The default is `5`.
#'
#' @param order A nonnegative integer giving the total order of the local
#'   polynomial. The default `1` corresponds to local linear regression.
#'
#' @param batch_size A positive integer giving the maximum number of prediction
#'   points processed in one batch. Smaller values reduce peak memory usage but
#'   may increase computation time. The default is `1000`.
#'
#' @param points A positive integer giving the number of points in the
#'   evaluation grids for `q0` and `q1`. The default is `100`.
#'
#' @param bootstrap_iter A positive integer giving the number of multiplier
#'   bootstrap draws. The default is `100`.
#'
#' @param parallel A logical value indicating whether the bootstrap draws should
#'   be evaluated using a parallel socket cluster. The default is `FALSE`.
#'
#' @param n_cores An optional positive integer giving the number of parallel
#'   workers. When `parallel = TRUE` and `n_cores = NULL`, the function uses one
#'   fewer than the number of detected physical cores, subject to a minimum of
#'   one and a maximum equal to `bootstrap_iter`.
#'
#' @details
#' For bootstrap draw \eqn{b}, the effective observation weight is
#'
#' \deqn{w_i e_{ib},}
#'
#' where \eqn{w_i} is the original observation weight and \eqn{e_{ib}} is an
#' independent standard exponential multiplier.
#'
#' The function stores draws of the estimated `q0` and `q1` curves and draws of
#' the imputed conditional mean potential outcomes. It does not itself construct
#' confidence intervals; the stored draws are consumed by
#' [plot_q_results()] and [counterfactual_policy_effect()].
#'
#' Set the random seed before calling this function when reproducible bootstrap
#' multipliers are required.
#'
#' @return A named list containing:
#'
#' - `q0_grid` and `q1_grid`: evaluation grids.
#' - `q0_draws` and `q1_draws`: matrices of bootstrap curve evaluations, with
#'   one bootstrap draw per column.
#' - `y0_draws` and `y1_draws`: matrices of bootstrap conditional mean
#'   potential-outcome estimates, with one draw per column.
#' - `multipliers`: matrix of standard exponential bootstrap multipliers.
#' - `complete_cases`: complete-case indicator for the original inputs.
#' - `band0` and `band1`: bandwidths inherited from the original fit.
#' - `bootstrap_iter`: number of bootstrap draws.
#' - `points`: number of curve-evaluation points.
#'
#' @seealso
#' [RDD_extrapolate_CV_band()], [plot_q_results()],
#' [counterfactual_policy_effect()]
#'
#' @export
RDD_extrapolate_bootstrap <- function(
    result,
    Y,
    X,
    D,
    kernel = "gaussian",
    weights = rep(1, length(Y)),
    num_folds = 5,
    order = 1,
    batch_size = 1000L,
    points = 100L,
    bootstrap_iter = 100L,
    parallel = FALSE,
    n_cores = NULL
) {
    Y <- as.numeric(Y)
    X <- as.matrix(X)
    D <- as.numeric(D)

    if (is.null(weights)) {
        weights <- rep(1, length(Y))
    }

    weights <- as.numeric(weights)

    complete_cases <- complete.cases(
        cbind(Y, X, D, weights)
    )

    Y <- Y[complete_cases]
    X <- X[complete_cases, , drop = FALSE]
    D <- D[complete_cases]
    weights <- weights[complete_cases]

    n <- length(Y)

    points <- as.integer(points)
    bootstrap_iter <- as.integer(bootstrap_iter)


    q0_training_data <- as.matrix(result$plotting1)
    q1_training_data <- as.matrix(result$plotting0)

    q0_grid <- seq(
        from = min(q0_training_data[, 1L], na.rm = TRUE),
        to = max(q0_training_data[, 1L], na.rm = TRUE),
        length.out = points
    )

    q1_grid <- seq(
        from = min(q1_training_data[, 1L], na.rm = TRUE),
        to = max(q1_training_data[, 1L], na.rm = TRUE),
        length.out = points
    )

    bootstrap_bands <- c(result$band0, result$band1)

    multipliers <- matrix(
        stats::rexp(
            n = n * bootstrap_iter,
            rate = 1
        ),
        nrow = n,
        ncol = bootstrap_iter
    )

    q0_draws <- matrix(NA_real_, nrow = length(q0_grid), ncol = bootstrap_iter)
    q1_draws <- matrix(NA_real_, nrow = length(q1_grid), ncol = bootstrap_iter)

    y0_draws <- matrix(NA_real_, nrow = n, ncol = bootstrap_iter)
    y1_draws <- matrix(NA_real_, nrow = n, ncol = bootstrap_iter)


        if (parallel) {
            if (is.null(n_cores)) {
                n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
            }
            n_cores <- min(as.integer(n_cores), bootstrap_iter)

            cl <- parallel::makeCluster(n_cores)

            bootstrap_results <- tryCatch(
                {

                    parallel::clusterExport(
                        cl,
                        varlist = c(
                            "Y",
                            "X",
                            "D",
                            "kernel",
                            "bootstrap_bands",
                            "weights",
                            "multipliers",
                            "num_folds",
                            "order",
                            "batch_size",
                            "q0_grid",
                            "q1_grid",
                            "RDD_extrapolate_CV_band",
                            "local_poly3",
                            "predict_in_batches"
                        ),
                        envir = environment()
                    )

                    parallel::parLapply(
                        cl,
                        seq_len(bootstrap_iter),
                        function(b) {

                            bootstrap_weights <-
                                weights * multipliers[, b]

                            bootstrap_fit <-
                                RDD_extrapolate_CV_band(
                                    Y = Y,
                                    X = X,
                                    D = D,
                                    kernel = kernel,
                                    bands = bootstrap_bands,
                                    weights = bootstrap_weights,
                                    num_folds = num_folds,
                                    order = order,
                                    batch_size = batch_size
                                )

                            list(
                                q0 = as.numeric(bootstrap_fit$q0(q0_grid)),
                                q1 = as.numeric(bootstrap_fit$q1(q1_grid)),
                                y0 = as.numeric(bootstrap_fit$y0),
                                y1 = as.numeric(bootstrap_fit$y1)
                            )
                        }
                    )
                },
                finally = {
                    parallel::stopCluster(cl)
                }
            )

            q0_draws <- do.call(cbind, lapply(bootstrap_results, function(x) x$q0))
            q1_draws <- do.call(cbind, lapply(bootstrap_results, function(x) x$q1))
            y0_draws <- do.call(cbind, lapply(bootstrap_results, function(x) x$y0))
            y1_draws <- do.call(cbind, lapply(bootstrap_results, function(x) x$y1))


            message("Completed ", bootstrap_iter, " bootstrap draws using ", n_cores, " cores")

        } else {

            for (b in seq_len(bootstrap_iter)) {

                bootstrap_weights <-
                    weights * multipliers[, b]

                bootstrap_fit <-
                    RDD_extrapolate_CV_band(
                        Y = Y,
                        X = X,
                        D = D,
                        kernel = kernel,
                        bands = bootstrap_bands,
                        weights = bootstrap_weights,
                        num_folds = num_folds,
                        order = order,
                        batch_size = batch_size
                    )

                q0_draws[, b] <- as.numeric(bootstrap_fit$q0(q0_grid))
                q1_draws[, b] <- as.numeric(bootstrap_fit$q1(q1_grid))

                y0_draws[, b] <- as.numeric(bootstrap_fit$y0)
                y1_draws[, b] <- as.numeric(bootstrap_fit$y1)

                message("Completed bootstrap draw ", b, " of ", bootstrap_iter)

            }
        }

        return(
        list(
            q0_grid = q0_grid,
            q1_grid = q1_grid,

            q0_draws = q0_draws,
            q1_draws = q1_draws,

            y0_draws = y0_draws,
            y1_draws = y1_draws,

            multipliers = multipliers,

            complete_cases = complete_cases,

            band0 = result$band0,
            band1 = result$band1,

            bootstrap_iter = bootstrap_iter,
            points = points
        )
    )
}
