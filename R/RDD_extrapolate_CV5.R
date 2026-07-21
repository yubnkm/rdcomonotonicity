RDD_extrapolate_CV5 <- function(
    Y,
    X,
    D,
    kernel = "gaussian",
    bands = seq(0.01, 1, length.out = 100),
    weights = rep(1, length(Y)),
    num_folds = 5,
    order = 1
) {

    Y <- as.numeric(Y)
    D <- as.numeric(D)
    weights <- as.numeric(weights)
    X <- as.matrix(X)
    p <- ncol(X)

    # Choose kernel function
    kernel <- match.arg(
        kernel,
        choices = c("gaussian", "uniform", "triangular")
        )
    if (kernel == "gaussian") {
        quart <- 2 * 0.67448
    } else if (kernel == "uniform") {
        quart <- 2 * 0.5
    } else if (kernel == "triangular") {
        quart <- 2 * (1 - 1/sqrt(2))
    }
    

    # Remove missing entries
    not_missing <- complete.cases(cbind(Y, X, D, weights))
    Y <- Y[not_missing]
    X <- X[not_missing, , drop = FALSE]
    D <- D[not_missing]
    weights <- weights[not_missing]

    # Split into d and 1-d groups
    X0 <- X[D == 0, , drop = FALSE]
    X1 <- X[D == 1, , drop = FALSE]
    Y0 <- Y[D == 0]
    Y1 <- Y[D == 1]
    weights0 <- weights[D == 0]
    weights1 <- weights[D == 1]

    fit0 <- local_poly3(
        Y = Y0,
        X = X0,
        kernel = kernel,
        bands = bands,
        weights = weights0,
        num_folds = num_folds,
        order = order
    )

    fit1 <- local_poly3(
        Y = Y1,
        X = X1,
        kernel = kernel,
        bands = bands,
        weights = weights1,
        num_folds = num_folds,
        order = order
    )

    g0_raw <- fit0$predict_func
    g1_raw <- fit1$predict_func
    band0 <- fit0$optimal_band
    band1 <- fit1$optimal_band

    # Construct q1
    nn1 <- FNN::get.knnx(data = X0, query = X1, k = 1)
    W1 <- nn1$nn.dist < quart * band0
    g0 <- g0_raw(X1, X0[nn1$nn.index, , drop = FALSE])

    y1_max <- max(g0[W1])
    y1_min <- min(g0[W1])

    q1_fit <- local_poly3(
        Y = Y1,
        X = matrix(g0, ncol = 1L),
        kernel = kernel,
        bands = bands,
        weights = weights1 * as.numeric(W1),
        num_folds = num_folds,
        order = order
    )
    q1_raw <- q1_fit$predict_func
    q1 <- function(x) {
        x <- matrix(as.numeric(x), ncol = 1L)
        as.numeric(q1_raw(x,x))
    }
    g0_func <- function(x) {
        x <- matrix(as.numeric(x), ncol = p)
        as.numeric(g0_raw(x,x))
    }
    g1_imputed_func <- function(x) q1(g0_func(x))

    # Construct q0
    nn0 <- FNN::get.knnx(data = X1, query = X0, k = 1)
    W0 <- nn0$nn.dist < quart * band1
    g1 <- g1_raw(X0, X1[nn0$nn.index, , drop = FALSE])

    y0_max <- max(g1[W0])
    y0_min <- min(g1[W0])

    q0_fit <- local_poly3(
        Y = Y0,
        X = matrix(g1, ncol = 1L),
        kernel = kernel,
        bands = bands,
        weights = weights0 * as.numeric(W0),
        num_folds = num_folds,
        order = order
    )
    q0_raw <- q0_fit$predict_func
    q0 <- function(x) {
        x <- matrix(as.numeric(x), ncol = 1L)
        as.numeric(q0_raw(x,x))
    }
    g1_func <- function(x) {
        x <- matrix(as.numeric(x), ncol = p)
        as.numeric(g1_raw(x,x))
    }
    g0_imputed_func <- function(x) q0(g1_func(x))

    plotting1 <- cbind(g1[W0], Y0[W0])
    plotting0 <- cbind(g0[W1], Y1[W1])

    y0 <- rep(NA_real_, length(Y))
    y1 <- rep(NA_real_, length(Y))

    y0[D == 0] <- as.numeric(g0_raw(X0, X0))
    y1[D == 1] <- as.numeric(g1_raw(X1, X1))

    # Indicator for observations whose counterfactual mean is supported (S)
    # S <- D * as.numeric( (y0_min < y1) & (y1 < y0_max) ) +
    #     (1-D) * as.numeric( (y1_min < y0) & (y0 < y1_max) )
    S <- integer(length(Y))
    S[D == 1] <- as.integer(y1[D == 1] > y0_min & y1[D == 1] < y0_max)
    S[D == 0] <- as.integer(y0[D == 0] > y1_min & y0[D == 0] < y1_max)
    
    y1[S == 1 & D == 0] <- q1(y0[S == 1 & D == 0])
    y0[S == 1 & D == 1] <- q0(y1[S == 1 & D == 1])

    return(
        list(
            g1_imputed_func = g1_imputed_func,
            g0_imputed_func = g0_imputed_func,
            g1_func = g1_func,
            g0_func = g0_func,
            q1 = q1,
            q0 = q0,
            plotting1 = plotting1,
            plotting0 = plotting0,
            y0 = y0,
            y1 = y1,
            S = S,
            band0 = band0,
            band1 = band1,
            y0_min = y0_min,
            y0_max = y0_max,
            y1_min = y1_min,
            y1_max = y1_max
        )
    )
}