#' Extrapolate Potential Outcomes in a Sharp RDD Using Comonotonicity
#'
#' The function first estimates conditional mean outcomes separately for treated
#' and untreated observations using multivariate local polynomial regression. It
#' then identifies observations near the treatment frontier using nearest
#' neighbors and estimates the mappings between treated and untreated conditional
#' mean potential outcomes.
#'
#' @param Y A numeric vector of observed outcomes. Its length must equal the number of rows of `X`.

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
#' @param bands A numeric vector specifying bandwidths.
#'
#'   When `length(bands) == 2`, the first and second values are used as the
#'   first-stage bandwidths for the untreated and treated regressions,
#'   respectively.
#'
#'   When `length(bands) > 2`, the values are treated as candidate bandwidths
#'   and bandwidth selection is performed by cross-validation.
#'
#'   A length-one vector is currently not accepted.
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
#'   may increase computation time.
#'
#' @return A named list containing:
#'
#' - `g1_imputed_func`: a function that imputes the conditional mean treated
#'   potential outcome from the estimated untreated conditional mean.
#'
#' - `g0_imputed_func`: a function that imputes the conditional mean untreated
#'   potential outcome from the estimated treated conditional mean.
#'
#' - `g1_func`: the estimated conditional mean outcome function fitted using the
#'   treated sample.
#'
#' - `g0_func`: the estimated conditional mean outcome function fitted using the
#'   untreated sample.
#'
#' - `q1`: the estimated mapping from untreated to treated conditional mean
#'   potential outcomes.
#'
#' - `q0`: the estimated mapping from treated to untreated conditional mean
#'   potential outcomes.
#'
#' - `plotting1` and `plotting0`: two-column matrices containing the generated
#'   regressors and observed outcomes used in the second-stage local polynomial
#'   regressions near the treatment frontier.
#'
#' - `y0`: estimated conditional mean untreated potential outcomes for the
#'   complete-case observations.
#'
#' - `y1`: estimated conditional mean treated potential outcomes for the
#'   complete-case observations.
#'
#' - `S`: an integer support indicator. A value of `1` indicates that the
#'   observation's factual conditional mean lies within the estimated domain of
#'   the relevant counterfactual mapping.
#'
#' - `band0` and `band1`: the first-stage bandwidths used for the untreated and
#'   treated regressions.
#'
#' - `y0_min`, `y0_max`, `y1_min`, and `y1_max`: estimated endpoints of the
#'   second-stage support regions.
#'
#' @export

RDD_extrapolate_CV_band <- function(
    Y,
    X,
    D,
    kernel = "gaussian",
    bands = seq(0.01, 1, length.out = 100),
    weights = rep(1, length(Y)),
    num_folds = 5,
    order = 1,
    batch_size = 1000L
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

    # Find unique rows and construct a mapping back to the original rows
    unique_rows_with_map <- function(x) {
        x <- as.matrix(x)
        x_unique <- unique(x)
        locs <- which(!duplicated(x))
        # identifier for each row
        row_key <- apply(x, MARGIN = 1L, FUN = function(z) paste(format(z, digits = 17), collapse = "\r"))
        inds <- match(row_key, row_key[locs])
        return(list(
            unique = x_unique,
            locs = locs,
            inds = inds
        ))
    }    

    # Choose bandwidth
    bands <- as.numeric(bands)

    if (length(bands) > 2L) {
        bands0 <- bands
        bands1 <- bands
    } else if (length(bands) == 2L) {
        bands0 <- bands[1L]
        bands1 <- bands[2L]
    } else {
        stop("bands must contain either two fixed bandwidths or more than two cross-validation candidates")
    }
    
    fit0 <- local_poly3(
        Y = Y0,
        X = X0,
        kernel = kernel,
        bands = bands0,
        weights = weights0,
        num_folds = num_folds,
        order = order
    )

    fit1 <- local_poly3(
        Y = Y1,
        X = X1,
        kernel = kernel,
        bands = bands1,
        weights = weights1,
        num_folds = num_folds,
        order = order
    )

    g0_raw <- fit0$predict_func
    g1_raw <- fit1$predict_func

    band0 <- fit0$optimal_band
    band1 <- fit1$optimal_band

    # Construct q1
    # Neareast nbh of X1
    nn1 <- FNN::get.knnx(data = X0, query = X1, k = 1)
    distance1 <- as.numeric(nn1$nn.dist[,1L]) # X1 & NN(X1) distance
    nnind1 <- as.integer(nn1$nn.index[, 1L])  # NN(X1) indices (in X0) 
    W1 <- distance1 < quart * band0
    # g0 <- g0_raw(X1, X0[nnind1, , drop = FALSE])  
    # Speed things up if there are repeated X values
    X1_sub <- X1[W1, , drop = FALSE]
    nnind1_sub <- nnind1[W1]

    X1_sub_map <- unique_rows_with_map(X1_sub)
    X1_sub_unique <- X1_sub_map$unique
    locs1_sub <- X1_sub_map$locs
    inds1_sub <- X1_sub_map$inds
    nnind1_sub_unique <- nnind1_sub[locs1_sub]

    # g0_unique <- g0_raw(X1_sub_unique, X0[nnind1_sub_unique, , drop = FALSE])
    # quick fix to avoid memory problems, do things in batches:
    g0_unique <- predict_in_batches(
        FUN = g0_raw, 
        x_fit = X1_sub_unique,
        x_center = X0[nnind1_sub_unique, , drop = FALSE],
        batch_size = batch_size
        )
    
    g0 <- g0_unique[inds1_sub]

    y1_max <- max(g0)
    y1_min <- min(g0)

    q1_fit <- local_poly3(
        Y = Y1[W1],
        X = matrix(g0, ncol = 1L),
        kernel = kernel,
        bands = bands,
        weights = weights1[W1],
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
    distance0 <- as.numeric(nn0$nn.dist[, 1L])
    nnind0 <- as.integer(nn0$nn.index[, 1L])
    W0 <- distance0 < quart * band1
    # g1 <- g1_raw(X0, X1[nnind0, , drop = FALSE])
    # Speed things up if there are repeated X values
    X0_sub <- X0[W0, , drop = FALSE]
    nnind0_sub <- nnind0[W0]

    X0_sub_map <- unique_rows_with_map(X0_sub)
    X0_sub_unique <- X0_sub_map$unique
    locs0_sub <- X0_sub_map$locs
    inds0_sub <- X0_sub_map$inds
    nnind0_sub_unique <- nnind0_sub[locs0_sub]

    # g1_unique <- g1_raw(X0_sub_unique, X1[nnind0_sub_unique, , drop = FALSE])
    # quick fix to avoid memory problems, do things in batches:
    g1_unique <- predict_in_batches(
        FUN = g1_raw, 
        x_fit = X0_sub_unique,
        x_center = X1[nnind0_sub_unique, , drop = FALSE],
        batch_size = batch_size
        )
    g1 <- g1_unique[inds0_sub]

    y0_max <- max(g1)
    y0_min <- min(g1)

    q0_fit <- local_poly3(
        Y = Y0[W0],
        X = matrix(g1, ncol = 1L),
        kernel = kernel,
        bands = bands,
        weights = weights0[W0],
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

    plotting1 <- cbind(g1, Y0[W0])
    plotting0 <- cbind(g0, Y1[W1])


    # Fitted values
    # y0 <- rep(NA_real_, length(Y))
    # y1 <- rep(NA_real_, length(Y))

    # y0[D == 0] <- as.numeric(g0_raw(X0, X0))
    # y1[D == 1] <- as.numeric(g1_raw(X1, X1))

    # Speed things up if there are repeated X values
    X_map <- unique_rows_with_map(X)
    X_unique <- X_map$unique
    inds_X <- X_map$inds

    # y0_unique <- g0_func(X_unique)
    # y1_unique <- g1_func(X_unique)
    # quick fix to avoid memory problems, do things in batches:
    y0_unique <- predict_in_batches(
        FUN = g0_func,
        x_fit = X_unique,
        batch_size = batch_size
    )
    y1_unique <- predict_in_batches(
        FUN = g1_func,
        x_fit = X_unique,
        batch_size = batch_size
    )
    y0 <- as.numeric(y0_unique)[inds_X]
    y1 <- as.numeric(y1_unique)[inds_X]
    

    # Indicator for observations whose counterfactual mean is supported (S)
    # S <- D * as.numeric( (y0_min < y1) & (y1 < y0_max) ) +
    #     (1-D) * as.numeric( (y1_min < y0) & (y0 < y1_max) )
    S <- integer(length(Y))
    S[D == 1] <- as.integer(y1[D == 1] > y0_min & y1[D == 1] < y0_max)
    S[D == 0] <- as.integer(y0[D == 0] > y1_min & y0[D == 0] < y1_max)
    
    # Counterfactual imputations
    # y1[S == 1 & D == 0] <- q1(y0[S == 1 & D == 0])
    # y0[S == 1 & D == 1] <- q0(y1[S == 1 & D == 1])

    # Speed things up if there are repeated X1 and X0 values
    X1_map <- unique_rows_with_map(X1)
    X1_unique <- X1_map$unique
    inds_X1 <- X1_map$inds
    X0_map <- unique_rows_with_map(X0)
    X0_unique <- X0_map$unique
    inds_X0 <- X0_map$inds

    # y0_imputed_unique <- g0_imputed_func(X1_unique)
    # y1_imputed_unique <- g1_imputed_func(X0_unique)
    # quick fix to avoid memory problems, do things in batches:
    y0_imputed_unique <- predict_in_batches(
        FUN = g0_imputed_func,
        x_fit = X1_unique,
        batch_size = batch_size
    )
    y1_imputed_unique <- predict_in_batches(
        FUN = g1_imputed_func,
        x_fit = X0_unique,
        batch_size = batch_size
    )
    
    y0_imputed <- as.numeric(y0_imputed_unique)[inds_X1]
    y1_imputed <- as.numeric(y1_imputed_unique)[inds_X0]

    y0[D == 1] <- y0_imputed
    y1[D == 0] <- y1_imputed

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