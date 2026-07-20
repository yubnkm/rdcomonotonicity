local_poly3 <- function(
    Y,
    X,
    kernel = "gaussian",
    bands = seq(0.01, 1, length.out = 100),
    weights = NULL,
    num_folds = 5,
    order = 1
) {

    Y <- as.numeric(Y)
    X <- as.matrix(X)

    if (is.null(weights)) {
    weights <- rep(1, length(Y))
    }
    weights <- as.numeric(weights)

    # Choose kernel function
    kernel <- match.arg(
        kernel,
        choices = c("gaussian", "uniform", "triangular")
        )
    if (kernel == "gaussian") {
    kernel_function <- function(x, bandwidth) {
        stats::dnorm(x, mean = 0, sd = bandwidth)
    }} else if (kernel == "uniform") {
    kernel_function <- function(x, bandwidth) {
        as.numeric(abs(x) <= bandwidth)
    }} else if (kernel == "triangular") {
    kernel_function <- function(x, bandwidth) {
        as.numeric(abs(x / bandwidth) <= 1) *
        (1 - abs(x / bandwidth))
    }}

    # Removing missing entries
    valid_indices <- complete.cases(Y, X, weights)
    Y <- Y[valid_indices]
    X <- X[valid_indices, , drop = FALSE]
    weights <- weights[valid_indices]

    # Shuffle data for cross-validation
    num_samples <- length(Y)
    num_features <- ncol(X) 
    random_indices <- sample(seq_len(num_samples))
    X <- X[random_indices, , drop = FALSE]
    Y <- Y[random_indices]
    weights <- weights[random_indices]

    # Create polynomial basis functions
    order_indices <- 0:order
    exponent_list <- rep( list(order_indices), num_features )
    ind_temp <- as.matrix( do.call(expand.grid, exponent_list) )
    ind_temp <- ind_temp[rowSums(ind_temp) <= order, , drop = FALSE]

    basis_function <- function(x) {
        basis_result <- matrix(1, nrow = nrow(x), ncol = nrow(ind_temp))

        for(basis_idx in seq_len(nrow(ind_temp))){
            for(feature_idx in seq_len(num_features)){
                basis_result[, basis_idx] <- basis_result[, basis_idx] * x[, feature_idx] ^ ind_temp[basis_idx, feature_idx]
            }
        }
        basis_result
    }
    basis_matrix <- basis_function(X)
    num_basis <- ncol(basis_matrix)

    # Cross-validation
    if(length(bands) > 1){

        fold_size <- floor(num_samples/num_folds)
        mse_values <- rep(NA_real_, length(bands))

        for (band_idx in seq_along(bands)){
            bandwidth <- bands[band_idx]
            sse_values <- rep(NA_real_, num_folds)

            for (fold_idx in seq_len(num_folds)){
                fold_indices <- seq.int((fold_idx -1)*fold_size + 1, fold_idx*fold_size)
                comp_indices <- setdiff(seq_len(num_samples), fold_indices)

                X_comp <- X[comp_indices, , drop = FALSE]
                basis_comp <- basis_matrix[comp_indices, , drop = FALSE]
                Y_comp <- Y[comp_indices]
                weights_comp <- weights[comp_indices]

                # Compute kernel weights
                weight_function <- function(x){
                    weights_comp * kernel_function(sqrt(rowSums(sweep(X_comp, 2, x, "-")^2)), bandwidth)
                    }
                XKX_function <- function(x){
                    crossprod(basis_comp, basis_comp  * weight_function(x))
                }
                XKY_function <- function(x){
                    crossprod(basis_comp, Y_comp * weight_function(x))
                }
                g_function <- function(x){
                    drop(basis_function(x) %*% solve(XKX_function(x), XKY_function(x)))
                }

                X_fold <- X[fold_indices, , drop = FALSE]
                Y_fold <- Y[fold_indices]
                weights_fold <- weights[fold_indices]

                fitted_fold <- vapply(seq_len(nrow(X_fold)), function(i) g_function(X_fold[i, , drop=FALSE]), numeric(1))
                sse_values[fold_idx] = sum(weights_fold * (Y_fold - fitted_fold)^2)
            }
            mse_values[band_idx] <- sum(sse_values) / num_samples
        }

        optimal_band_idx <- which.min(mse_values)
        optimal_band <- bands[optimal_band_idx]

    } else { optimal_band <- bands[1]}

    # Defind final prediction function
    weight_function <- function(x) {
        weights * kernel_function(sqrt(rowSums(sweep(X, 2, x, "-")^2)), optimal_band)
    }
    XKX_function <- function(x) {
        crossprod(basis_matrix, basis_matrix * weight_function(x))
    }
    XKY_function <- function(x) {
        crossprod(basis_matrix, Y * weight_function(x))
    }
    predict_func <- function(x_fit, x) {
        vapply(
            seq_len(nrow(x_fit)),
            function(i) {
                beta_i <- solve(XKX_function(x[i, , drop = FALSE]), XKY_function(x[i, , drop=FALSE]))
                drop(basis_function(x_fit[i, , drop=FALSE]) %*% beta_i)
            },
            numeric(1)
        )
    }

    return(
        list(
        predict_func = predict_func,
        optimal_band = optimal_band
        )
    )
}