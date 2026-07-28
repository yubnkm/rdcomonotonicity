plot_q_results <- function(
    result, 
    points = 100L, 
    show_points = FALSE,
    bootstrap = TRUE,
    bootstrap_iter = 100L,
    level = 0.90,
    Y = NULL,
    X = NULL,
    D = NULL,
    kernel = "gaussian",
    weights = NULL,
    num_folds = 5,
    order = 1,
    batch_size = 1000L,
    parallel = FALSE,
    n_cores = NULL
    ) {

    if (bootstrap) {
        if(is.null(Y) || is.null(X) || is.null(D) ) {
            stop("Y, X, and D must be supplied when bootstrap = TRUE")
        }

        Y <- as.numeric(Y)
        X <- as.matrix(X)
        D <- as.numeric(D)
        
        if(is.null(weights)) weights <- rep(1, length(Y))
        weights <- as.numeric(weights)

        complete_cases <- complete.cases(cbind(Y, X, D, weights))
        Y <- Y[complete_cases]
        X <- X[complete_cases, , drop = FALSE]
        D <- D[complete_cases]
        weights <- weights[complete_cases]
    }


    # Extract second-stage estimation data
    q0_training_data <- as.data.frame(result$plotting1)
    q1_training_data <- as.data.frame(result$plotting0)
    names(q0_training_data) <- c("EY1", "EY0")
    names(q1_training_data) <- c("EY0", "EY1")

    q0_grid <- seq(
        from = min(q0_training_data$EY1, na.rm = TRUE),
        to = max(q0_training_data$EY1, na.rm = TRUE),
        length.out = points
    )
    q0_fitted <- as.numeric(result$q0(q0_grid))
    q0_data <- data.frame(
        EY1 = q0_grid,
        EY0 = q0_fitted
    )

    q1_grid <- seq(
        from = min(q1_training_data$EY0, na.rm = TRUE),
        to = max(q1_training_data$EY0, na.rm = TRUE),
        length.out = points
    )
    q1_fitted <- as.numeric(result$q1(q1_grid))
    q1_data <- data.frame(
        EY0 = q1_grid,
        EY1 = q1_fitted
    )

    q0_band_data <- NULL
    q1_band_data <- NULL

    if (bootstrap){
        bootstrap_bands <- c(result$band0, result$band1)
        n <- length(Y)

        multipliers <- matrix(
        stats::rexp(n = n*bootstrap_iter, rate = 1),
        nrow = n,
        ncol = bootstrap_iter
        )

        q0_draws <- matrix(NA_real_, nrow = length(q0_grid), ncol = bootstrap_iter)
        q1_draws <- matrix(NA_real_, nrow = length(q1_grid), ncol = bootstrap_iter)
        
        # Parallel for loop
        if(parallel){
            if(is.null(n_cores)) {
                n_cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
            }
            n_cores <- min(as.integer(n_cores), as.integer(bootstrap_iter))

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

                    parallel::clusterSetRNGStream(cl)

                    parallel::parLapply(
                        cl,
                        seq_len(bootstrap_iter),
                        function(b) {

                            bootstrap_weights <- weights * multipliers[, b]

                            bootstrap_fit <- RDD_extrapolate_CV_band(
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
                                q0 = as.numeric(
                                    bootstrap_fit$q0(q0_grid)
                                ),
                                q1 = as.numeric(
                                    bootstrap_fit$q1(q1_grid)
                                )
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

            message("Completed ", bootstrap_iter, " bootstrap draws using ", n_cores, " cores")

        # For loop
        } else { 

            for(b in seq_len(bootstrap_iter)) {

                bootstrap_weights <- weights * multipliers[, b]

                bootstrap_fit <- RDD_extrapolate_CV_band(
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

                message("Completed bootstrap draw ", b, " of ", bootstrap_iter)
            }

        }

        q0_dev <- abs(sweep(q0_draws, MARGIN = 1L, STATS = q0_fitted, FUN = "-"))
        q1_dev <- abs(sweep(q1_draws, MARGIN = 1L, STATS = q1_fitted, FUN = "-"))
        q0_rad <- apply(q0_dev, MARGIN = 1L, FUN = stats::quantile, probs = level, names = FALSE, na.rm=TRUE)
        q1_rad <- apply(q1_dev, MARGIN = 1L, FUN = stats::quantile, probs = level, names = FALSE, na.rm=TRUE)
        
        q0_band_data <- data.frame(
            EY1 = q0_grid,
            lower = q0_fitted - q0_rad,
            upper = q0_fitted + q0_rad
        )

        q1_band_data <- data.frame(
            EY0 = q1_grid,
            lower = q1_fitted - q1_rad,
            upper = q1_fitted + q1_rad
        )

    }

    # Plot q0
    q0_plot <- ggplot2::ggplot(
        q0_data,
        ggplot2::aes(x = EY1, y = EY0)
    ) +
        ggplot2::geom_line(color = "black") +
        ggplot2::geom_abline(
            intercept = 0, 
            slope = 1, 
            linetype = "dashed", 
            color = "gray40"
        ) +
        ggplot2::geom_vline(
            xintercept = c(min(q0_training_data$EY1, na.rm = TRUE), max(q0_training_data$EY1, na.rm = TRUE)),
            linetype = "dotted",
            color = "gray40"
        ) +
        ggplot2::labs(x = "EY1", y = "EY0") +
        ggplot2::theme_classic()

    if (show_points) {
        q0_plot <- q0_plot +
            ggplot2::geom_point(
                data = q0_training_data,
                ggplot2::aes(x = EY1, y = EY0),
                inherit.aes = FALSE,
                shape = 4,
                size = 1,
                stroke = 0.5,
                alpha = 0.4
            )
    }

    if (bootstrap) {
        q0_plot <- q0_plot +
            ggplot2::geom_ribbon(
                data = q0_band_data,
                ggplot2::aes(x = EY1, ymin = lower, ymax = upper),
                inherit.aes = FALSE,
                fill = "#82808074",
                alpha = 0.15
            )
    }

    # Plot q1
    q1_plot <- ggplot2::ggplot(
        q1_data,
        ggplot2::aes(x = EY0, y = EY1)
    ) +
        ggplot2::geom_line(color = "black") +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
        ggplot2::geom_vline(
            xintercept = c(min(q1_training_data$EY0, na.rm = TRUE), max(q1_training_data$EY0, na.rm = TRUE)),
            linetype = "dotted",
            color = "gray40"
        ) +
        ggplot2::labs(x = "EY0", y = "EY1") +
        ggplot2::theme_classic()

    if (show_points) {
        q1_plot <- q1_plot +
            ggplot2::geom_point(
                data = q1_training_data,
                ggplot2::aes(x = EY0, y = EY1),
                inherit.aes = FALSE,
                shape = 4,
                size = 1,
                stroke = 0.5,
                alpha = 0.4
            )
    }

    if (bootstrap) {
        q1_plot <- q1_plot +
            ggplot2::geom_ribbon(
                data = q1_band_data,
                ggplot2::aes(x = EY0, ymin = lower, ymax = upper),
                inherit.aes = FALSE,
                fill = "#82808074",
                alpha = 0.15
            )
    }

    # Compare q1 and q0
    q1_comp_data <- data.frame(
        EY0 = q1_data$EY0,
        EY1 = q1_data$EY1,
        func = "q1"
    )
    q0_inv_data <- data.frame(
        EY0 = q0_data$EY0,
        EY1 = q0_data$EY1,
        func = "q0_inverse"
    )
    
    comp_data <- rbind(q1_comp_data, q0_inv_data)

    comp_x_min <- max(min(q1_comp_data$EY0), min(q0_inv_data$EY0))
    comp_x_max <- min(max(q1_comp_data$EY0), max(q0_inv_data$EY0))

    comp_plot <- ggplot2::ggplot(
        comp_data,
        ggplot2::aes(x = EY0, y = EY1, color = func)
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    ggplot2::scale_color_manual(
        name = NULL,
        breaks = c("q1", "q0_inverse"),
        values = c(q1 = "blue", q0_inverse = "red"),
        labels = expression(hat(q)[1](y), hat(q)[0]^{-1}(y))
    ) +
    ggplot2::coord_cartesian(xlim = c(comp_x_min, comp_x_max)) +
    ggplot2::labs(
        x = expression(E[Y(0) ~ "|" ~ X == x]),
        y = expression(E[Y(1) ~ "|" ~ X == x])
    ) +
    ggplot2::theme_classic(base_size = 14) 

    return(
        list(
            q0_plot = q0_plot,
            q1_plot = q1_plot,
            comp_plot = comp_plot
        )
    )

}
