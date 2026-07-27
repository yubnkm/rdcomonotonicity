plot_q_results <- function(result, points = 100L, show_points = FALSE) {

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

    # Plot q0
    q0_plot <- ggplot2::ggplot(
        q0_data,
        ggplot2::aes(x = EY1, y = EY0)
    ) +
        ggplot2::geom_line(color = "black") +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
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
        labels = expression(hat(q)[1](y), hat(q)^{-1}(y))
    ) +
    ggplot2::coord_cartesian(xlim = c(comp_x_min, comp_x_max)) +
    ggplot2::labs(
        x = expression(E[Y(0) ~ "|" ~ X == x]),
        y = expression(E[Y(1) ~ "|" ~ X == x])
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme_classic()

    return(
        list(
            q0_plot = q0_plot,
            q1_plot = q1_plot,
            comp_plot = comp_plot
        )
    )

}