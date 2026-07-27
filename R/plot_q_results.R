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

    q0_plot
}