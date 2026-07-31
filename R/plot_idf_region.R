#' Plot the Region of Identified Conditional Average Treatment Effects
#'
#' Produces a scatter plot showing which observations have identified
#' conditional average treatment effects under the estimated comonotonicity
#' model. The plot is available only when `X` contains exactly two covariates.
#'
#' @param X A numeric matrix or data frame of assignment variables or covariates.
#'   Rows correspond to observations and columns correspond to covariates.
#'
#' @param result A fitted object returned by
#'   [RDD_extrapolate_CV_band()]. It must contain the support indicator `S`.
#'
#' @param threshold_function A function that accepts `X` and returns one
#'   treatment indicator per row. Returned values should be zero or one, or
#'   logical values that can be converted to zero and one. This function is used
#'   only to distinguish the treated and untreated regions in the plot.
#'
#' @return A `ggplot` object.
#'

#' @export
plot_idf_region <- function(
    X,
    result,
    threshold_function
){
    X <- as.matrix(X)
    if(ncol(X) != 2) stop("X must contain exactly two variables.")

    S <- as.integer(result$S)

    treatment <- as.integer(threshold_function(X))

    x_names <- colnames(X)
    if(is.null(x_names)) x_names <- c("X1", "X2")

    plot_data <- data.frame(
        x1 = X[, 1L],
        x2 = X[, 2L], 
        treatment = factor(
            treatment, 
            levels = c(0, 1),
            labels = c("Untreated", "Treated")
        ),
        identified = factor(
            S,
            levels = c(0, 1),
            labels = c("Not identified", "Identified")
        )
    )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = x1,
            y = x2,
            color = treatment, 
            shape = identified
        )
    ) +
        ggplot2::geom_point(
            ggplot2::aes(size = identified),
            stroke = 0.6,
            alpha = 0.8
        ) +
        ggplot2::scale_shape_manual(
            values = c(
                "Not identified" = 1,
                "Identified" = 4
            )
        ) +
        ggplot2::scale_size_manual(
            values = c(
                "Not identified" = 0.8,
                "Identified" = 2
            ),
            guide = "none"
        ) +
        ggplot2::labs(
            x = x_names[1L],
            y = x_names[2L],
            color = "Treatment region",
            shape = "CATE"
        ) +
        ggplot2::theme_classic(base_size = 14) +
        ggplot2::theme(
            legend.position = "bottom"
        )





        
}