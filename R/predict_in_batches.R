#' Evaluate Predictions in Batches
#'
#' Evaluates a prediction function in smaller batches to reduce memory usage.
#'
#' This is an internal computational helper used by
#' [RDD_extrapolate_CV_band()].
#' 
#'
#' @param x_fit A numeric matrix containing the points at which predictions are evaluated.
#' @param x_center An optional numeric matrix containing the local-regression centers corresponding row-by-row to `x_fit`.
#' @param batch_size A positive integer specifying the maximum number of rows evaluated in each batch.
#' 
#' @return A numeric vector containing one prediction for each row of `x_fit`.
#'
#' @keywords internal
#' 
#' 
predict_in_batches <- function(FUN, x_fit, x_center = NULL, batch_size) {
    x_fit <- as.matrix(x_fit)
    num_predictions <- nrow(x_fit)
    batch_size <- as.integer(batch_size)

    predictions <- rep(NA_real_, num_predictions)

    batch_starts <- seq.int(from = 1L, to = num_predictions, by = batch_size)

    for(batch_start in batch_starts) {
        batch_end <- min(batch_start + batch_size -1L, num_predictions)
    

    batch_indices <- batch_start:batch_end

    if (is.null(x_center)) {
            batch_predictions <- FUN(
                x_fit[batch_indices, , drop = FALSE]
            )
        } else {
            batch_predictions <- FUN(
                x_fit[batch_indices, , drop = FALSE],
                x_center[batch_indices, , drop = FALSE]
            )
        }
    
    batch_predictions <- as.numeric(batch_predictions)

    predictions[batch_indices] <- batch_predictions
    }
    return(predictions)
}
