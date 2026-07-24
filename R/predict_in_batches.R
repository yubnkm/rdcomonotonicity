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
