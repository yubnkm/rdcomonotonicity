make_mock_result <- function() {
    list(
        y0 = c(0.2, 0.4, 0.6, 0.8),
        y1 = c(0.3, 0.5, 0.7, 0.9),
        S = c(1L, 1L, 0L, 1L),

        q0 = function(y) y - 0.1,
        q1 = function(y) y + 0.1,

        plotting1 = cbind(
            EY1 = c(0.3, 0.5, 0.7),
            EY0 = c(0.2, 0.4, 0.6)
        ),

        plotting0 = cbind(
            EY0 = c(0.2, 0.4, 0.6),
            EY1 = c(0.3, 0.5, 0.7)
        ),

        band0 = 0.4,
        band1 = 0.4
    )
}


make_mock_bootstrap_result <- function() {
    y0 <- c(0.2, 0.4, 0.6, 0.8)
    y1 <- c(0.3, 0.5, 0.7, 0.9)

    list(
        q0_grid = c(0.3, 0.5, 0.7),
        q1_grid = c(0.2, 0.4, 0.6),

        q0_draws = matrix(
            rep(c(0.2, 0.4, 0.6), 3L),
            nrow = 3L,
            ncol = 3L
        ),

        q1_draws = matrix(
            rep(c(0.3, 0.5, 0.7), 3L),
            nrow = 3L,
            ncol = 3L
        ),

        y0_draws = matrix(
            rep(y0, 3L),
            nrow = 4L,
            ncol = 3L
        ),

        y1_draws = matrix(
            rep(y1, 3L),
            nrow = 4L,
            ncol = 3L
        ),

        multipliers = matrix(
            1,
            nrow = 4L,
            ncol = 3L
        ),

        complete_cases = rep(TRUE, 4L),
        bootstrap_iter = 3L,
        points = 3L
    )
}


make_small_dgp <- function(include_x3 = FALSE) {
    grid <- seq(0.05, 0.95, length.out = 16L)

    X <- as.matrix(
        expand.grid(
            X1 = grid,
            X2 = grid
        )
    )

    if (include_x3) {
        X3 <- ((seq_len(nrow(X)) * 37L) %% nrow(X)) / nrow(X)
        X <- cbind(X, X3 = X3)
    }

    index <- sin(X[, 1L]) + 0.3 * sin(X[, 2L])

    if (include_x3) {
        index <- index + 0.1 * X[, 3L]
    }

    mu1 <- index
    mu0 <- 0.8 - 0.8 * cos(mu1)

    D <- as.integer(X[, 2L] < 0.7 - 0.4 * X[, 1L])

    set.seed(123)
    Y <- ifelse(D == 1L, mu1, mu0) +
        stats::rnorm(nrow(X), sd = 0.01)

    list(
        Y = Y,
        X = X,
        D = D
    )
}