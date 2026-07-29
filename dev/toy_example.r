devtools::load_all(".")
library(ggplot2)

set.seed(1)

g1_true <- function(x1, x2) sin(x1) + 0.3 * sin(x2)
g0_true <- function(x1,x2) 0.8 - 0.8 * cos(g1_true(x1, x2))

frontier <- function(x1) 0.7 - 0.4 * x1
d_func <- function(x1, x2) as.integer(x2 < frontier(x1))

q0_true <- function(y1) 0.8 - 0.8 * cos(y1)
q1_true <- function(y0) acos(pmax(-1, pmin(1, 1 - y0 / 0.8)))
  # Numerical truncation ensures the argument remains in [-1, 1].

simul_dgp <- function(n=500, noise_sd = 0.02) {

    X <- matrix( runif(2*n), nrow = n, ncol =2 )
    colnames(X) <- c("X1", "X2")

    mu1 <- g1_true(X[, 1], X[, 2])
    mu0 <- g0_true(X[, 1], X[, 2])

    D <- d_func(X[, 1], X[, 2])

    factual_mean <- ifelse(D == 1, mu1, mu0)
    Y <- factual_mean + rnorm(n, mean = 0, sd = noise_sd)

    list(
        Y = Y,
        X = X,
        D = D,
        mu0 = mu0,
        mu1 = mu1,
        factual_mean = factual_mean
    )
}

dat <- simul_dgp(n = 1000, noise_sd = 0.02)

result <- RDD_extrapolate_CV_band(
  Y = dat$Y,
  X = dat$X,
  D = dat$D,
  kernel = "gaussian",
  bands = seq(0.2,0.6, length.out = 5),
  num_folds = 5,
  order = 1
)

bootstrap_result <- RDD_extrapolate_bootstrap(
    result,
    Y = dat$Y,
    X = dat$X,
    D = dat$D,
    kernel = "gaussian",
    num_folds = 5,
    order = 1,
    points = 100L,
    bootstrap_iter = 100L,
    parallel = TRUE
)

plot <- plot_q_results(
    result, 
    points = 100L, 
    show_points = FALSE,
    bootstrap_result = bootstrap_result ,
    level = 0.90
)

plot$q0_plot
plot$q1_plot
plot$comp_plot
