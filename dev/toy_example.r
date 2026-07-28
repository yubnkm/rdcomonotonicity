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

## 
Y <- dat$Y
X <- dat$X
D <- dat$D
weights <- rep(1, length(Y))
kernel <- "gaussian"
bands <- seq(0.2, 0.6, length.out = 5)
num_folds <- 5
order <- 1
batch_size <- 1000L

bootstrap_iter <- 10L
level <- 0.9
points <- 100L

## batch function test
# set.seed(123)
# result_large_batch <- RDD_extrapolate_CV_band(
#   Y = dat$Y,
#   X = dat$X,
#   D = dat$D,
#   kernel = "gaussian",
#   bands = seq(0.2,0.6, length.out = 5),
#   num_folds = 5,
#   order = 1,
#   batch_size = 1000000000L
# )
# set.seed(123)

# result_small_batch <- RDD_extrapolate_CV_band(
#   Y = dat$Y,
#   X = dat$X,
#   D = dat$D,
#   kernel = "gaussian",
#   bands = seq(0.2,0.6, length.out = 5),
#   num_folds = 5,
#   order = 1,
#     batch_size = 25L
# )
# all.equal(
#     result_large_batch$y0,
#     result_small_batch$y0,
#     tolerance = 1e-12
# )

# all.equal(
#     result_large_batch$y1,
#     result_small_batch$y1,
#     tolerance = 1e-12
# )

# all.equal(
#     result_large_batch$S,
#     result_small_batch$S
# )

# all.equal(
#     result_large_batch$plotting0,
#     result_small_batch$plotting0,
#     tolerance = 1e-12
# )

# all.equal(
#     result_large_batch$plotting1,
#     result_small_batch$plotting1,
#     tolerance = 1e-12
# )

result <- RDD_extrapolate_CV_band(
  Y = dat$Y,
  X = dat$X,
  D = dat$D,
  kernel = "gaussian",
  bands = seq(0.2,0.6, length.out = 5),
  num_folds = 5,
  order = 1
)

plot <- plot_q_results(
    result, 
    points = 100L, 
    show_points = TRUE,
    bootstrap = TRUE,
    bootstrap_iter = 10L,
    level = 0.90,
    Y = dat$Y,
    X = dat$X,
    D = dat$D,
    kernel = "gaussian",
    num_folds = 5,
    order = 1
    )
plot$q1_plot
plot$q0_plot
plot$comp_plot


## Plotting figure 3.1 (b)
# q0 <- result$q0
# plotting1 <- as.data.frame(result$plotting1)

# names(plotting1) <- c("g1_tilde", "Y0")

# grid <- seq(0, 1, length.out = 1000)

# true_curve <- data.frame(
#     y1 = g1_true(grid, frontier(grid)),
#     y0 = g0_true(grid, frontier(grid))
# )
# true_curve <- true_curve[order(true_curve$y1), ]

# gsgrd <- seq(from = result$y0_min, to = result$y0_max, length.out = 100)

# fit <- vapply(
#   gsgrd,
#   function(y) {
#     as.numeric(q0(matrix(y, nrow = 1, ncol = 1)))[1]
#   },
#   numeric(1)
# )

# estimated_curve <- data.frame(y1 = gsgrd, y0 = fit)

# fig9 <- ggplot() +
#     # True q0(y)
#     geom_line(
#         data = true_curve,
#         aes(x = y1, y = y0,linetype = "true"),
#         color = "blue",
#         linewidth = 0.9
#     ) + 
#     # Estimated q0(y)
#     geom_line(
#         data = estimated_curve,
#         aes(x = y1, y= y0, linetype = "estimated"),
#         color = "black",
#         linewidth = 0.9

#     ) +
#     # Domain of q0 (y0_min and y0_max)
#     geom_vline(
#         xintercept = c(result$y0_min, result$y0_max),
#         linetype = "dashed",
#         color = "black",
#         linewidth = 0.5
#     ) +
#     # observation scatter plots
#     geom_point(
#         data = plotting1,
#         aes(x = g1_tilde, y = Y0),
#         shape = 4,
#         size = 1.8,
#         stroke = 1,
#         color = "#828181"
#     ) +
#     scale_linetype_manual(
#         name = NULL,
#         breaks = c("true", "estimated"),
#         values = c(true = "dashed", estimated = "solid"),
#         labels = expression(q[0](y), hat(q)[0](y))
#     ) +
#     labs(
#         x = expression(E[Y(1) ~ "|" ~ X == x]),
#         y = expression(E[Y(0) ~ "|" ~ X == x])
#     ) +
#     theme_classic(base_size = 14) +

#     theme(
#         axis.line = element_line(
#         color = "black",
#         linewidth = 0.5
#     ),
#     axis.ticks = element_line(
#         color = "black",
#         linewidth = 0.5
#     ),
#         legend.position = "top",
#         legend.justification = "left",
#         legend.box.just = "left",
#         legend.background = element_blank(),
#         legend.key = element_blank()
#   )
# fig9
