devtools::load_all(".")
library(ggplot2)

set.seed(1)

data_path <- file.path("dev","ssextract_karthik.csv")
output_dir <- file.path("dev","figures", "summer_school")

level <- 0.9
bootstrap_iter <- 100
points <- 100
batch_size <- 1000L

data <- utils::read.csv(data_path)

# Imbens and Wager only use this subsample:
imbens_sampe <- with(
    data,
    mdcut01 >= -40 & mdcut01 <= 40
    & rdcut01 >= -40 & rdcut01 <= 40
)

data <- data[imbens_sampe, , drop = FALSE]

# Treatment: satisfies cutoff
D <- as.integer(!(data$mdcut01 > 0 & data$rdcut01 > 0))

# Covariates
X_raw <- as.matrix(data[, c("mdcut01", "rdcut01")])
colnames(X_raw) <- c("math_score", "reading_score")

min_X <- apply(X_raw, 2L, min)
max_X <- apply(X_raw, 2L, max)

cutoff <- -min_X / (max_X - min_X)

# normalize the covariates
X <- sweep(X_raw, 2L, min_X, FUN = "-")
X <- sweep(X, 2L, max_X - min_X, FUN = "/")
colnames(X) <- c("math_score", "reading_score")

# Outcomes
Y_math <- as.numeric(data$zmscr02)
Y_reading <- as.numeric(data$zrscr02)

n <- length(Y_math)

## Estimation
# bands <- seq(0.2, 0.5, length.out = 10)
bands <- c(0.2, 0.2)

fit_math <- RDD_extrapolate_CV_band(
    Y = Y_math,
    X = X,
    D = D,
    kernel = "gaussian",
    bands = bands,
    batch_size = batch_size
)

fit_read <- RDD_extrapolate_CV_band(
    Y = Y_reading,
    X = X,
    D = D,
    kernel = "gaussian",
    bands = bands,
    batch_size = batch_size
)

bootstrap_math <- RDD_extrapolate_bootstrap(
    fit_math,
    Y = Y_math,
    X = X,
    D = D,
    kernel = "gaussian",
    num_folds = 5,
    order = 1,
    points = 100L,
    bootstrap_iter = bootstrap_iter,
    parallel = TRUE
)

bootstrap_read <- RDD_extrapolate_bootstrap(
    fit_read,
    Y = Y_reading,
    X = X,
    D = D,
    kernel = "gaussian",
    num_folds = 5,
    order = 1,
    points = 100L,
    bootstrap_iter = bootstrap_iter,
    parallel = TRUE
)

# saveRDS(
#   list(
#     fit_math = fit_math,
#     fit_read = fit_read,
#     bootstrap_math = bootstrap_math,
#     bootstrap_read = bootstrap_read
#   ),
#   file = file.path("dev", "summerschool_results.rds")
# )

# Plotting q
plot_math <- plot_q_results(
    fit_math, 
    points = 100L, 
    show_points = FALSE,
    bootstrap_result = bootstrap_math,
    level = 0.90
)
ggsave(paste0(output_dir,"/math_q0.png"), plot = plot_math$q0_plot)
ggsave(paste0(output_dir,"/math_q1.png"), plot = plot_math$q1_plot)
ggsave(paste0(output_dir,"/math_q_comp.png"), plot = plot_math$comp_plot)

plot_read <- plot_q_results(
    fit_read, 
    points = 100L, 
    show_points = FALSE,
    bootstrap_result = bootstrap_read,
    level = 0.90
)
ggsave(paste0(output_dir,"/read_q0.png"), plot = plot_read$q0_plot)
ggsave(paste0(output_dir,"/read_q1.png"), plot = plot_read$q1_plot)
ggsave(paste0(output_dir,"/read_q_comp.png"), plot = plot_read$comp_plot)

# Figure 5.4
threshold_function <- function(x) {
    as.integer(!(x[,1L] > cutoff[1L] & x[,2L] > cutoff[2L]))
}
plot_idf_region(
    X = X,
    result = fit_math,
    threshold_function = threshold_function 
)
plot_idf_region(
    X = X,
    result = fit_read,
    threshold_function = threshold_function 
)

## Counterfactual policy effects
# math
math_cutoffs <- sort(unique(X[X[, 1L] >= cutoff[1L], 1L]))

math_policy_results <- lapply(
    math_cutoffs,
    function(math_cutoff) {
        math_policy <- function(x) {
            as.numeric(!(x[, 1L] > math_cutoff & x[, 2L] > cutoff[2L]))
        }
        counterfactual_policy_effect(
            result = fit_math,
            bootstrap_result = bootstrap_math,
            policy = math_policy,
            Y = Y_math,
            X = X,
            D = D,
            level = level
        )
    }
)

math_counterfactual <- data.frame(
    cutoff = math_cutoffs,
    est = vapply(math_policy_results, function(x) x$estimate, numeric(1L)),
    conf_low = vapply(math_policy_results, function(x) x$conf_low, numeric(1L)),
    conf_high = vapply(math_policy_results, function(x) x$conf_high, numeric(1L)),
    num_affected = vapply(math_policy_results, function(x) x$num_affected, numeric(1L))
)

# reading
read_cutoffs <- sort(unique(X[X[, 2L] >= cutoff[2L], 2L]))

read_policy_results <- lapply(
    read_cutoffs,
    function(read_cutoff) {
        read_policy <- function(x) {
            as.numeric(!(x[, 1L] > cutoff[1L] & x[, 2L] > read_cutoff))
        }
        counterfactual_policy_effect(
            result = fit_read,
            bootstrap_result = bootstrap_read,
            policy = read_policy,
            Y = Y_reading,
            X = X,
            D = D,
            level = level
        )
    }
)

read_counterfactual <- data.frame(
    cutoff = read_cutoffs,
    est = vapply(read_policy_results, function(x) x$estimate, numeric(1L)),
    conf_low = vapply(read_policy_results, function(x) x$conf_low, numeric(1L)),
    conf_high = vapply(read_policy_results, function(x) x$conf_high, numeric(1L)),
    num_affected = vapply(read_policy_results, function(x) x$num_affected, numeric(1L))
)


\library(patchwork)

# Scale the number affected onto the left y-axis
math_scale <- 0.025 / 15000
read_scale <- 0.025 / 10000


# Figure 5.5(a): Math cutoff
p_math_counterfactual <- ggplot(
    math_counterfactual,
    aes(x = cutoff)
) +
    geom_line(
        aes(y = conf_low),
        color = "green3",
        linewidth = 0.5
    ) +
    geom_line(
        aes(y = conf_high),
        color = "green3",
        linewidth = 0.5
    ) +
    geom_line(
        aes(y = est),
        color = "blue",
        linewidth = 0.9
    ) +
    geom_line(
        aes(y = num_affected * math_scale),
        color = "red",
        linetype = "dashed",
        linewidth = 0.9
    ) +
    scale_x_continuous(
        limits = c(cutoff[1L], 1),
        expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
        name = expression(theta[0]),
        breaks = seq(0, 0.025, by = 0.005),
        sec.axis = sec_axis(
            transform = ~ . / math_scale,
            name = "Number affected",
            breaks = seq(0, 15000, by = 3000)
        )
    ) +
    coord_cartesian(
        ylim = c(0, 0.025)
    ) +
    labs(
        x = "Math Score Threshold",
        title = "(a)"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.y.left = element_text(color = "blue"),
        axis.text.y.left = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red")
    )


# Figure 5.5(b): Reading cutoff
p_read_counterfactual <- ggplot(
    read_counterfactual,
    aes(x = cutoff)
) +
    geom_line(
        aes(y = conf_low),
        color = "green3",
        linewidth = 0.5
    ) +
    geom_line(
        aes(y = conf_high),
        color = "green3",
        linewidth = 0.5
    ) +
    geom_line(
        aes(y = est),
        color = "blue",
        linewidth = 0.9
    ) +
    geom_line(
        aes(y = num_affected * read_scale),
        color = "red",
        linetype = "dashed",
        linewidth = 0.9
    ) +
    scale_x_continuous(
        limits = c(cutoff[2L], 1),
        expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
        name = expression(theta[0]),
        breaks = seq(0, 0.025, by = 0.005),
        sec.axis = sec_axis(
            transform = ~ . / read_scale,
            name = "Number affected",
            breaks = seq(0, 10000, by = 2000)
        )
    ) +
    coord_cartesian(
        ylim = c(0, 0.025)
    ) +
    labs(
        x = "Reading Score Threshold",
        title = "(b)"
    ) +
    theme_classic(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.y.left = element_text(color = "blue"),
        axis.text.y.left = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "red"),
        axis.text.y.right = element_text(color = "red")
    )


figure_5_5 <- p_math_counterfactual + p_read_counterfactual
figure_5_5
ggsave(paste0(output_dir,"/Figure5_5.png"), plot = figure_5_5)

