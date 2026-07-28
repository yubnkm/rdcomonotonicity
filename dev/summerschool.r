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
D <- as.integer(!(data$mdcut01 > 0 & data$rdcut01 > 01))

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

plot_math <- plot_q_results(fit_math, points = 100L, show_points = FALSE, bootstrap = FALSE)
ggsave(paste0(output_dir,"/math_q0.png"), plot = plot_math$q0_plot)
ggsave(paste0(output_dir,"/math_q1.png"), plot = plot_math$q1_plot)
ggsave(paste0(output_dir,"/math_q_comp.png"), plot = plot_math$comp_plot)

plot_read <- plot_q_results(fit_read, points = 100L, show_points = FALSE, bootstrap = FALSE)
ggsave(paste0(output_dir,"/read_q0.png"), plot = plot_read$q0_plot)
ggsave(paste0(output_dir,"/read_q1.png"), plot = plot_read$q1_plot)
ggsave(paste0(output_dir,"/read_q_comp.png"), plot = plot_read$comp_plot)

plot_math_bs <- plot_q_results(
    fit_math, 
    points = 100L, 
    show_points = FALSE, 
    bootstrap = TRUE,
    bootstrap_iter = 7L,
    level = 0.90,
    Y = Y_math,
    X = X,
    D = D,
    kernel = "gaussian",
    batch_size = batch_size,
    parallel = TRUE
    )
plot_math_bs$q0_plot
plot_math_bs$q1_plot
ggsave(paste0(output_dir,"/math_q0_bs.png"), plot = plot_math_bs$q0_plot)
ggsave(paste0(output_dir,"/math_q1_bs.png"), plot = plot_math_bs$q1_plot)

plot_read_bs <- plot_q_results(
    fit_read, 
    points = 100L, 
    show_points = FALSE, 
    bootstrap = TRUE,
    bootstrap_iter = 100L,
    level = 0.90,
    Y = Y_reading,
    X = X,
    D = D,
    kernel = "gaussian",
    batch_size = batch_size
    )
ggsave(paste0(output_dir,"/read_q0_bs.png"), plot = plot_read_bs$q0_plot)
ggsave(paste0(output_dir,"/read_q1_bs.png"), plot = plot_read_bs$q1_plot)


