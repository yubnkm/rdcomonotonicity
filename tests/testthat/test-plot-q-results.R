test_that("plot_q_results returns all three plots", {
    result <- make_mock_result()

    plots <- plot_q_results(
        result = result,
        points = 20L,
        show_points = FALSE
    )

    expect_named(
        plots,
        c("q0_plot", "q1_plot", "comp_plot"),
        ignore.order = FALSE
    )

    expect_s3_class(plots$q0_plot, "ggplot")
    expect_s3_class(plots$q1_plot, "ggplot")
    expect_s3_class(plots$comp_plot, "ggplot")
})


test_that("plot_q_results accepts bootstrap results", {
    result <- make_mock_result()
    bootstrap_result <- make_mock_bootstrap_result()

    plots <- plot_q_results(
        result = result,
        bootstrap_result = bootstrap_result,
        level = 0.9
    )

    expect_s3_class(plots$q0_plot, "ggplot")
    expect_s3_class(plots$q1_plot, "ggplot")
    expect_s3_class(plots$comp_plot, "ggplot")
})