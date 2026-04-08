# Joshua Alley
# Simulation: robustness to non-normal group effects
# Fix k=5 (32 groups), vary tail weight of group-effect distribution


# parameters
k <- 5
n <- 2000
n_groups <- 2^k
mu_tau <- 0.2
sigma_tau <- 0.3
sigma_y <- 1
group_cols <- paste0("x", 1:k)
beta_controls <- c(0.4, -0.3, 0.25)

# distribution scenarios: normal and t with varying df
# scale t draws to match the same mean and variance as normal case
dist_scenarios <- list(
  normal = list(label = "Normal", df = Inf),
  t5     = list(label = "t(5)", df = 5),
  t3     = list(label = "t(3)", df = 3)
)

draw_group_effects <- function(n_groups, mu, sigma, df) {
  if (is.infinite(df)) {
    # normal
    rnorm(n_groups, mu, sigma)
  } else {
    # t distribution scaled to match mean and variance
    # var of t(df) = df/(df-2), so scale by sigma / sqrt(df/(df-2))
    #raw <- rt(n_groups, df = df)
    #scale_factor <- sigma / sqrt(df / (df - 2))
    #mu + raw # * scale_factor
    rstudent_t(n = n_groups, df = df, mu = mu, sigma = sigma)
  }
}


# storage
results_summary <- vector("list", length(dist_scenarios))
results_group <- vector("list", length(dist_scenarios))


for (s in seq_along(dist_scenarios)) {

  scenario <- names(dist_scenarios)[s]
  sc <- dist_scenarios[[s]]
  cat("\n---", sc$label, "group effects ---\n")

  # generate data
  sim_data <- data.frame(obs = 1:n)
  for (j in 1:k) {
    sim_data[[paste0("x", j)]] <- rbinom(n, 1, 0.5)
  }
  sim_data$treat <- rbinom(n, 1, 0.5)

  # control variables
  sim_data$z1 <- rnorm(n)
  sim_data$z2 <- rbinom(n, 1, 0.5)
  sim_data$z3 <- rnorm(n)

  # group from combination of x variables
  sim_data$grp <- apply(sim_data[, group_cols], 1, paste, collapse = "_")

  # true group treatment effects from specified distribution
  unique_groups <- sort(unique(sim_data$grp))
  true_effects <- data.frame(
    grp = unique_groups,
    tau = draw_group_effects(length(unique_groups), mu_tau, sigma_tau, sc$df)
  )
  sim_data <- left_join(sim_data, true_effects, by = "grp")

  # flag outlier groups (beyond 2 SD from mean under normal)
  true_effects$outlier <- abs(true_effects$tau - mu_tau) > 2 * sigma_tau

  # outcome
  control_effect <- sim_data$z1 * beta_controls[1] +
                    sim_data$z2 * beta_controls[2] +
                    sim_data$z3 * beta_controls[3]
  sim_data$y <- control_effect +
    sim_data$tau * sim_data$treat +
    rnorm(n, 0, sigma_y)

  # formulas
  x_add <- paste(group_cols, collapse = " + ")
  x_sat <- paste(group_cols, collapse = " * ")

  # prediction data
  pred_data <- sim_data %>%
    distinct(grp, across(all_of(group_cols))) %>%
    mutate(treat = 1, z1 = 0, z2 = 0, z3 = 0)


  # 1. OLS with fully crossed interactions
  cat("  OLS...")
  f_ols <- as.formula(paste("y ~ treat * (", x_sat, ") + z1 + z2 + z3"))
  fit_ols <- lm(f_ols, data = sim_data)

  est_ols <- slopes(fit_ols, variables = "treat",
                    newdata = pred_data) %>%
    select(grp, estimate) %>%
    rename(est_ols = estimate)
  cat(" done\n")


  # 2. Lambda model: main effects in heterogeneity equation
  cat("  Lambda model...")
  f_lambda <- bf(
    y ~ lambda * treat + controls,
    as.formula(paste("lambda ~", x_add, "+ (1 | grp)")),
    as.formula(paste("controls ~", x_add, "+ z1 + z2 + z3")),
    nl = TRUE
  )

  lambda_prior <- c(
    prior(normal(0, 1), nlpar = "lambda"),
    prior(normal(0, 1), nlpar = "controls")
  )

  fit_lambda <- brm(f_lambda,
                     data = sim_data,
                     prior = lambda_prior,
                     family = gaussian(),
                     cores = 4,
                     backend = "cmdstanr",
                     refresh = 0, silent = 2)

  lambda_draws <- posterior_epred(fit_lambda, nlpar = "lambda",
                                  newdata = pred_data)
  lambda_med <- apply(lambda_draws, 2, median)

  est_lambda <- pred_data %>%
    mutate(est_lambda = lambda_med) %>%
    select(grp, est_lambda)
  cat(" done\n")


  # evaluate
  eval_data <- true_effects %>%
    left_join(est_ols, by = "grp") %>%
    left_join(est_lambda, by = "grp") %>%
    mutate(
      bias_ols = est_ols - tau,
      bias_lambda = est_lambda - tau,
      scenario = sc$label
    )

  results_group[[s]] <- eval_data

  # separate RMSE for outlier and non-outlier groups
  eval_outlier <- eval_data %>% filter(outlier)
  eval_regular <- eval_data %>% filter(!outlier)

  results_summary[[s]] <- data.frame(
    scenario = sc$label,
    df = sc$df,
    n_outliers = sum(eval_data$outlier),
    rmse_ols = sqrt(mean(eval_data$bias_ols^2)),
    rmse_lambda = sqrt(mean(eval_data$bias_lambda^2)),
    rmse_ols_outlier = if (nrow(eval_outlier) > 0)
      sqrt(mean(eval_outlier$bias_ols^2)) else NA_real_,
    rmse_lambda_outlier = if (nrow(eval_outlier) > 0)
      sqrt(mean(eval_outlier$bias_lambda^2)) else NA_real_,
    rmse_ols_regular = if (nrow(eval_regular) > 0)
      sqrt(mean(eval_regular$bias_ols^2)) else NA_real_,
    rmse_lambda_regular = if (nrow(eval_regular) > 0)
      sqrt(mean(eval_regular$bias_lambda^2)) else NA_real_
  )

  cat("  RMSE - OLS:", round(results_summary[[s]]$rmse_ols, 3),
      "| Lambda:", round(results_summary[[s]]$rmse_lambda, 3), "\n")
  cat("  Outlier groups:", results_summary[[s]]$n_outliers, "\n")
}



### compile and plot
summary_all <- bind_rows(results_summary)
group_all <- bind_rows(results_group)

# overall RMSE comparison
rmse_long <- summary_all %>%
  select(scenario, rmse_ols, rmse_lambda) %>%
  pivot_longer(
    cols = starts_with("rmse"),
    names_to = "model",
    values_to = "rmse",
    names_prefix = "rmse_"
  ) %>%
  mutate(
    model = case_when(
      model == "ols" ~ "OLS",
      model == "lambda" ~ "Hierarchical"
    ),
    model = factor(model, levels = c("OLS", "Hierarchical")),
    scenario = factor(scenario, levels = c("Normal", "t(5)", "t(3)"))
  )

ht_rmse <- ggplot(rmse_long, aes(x = scenario, y = rmse,
                       color = model, group = model)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_grey(start = 0.5, end = 0.1) +
  labs(
    title = "RMSE by Distribution",
    x = "Effect Distribution",
    y = "RMSE",
    color = "Model"
  ) +
  theme(legend.position = "bottom")
ht_rmse

# scatter: true vs estimated by distribution
group_long <- group_all %>%
  pivot_longer(
    cols = c(est_ols, est_lambda),
    names_to = "model",
    values_to = "estimate",
    names_prefix = "est_"
  ) %>%
  mutate(
    model = case_when(
      model == "ols" ~ "OLS",
      model == "lambda" ~ "Hierarchical"
    ),
    model = factor(model, levels = c("OLS", "Hierarchical")),
    scenario = factor(scenario, levels = c("Normal", "t(5)", "t(3)"))
  )

ht_scatter <- ggplot(group_long, aes(x = tau, y = estimate)) +
  facet_grid(model ~ scenario) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(shape = outlier), alpha = 0.6) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("Regular", "Outlier"),
                     name = "Group Type") +
  labs(
    title = "Estimated vs. True Effects",
    x = "True",
    y = "Estimated"
  ) +
  theme(legend.position = "bottom")

ht_rmse + ht_scatter + plot_annotation(title = "Student-T vs Normal Effect Distributions",     
                                       subtitle = paste0(n_groups, " groups, n = ", n),)
ggsave("figures/sim-heavy-tails.png", height = 6, width = 8)


# bias by outlier status across distributions
bias_summary <- group_all %>%
  pivot_longer(
    cols = c(bias_ols, bias_lambda),
    names_to = "model",
    values_to = "bias",
    names_prefix = "bias_"
  ) %>%
  mutate(
    model = case_when(
      model == "ols" ~ "OLS",
      model == "lambda" ~ "Hierarchical"
    ),
    model = factor(model, levels = c("OLS", "Hierarchical")),
    scenario = factor(scenario, levels = c("Normal", "t(5)", "t(3)")),
    group_type = ifelse(outlier, "Outlier", "Regular")
  )

ggplot(bias_summary, aes(x = abs(tau - mu_tau), y = abs(bias),
                          color = model)) +
  facet_wrap(~ scenario) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_grey(start = 0.5, end = 0.1) +
  labs(
    title = "Absolute Bias by Distance from Mean Effect",
    subtitle = paste0(n_groups, " groups, n = ", n),
    x = "Distance of True Effect from Mean",
    y = "|Bias|",
    color = "Model"
  ) +
  theme(legend.position = "bottom")
