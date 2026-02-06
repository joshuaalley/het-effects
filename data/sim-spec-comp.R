# Joshua Alley
# Simulation: compare OLS, interaction, and lambda model specifications
# Vary number of groups to illustrate regularization differences


# simulation parameters
num_vars <- c(3, 4, 5, 6, 7)  # binary grouping vars -> 2^k groups
n <- 2000
mu_tau <- 0.2                  # mean group treatment effect
sigma_tau <- 0.3               # SD of group treatment effects
sigma_y <- 1                   # outcome noise

# control variable coefficients (affect baseline only, not treatment)
beta_controls <- c(0.4, -0.3, 0.25)

# storage
results_summary <- vector("list", length(num_vars))
results_group <- vector("list", length(num_vars))


for(v in seq_along(num_vars)){

  k <- num_vars[v]
  n_groups <- 2^k
  cat("\n--- k =", k, ":", n_groups, "groups",
      "(~", round(n / n_groups), "obs per group) ---\n")

  # generate data
  sim_data <- data.frame(obs = 1:n)
  for(j in 1:k){
    sim_data[[paste0("x", j)]] <- rbinom(n, 1, 0.5)
  }
  sim_data$treat <- rbinom(n, 1, 0.5)

  # control variables (affect baseline only, not treatment)
  sim_data$z1 <- rnorm(n)
  sim_data$z2 <- rbinom(n, 1, 0.5)
  sim_data$z3 <- rnorm(n)

  # group from combination of x variables
  group_cols <- paste0("x", 1:k)
  sim_data$grp <- apply(sim_data[, group_cols], 1, paste, collapse = "_")

  # true group treatment effects
  unique_groups <- sort(unique(sim_data$grp))
  true_effects <- data.frame(
    grp = unique_groups,
    tau = rnorm(length(unique_groups), mu_tau, sigma_tau)
  )
  sim_data <- left_join(sim_data, true_effects, by = "grp")

  # outcome: controls affect baseline, x's define treatment heterogeneity
  control_effect <- sim_data$z1 * beta_controls[1] +
                    sim_data$z2 * beta_controls[2] +
                    sim_data$z3 * beta_controls[3]
  sim_data$y <- control_effect +
    sim_data$tau * sim_data$treat +
    rnorm(n, 0, sigma_y)


  # formulas
  x_add <- paste(group_cols, collapse = " + ")
  x_sat <- paste(group_cols, collapse = " * ")

  # one row per group for evaluation
  pred_data <- sim_data %>%
    distinct(grp, across(all_of(group_cols))) %>%
    mutate(treat = 1, z1 = 0, z2 = 0, z3 = 0)


  # 1. OLS 
  cat("  OLS...")
  f_ols <- as.formula(paste("y ~ treat * (", x_sat, ") + z1 + z2 + z3"))
  fit_ols <- lm(f_ols, data = sim_data)

  est_ols <- slopes(fit_ols, variables = "treat",
                    newdata = pred_data) %>%
    select(grp, estimate) %>%
    rename(est_ols = estimate)
  cat(" done\n")


  # 2. Interaction model with varying slopes
  cat("  Interaction model...")
  f_inter <- bf(as.formula(
    paste("y ~ (1 + treat || grp) + treat * (", x_sat, ") + z1 + z2 + z3")
  ))

  inter_prior <- c(
    prior(normal(0, 1), class = "sigma"),
    prior(normal(0, 1), class = "b")
  )

  fit_inter <- brm(f_inter,
                    data = sim_data,
                    prior = inter_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 0, silent = 2)

  est_inter <- slopes(fit_inter, variables = "treat",
                      newdata = pred_data) %>%
    select(grp, estimate) %>%
    rename(est_inter = estimate)
  cat(" done\n")


  # 3. Lambda model
  # x variables modify treatment (lambda), z variables are pure controls
  cat("  Lambda model...")
  f_lambda <- bf(
    y ~ lambda * treat + controls,
    as.formula(paste("lambda ~", x_sat, "+ (1 | grp)")),
    controls ~ z1 + z2 + z3,
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

  # extract group-level lambda (= treatment effect) estimates
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
    left_join(est_inter, by = "grp") %>%
    left_join(est_lambda, by = "grp") %>%
    mutate(
      bias_ols = est_ols - tau,
      bias_inter = est_inter - tau,
      bias_lambda = est_lambda - tau,
      k = k,
      num_groups = n_groups
    )

  results_group[[v]] <- eval_data

  results_summary[[v]] <- data.frame(
    k = k,
    num_groups = n_groups,
    obs_per_group = round(n / n_groups),
    rmse_ols = sqrt(mean(eval_data$bias_ols^2)),
    rmse_inter = sqrt(mean(eval_data$bias_inter^2)),
    rmse_lambda = sqrt(mean(eval_data$bias_lambda^2)),
    mean_bias_ols = mean(eval_data$bias_ols),
    mean_bias_inter = mean(eval_data$bias_inter),
    mean_bias_lambda = mean(eval_data$bias_lambda)
  )

  cat("  RMSE - OLS:", round(results_summary[[v]]$rmse_ols, 3),
      "| Interaction:", round(results_summary[[v]]$rmse_inter, 3),
      "| Lambda:", round(results_summary[[v]]$rmse_lambda, 3), "\n")
}



### compile results and plot
summary_all <- bind_rows(results_summary)

# RMSE by number of groups
rmse_long <- summary_all %>%
  select(num_groups, obs_per_group, starts_with("rmse")) %>%
  pivot_longer(
    cols = starts_with("rmse"),
    names_to = "model",
    values_to = "rmse",
    names_prefix = "rmse_"
  ) %>%
  mutate(
    model = case_when(
      model == "ols" ~ "OLS",
      model == "inter" ~ "Interaction",
      model == "lambda" ~ "Hierarchical"
    ),
    model = factor(model, levels = c("OLS", "Interaction", "Hierarchical"))
  )  %>%
  filter(model != "Interaction")

ggplot(rmse_long, aes(x = factor(num_groups), y = rmse,
                       color = model, group = model)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_grey(start = 0.5, end = 0.1) +
  labs(
    title = "Equal Group Size: RMSE of Group Treatment Effect Estimates",
    subtitle = paste0("n = ", n), # , ", treatment effect SD = ", sigma_tau
    x = "Number of Groups",
    y = "RMSE",
    color = "Model"
  ) +
  theme(legend.position = "bottom")
ggsave("figures/sim-spec-rmse.png", height = 6, width = 8)


# mean bias by number of groups
bias_long <- summary_all %>%
  select(num_groups, obs_per_group, starts_with("mean_bias")) %>%
  pivot_longer(
    cols = starts_with("mean_bias"),
    names_to = "model",
    values_to = "bias",
    names_prefix = "mean_bias_"
  ) %>%
  mutate(
    model = case_when(
      model == "ols" ~ "OLS",
      model == "inter" ~ "Interaction",
      model == "lambda" ~ "Lambda"
    ),
    model = factor(model, levels = c("OLS", "Interaction", "Lambda"))
  )

ggplot(bias_long, aes(x = factor(num_groups), y = bias,
                       color = model, group = model)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_grey(start = 0.1, end = 0.6) +
  labs(
    title = "Mean Bias of Group Treatment Effect Estimates",
    subtitle = paste0("n = ", n, ", treatment effect SD = ", sigma_tau),
    x = "Number of Groups",
    y = "Mean Bias",
    color = "Model"
  ) +
  theme(legend.position = "bottom")


# estimated vs true treatment effects
group_all <- bind_rows(results_group) %>%
  pivot_longer(
    cols = starts_with("est_"),
    names_to = "model",
    values_to = "estimate",
    names_prefix = "est_"
  ) %>%
  mutate(
    model = case_when(
      model == "ols" ~ "OLS",
      model == "inter" ~ "Interaction",
      model == "lambda" ~ "Hierarchical"
    ),
    model = factor(model, levels = c("OLS", "Interaction", "Hierarchical")),
    group_size = factor(paste0(num_groups, " Groups"),
     ordered = TRUE, 
    levels = c("8 Groups", "16 Groups", "32 Groups", "64 Groups", "128 Groups"))
  ) %>%
  filter(model != "Interaction")

ggplot(group_all, aes(x = tau, y = estimate)) +
  facet_grid(model ~ group_size) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  labs(
    title = "Equal Group Size: Estimated vs. True Effects",
    x = "True",
    y = "Estimated",
  ) +
  theme(legend.position = "bottom")
ggsave("figures/sim-spec-scatter.png", height = 6, width = 8)
