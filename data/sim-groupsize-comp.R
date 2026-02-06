# Joshua Alley
# Simulation: effect of group size variation on model comparison
# Fix number of groups (k=6, 64 groups), vary balance of group sizes


# parameters
k <- 6
n <- 2000
mu_tau <- 0.2
sigma_tau <- 0.3
sigma_y <- 1
group_cols <- paste0("x", 1:k)

# true group effects (fixed across scenarios)
all_groups <- expand.grid(lapply(1:k, function(j) 0:1))
colnames(all_groups) <- group_cols
all_groups$grp <- apply(all_groups[, group_cols], 1, paste, collapse = "_")
all_groups$tau <- rnorm(nrow(all_groups), mu_tau, sigma_tau)

# control variable coefficients (affect baseline only, not treatment effect)
beta_controls <- c(0.4, -0.3, 0.25)

# group size scenarios: Bernoulli probs for each x variable
# unequal probs create groups with very different sizes
prob_scenarios <- list(
  balanced = rep(0.5, k),
  mild     = c(0.35, 0.4, 0.5, 0.5, 0.6, 0.65),
  moderate = c(0.25, 0.35, 0.5, 0.5, 0.65, 0.7),
  strong   = c(0.2, 0.3, 0.5, 0.5, 0.675, 0.75)
)


# storage
results_summary <- vector("list", length(prob_scenarios))
results_group <- vector("list", length(prob_scenarios))


for(s in seq_along(prob_scenarios)){

  scenario <- names(prob_scenarios)[s]
  probs <- prob_scenarios[[s]]
  cat("\n---", scenario, "imbalance: probs =",
      paste(probs, collapse = ", "), "---\n")

  # generate data, ensuring all 64 groups exist (>= 1 obs each)
  sim_data <- data.frame(obs = 1:n)
  n_possible_groups <- 2^k
  repeat {
    for(j in 1:k){
      sim_data[[group_cols[j]]] <- rbinom(n, 1, probs[j])
    }
    sim_data$grp <- apply(sim_data[, group_cols], 1,
                          paste, collapse = "_")
    # check ALL possible groups are present
    if(length(unique(sim_data$grp)) == n_possible_groups) break
    cat("  Regenerating data (missing groups)...\n")
  }
  sim_data$treat <- rbinom(n, 1, 0.5)

  # control variables (affect baseline only, not treatment)
  sim_data$z1 <- rnorm(n)
  sim_data$z2 <- rbinom(n, 1, 0.5)
  sim_data$z3 <- rnorm(n)

  # merge true treatment effects
  sim_data <- sim_data %>%
    select(-any_of("tau")) %>%
    left_join(select(all_groups, grp, tau), by = "grp")

  # group sizes
  grp_sizes <- sim_data %>% count(grp, name = "n_obs")
  cat("  Min:", min(grp_sizes$n_obs),
      "| Median:", median(grp_sizes$n_obs),
      "| Max:", max(grp_sizes$n_obs),
      "| CV:", round(sd(grp_sizes$n_obs) / mean(grp_sizes$n_obs), 2),
      "\n")

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

  pred_data <- sim_data %>%
    distinct(grp, across(all_of(group_cols))) %>%
    mutate(treat = 1, z1 = 0, z2 = 0, z3 = 0)


  # 1. OLS 
  cat("  OLS...")
  f_ols <- as.formula(paste("y ~ treat * (", x_add, ") + z1 + z2 + z3"))
  fit_ols <- lm(f_ols, data = sim_data)

  est_ols <- slopes(fit_ols, variables = "treat",
                    newdata = pred_data) %>%
    select(grp, estimate) %>%
    rename(est_ols = estimate)
  cat(" done\n")


  # 2. Interaction model with varying slopes
  cat("  Interaction model...")
  f_inter <- bf(as.formula(
    paste("y ~ (1 + treat || grp) + treat * (", x_add, ") + z1 + z2 + z3")
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
    as.formula(paste("lambda ~", x_add, "+ (1 | grp)")),
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

  lambda_draws <- posterior_epred(fit_lambda, nlpar = "lambda",
                                  newdata = pred_data)
  lambda_med <- apply(lambda_draws, 2, median)

  est_lambda <- pred_data %>%
    mutate(est_lambda = lambda_med) %>%
    select(grp, est_lambda)
  cat(" done\n")


  # evaluate
  eval_data <- select(all_groups, grp, tau) %>%
    left_join(est_ols, by = "grp") %>%
    left_join(est_inter, by = "grp") %>%
    left_join(est_lambda, by = "grp") %>%
    left_join(grp_sizes, by = "grp") %>%
    mutate(
      bias_ols = est_ols - tau,
      bias_inter = est_inter - tau,
      bias_lambda = est_lambda - tau,
      scenario = scenario
    )

  results_group[[s]] <- eval_data

  results_summary[[s]] <- data.frame(
    scenario = scenario,
    min_obs = min(grp_sizes$n_obs),
    median_obs = median(grp_sizes$n_obs),
    max_obs = max(grp_sizes$n_obs),
    cv_obs = sd(grp_sizes$n_obs) / mean(grp_sizes$n_obs),
    rmse_ols = sqrt(mean(eval_data$bias_ols^2, na.rm = TRUE)),
    rmse_inter = sqrt(mean(eval_data$bias_inter^2, na.rm = TRUE)),
    rmse_lambda = sqrt(mean(eval_data$bias_lambda^2, na.rm = TRUE)),
    mean_bias_ols = mean(eval_data$bias_ols, na.rm = TRUE),
    mean_bias_inter = mean(eval_data$bias_inter, na.rm = TRUE),
    mean_bias_lambda = mean(eval_data$bias_lambda, na.rm = TRUE)
  )

  cat("  RMSE - OLS:", round(results_summary[[s]]$rmse_ols, 3),
      "| Interaction:", round(results_summary[[s]]$rmse_inter, 3),
      "| Lambda:", round(results_summary[[s]]$rmse_lambda, 3), "\n")
}



### compile and plot
summary_all <- bind_rows(results_summary) %>%
  mutate(
    scenario = case_when(
             scenario ==  "balanced" ~ "Balanced: (0.5 to 0.5)",
             scenario == "mild" ~ "Mild: (0.35 to 0.65)",
             scenario == "moderate" ~ "Modest: (0.25, to 0.7)",
              scenario == "strong"   ~ "Strong: (0.2, 0.75)"
  )
)


# RMSE by imbalance scenario
rmse_long <- summary_all %>%
  select(scenario, cv_obs, starts_with("rmse")) %>%
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

ggplot(rmse_long, aes(x = scenario, y = rmse,
                       color = model, group = model)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_grey(start = 0.1, end = 0.6) +
  labs(
    title = "RMSE by Group Size Imbalance",
    subtitle = paste0("64 groups"),
    x = "Group Size Imbalance",
    y = "RMSE",
    color = "Model"
  ) +
  theme(legend.position = "bottom")
ggsave("figures/sim-groupsize-rmse.png", height = 6, width = 8)


# mean bias by imbalance
bias_long <- summary_all %>%
  select(scenario, starts_with("mean_bias")) %>%
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

ggplot(bias_long, aes(x = scenario, y = bias,
                       color = model, group = model)) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_grey(start = 0.1, end = 0.6) +
  labs(
    title = "Mean Bias by Group Size Imbalance",
    subtitle = paste0("k = ", k, " (32 groups), n = ", n),
    x = "Group Size Imbalance",
    y = "Mean Bias",
    color = "Model"
  ) +
  theme(legend.position = "bottom")

# absolute bias by group size within each scenario
group_all <- bind_rows(results_group) %>%
  pivot_longer(
    cols = starts_with("bias_"),
    names_to = "model",
    values_to = "bias",
    names_prefix = "bias_"
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

ggplot(group_all, aes(x = n_obs, y = abs(bias))) +
  facet_wrap(model ~ scenario, scales = "free_x",
      ncol = 4) +
  geom_hline(yintercept = 0) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "darkgrey") +
  scale_color_grey(start = 0.1, end = 0.6) +
  labs(
    title = "Absolute Bias by Group Size",
    x = "Observations in Group",
    y = "|Bias|"
  ) +
  theme(legend.position = "bottom")
ggsave("figures/sim-groupsize-bias.png", height = 6, width = 8)
