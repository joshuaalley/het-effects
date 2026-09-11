# Joshua Alley
# Simulation: hierarchical models against regularized machine learning
# see sim-additive-comp.R for the OLS comparisons and the data generating
# processes, which this script reuses
#
# R3 asks for the comparisons the paper spends several pages describing as the
# alternative: the LASSO on interactions (Ratkovic and Tingley 2017; Blackwell
# and Olson 2022) and BART or the causal forest. The point of the request is
# that beating unregularized OLS shows little, since the interesting question
# is how partial pooling compares to other ways of regularizing.
#
# Estimators, all recovering the same quantity - the treatment effect in each
# of the 2^k groups:
#   ols_sat   fully crossed OLS, the published comparison
#   hier      additive systematic part plus a pooled varying slope (lme4)
#   lasso     LASSO over the saturated interaction expansion (glmnet)
#   forest    causal forest (grf)
#   bart      Bayesian additive regression trees (bartCause)
#
# Fewer replications than sim-additive-comp.R because the forest and BART fits
# are several orders of magnitude slower than a mixed model.


library(lme4)
library(glmnet)
library(grf)
library(bartCause)
library(future.apply)

ML_REPS <- 50
ML_VARS <- c(3, 5, 7)
ML_DGPS <- c("normal", "systematic", "sparse", "additive_only")

# the data generating processes live in sim-additive-comp.R; source everything
# above its own run block so this script cannot trigger that simulation
sim_src <- readLines("data/sim-additive-comp.R")
sim_stop <- grep("^plan\\(multisession", sim_src)[1]
eval(parse(text = paste(sim_src[1:(sim_stop - 1)], collapse = "\n")))

# defined below the cut point in that script, so restated here
mc_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))


### estimators
# Each returns one treatment effect per group, in the order of s$grid.

est_ols_sat <- function(s, p1, p0) {
  x_sat <- paste(s$group_cols, collapse = " * ")
  fit <- lm(as.formula(paste("y ~ treat * (", x_sat, ") + z1 + z2 + z3")),
            data = s$d)
  as.vector(predict(fit, p1) - predict(fit, p0))
}

est_hier <- function(s, p1, p0) {
  x_add <- paste(s$group_cols, collapse = " + ")
  fit <- suppressMessages(suppressWarnings(lme4::lmer(
    as.formula(paste("y ~ treat * (", x_add,
                     ") + z1 + z2 + z3 + (0 + treat | grp)")),
    data = s$d, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))))
  as.vector(predict(fit, p1, allow.new.levels = TRUE) -
              predict(fit, p0, allow.new.levels = TRUE))
}

# LASSO over the saturated expansion. Penalizing the interaction basis is the
# regularizer that most directly competes with partial pooling: both shrink
# group-specific departures, one to zero and one toward the mean.
est_lasso <- function(s, p1, p0) {
  x_sat <- paste(s$group_cols, collapse = " * ")
  f <- as.formula(paste("~ treat * (", x_sat, ") + z1 + z2 + z3"))
  X <- model.matrix(f, data = s$d)[, -1, drop = FALSE]
  fit <- glmnet::cv.glmnet(X, s$d$y, alpha = 1, nfolds = 5)
  X1 <- model.matrix(f, data = p1)[, -1, drop = FALSE]
  X0 <- model.matrix(f, data = p0)[, -1, drop = FALSE]
  as.vector(predict(fit, X1, s = "lambda.min") -
              predict(fit, X0, s = "lambda.min"))
}

# Causal forest. Estimates a conditional average treatment effect for each
# respondent, which is averaged within group to match the estimand.
est_forest <- function(s, p1, p0) {
  X <- as.matrix(s$d[, c(s$group_cols, "z1", "z2", "z3")])
  cf <- grf::causal_forest(X, s$d$y, s$d$treat, num.trees = 1000, seed = 1)
  tau_hat <- predict(cf)$predictions
  as.vector(tapply(tau_hat, s$d$grp, mean)[s$grid$grp])
}

est_bart <- function(s, p1, p0) {
  X <- as.matrix(s$d[, c(s$group_cols, "z1", "z2", "z3")])
  fit <- bartCause::bartc(response = s$d$y, treatment = s$d$treat, confounders = X,
               n.samples = 400, n.burn = 200, n.chains = 2, verbose = FALSE)
  tau_hat <- apply(bartCause::extract(fit, type = "icate"), 2, mean)
  as.vector(tapply(tau_hat, s$d$grp, mean)[s$grid$grp])
}

ESTIMATORS <- list(ols_sat = est_ols_sat, hier = est_hier, lasso = est_lasso,
                   forest = est_forest, bart = est_bart)

ML_LABS <- c(
  ols_sat = "OLS, fully crossed",
  hier = "Hierarchical (partial pooling)",
  lasso = "LASSO on interactions",
  forest = "Causal forest",
  bart = "BART"
)


sim_ml_once <- function(rep_id, k, dgp) {
  s <- make_data(rep_id, k, dgp)
  p1 <- s$pred
  p1$treat <- 1
  p0 <- s$pred
  p0$treat <- 0

  bind_rows(lapply(names(ESTIMATORS), function(m) {
    # record why a fit failed rather than silently returning NA; a whole
    # estimator dropping out of the results is otherwise invisible
    err <- NA_character_
    est <- tryCatch(ESTIMATORS[[m]](s, p1, p0), error = function(e) {
      err <<- conditionMessage(e)
      rep(NA_real_, nrow(s$grid))
    })
    data.frame(dgp = dgp, k = k, num_groups = 2^k, rep = rep_id, model = m,
               grp = s$grid$grp, tau = s$grid$tau, est = est, error = err)
  }))
}


### run
plan(multisession, workers = max(1, parallel::detectCores() - 1))

ml_runs <- expand.grid(rep_id = 1:ML_REPS, k = ML_VARS, dgp = ML_DGPS,
                       stringsAsFactors = FALSE)

t0 <- Sys.time()
ml_raw <- bind_rows(future_lapply(seq_len(nrow(ml_runs)), function(i) {
  sim_ml_once(ml_runs$rep_id[i], ml_runs$k[i], ml_runs$dgp[i])
}, future.seed = TRUE))
message("machine learning comparison elapsed: ", format(Sys.time() - t0))


### summarise
fails <- ml_raw %>% filter(!is.na(error)) %>% distinct(model, error)
if (nrow(fails)) {
  cat("\n!!! some fits failed\n")
  print(as.data.frame(fails %>% mutate(error = substr(error, 1, 90))),
        row.names = FALSE)
}

ml_per_rep <- ml_raw %>%
  filter(!is.na(est)) %>%
  group_by(dgp, k, num_groups, rep, model) %>%
  summarize(rmse = sqrt(mean((est - tau)^2)),
            bias = mean(est - tau), .groups = "drop")

ml_summary <- ml_per_rep %>%
  group_by(dgp, k, num_groups, model) %>%
  summarize(rmse_se = mc_se(rmse), rmse = mean(rmse),
            bias = mean(bias), n_reps = n(), .groups = "drop") %>%
  mutate(rmse_lo = rmse - 1.96 * rmse_se, rmse_hi = rmse + 1.96 * rmse_se)

cat("\n=========== RMSE by estimator (Monte Carlo mean, 95% CI) ===========\n")
print(as.data.frame(
  ml_summary %>%
    transmute(dgp, groups = num_groups, model,
              rmse = sprintf("%.3f [%.3f, %.3f]", rmse, rmse_lo, rmse_hi),
              bias = sprintf("%+.3f", bias), n_reps) %>%
    arrange(dgp, groups, model)
), row.names = FALSE)

# paired against the hierarchical model, within replication
ml_paired <- ml_per_rep %>%
  select(dgp, k, num_groups, rep, model, rmse) %>%
  pivot_wider(names_from = model, values_from = rmse) %>%
  pivot_longer(cols = c(ols_sat, lasso, forest, bart),
               names_to = "model", values_to = "rmse") %>%
  filter(!is.na(rmse), !is.na(hier)) %>%
  group_by(dgp, num_groups, model) %>%
  summarize(d_se = mc_se(rmse - hier), d = mean(rmse - hier),
            win_hier = mean(rmse > hier), .groups = "drop") %>%
  mutate(lo = d - 1.96 * d_se, hi = d + 1.96 * d_se)

cat("\n====== Paired RMSE difference from the hierarchical model ======\n")
cat("Positive means the hierarchical model is more accurate.\n\n")
print(as.data.frame(
  ml_paired %>%
    transmute(dgp, groups = num_groups, model,
              difference = sprintf("%+.3f [%+.3f, %+.3f]", d, lo, hi),
              hier_wins = sprintf("%.2f", win_hier)) %>%
    arrange(dgp, groups, model)
), row.names = FALSE)


### figure
p_ml <- ml_summary %>%
  mutate(dgp = factor(DGP_LABS[dgp], levels = unname(DGP_LABS)),
         model = factor(ML_LABS[model], levels = unname(ML_LABS))) %>%
  ggplot(aes(x = factor(num_groups), y = rmse, group = model)) +
  facet_wrap(~dgp, nrow = 2) +
  geom_ribbon(aes(ymin = rmse_lo, ymax = rmse_hi), alpha = .18) +
  geom_line(aes(linetype = model)) +
  geom_point(aes(shape = model), size = 2.2, fill = "white") +
  scale_shape_manual(values = c(21, 24, 22, 23, 25)) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed",
                                   "dotdash", "longdash")) +
  labs(
    x = "Number of Groups", y = "RMSE of group treatment effects",
    linetype = NULL, shape = NULL,
    title = "Partial pooling against other ways of regularizing",
    subtitle = sprintf(paste("n = 2,000; %d Monte Carlo replications;",
                             "bands give 95%% Monte Carlo intervals"), ML_REPS)
  ) +
  theme(legend.position = "bottom")

ggsave("figures/sim-ml-comp.png", p_ml, height = 7.5, width = 9)
