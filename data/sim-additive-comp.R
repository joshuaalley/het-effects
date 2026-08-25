# Joshua Alley
# Where does the hierarchical advantage actually come from?
#
# Models compared, all estimating the group treatment effect lambda_g:
#   ols_sat  y ~ treat * (x1 * x2 * ... * xk)   fully crossed (published race)
#   ols_add  y ~ treat * (x1 + x2 + ... + xk)   additive
#   lmer     ols_add + (0 + treat | grp)        partial pooling
#
# The Monte Carlo loop uses lme4 so that 200 replications across six data
# generating processes is affordable. REML corresponds to no prior at all, so
# it regularizes LESS than any half-t or half-normal prior --- if the additive
# OLS loses to lme4 it loses to the brms model a fortiori. Set RUN_BRMS_CHECK
# below to confirm that on a smaller grid; it takes roughly two hours.


RUN_BRMS_CHECK <- TRUE # slow: validates lme4 against brms and varies the prior

N_REPS <- 200
NUM_VARS <- c(3, 4, 5, 6, 7) # binary grouping vars -> 2^k groups
N <- 2000
MU_TAU <- 0.2
SIGMA_Y <- 1
BETA_CONTROLS <- c(0.4, -0.3, 0.25)
SIGMA_TAU <- 0.3 # target marginal SD of group effects, held roughly
                 # constant across DGPs so RMSE stays comparable

# NOTE: make_data() seeds on match(dgp, DGPS), so a DGP's position in this
# vector determines its random draws. Append new DGPs at the end rather than
# inserting them, or every result for the DGPs after the insertion point
# changes. Display order is set by DGP_LABS, not by this vector.
DGPS <- c("normal", "systematic", "sparse", "bimodal", "size_corr",
          "student_t", "imbalanced", "additive_only")

DGP_LABS <- c(
  normal      = "Normal",
  systematic  = "Systematic + deviation",
  sparse      = "Sparse\n(most groups identical)",
  bimodal     = "Bimodal",
  student_t   = "Student-t\n(heavy tails)",
  imbalanced  = "Unequal group sizes",
  size_corr   = "Unequal sizes, effects\ncorrelated with size"
)
MODEL_LABS <- c(
  ols_sat = "OLS, fully crossed",
  ols_add = "OLS, additive interaction",
  lmer    = "Hierarchical (partial pooling)"
)

GAMMA_POOL <- c(0.20, -0.15, 0.12, -0.10, 0.08, -0.06, 0.05)
# Assignment probabilities for the two unbalanced DGPs. Both use the same
# range, so the only thing separating "imbalanced" from "size_corr" is whether
# the effects happen to correlate with how big a group is.
IMBALANCE_RANGE <- c(0.35, 0.65)
SPARSE_SHARE <- 0.15 # fraction of groups that deviate at all
BIMODAL_GAP <- 0.6   # distance between the two modes


# ---------------------------------------------------------------------------
# data generating processes
# ---------------------------------------------------------------------------

# Assignment probabilities for the modifiers. Only size_corr uses unequal
# probabilities, which is what makes group sizes --- and therefore the
# size-effect correlation --- vary.
assign_probs <- function(k, dgp) {
  if (!dgp %in% c("imbalanced", "size_corr")) return(rep(0.5, k))
  # Unequal assignment probabilities are what make group sizes uneven. At
  # k = 7 this range already produces severe imbalance: expected cell sizes
  # run from about 1 observation to about 100. This generalises the
  # prob_scenarios grid in sim-groupsize-comp.R, which fixed k = 6 and stepped
  # through balanced / mild / moderate / strong; widen IMBALANCE_RANGE to
  # recover the stronger settings.
  seq(IMBALANCE_RANGE[1], IMBALANCE_RANGE[2], length.out = k)
}

# True group treatment effects. grid holds one row per possible group.
draw_effects <- function(grid, group_cols, k, dgp) {
  n_groups <- nrow(grid)
  X <- as.matrix(grid[, group_cols, drop = FALSE])

  if (dgp == "normal") {
    MU_TAU + rnorm(n_groups, 0, SIGMA_TAU)

  } else if (dgp == "systematic") {
    # part of the heterogeneity really is additive in the modifiers; the rest
    # is idiosyncratic. Additive model should be closer
    systematic <- as.vector(X %*% GAMMA_POOL[1:k])
    systematic <- systematic - mean(systematic)
    scale_sys <- SIGMA_TAU / sqrt(2) / max(sd(systematic), 1e-8)
    MU_TAU + systematic * scale_sys + rnorm(n_groups, 0, SIGMA_TAU / sqrt(2))

  } else if (dgp == "sparse") {
    # most groups share one common effect, a minority genuinely differ.
    deviates <- rbinom(n_groups, 1, SPARSE_SHARE)
    spike <- SIGMA_TAU / sqrt(max(SPARSE_SHARE, 1e-8))
    MU_TAU + deviates * rnorm(n_groups, 0, spike)

  } else if (dgp == "bimodal") {
    # two clusters of effects. A normal prior centred on the grand mean 
    # xenters where there are few true effects
    hi <- rbinom(n_groups, 1, 0.5)
    within <- sqrt(max(SIGMA_TAU^2 - (BIMODAL_GAP / 2)^2, 0.01))
    MU_TAU + ifelse(hi, BIMODAL_GAP / 2, -BIMODAL_GAP / 2) +
      rnorm(n_groups, 0, within)

  } else if (dgp == "additive_only") {
    # The adversarial case. Effects are EXACTLY additive in the modifiers with
    # no idiosyncratic group component, so the additive OLS is correctly
    # specified and the hierarchical model is paying to estimate a sigma_theta
    # that is genuinely zero. Every other DGP here gives the group effects a
    # component the additive model structurally cannot represent, which is a
    # race it cannot win; this one it can. The hierarchical model does lose
    # here, by roughly 3% of RMSE --- small against the 20-100% it gains when
    # pooling is warranted, but a real scope condition rather than a rounding
    # error, and it belongs in the paper.
    sysc <- as.vector(X %*% GAMMA_POOL[1:k])
    sysc <- sysc - mean(sysc)
    MU_TAU + sysc * (SIGMA_TAU / max(sd(sysc), 1e-8))

  } else if (dgp == "imbalanced") {
    # Group sizes are very uneven, but the effects themselves are drawn the
    # same way as under "normal". This isolates imbalance per se: partial
    # pooling should shrink the small, noisily estimated groups hardest while
    # leaving the large ones alone, and a saturated model has nothing to fall
    # back on when a cell holds a handful of observations.
    MU_TAU + rnorm(n_groups, 0, SIGMA_TAU)

  } else if (dgp == "size_corr") {
    # larger groups have larger effects, so the groups with the most data pull
    # the grand mean toward themselves and the small groups get shrunk a different mean
    p <- assign_probs(k, dgp)
    log_size <- as.vector(X %*% log(p) + (1 - X) %*% log(1 - p))
    z <- scale(log_size)[, 1]
    MU_TAU + SIGMA_TAU * (0.8 * z + sqrt(1 - 0.8^2) * rnorm(n_groups))

  } else if (dgp == "student_t") {
    # heavy tails: true outliers from a t-distribution
    df <- 3
    MU_TAU + SIGMA_TAU * rt(n_groups, df) / sqrt(df / (df - 2))

  } else {
    stop("unknown dgp: ", dgp)
  }
}


# ---------------------------------------------------------------------------
# estimation
# ---------------------------------------------------------------------------

# delta-method SE for fitted(treat = 1) - fitted(treat = 0) at each grid row
se_effect_lm <- function(fit, p1, p0) {
  tt <- delete.response(terms(fit))
  D <- model.matrix(tt, data = p1, xlev = fit$xlevels) -
    model.matrix(tt, data = p0, xlev = fit$xlevels)
  V <- vcov(fit)
  keep <- !is.na(coef(fit))
  sqrt(pmax(rowSums((D[, keep, drop = FALSE] %*% V[keep, keep]) *
    D[, keep, drop = FALSE]), 0))
}

# For the mixed model the effect is (fixed contrast) + theta_g, so the variance
# is the delta-method variance of the fixed part plus the conditional variance
# of the varying slope. Ignoring the covariance between the two is the usual
# approximation; the brms block reports exact posterior intervals.
se_effect_mer <- function(fit, p1, p0, grps) {
  tt <- delete.response(terms(fit, fixed.only = TRUE))
  D <- model.matrix(tt, data = p1) - model.matrix(tt, data = p0)
  v_fix <- pmax(rowSums((D %*% as.matrix(vcov(fit))) * D), 0)
  re <- ranef(fit, condVar = TRUE)$grp
  v_re <- setNames(as.vector(attr(re, "postVar")[1, 1, ]), rownames(re))
  sqrt(v_fix + ifelse(is.na(v_re[grps]), 0, v_re[grps]))
}


make_data <- function(rep_id, k, dgp) {
  set.seed(as.integer(1e6 * match(dgp, DGPS) + 1e4 * k + rep_id))

  group_cols <- paste0("x", 1:k)
  p <- assign_probs(k, dgp)

  d <- data.frame(obs = 1:N)
  for (j in 1:k) d[[paste0("x", j)]] <- rbinom(N, 1, p[j])
  d$treat <- rbinom(N, 1, 0.5)
  d$z1 <- rnorm(N)
  d$z2 <- rbinom(N, 1, 0.5)
  d$z3 <- rnorm(N)
  d$grp <- apply(d[, group_cols, drop = FALSE], 1, paste, collapse = "_")

  grid <- expand.grid(rep(list(0:1), k))
  names(grid) <- group_cols
  grid$grp <- apply(grid[, group_cols, drop = FALSE], 1, paste, collapse = "_")
  grid$tau <- draw_effects(grid, group_cols, k, dgp)

  d <- left_join(d, grid[, c("grp", "tau")], by = "grp")
  d$y <- d$z1 * BETA_CONTROLS[1] + d$z2 * BETA_CONTROLS[2] +
    d$z3 * BETA_CONTROLS[3] + d$tau * d$treat + rnorm(N, 0, SIGMA_Y)

  pred <- grid
  pred$z1 <- 0
  pred$z2 <- 0
  pred$z3 <- 0
  list(d = d, grid = grid, pred = pred, group_cols = group_cols)
}


sim_once <- function(rep_id, k, dgp) {
  s <- make_data(rep_id, k, dgp)
  d <- s$d
  x_add <- paste(s$group_cols, collapse = " + ")
  x_sat <- paste(s$group_cols, collapse = " * ")

  p1 <- s$pred
  p1$treat <- 1
  p0 <- s$pred
  p0$treat <- 0

  est <- list()
  se <- list()

  # 1. fully crossed OLS: the published comparison
  fit_sat <- lm(as.formula(paste("y ~ treat * (", x_sat, ") + z1 + z2 + z3")),
                data = d)
  est$ols_sat <- as.vector(predict(fit_sat, p1) - predict(fit_sat, p0))
  se$ols_sat <- se_effect_lm(fit_sat, p1, p0)

  # 2. additive interaction OLS
  fit_add <- lm(as.formula(paste("y ~ treat * (", x_add, ") + z1 + z2 + z3")),
                data = d)
  est$ols_add <- as.vector(predict(fit_add, p1) - predict(fit_add, p0))
  se$ols_add <- se_effect_lm(fit_add, p1, p0)

  # 3. additive systematic part plus a partially pooled varying slope. This
  # matches the brms non-linear model in the paper: lambda_g has an additive
  # systematic component plus a group varying intercept, and there is no
  # separate varying intercept on the outcome.
  fit_mer <- suppressMessages(suppressWarnings(lmer(
    as.formula(paste("y ~ treat * (", x_add,
                     ") + z1 + z2 + z3 + (0 + treat | grp)")),
    data = d, REML = TRUE, control = lmerControl(calc.derivs = FALSE)
  )))
  est$lmer <- as.vector(predict(fit_mer, p1, allow.new.levels = TRUE) -
    predict(fit_mer, p0, allow.new.levels = TRUE))
  # No standard error is recorded for the mixed model, deliberately. The
  # prediction error of lambda_g is x_g'(betahat - beta) + (thetahat_g -
  # theta_g), and lme4 does not expose the covariance between those two terms.
  # It is strongly negative, because the varying slope absorbs error in the
  # fixed part: adding the two variances gives intervals covering 0.99 against
  # a nominal 0.95 at eight groups, while condVar alone covers 0.84. Sign and
  # exaggeration errors for this model are computed in the brms block below,
  # from posterior intervals whose coverage is 0.94 to 0.97.
  se$lmer <- rep(NA_real_, nrow(s$grid))

  bind_rows(lapply(names(est), function(m) {
    data.frame(
      dgp = dgp, k = k, num_groups = 2^k, rep = rep_id, model = m,
      grp = s$grid$grp, tau = s$grid$tau, est = est[[m]], se = se[[m]],
      singular = if (m == "lmer") isSingular(fit_mer) else NA,
      sigma_theta = if (m == "lmer") {
        as.data.frame(VarCorr(fit_mer))$sdcor[1]
      } else {
        NA_real_
      }
    )
  }))
}


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

plan(multisession, workers = max(1, parallel::detectCores() - 1))

runs <- expand.grid(rep_id = 1:N_REPS, k = NUM_VARS, dgp = DGPS,
                    stringsAsFactors = FALSE)

t0 <- Sys.time()
raw <- bind_rows(future_lapply(seq_len(nrow(runs)), function(i) {
  sim_once(runs$rep_id[i], runs$k[i], runs$dgp[i])
}, future.seed = TRUE))
message("simulation elapsed: ", format(Sys.time() - t0))


# ---------------------------------------------------------------------------
# summarise
# ---------------------------------------------------------------------------

# A saturated OLS cannot estimate a group with no treated or no control
# observations. Dropping those groups only for OLS would flatter it, so
# accuracy is computed on groups where every model returns an estimate and the
# failures are counted separately.
raw <- raw %>%
  group_by(dgp, k, rep, grp) %>%
  mutate(common = all(!is.na(est))) %>%
  ungroup()

nonestim <- raw %>%
  filter(model == "ols_sat") %>%
  group_by(dgp, k) %>%
  summarize(nonestim_sat = sum(is.na(est)) / n_distinct(rep), .groups = "drop")

ok <- raw %>% filter(common)

mc_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# RMSE is a per-replication quantity, so the Monte Carlo SE comes from the
# spread across replications.
per_rep <- ok %>%
  group_by(dgp, k, num_groups, rep, model) %>%
  summarize(rmse = sqrt(mean((est - tau)^2)), bias = mean(est - tau),
            .groups = "drop")

rmse_tab <- per_rep %>%
  group_by(dgp, k, num_groups, model) %>%
  summarize(rmse_se = mc_se(rmse), rmse = mean(rmse),
            bias_se = mc_se(bias), bias = mean(bias),
            n_reps = n(), .groups = "drop") %>%
  mutate(rmse_lo = rmse - 1.96 * rmse_se, rmse_hi = rmse + 1.96 * rmse_se)

# Type S and type M are conditional on being declared significant, so they are
# pooled over all group-by-replication rows in a cell rather than averaged
# across per-replication means. Averaging per replication would silently drop
# replications in which nothing was significant, conditioning the result on
# significance and biasing type S upward for the most
# conservative estimator.
# Only the two OLS models appear here. Their delta-method intervals are
# correctly calibrated --- measured coverage 0.947 to 0.958 for the saturated
# model --- so sign and exaggeration errors computed from them are honest. The
# additive model's coverage of about 0.55 is NOT a standard error problem but
# genuine misspecification bias: it cannot represent group deviations, so its
# error is mostly bias that a sampling variance never sees. That is precisely
# the overconfidence the type S column reports. The hierarchical model has no
# usable frequentist interval here and is handled in the brms block.
sm_tab <- ok %>%
  filter(!is.na(se)) %>%
  mutate(sig = abs(est) > 1.96 * se) %>%
  group_by(dgp, k, num_groups, model) %>%
  summarize(
    n_sig = sum(sig) / n_distinct(rep),
    frac_sig = mean(sig),
    type_s = sum(sign(est[sig]) != sign(tau[sig])) / sum(sig),
    # Exaggeration ratio. With a single true effect this is E[|est| | sig] /
    # |tau|; with true effects varying by group the natural aggregate is the
    # ratio of means, which is what is reported. The MEAN of the per-group
    # ratio is useless here --- groups whose true effect sits near zero send it
    # to absurd values --- so the median of that ratio is carried alongside as
    # a robustness check. The two agree except under the sparse DGP, where most
    # groups share one modest true effect and the median is the larger number.
    type_m = mean(abs(est[sig])) / mean(abs(tau[sig])),
    type_m_median = median(abs(est[sig]) / abs(tau[sig])),
    # how strongly does declaring significance select on the true effect?
    sel_ratio = mean(abs(tau[sig])) / mean(abs(tau)),
    .groups = "drop"
  )

# Paired differences remove the shared data draw, which is the right way to ask
# whether an RMSE gap is real.
paired <- per_rep %>%
  select(dgp, k, num_groups, rep, model, rmse) %>%
  pivot_wider(names_from = model, values_from = rmse) %>%
  mutate(d_sat = ols_sat - lmer, d_add = ols_add - lmer) %>%
  group_by(dgp, k, num_groups) %>%
  summarize(
    diff_sat_se = mc_se(d_sat), diff_sat = mean(d_sat),
    diff_add_se = mc_se(d_add), diff_add = mean(d_add),
    win_add = mean(d_add > 0), win_sat = mean(d_sat > 0),
    .groups = "drop"
  ) %>%
  mutate(
    sat_lo = diff_sat - 1.96 * diff_sat_se, sat_hi = diff_sat + 1.96 * diff_sat_se,
    add_lo = diff_add - 1.96 * diff_add_se, add_hi = diff_add + 1.96 * diff_add_se,
    # What fraction of the edge over the published straw man survives once the
    # competitor simply stops saturating? Undefined where the saturated model
    # is not actually losing.
    share_pooling = ifelse(sat_lo > 0, diff_add / diff_sat, NA_real_)
  )


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

lab_dgp <- function(x) factor(DGP_LABS[x], levels = unname(DGP_LABS))
lab_mod <- function(x) factor(MODEL_LABS[x], levels = unname(MODEL_LABS))

SHAPES <- c(21, 22, 24)
LINES <- c("dotted", "dashed", "solid")

p_rmse <- rmse_tab %>%
  mutate(dgp = lab_dgp(dgp), model = lab_mod(model)) %>%
  ggplot(aes(x = factor(num_groups), y = rmse, group = model)) +
  facet_wrap(~dgp, nrow = 2) +
  geom_ribbon(aes(ymin = rmse_lo, ymax = rmse_hi), alpha = 0.18) +
  geom_line(aes(linetype = model)) +
  geom_point(aes(shape = model), size = 2.2, fill = "white") +
  scale_shape_manual(values = SHAPES) +
  scale_linetype_manual(values = LINES) +
  labs(
    x = "Number of Groups", y = "RMSE of group treatment effects",
    linetype = NULL, shape = NULL,
    title = "Accuracy of group treatment effect estimates",
    subtitle = sprintf(paste("n = %s; %d Monte Carlo replications"),
                       format(N, big.mark = ","), N_REPS)
  ) +
  theme(legend.position = "bottom")
p_rmse
ggsave("figures/sim-additive-rmse.png", p_rmse, height = 7.5, width = 11)


lab_sat <- "vs. fully crossed OLS"
lab_add <- "vs. additive OLS"

p_decomp <- bind_rows(
  paired %>% transmute(dgp, num_groups, source = lab_sat,
                       est = diff_sat, lo = sat_lo, hi = sat_hi),
  paired %>% transmute(dgp, num_groups, source = lab_add,
                       est = diff_add, lo = add_lo, hi = add_hi)
) %>%
  mutate(dgp = lab_dgp(dgp),
         source = factor(source, levels = c(lab_sat, lab_add))) %>%
  ggplot(aes(x = factor(num_groups), y = est, group = source)) +
  facet_wrap(~dgp, nrow = 2) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18) +
  geom_line(aes(linetype = source)) +
  geom_point(aes(shape = source), size = 2.2, fill = "white") +
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(
    x = "Number of Groups", y = "RMSE reduction from hierarchical model",
    linetype = NULL, shape = NULL,
    title = "Decomposing the hierarchical advantage",
    subtitle = paste("Paired within replication. The distance above zero of",
                     "the solid line is the part of the advantage that",
                     "partial\npooling contributes once the competitor stops",
                     "saturating the interaction.")
  ) +
  theme(legend.position = "bottom")
p_decomp


# The type S / type M figure is built in the brms block below rather than here.
# Two of the three models have calibrated frequentist intervals and one does
# not, so drawing them on the same axes from this tier would compare an exact
# interval against an approximation that covers 0.99.


# --- shrinkage: what partial pooling is actually doing --------------------
# For readers without partial-pooling intuition this is the figure that
# explains the mechanism, and it is worth more than another accuracy curve.
# Each point is one group in one replication: its true effect against the
# estimate. The dashed diagonal is perfect recovery. The solid line is the
# fitted slope of estimate on truth, which IS the shrinkage factor --- a slope
# of one means no pooling, a slope of zero means the model has given every
# group the same answer.
SHRINK_DGP <- "normal"
SHRINK_GROUPS <- c(8, 128)
SHRINK_REPS <- 15 # enough points to see the cloud, few enough to plot

shrink <- ok %>%
  filter(dgp == SHRINK_DGP, num_groups %in% SHRINK_GROUPS, rep <= SHRINK_REPS) %>%
  mutate(model = lab_mod(model),
         panel = factor(paste(num_groups, "groups"),
                        levels = paste(SHRINK_GROUPS, "groups")))

# slope of estimate on truth, per panel: the attenuation the reader can see
slopes_tab <- shrink %>%
  group_by(panel, model) %>%
  summarize(slope = coef(lm(est ~ tau))[2], .groups = "drop") %>%
  mutate(lab = sprintf("slope = %.2f", slope))

p_shrink <- ggplot(shrink, aes(x = tau, y = est)) +
  facet_grid(panel ~ model) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40") +
  geom_point(alpha = 0.12, size = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linewidth = 0.6) +
  geom_text(data = slopes_tab, aes(label = lab), x = -Inf, y = Inf,
            hjust = -0.12, vjust = 1.6, size = 3.4, inherit.aes = FALSE) +
  coord_fixed(xlim = c(-0.8, 1.2), ylim = c(-0.8, 1.2)) +
  labs(
    x = "True group treatment effect",
    y = "Estimated group treatment effect",
    title = "What partial pooling does to individual group estimates",
    subtitle = paste(
      "Dashed line is perfect recovery; the solid line is the fitted slope.",
      "The fully crossed model scatters\nwidely around the truth, the additive",
      "model flattens toward a single common answer, and partial\npooling",
      "trades a little attenuation for a large reduction in scatter.")
  )

p_shrink
ggsave("figures/sim-additive-shrinkage.png", p_shrink, height = 7, width = 10)


# ---------------------------------------------------------------------------
# printed results
# ---------------------------------------------------------------------------

fmt_ci <- function(m, lo, hi) sprintf("%.3f [%.3f, %.3f]", m, lo, hi)

cat("\n\n=================== RMSE (Monte Carlo mean, 95% CI) ===================\n")
print(as.data.frame(
  rmse_tab %>%
    left_join(nonestim, by = c("dgp", "k")) %>%
    transmute(dgp, groups = num_groups, obs_per_group = round(N / num_groups),
              model, rmse = fmt_ci(rmse, rmse_lo, rmse_hi),
              bias = sprintf("%+.3f", bias),
              sat_nonestimable = round(nonestim_sat, 1)) %>%
    arrange(dgp, groups, model)
), row.names = FALSE)

cat("\n\n=========== Paired RMSE difference vs. hierarchical model ===========\n")
cat("Positive = hierarchical model is more accurate. 'share_pooling' is the\n")
cat("fraction of the advantage over the published comparison that survives\n")
cat("once the competitor stops saturating the interaction.\n\n")
print(as.data.frame(
  paired %>%
    transmute(dgp, groups = num_groups,
              vs_fully_crossed = fmt_ci(diff_sat, sat_lo, sat_hi),
              vs_additive = fmt_ci(diff_add, add_lo, add_hi),
              win_rate_vs_additive = sprintf("%.3f", win_add),
              share_pooling = sprintf("%.2f", share_pooling)) %>%
    arrange(dgp, groups)
), row.names = FALSE)

cat("\n\n================= Type S and type M errors =================\n")
cat("Pooled over all group-by-replication rows declared significant, rather\n")
cat("than averaged across per-replication means: averaging would drop the\n")
cat("replications in which nothing was significant, conditioning the result\n")
cat("on the model having found something.\n")
cat("type_m is mean|estimate| / mean|true| among significant estimates;\n")
cat("type_m_med is the median of the per-group ratio. sel_ratio > 1 means\n")
cat("significance selects on groups whose true effects are genuinely larger.\n\n")
print(as.data.frame(
  sm_tab %>%
    transmute(dgp, groups = num_groups, model,
              n_sig = sprintf("%.1f", n_sig),
              type_s = sprintf("%.3f", type_s),
              type_m = sprintf("%.2f", type_m),
              type_m_med = sprintf("%.2f", type_m_median),
              sel_ratio = sprintf("%.2f", sel_ratio)) %>%
    arrange(dgp, groups, model)
), row.names = FALSE)

cat("\n\n=========== Mixed model fitting: singular fits and sigma_theta ===========\n")
print(as.data.frame(
  ok %>%
    filter(model == "lmer") %>%
    group_by(dgp, groups = num_groups) %>%
    summarize(singular_rate = sprintf("%.3f", mean(singular)),
              sigma_theta = sprintf("%.3f", mean(sigma_theta)),
              .groups = "drop")
), row.names = FALSE)


# ---------------------------------------------------------------------------
# optional: validate lme4 against the brms model the paper actually fits, and
# vary the prior on sigma_theta
# ---------------------------------------------------------------------------
#
# Two things need checking before the lme4 result can be reported as a claim
# about the paper's model.
#
# First, does lme4 behave like the brms non-linear model? Second, which prior
# on sigma_theta is actually in force? The appendix table reports half-N(0, 1),
# but none of the scripts in this project set a prior on class = "sd" --- only
# on the two nlpar coefficient blocks --- so brms falls back to its default,
# student_t(3, 0, 2.5). That is a substantially weaker prior than the one
# documented. Both are fitted here so the difference is reported rather than
# discovered by a referee. Reviewer 3's major comment 4 asks for exactly this.

if (RUN_BRMS_CHECK) {
  library(brms)

  BRMS_REPS <- 12
  # k = 3 gives only 8 groups, where sigma_theta is least well identified and
  # the hyperprior has the most room to matter; k = 7 is the many-small-groups
  # case. Those two bracket the range where a prior could plausibly bite.
  BRMS_K <- c(3, 5, 7)
  # normal and systematic were checked first and the prior was irrelevant in
  # both (RMSE differences within 0.001). Sparse and bimodal are the cases
  # where sigma_theta is poorly identified --- lme4 returns singular fits in a
  # quarter of sparse replications at 8 groups --- so they are where the two
  # priors should diverge if they ever do.
  # Type S and type M for all three models are computed here, so this needs to
  # span the DGPs that appear in the paper, not just the two used for the prior
  # check. Budget roughly 12 x 3 x length(BRMS_DGPS) x 2 model fits.
  BRMS_DGPS <- c("normal", "sparse", "bimodal", "additive_only")
  # The prior sensitivity question is settled: across normal, systematic,
  # sparse and bimodal, at 8 and 128 groups, swapping brms's default
  # student_t(3, 0, 2.5) for the half-N(0, 1) the appendix documents moved RMSE
  # by at most 0.001. Running both priors here would double a job that is
  # already several hours, so the second is off by default. Set
  # BRMS_BOTH_PRIORS to TRUE to reproduce the sensitivity check.
  BRMS_BOTH_PRIORS <- FALSE

  PRIOR_SETS <- list(
    default_half_t = NULL, # brms default: student_t(3, 0, 2.5)
    half_normal_1 = prior(normal(0, 1), class = "sd", nlpar = "lambda")
  )
  if (!BRMS_BOTH_PRIORS) PRIOR_SETS <- PRIOR_SETS["default_half_t"]

  brms_one <- function(rep_id, k, dgp) {
    s <- make_data(rep_id, k, dgp)
    x_add <- paste(s$group_cols, collapse = " + ")
    pred <- s$pred
    pred$treat <- 1 # brms validates the whole frame even for one nlpar

    f <- bf(y ~ lambda * treat + controls,
            as.formula(paste("lambda ~", x_add, "+ (1 | grp)")),
            as.formula(paste("controls ~", x_add, "+ z1 + z2 + z3")),
            nl = TRUE)
    base <- c(prior(normal(0, 1), nlpar = "lambda"),
              prior(normal(0, 1), nlpar = "controls"))

    # The two OLS competitors are refit on the same data so that sign and
    # exaggeration errors for all three models come from identical
    # replications, rather than being spliced across tiers with different
    # replication counts.
    p1 <- s$pred; p1$treat <- 1
    p0 <- s$pred; p0$treat <- 0
    x_sat <- paste(s$group_cols, collapse = " * ")
    ols <- bind_rows(lapply(c("ols_sat", "ols_add"), function(m) {
      rhs <- if (m == "ols_sat") x_sat else x_add
      fo <- lm(as.formula(paste("y ~ treat * (", rhs, ") + z1 + z2 + z3")),
               data = s$d)
      e <- as.vector(predict(fo, p1) - predict(fo, p0))
      sd_e <- se_effect_lm(fo, p1, p0)
      data.frame(dgp = dgp, k = k, num_groups = 2^k, rep = rep_id,
                 prior = "n_a", model = m, tau = s$grid$tau, est = e,
                 lo = e - 1.96 * sd_e, hi = e + 1.96 * sd_e,
                 max_rhat = NA_real_, min_ess = NA_real_, divergent = NA_integer_)
    }))

    bind_rows(ols, bind_rows(lapply(names(PRIOR_SETS), function(pn) {
      pr <- base
      if (!is.null(PRIOR_SETS[[pn]])) pr <- pr + PRIOR_SETS[[pn]]
      # cores = chains so the two chains run concurrently. The outer loop stays
      # sequential deliberately: parallelising it would have several workers
      # racing on the same cmdstanr compiled-model cache.
      fit <- brm(f, data = s$d, prior = pr, family = gaussian(),
                 chains = 2, iter = 2000, cores = 2,
                 backend = "cmdstanr", refresh = 0, silent = 2)
      draws <- posterior_epred(fit, nlpar = "lambda", newdata = pred)
      np <- brms::nuts_params(fit)
      data.frame(
        dgp = dgp, k = k, num_groups = 2^k, rep = rep_id, prior = pn,
        model = "hierarchical",
        tau = s$grid$tau, est = apply(draws, 2, median),
        lo = apply(draws, 2, quantile, 0.025),
        hi = apply(draws, 2, quantile, 0.975),
        max_rhat = max(brms::rhat(fit), na.rm = TRUE),
        min_ess = min(summary(fit)$fixed[, "Bulk_ESS"], na.rm = TRUE),
        divergent = sum(np$Value[np$Parameter == "divergent__"])
      )
    })))
  }

  bruns <- expand.grid(rep_id = 1:BRMS_REPS, k = BRMS_K, dgp = BRMS_DGPS,
                       stringsAsFactors = FALSE)

  # This loop runs for hours, so completed iterations are cached to disk and a
  # restart picks up where it left off. The cache is keyed on the settings that
  # affect the answer, so changing the DGPs, the group sizes, or the priors
  # starts a fresh cache rather than silently mixing incompatible runs.
  cache_key <- paste0(paste(sort(BRMS_DGPS), collapse = "-"), "_k",
                      paste(BRMS_K, collapse = "-"), "_r", BRMS_REPS)
  cache <- file.path("data", paste0("brms-check-cache_", cache_key, ".rds"))
  bres_list <- if (file.exists(cache)) readRDS(cache) else list()
  message("brms cache: ", length(bres_list), " of ", nrow(bruns),
          " iterations already done")

  t1 <- Sys.time()
  for (i in seq_len(nrow(bruns))) {
    key <- paste(bruns$dgp[i], bruns$k[i], bruns$rep_id[i], sep = "_")
    if (!is.null(bres_list[[key]])) next
    bres_list[[key]] <- brms_one(bruns$rep_id[i], bruns$k[i], bruns$dgp[i])
    saveRDS(bres_list, cache)
    if (i %% 5 == 0) {
      message("  brms ", i, " / ", nrow(bruns), "  (",
              format(Sys.time() - t1), ")")
    }
  }
  bres <- bind_rows(bres_list)
  message("brms check elapsed: ", format(Sys.time() - t1))

  # brms vs the lme4 run above, on the same (dgp, k, rep) draws
  lme_ref <- per_rep %>%
    filter(model == "lmer", rep <= BRMS_REPS, k %in% BRMS_K,
           dgp %in% BRMS_DGPS) %>%
    select(dgp, k, rep, rmse_lmer = rmse)

  cmp <- bres %>%
    filter(model == "hierarchical") %>%
    group_by(dgp, k, num_groups, rep, prior) %>%
    summarize(rmse = sqrt(mean((est - tau)^2)),
              covered = mean(tau >= lo & tau <= hi), .groups = "drop") %>%
    pivot_wider(names_from = prior, values_from = c(rmse, covered)) %>%
    left_join(lme_ref, by = c("dgp", "k", "rep"))

  cat("\n\n====== brms validation: RMSE against lme4 ======\n")
  cat("lme4 regularizes less than any proper prior, so where it loses to brms\n")
  cat("the Monte Carlo tier is understating the Bayesian model rather than\n")
  cat("flattering it.\n\n")
  print(as.data.frame(
    cmp %>%
      group_by(dgp, groups = num_groups) %>%
      summarize(
        lme4 = sprintf("%.3f", mean(rmse_lmer)),
        brms = sprintf("%.3f", mean(rmse_default_half_t)),
        coverage = sprintf("%.3f", mean(covered_default_half_t)),
        .groups = "drop")
  ), row.names = FALSE)

  if (BRMS_BOTH_PRIORS) {
    cat("\n====== Prior sensitivity: code's prior minus documented prior ======\n")
    print(as.data.frame(
      cmp %>%
        group_by(dgp, groups = num_groups) %>%
        summarize(
          d_m = mean(rmse_default_half_t - rmse_half_normal_1),
          d_se = mc_se(rmse_default_half_t - rmse_half_normal_1),
          .groups = "drop") %>%
        transmute(dgp, groups,
                  prior_gap = fmt_ci(d_m, d_m - 1.96 * d_se, d_m + 1.96 * d_se))
    ), row.names = FALSE)
  }

  # ---- type S and type M, all three models, identical replications --------
  # An estimate counts as significant when its interval excludes zero: the 95%
  # posterior interval for the hierarchical model, the 95% delta-method
  # interval for the two OLS models. Pooled over group-by-replication rows
  # rather than averaged per replication, so replications that find nothing are
  # not silently dropped.
  #
  # This lives here rather than in the lme4 tier because the mixed model has no
  # usable frequentist interval for lambda_g: adding the fixed-effect variance
  # to the varying slope's conditional variance ignores a strongly negative
  # covariance and covers 0.99 against a nominal 0.95, while the conditional
  # variance alone covers 0.84. Posterior intervals are exact.
  sm_brms <- bres %>%
    filter(prior %in% c("n_a", "default_half_t")) %>%
    mutate(
      model = factor(
        dplyr::recode(model,
                      ols_sat = "OLS, fully crossed",
                      ols_add = "OLS, additive interaction",
                      hierarchical = "Hierarchical (partial pooling)"),
        levels = unname(MODEL_LABS)),
      sig = (lo > 0) | (hi < 0)
    ) %>%
    group_by(dgp, num_groups, model) %>%
    summarize(
      coverage = mean(tau >= lo & tau <= hi),
      n_sig = sum(sig) / n_distinct(rep),
      type_s = sum(sign(est[sig]) != sign(tau[sig])) / sum(sig),
      type_m = mean(abs(est[sig])) / mean(abs(tau[sig])),
      type_m_median = median(abs(est[sig]) / abs(tau[sig])),
      .groups = "drop")

  cat("\n====== Type S and type M, from calibrated intervals ======\n")
  cat("coverage checks the intervals themselves and should sit near 0.95.\n")
  cat("The additive model's low coverage is misspecification bias, not a bad\n")
  cat("standard error: its error is mostly bias, which a sampling variance\n")
  cat("never sees. That overconfidence is what the type S column reports.\n\n")
  print(as.data.frame(sm_brms %>%
    transmute(dgp, groups = num_groups, model,
              coverage = sprintf("%.3f", coverage),
              n_sig = sprintf("%.1f", n_sig),
              type_s = sprintf("%.3f", type_s),
              type_m = sprintf("%.2f", type_m),
              type_m_med = sprintf("%.2f", type_m_median))), row.names = FALSE)

  p_err <- sm_brms %>%
    select(dgp, num_groups, model, type_s, type_m) %>%
    pivot_longer(c(type_s, type_m), names_to = "metric", values_to = "value") %>%
    mutate(dgp = lab_dgp(dgp),
           metric = factor(ifelse(metric == "type_s", "Type S: wrong sign",
                                  "Type M: exaggeration ratio"),
                           levels = c("Type S: wrong sign",
                                      "Type M: exaggeration ratio"))) %>%
    ggplot(aes(x = factor(num_groups), y = value, group = model)) +
    facet_grid(metric ~ dgp, scales = "free_y") +
    geom_line(aes(linetype = model)) +
    geom_point(aes(shape = model), size = 2, fill = "white") +
    scale_shape_manual(values = SHAPES) +
    scale_linetype_manual(values = LINES) +
    labs(x = "Number of Groups", y = NULL, linetype = NULL, shape = NULL,
         title = "Errors of inference among subgroup effects declared significant",
         subtitle = paste("Intervals excluding zero: posterior for the",
                          "hierarchical model, delta method for OLS.")) +
    theme(legend.position = "bottom")
  p_err
  ggsave("figures/sim-additive-errors.png", p_err, height = 6, width = 10)

  cat("\n====== Sampler diagnostics ======\n")
  print(as.data.frame(
    bres %>%
      filter(model == "hierarchical") %>%
      group_by(dgp, groups = num_groups, prior) %>%
      summarize(max_rhat = sprintf("%.4f", max(max_rhat)),
                min_ess = round(min(min_ess)),
                reps_with_divergences = n_distinct(rep[divergent > 0]),
                .groups = "drop")
  ), row.names = FALSE)
}
