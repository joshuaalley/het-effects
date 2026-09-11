# Joshua Alley
# Simulation: continuous modifiers and the cost of binning
# see sim-additive-comp.R for the binary-modifier simulations
#
# Grouping a continuous modifier means choosing how many bins and where to cut,
# and both are researcher degrees of freedom the paper otherwise leaves open.
# This asks what those choices cost.
#
# The estimand has to change. With binary modifiers there is a finite set of
# groups each with a true effect. With a continuous modifier there are no
# groups until the analyst makes them, and every binning scheme induces a
# different set of "true group effects", so comparing schemes against their own
# induced estimands is circular. Accuracy is therefore measured against the
# true conditional effect at each respondent, tau(x_i), which is invariant to
# the binning choice.
#
# Estimators:
#   ols_bin      saturated OLS on the bins
#   hier_exch    partial pooling across bins, treated as exchangeable - the
#                paper's model applied to a binned modifier
#   hier_smooth  partial pooling that respects the ordering of the bins, via a
#                penalized spline over the bin index
#   ols_cont     x kept continuous and interacted with treatment linearly
#   gam_cont     x kept continuous with a smooth varying slope - no binning
#
# The last two never bin, and are the honest benchmark: if the true effect
# varies smoothly, a model that keeps the modifier continuous should win, and
# the paper should say so rather than let a reader discover it.


library(lme4)
library(mgcv)
library(future.apply)

# brms also exports s(), and setup-script.R has both attached
conflict_prefer("s", "mgcv")

BIN_REPS <- 100
BIN_COUNTS <- c(2, 3, 4, 5, 8)
BIN_RULES <- c("equal_width", "quantile", "oracle")
BIN_SHAPES <- c("linear", "smooth", "threshold")
BIN_N <- 2000
BIN_SIGMA_Y <- 1
BIN_BETA <- c(0.4, -0.3, 0.25)

mc_se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

SHAPE_LABS <- c(
  linear = "Linear",
  smooth = "Smooth (logistic)",
  threshold = "Threshold"
)
RULE_LABS <- c(
  equal_width = "Equal width",
  quantile = "Quantile",
  oracle = "Oracle cuts"
)
BIN_MODEL_LABS <- c(
  ols_bin = "OLS on bins",
  hier_exch = "Hierarchical, exchangeable bins",
  hier_smooth = "Hierarchical, ordered bins",
  ols_cont = "No binning (linear)",
  gam_cont = "No binning (smooth)"
)


### the true conditional effect
# Each shape is scaled to roughly the same marginal spread so that root mean
# squared errors are comparable, matching the convention in the binary
# simulations.
true_tau <- function(x, shape) {
  if (shape == "linear") {
    0.2 + 0.6 * (x - 0.5)
  } else if (shape == "smooth") {
    0.2 + 0.6 * (plogis(8 * (x - 0.5)) - 0.5)
  } else if (shape == "threshold") {
    # a jump at 0.5: bin edges that straddle it cannot represent the effect,
    # however many bins are used
    0.2 + ifelse(x > 0.5, 0.3, -0.3)
  } else {
    stop("unknown shape: ", shape)
  }
}


make_bin_data <- function(rep_id, shape) {
  set.seed(as.integer(3e6 + 1e4 * match(shape, BIN_SHAPES) + rep_id))

  d <- data.frame(
    x = runif(BIN_N),
    treat = rbinom(BIN_N, 1, 0.5),
    z1 = rnorm(BIN_N),
    z2 = rbinom(BIN_N, 1, 0.5),
    z3 = rnorm(BIN_N)
  )
  d$tau <- true_tau(d$x, shape)
  d$y <- d$z1 * BIN_BETA[1] + d$z2 * BIN_BETA[2] + d$z3 * BIN_BETA[3] +
    d$tau * d$treat + rnorm(BIN_N, 0, BIN_SIGMA_Y)
  d
}


# Cut points. The oracle rule cuts on quantiles of the true effect rather than
# of x, which is infeasible in practice and included only as an upper bound on
# what better guidance could buy.
make_bins <- function(d, n_bins, rule, shape) {
  brk <- if (rule == "equal_width") {
    seq(0, 1, length.out = n_bins + 1)
  } else if (rule == "quantile") {
    quantile(d$x, probs = seq(0, 1, length.out = n_bins + 1))
  } else {
    x_at <- sort(d$x)
    tau_at <- true_tau(x_at, shape)
    qs <- quantile(tau_at, probs = seq(0, 1, length.out = n_bins + 1))
    # translate effect quantiles back to cut points on x
    cuts <- sapply(qs, function(q) x_at[which.min(abs(tau_at - q))])
    unique(c(0, sort(cuts[-c(1, length(cuts))]), 1))
  }
  brk <- unique(brk)
  brk[1] <- -Inf
  brk[length(brk)] <- Inf
  bin <- cut(d$x, breaks = brk, labels = FALSE, include.lowest = TRUE)
  d$bin <- factor(bin)
  d$bin_idx <- as.numeric(bin)
  # midpoint of each bin, the systematic part of the heterogeneity equation
  d$bin_mid <- ave(d$x, d$bin, FUN = mean)
  d
}


### estimators
# each returns one estimated effect per respondent, compared against d$tau

eff_from_fit <- function(fit, d, ...) {
  d1 <- d
  d1$treat <- 1
  d0 <- d
  d0$treat <- 0
  as.vector(predict(fit, newdata = d1, ...) - predict(fit, newdata = d0, ...))
}

est_ols_bin <- function(d) {
  fit <- lm(y ~ treat * bin + z1 + z2 + z3, data = d)
  eff_from_fit(fit, d)
}

est_hier_exch <- function(d) {
  fit <- suppressMessages(suppressWarnings(lme4::lmer(
    y ~ treat * bin_mid + z1 + z2 + z3 + (0 + treat | bin),
    data = d, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))))
  eff_from_fit(fit, d, allow.new.levels = TRUE)
}

# Partial pooling that knows the bins are ordered. A penalized spline over the
# bin index is a smoothing prior on the group effects: adjacent bins are pulled
# toward each other rather than all bins toward the grand mean. Needs at least
# a few bins to have a basis, so it is not fitted for the coarsest grids.
est_hier_smooth <- function(d) {
  n_bins <- length(unique(d$bin_idx))
  if (n_bins < 4) return(rep(NA_real_, nrow(d)))
  k <- max(3, min(n_bins, 10))
  fit <- mgcv::gam(
    y ~ treat + s(bin_idx, by = treat, bs = "ps", k = k) +
      s(bin_idx, bs = "ps", k = k) + z1 + z2 + z3,
    data = d, method = "REML")
  eff_from_fit(fit, d)
}

est_ols_cont <- function(d) {
  fit <- lm(y ~ treat * x + z1 + z2 + z3, data = d)
  eff_from_fit(fit, d)
}

est_gam_cont <- function(d) {
  fit <- mgcv::gam(
    y ~ treat + s(x, by = treat) + s(x) + z1 + z2 + z3,
    data = d, method = "REML")
  eff_from_fit(fit, d)
}

# The estimators are called by name rather than looked up in a list. future's
# global detection reads the body of the function it is given and does not
# traverse closures stored in a list, so a list-based dispatch silently fails
# to export the helpers to the workers.
sim_bin_once <- function(rep_id, shape, n_bins, rule) {
  d <- make_bin_data(rep_id, shape)
  d <- make_bins(d, n_bins, rule, shape)

  score <- function(m, est) {
    data.frame(shape = shape, n_bins = n_bins, rule = rule, rep = rep_id,
               model = m,
               rmse = sqrt(mean((est - d$tau)^2)),
               bias = mean(est - d$tau))
  }
  safe <- function(f) tryCatch(f(d), error = function(e) rep(NA_real_, nrow(d)))

  bind_rows(
    score("ols_bin", safe(est_ols_bin)),
    score("hier_exch", safe(est_hier_exch)),
    score("hier_smooth", safe(est_hier_smooth)),
    score("ols_cont", safe(est_ols_cont)),
    score("gam_cont", safe(est_gam_cont))
  )
}


### run
plan(multisession, workers = max(1, parallel::detectCores() - 1))

bin_runs <- expand.grid(rep_id = 1:BIN_REPS, n_bins = BIN_COUNTS,
                        rule = BIN_RULES, shape = BIN_SHAPES,
                        stringsAsFactors = FALSE)

t0 <- Sys.time()
bin_raw <- bind_rows(future_lapply(seq_len(nrow(bin_runs)), function(i) {
  sim_bin_once(bin_runs$rep_id[i], bin_runs$shape[i],
               bin_runs$n_bins[i], bin_runs$rule[i])
}, future.seed = TRUE))
message("binning simulation elapsed: ", format(Sys.time() - t0))


### summarise
n_missing <- bin_raw %>% filter(is.na(rmse)) %>% count(model, n_bins)
if (nrow(n_missing)) {
  cat("\nconditions with no estimate",
      "(the ordered-bin model needs at least four bins)\n")
  print(as.data.frame(n_missing), row.names = FALSE)
}

bin_summary <- bin_raw %>%
  filter(!is.na(rmse)) %>%
  group_by(shape, n_bins, rule, model) %>%
  summarize(rmse_se = mc_se(rmse), rmse = mean(rmse),
            bias = mean(bias), n_reps = n(), .groups = "drop") %>%
  mutate(rmse_lo = rmse - 1.96 * rmse_se, rmse_hi = rmse + 1.96 * rmse_se)

cat("\n=========== RMSE against the true conditional effect ===========\n")
cat("The continuous estimators do not bin, so their accuracy does not vary\n")
cat("with the number of bins; they are repeated across rows as a benchmark.\n\n")
print(as.data.frame(
  bin_summary %>%
    transmute(shape, bins = n_bins, rule, model,
              rmse = sprintf("%.3f [%.3f, %.3f]", rmse, rmse_lo, rmse_hi)) %>%
    arrange(shape, rule, bins, model)
), row.names = FALSE)

cat("\n====== what does the ordering of the bins buy? ======\n")
cat("Exchangeable pooling treats bin 1 and bin 5 as equally similar to bin 2.\n")
cat("Positive means the ordered version is more accurate.\n\n")
print(as.data.frame(
  bin_raw %>%
    filter(!is.na(rmse), model %in% c("hier_exch", "hier_smooth")) %>%
    select(shape, n_bins, rule, rep, model, rmse) %>%
    pivot_wider(names_from = model, values_from = rmse) %>%
    filter(!is.na(hier_smooth), !is.na(hier_exch)) %>%
    group_by(shape, bins = n_bins, rule) %>%
    summarize(d_se = mc_se(hier_exch - hier_smooth),
              d = mean(hier_exch - hier_smooth), .groups = "drop") %>%
    transmute(shape, bins, rule,
              difference = sprintf("%+.3f [%+.3f, %+.3f]",
                                   d, d - 1.96 * d_se, d + 1.96 * d_se)) %>%
    arrange(shape, rule, bins)
), row.names = FALSE)

cat("\n====== best bin count for each shape and rule ======\n")
print(as.data.frame(
  bin_summary %>%
    filter(model %in% c("hier_exch", "hier_smooth", "ols_bin")) %>%
    group_by(shape, rule, model) %>%
    slice_min(rmse, n = 1) %>%
    transmute(shape, rule, model, best_bins = n_bins,
              rmse = sprintf("%.3f", rmse)) %>%
    arrange(shape, rule, model)
), row.names = FALSE)


### figure
p_bin <- bin_summary %>%
  filter(model != "ols_cont") %>%
  mutate(shape = factor(SHAPE_LABS[shape], levels = unname(SHAPE_LABS)),
         rule = factor(RULE_LABS[rule], levels = unname(RULE_LABS)),
         model = factor(BIN_MODEL_LABS[model],
                        levels = unname(BIN_MODEL_LABS))) %>%
  ggplot(aes(x = factor(n_bins), y = rmse, group = model)) +
  facet_grid(rule ~ shape) +
  geom_ribbon(aes(ymin = rmse_lo, ymax = rmse_hi), alpha = .18) +
  geom_line(aes(linetype = model)) +
  geom_point(aes(shape = model), size = 2, fill = "white") +
  scale_shape_manual(values = c(21, 24, 22, 23)) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed", "dotdash")) +
  labs(
    x = "Number of Bins", y = "RMSE against the true conditional effect",
    linetype = NULL, shape = NULL,
    title = "What binning a continuous modifier costs",
    subtitle = sprintf(paste("n = %s; %d Monte Carlo replications. The smooth",
                             "varying slope never bins and is\nrepeated across",
                             "bin counts as a benchmark."),
                       format(BIN_N, big.mark = ","), BIN_REPS)
  ) +
  theme(legend.position = "bottom")

ggsave("figures/sim-binning-comp.png", p_bin, height = 8.5, width = 10)
