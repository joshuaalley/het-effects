# Joshua Alley
# Barber and Pope (2019): four modifiers of a Trump cue, estimated together
#
# Barber and Pope ask how political knowledge, partisanship, Trump approval and
# ideology modify the effect of a Trump cue on policy support. Each modifier
# gets its own interaction model, with the others entered additively. Three of
# the four correlate between .62 and .68, so the separate models cannot say
# which one carries the moderation.
#
# This script fits the separate models, one model with every modifier, and the
# hierarchical model, then compares them in the full sample and out of sample.
#
# Every respondent answers ten policy items, so the models here include a
# respondent intercept and item effects, and the least squares comparisons
# cluster by respondent. The published table treats the 9,000 rows as
# independent, which matters for the liberal cue and not the conservative one.
#
# data: Ideology_Trump.dta from doi:10.7910/DVN/38BFML
# Assumes data/setup-script.R has been sourced.

N_SPLITS <- 200
dir.create("data/barber-pope-rep/fits", showWarnings = FALSE)


### data
bp_raw <- readstata13::read.dta13("data/barber-pope-rep/Ideology_Trump.dta")

# the three arms: Trump takes the liberal position, the conservative position,
# or the policy appears without a cue
bp_arms <- bp_raw %>% filter(libtrump == 1 | contrump == 1 | self == 1)

bp <- bp_arms %>%
  filter(pid7 %in% 1:7) %>%
  mutate(
    party = case_when(pid7 >= 5 ~ "Republican",
                      pid7 <= 3 ~ "Democrat",
                      .default = "Independent"),
    approval = case_when(trump_approve >= 4 ~ "Approves",
                         trump_approve == 3 ~ "Neither",
                         .default = "Disapproves"),
    knows = ifelse(knowledge >= 5, "High knowledge", "Low knowledge"),
    # groups for partial pooling: party by Trump approval by knowledge
    grp = paste(party, approval, knows, sep = "_"),
    grp_coarse = paste(party, approval, sep = "_")
  ) %>%
  drop_na(Support, knowledge, pid7, trump_approve, ideo5b, race_white)

cat(sprintf("\n%d rows, %d respondents, %d items, %d groups (%d in the coarser grouping)\n",
            nrow(bp), length(unique(bp$caseid)), length(unique(bp$Question)),
            length(unique(bp$grp)), length(unique(bp$grp_coarse))))
cat("\ncorrelations among the modifiers, one row per respondent\n")
print(round(cor(bp[!duplicated(bp$caseid), c("knowledge", "pid7", "trump_approve", "ideo5b")]), 2))
cat("\nrespondents per group\n")
print(summary(as.numeric(table(bp$grp[!duplicated(bp$caseid)]))))


### tools
bp_terms <- c("knowledge", "pid7", "trump_approve", "ideo5b")
bp_labels <- c(knowledge = "Political knowledge", pid7 = "Party identification",
               trump_approve = "Trump approval", ideo5b = "Ideology")

# interaction terms for one cue, with respondent-clustered standard errors
int_rows <- function(fit, cue, vc = NULL, label = "") {
  b <- coef(fit)
  if (is.null(vc)) vc <- sandwich::vcovCL(fit, cluster = ~caseid)
  pat <- paste0("(^", cue, ":)|(:", cue, "$)")
  nm <- names(b)[grepl(pat, names(b)) & !is.na(b)]
  se <- sqrt(diag(vc)[nm])
  tibble(model = label, cue = cue, term = gsub(pat, "", nm), est = unname(b[nm]),
         se = unname(se), p = 2 * pnorm(-abs(unname(b[nm]) / unname(se))))
}

# posterior draws of the cue effect for every row, as a difference in expected
# values, so the Gaussian and Bernoulli models sit on one scale
cue_draws <- function(model, cue, data = bp, n_draws = 1000) {
  ids <- sort(sample.int(brms::ndraws(model), min(n_draws, brms::ndraws(model))))
  d1 <- data
  d1$libtrump <- 0
  d1$contrump <- 0
  d1[[cue]] <- 1
  d0 <- data
  d0$libtrump <- 0
  d0$contrump <- 0
  posterior_epred(model, newdata = d1, draw_ids = ids) -
    posterior_epred(model, newdata = d0, draw_ids = ids)
}

# median and interval of the average effect within each set of rows
effect_by <- function(draws, index) {
  tibble(group = names(index),
         n = lengths(index),
         est = vapply(index, function(i) median(rowMeans(draws[, i, drop = FALSE])), numeric(1)),
         lo = vapply(index, function(i) unname(quantile(rowMeans(draws[, i, drop = FALSE]), .025)), numeric(1)),
         hi = vapply(index, function(i) unname(quantile(rowMeans(draws[, i, drop = FALSE]), .975)), numeric(1)))
}

rows_by <- function(d, defs) lapply(defs, function(e) which(eval(e, d)))

# the cue effect within a set of rows, against the no-cue arm
raw_effect <- function(d, idx, cue) {
  vapply(idx, function(i) {
    x <- d[i, ]
    treated <- x[[cue]] == 1
    control <- x$libtrump == 0 & x$contrump == 0
    if (!any(treated) || !any(control)) return(NA_real_)
    mean(x$Support[treated]) - mean(x$Support[control])
  }, numeric(1))
}

bp_subgroups <- list(
  Democrats = quote(party == "Democrat"),
  Independents = quote(party == "Independent"),
  Republicans = quote(party == "Republican"),
  Disapproves = quote(trump_approve <= 2),
  Neither = quote(trump_approve == 3),
  Approves = quote(trump_approve >= 4),
  `Low knowledge` = quote(knowledge <= 4),
  `High knowledge` = quote(knowledge >= 5),
  Liberals = quote(ideo5b <= 2),
  Moderates = quote(ideo5b == 3),
  Conservatives = quote(ideo5b >= 4)
)


### what Barber and Pope report: one modifier per model
bp_published <- list(
  knowledge = lm(Support ~ libtrump*knowledge + contrump*knowledge + trump_approve +
                   ideo5b + republican + party_strength + race_white, data = bp_arms),
  party_strength = lm(Support ~ libtrump*party_strength + contrump*party_strength +
                        knowledge + ideo5b + trump_approve + republican + race_white,
                      data = bp_arms %>% filter(pid7 %in% 4:7)),
  trump_approve = lm(Support ~ libtrump*trump_approve + contrump*trump_approve +
                       knowledge + ideo5b + republican + party_strength + race_white,
                     data = bp_arms),
  ideo5b = lm(Support ~ libtrump*ideo5b + contrump*ideo5b + knowledge + trump_approve +
                republican + party_strength + race_white, data = bp_arms)
)

cat("\n### Table 1 as published, with the standard errors it reports and with clustering\n")
print(as.data.frame(bind_rows(lapply(names(bp_published), function(m) {
  fit <- bp_published[[m]]
  bind_rows(lapply(c("libtrump", "contrump"), function(cue) {
    clustered <- int_rows(fit, cue, label = m)
    plain <- int_rows(fit, cue, vc = vcov(fit), label = m)
    clustered %>% transmute(modifier = m, cue, est = round(est, 4),
                            p_published = round(plain$p, 3), p_clustered = round(p, 3))
  }))
})) %>% arrange(cue, modifier)), row.names = FALSE)


### the same modifiers, one at a time and then together
# item effects and clustered standard errors throughout, so the only thing that
# changes between them is whether the modifiers are estimated jointly
bp_sep <- lapply(setNames(bp_terms, bp_terms), function(m) {
  lm(as.formula(paste("Support ~ (libtrump + contrump) *", m, "+",
                      paste(setdiff(bp_terms, m), collapse = " + "),
                      "+ race_white + factor(Question)")), data = bp)
})
bp_one <- lm(as.formula(paste("Support ~ (libtrump + contrump) * (",
                              paste(bp_terms, collapse = " + "),
                              ") + race_white + factor(Question)")), data = bp)


### the hierarchical model
# The paper writes this as two linked equations, which in brms is the
# non-linear syntax:
#
#   bf(Support ~ lamlib * libtrump + lamcon * contrump + ctrl,
#      lamlib ~ knowledge + pid7 + trump_approve + ideo5b + (1 | grp),
#      lamcon ~ knowledge + pid7 + trump_approve + ideo5b + (1 | grp),
#      ctrl ~ knowledge + pid7 + trump_approve + ideo5b + race_white +
#        (1 | grp) + (1 | caseid) + (1 | Question),
#      nl = TRUE)
#
# Multiplying that out gives the formula below, which is the same model: the
# modifiers enter both equations, each group gets an intercept and a varying
# slope for each cue, and respondents and items get intercepts. With 845
# respondent intercepts the non-linear version samples for hours and this one
# for minutes, so the script fits this and the appendix prints the two-equation
# form.
#
# The group intercept belongs in the outcome equation for the reason the paper
# gives: without it, baseline differences across groups can only show up in the
# cue slopes. Leaving it out puts the liberal cue's group SD at 0.096 rather
# than 0.029, and leave-one-out prefers the model with it by 13.4 (SE 5.8).
formula_bp <- Support ~ (libtrump + contrump) * (knowledge + pid7 + trump_approve + ideo5b) +
  race_white + (1 | caseid) + (1 | Question) + (1 + libtrump + contrump || grp)

formula_coarse <- Support ~ (libtrump + contrump) * (knowledge + pid7 + trump_approve + ideo5b) +
  race_white + (1 | caseid) + (1 | Question) + (1 + libtrump + contrump || grp_coarse)

bp_prior <- prior(normal(0, 1), class = "b")

fit_bp <- function(formula, label, family = gaussian()) {
  brm(formula,
      data = bp,
      prior = bp_prior,
      family = family,
      cores = 4,
      backend = "cmdstanr",
      control = list(adapt_delta = .95),
      file = paste0("data/barber-pope-rep/fits/bp-", label),
      file_refit = "on_change",
      refresh = 0)
}

bp_hier <- fit_bp(formula_bp, "hier")
bp_bern <- fit_bp(formula_bp, "bernoulli", family = bernoulli())
bp_coarse <- fit_bp(formula_coarse, "coarse")


### which modifier carries the moderation?
bp_compare <- bind_rows(lapply(c("libtrump", "contrump"), function(this_cue) {
  sep <- bind_rows(lapply(bp_terms, function(m) {
    int_rows(bp_sep[[m]], this_cue, label = "separate") %>% filter(term == m)
  }))
  one <- int_rows(bp_one, this_cue, label = "one model") %>% filter(term %in% bp_terms)
  fe <- brms::fixef(bp_hier)
  hnm <- paste0(this_cue, ":", bp_terms)
  hier <- tibble(model = "hierarchical", cue = this_cue, term = bp_terms,
                 est = fe[hnm, "Estimate"], se = fe[hnm, "Est.Error"],
                 lo = fe[hnm, "Q2.5"], hi = fe[hnm, "Q97.5"])
  bind_rows(sep, one, hier)
}))

cat("\n### moderator estimates: separately, together, and hierarchically\n")
print(as.data.frame(bp_compare %>%
  transmute(cue, modifier = bp_labels[term], model,
            est = round(est, 4), se = round(se, 4), p = round(p, 3),
            ci = ifelse(is.na(lo), "", sprintf("[%.3f, %.3f]", lo, hi))) %>%
  arrange(cue, modifier, model)), row.names = FALSE)

# the three correlated modifiers, tested together
bp_vc <- sandwich::vcovCL(bp_one, cluster = ~caseid)
bp_nm <- names(coef(bp_one))
bp_wald <- function(terms) {
  b <- coef(bp_one)[terms]
  stat <- as.numeric(t(b) %*% solve(bp_vc[terms, terms]) %*% b)
  sprintf("chi2(%d) = %.2f, p = %.4f", length(terms), stat,
          pchisq(stat, length(terms), lower.tail = FALSE))
}
cat("\n### joint tests in the one model\n")
for (cue in c("libtrump", "contrump")) {
  pick <- function(m) bp_nm[grepl(paste0("(^", cue, ":", m, "$)|(^", m, ":", cue, "$)"), bp_nm)]
  cat(sprintf("%s: approval, ideology and party together %s; all four modifiers %s\n", cue,
              bp_wald(c(pick("trump_approve"), pick("ideo5b"), pick("pid7"))),
              bp_wald(c(pick("knowledge"), pick("trump_approve"), pick("ideo5b"), pick("pid7")))))
}


### how much variation is left after the modifiers?
sd_rows <- function(model, label, group = "grp") {
  dr <- as_draws_df(model)
  pars <- paste0("sd_", group, "__", c("libtrump", "contrump"))
  tibble(model = label, cue = c("liberal cue", "conservative cue"),
         sigma_theta = c(median(dr[[pars[1]]]), median(dr[[pars[2]]])),
         lo = c(unname(quantile(dr[[pars[1]]], .05)), unname(quantile(dr[[pars[2]]], .05))),
         hi = c(unname(quantile(dr[[pars[1]]], .95)), unname(quantile(dr[[pars[2]]], .95))))
}
cat("\n### group standard deviation of the cue effects\n")
print(as.data.frame(bind_rows(
  sd_rows(bp_hier, "main model"),
  sd_rows(bp_coarse, "coarser grouping", "grp_coarse")) %>%
    mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)


### diagnostics
diag_row <- function(model, label) {
  dr <- as_draws_df(model)
  np <- brms::nuts_params(model)
  sm <- posterior::summarise_draws(dr, rhat = posterior::rhat, ess_bulk = posterior::ess_bulk)
  sm <- sm[is.finite(sm$rhat), ]
  tibble(model = label, chains = posterior::nchains(dr), iter = posterior::niterations(dr),
         max_rhat = round(max(sm$rhat), 4), min_ess = round(min(sm$ess_bulk)),
         divergent = sum(np$Value[np$Parameter == "divergent__"]))
}
cat("\n### sampler diagnostics\n")
print(as.data.frame(bind_rows(
  diag_row(bp_hier, "Gaussian, party x approval x knowledge"),
  diag_row(bp_bern, "Bernoulli"),
  diag_row(bp_coarse, "coarser grouping"))), row.names = FALSE)


### cue effects by subgroup and by group
sg <- rows_by(bp, bp_subgroups)
grp_rows <- split(seq_len(nrow(bp)), bp$grp)

bp_effects <- bind_rows(lapply(c("libtrump", "contrump"), function(this_cue) {
  raw <- raw_effect(bp, sg, this_cue)
  effect_by(cue_draws(bp_hier, this_cue), sg) %>%
    mutate(cue = this_cue, raw = raw)
}))
cat("\n### cue effects by subgroup: difference in means and the hierarchical model\n")
print(as.data.frame(bp_effects %>%
  transmute(cue, subgroup = group, rows = n, raw = round(raw, 3),
            hierarchical = round(est, 3), ci = sprintf("[%.3f, %.3f]", lo, hi))),
  row.names = FALSE)

bp_groups <- bind_rows(lapply(c("libtrump", "contrump"), function(this_cue) {
  d1 <- bp
  d1$libtrump <- 0
  d1$contrump <- 0
  d1[[this_cue]] <- 1
  d0 <- bp
  d0$libtrump <- 0
  d0$contrump <- 0
  one_eff <- predict(bp_one, newdata = d1) - predict(bp_one, newdata = d0)
  raw <- raw_effect(bp, grp_rows, this_cue)
  effect_by(cue_draws(bp_hier, this_cue), grp_rows) %>%
    mutate(cue = this_cue,
           one_model = vapply(grp_rows, function(i) mean(one_eff[i]), numeric(1)),
           raw = raw)
}))
cat("\n### spread of the group effects\n")
print(as.data.frame(bp_groups %>% group_by(cue) %>%
  summarize(groups = n(), sd_raw = sd(raw, na.rm = TRUE), sd_one_model = sd(one_model),
            sd_hierarchical = sd(est), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)


### does the likelihood or the grouping change this?
bp_alt <- bind_rows(lapply(c("libtrump", "contrump"), function(this_cue) {
  main <- effect_by(cue_draws(bp_hier, this_cue), sg)$est
  bern <- effect_by(cue_draws(bp_bern, this_cue), sg)$est
  coarse <- effect_by(cue_draws(bp_coarse, this_cue), sg)$est
  tibble(cue = this_cue, max_diff_bernoulli = max(abs(main - bern)),
         cor_bernoulli = cor(main, bern),
         max_diff_coarse = max(abs(main - coarse)), cor_coarse = cor(main, coarse))
}))
cat("\n### subgroup effects under a Bernoulli likelihood and a coarser grouping\n")
print(as.data.frame(bp_alt %>% mutate(across(where(is.numeric), ~ round(.x, 4)))), row.names = FALSE)


### the two parameterizations agree
# bp-hier-nl.rds is the paper's two-equation syntax, kept from an earlier run
# because it samples for hours rather than minutes. It predates the group
# intercept in the outcome equation, so it matches the current model only if it
# carries sd_grp__ctrl_Intercept; the script says so rather than comparing two
# different models silently.
nl_file <- "data/barber-pope-rep/fits/bp-hier-nl.rds"
if (file.exists(nl_file)) {
  bp_nl <- readRDS(nl_file)
  if (!"sd_grp__ctrl_Intercept" %in% names(as_draws_df(bp_nl))) {
    cat("
### NOTE: the cached two-equation fit has no group intercept in the",
        "outcome equation,
### so it is the earlier specification. Comparison",
        "below is for that spec.
")
  }
  fe_nl <- brms::fixef(bp_nl)
  fe_lin <- brms::fixef(bp_hier)
  nl_nm <- c(paste0("lamlib_", bp_terms), paste0("lamcon_", bp_terms))
  lin_nm <- c(paste0("libtrump:", bp_terms), paste0("contrump:", bp_terms))
  cat("\n### two-equation syntax against the multiplied-out formula\n")
  print(as.data.frame(tibble(
    term = lin_nm,
    two_equation = round(fe_nl[nl_nm, "Estimate"], 4),
    multiplied_out = round(fe_lin[lin_nm, "Estimate"], 4))), row.names = FALSE)
  dr_nl <- as_draws_df(bp_nl)
  cat(sprintf("group SD, two-equation syntax: liberal %.4f, conservative %.4f\n",
              median(dr_nl[["sd_grp__lamlib_Intercept"]]),
              median(dr_nl[["sd_grp__lamcon_Intercept"]])))
}


### out of sample: estimate on half the respondents, score on the other half
# The estimand is the set of subgroup effects the paper reports. lme4 stands in
# for brms, as in the simulations, because this is 200 fits per cue.
row_effects <- function(fit, d) {
  d1 <- d
  d1$trt <- 1
  d0 <- d
  d0$trt <- 0
  suppressWarnings(predict(fit, newdata = d1) - predict(fit, newdata = d0))
}

split_once <- function(split_id, d) {
  ids <- unique(d$caseid)
  first <- match(ids, d$caseid)
  strata <- paste(d$grp[first], d$trt[first])
  half <- integer(length(ids))
  for (s in unique(strata)) {
    k <- which(strata == s)
    half[k] <- sample(rep(1:2, length.out = length(k)))
  }
  a <- d[d$caseid %in% ids[half == 1], ]
  b <- d[d$caseid %in% ids[half == 2], ]

  idx_a <- rows_by(a, bp_subgroups)
  idx_b <- rows_by(b, bp_subgroups)
  target <- vapply(idx_b, function(i) {
    x <- b[i, ]
    if (!any(x$trt == 1) || !any(x$trt == 0)) return(NA_real_)
    mean(x$Support[x$trt == 1]) - mean(x$Support[x$trt == 0])
  }, numeric(1))

  fits <- list(party = stats::lm(Support ~ trt * factor(party) + race_white + factor(Question),
                                 data = a))
  for (m in c("knowledge", "trump_approve", "ideo5b")) {
    fits[[m]] <- stats::lm(as.formula(paste("Support ~ trt *", m, "+",
                                            paste(setdiff(bp_terms, m), collapse = " + "),
                                            "+ race_white + factor(Question)")), data = a)
  }
  map_sub <- list(party = c("Democrats", "Independents", "Republicans"),
                  trump_approve = c("Disapproves", "Neither", "Approves"),
                  knowledge = c("Low knowledge", "High knowledge"),
                  ideo5b = c("Liberals", "Moderates", "Conservatives"))
  sep <- setNames(rep(NA_real_, length(idx_a)), names(idx_a))
  for (m in names(map_sub)) {
    e <- row_effects(fits[[m]], a)
    for (nm in map_sub[[m]]) sep[nm] <- mean(e[idx_a[[nm]]])
  }

  one <- stats::lm(as.formula(paste("Support ~ trt * (", paste(bp_terms, collapse = " + "),
                                    ") + race_white + factor(Question)")), data = a)
  hier <- suppressMessages(suppressWarnings(lme4::lmer(
    as.formula(paste("Support ~ trt * (", paste(bp_terms, collapse = " + "),
                     ") + race_white + factor(Question) + (1 | caseid) + (1 + trt || grp)")),
    data = a, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))))
  vc <- as.data.frame(lme4::VarCorr(hier))

  ok <- !is.na(target)
  err <- function(est) mean((est[ok] - target[ok])^2)
  avg <- function(e) vapply(idx_a, function(i) mean(e[i]), numeric(1))
  tibble(split = split_id, separate = err(sep), one_model = err(avg(row_effects(one, a))),
         hierarchical = err(avg(row_effects(hier, a))),
         sigma_theta = vc$sdcor[vc$var1 == "trt" & is.na(vc$var2)][1])
}

future::plan(future::multisession, workers = max(1, min(8, parallelly::availableCores() - 1)))
bp_split <- bind_rows(lapply(c("libtrump", "contrump"), function(this_cue) {
  d <- bp %>%
    filter(.data[[this_cue]] == 1 | (libtrump == 0 & contrump == 0)) %>%
    mutate(trt = .data[[this_cue]])
  bind_rows(future.apply::future_lapply(seq_len(N_SPLITS), split_once, d = d,
                                        future.seed = 12L)) %>%
    mutate(cue = this_cue)
}))
future::plan(future::sequential)

cat("\n### out of sample, by cue\n")
print(as.data.frame(bp_split %>%
  group_by(cue) %>%
  # summary names differ from the column names, so that later expressions
  # still see the per-split values rather than the summaries
  summarize(sep_error = mean(separate), one_model_error = mean(one_model),
            hier_error = mean(hierarchical),
            hier_minus_sep = mean(hierarchical - separate),
            mc_se = sd(hierarchical - separate) / sqrt(n()),
            hier_beats_sep = mean(hierarchical < separate),
            one_beats_sep = mean(one_model < separate),
            hier_beats_one = mean(hierarchical < one_model),
            sigma_theta_zero = mean(sigma_theta < 1e-6), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)


### figure
cue_labels <- c(libtrump = "Liberal cue", contrump = "Conservative cue")
model_levels <- c("separate", "one model", "hierarchical")

panel_a <- bp_compare %>%
  mutate(lo = ifelse(is.na(lo), est - 1.96 * se, lo),
         hi = ifelse(is.na(hi), est + 1.96 * se, hi),
         model = factor(model, levels = model_levels),
         modifier = factor(bp_labels[term], levels = rev(unname(bp_labels))),
         cue = factor(cue_labels[cue], levels = unname(cue_labels))) %>%
  ggplot(aes(x = est, y = modifier, shape = model, linetype = model)) +
  facet_wrap(~ cue) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .6) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .6),
                  fill = "white") +
  scale_shape_manual(values = c(21, 24, 16)) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(x = "Change in the cue effect per unit of the modifier", y = "",
       shape = "", linetype = "",
       title = "Interaction Effects: Barber and Pope 2019") +
  theme(legend.position = "bottom")
panel_a

panel_b <- bp_groups %>%
  filter(cue == "contrump") %>%
  mutate(group = reorder(gsub("_", ", ", group), est)) %>%
  ggplot(aes(y = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .6) +
  geom_point(aes(x = raw, shape = "Difference in means"), alpha = .7) +
  geom_point(aes(x = one_model, shape = "One model")) +
  geom_pointrange(aes(x = est, xmin = lo, xmax = hi, shape = "Hierarchical")) +
  scale_shape_manual(values = c("Difference in means" = 4, "One model" = 2,
                                "Hierarchical" = 16)) +
  labs(x = "Effect of the conservative cue", y = "", shape = "",
       title = "Conservative cue by group") +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 8))
panel_b

panel_a
#+ panel_b + patchwork::plot_layout(heights = c(1, 1.3))
ggsave("figures/bp-modifiers.png", height = 9, width = 8)

