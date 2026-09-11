# Joshua Alley
# Tomz and Weeks reanalysis - robustness checks
# see tomz-weeks-reanalysis.R for the main analysis
#
# Addresses four concerns about the application:
#  - is the shift from OLS pooling, or just covariate adjustment?
#  - does the result survive a Bernoulli likelihood? 
#  - how sensitive is it to the prior on the group SD? (R3)
#  - how often would the two models give a researcher opposite signs? (R2)


# load data: appendix data with all controls
tw_rep <- read_dta("data/tomz-weeks-rep/2017-04-YouGov-extracted.dta")
tw_rep <- sjlabelled::remove_all_labels(tw_rep)

# clean data- 0/1 for treatments
tw_rep <- tw_rep %>%
            mutate(
              high_newsint = ifelse(newsint == 1, 1, 0),
              rep = case_when(
                pid7 == 6 | pid7 == 7 ~ 1,
                .default = 0
              ),
              dem = ifelse(pid7 <= 2, 1, 0),
              force = ifelse(pref >= 4, 1, 0),
              white = ifelse(race == 1, 1, 0),
              male = ifelse(gender == 1, 1, 0),
              hawk = abs(6 - hawk),
              ed4 = case_when(educ ==  1 | educ == 2 ~ 1,
                              educ ==  3 | educ == 4 ~ 2,
                              educ ==  5 ~ 3,
                              educ ==  6 ~ 4),
              age = (2016-birthyr) / 10,
              alliance = recode(alliance, `1` = 0, `2` = 1),
              regime = recode(regime, `1` = 0, `2` = 1),
              stakes = recode(stakes, `1` = 1, `2` = 0),
              costs = recode(costs, `1` = 0, `2` = 1),
              region = as.integer(region),
              region_txt = case_when(
                region == 1 ~ "Africa",
                region == 2 ~ "Asia",
                region == 3 ~ "Eastern Europe",
                region == 4 ~ "South America"
              ),
              africa = ifelse(region_txt == "Africa", 1, 0),
              asia = ifelse(region_txt == "Asia", 1, 0),
              europe = ifelse(region_txt == "Eastern Europe", 1, 0),

              treat_group = paste(regime, stakes, costs, region_txt,
                            sep = "_"),
              # coarser grouping: drop region, leaving 8 cells instead of 32
              treat_group_coarse = paste(regime, stakes, costs,
                            sep = "_")
            )
tw_rep$obs_id <- 1:nrow(tw_rep)


# group treatment effects on the response scale, averaged over respondents
# within each cell. For the Gaussian models this equals lambda; for the
# Bernoulli model it is the marginal effect on the probability of supporting
# force, which is what makes the two comparable.
group_effects <- function(model, group_var = "treat_group") {
  d1 <- tw_rep
  d1$alliance <- 1
  d0 <- tw_rep
  d0$alliance <- 0

  diff_draws <- posterior_epred(model, newdata = d1) -
                posterior_epred(model, newdata = d0)

  grp <- tw_rep[[group_var]]
  cells <- sort(unique(grp))
  by_cell <- sapply(cells, function(g) rowMeans(diff_draws[, grp == g,
                                                           drop = FALSE]))

  tibble(
    group = cells,
    est_10 = apply(by_cell, 2, quantile, probs = .1),
    est_med = apply(by_cell, 2, median),
    est_90 = apply(by_cell, 2, quantile, probs = .9)
  )
}

# same quantity from an OLS fit, so the models are compared on one estimand
group_effects_ols <- function(model, group_var = "treat_group") {
  d1 <- tw_rep
  d1$alliance <- 1
  d0 <- tw_rep
  d0$alliance <- 0

  diff <- predict(model, newdata = d1) - predict(model, newdata = d0)
  tibble(group = tw_rep[[group_var]], diff = diff) %>%
    group_by(group) %>%
    summarize(est_med = mean(diff), .groups = "drop")
}


### the two OLS competitors, same modifiers as the hierarchical model
# R3 asks whether the shift comes from partial pooling or from putting several
# correlated modifiers in one equation. The additive model has the same
# covariates and no pooling, so the gap between it and the saturated model is
# specification and the gap between it and the hierarchical model is pooling.
lm_sat <- lm(force ~
               alliance*(regime*stakes*costs*region_txt) +
               white + male + hawk + intl +
               dem + rep + age + ed4,
             data = tw_rep)

lm_add <- lm(force ~
               alliance*(regime + stakes + costs + region_txt) +
               white + male + hawk + intl +
               dem + rep + age + ed4,
             data = tw_rep)


### hierarchical models: the paper's prior, and two alternatives
formula_het_treat <- bf(
  force ~ lambda*alliance + controls,

  lambda ~ regime + stakes + costs + region_txt + (1|treat_group),

  controls ~ white + male + hawk + intl +
    dem + rep + age + ed4,

  nl = TRUE
)

het_treat_prior <- c(
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, 1), nlpar = "controls")
)

# the appendix documents half-normal(0, 1) on the group SD, but no prior was
# ever set on class "sd", so brms used its default student_t(3, 0, 2.5).
# Fitting both settles which one the results depend on.
prior_default <- het_treat_prior
prior_hn1 <- het_treat_prior + prior(normal(0, 1), class = "sd", nlpar = "lambda")
prior_hn5 <- het_treat_prior + prior(normal(0, .5), class = "sd", nlpar = "lambda")

# file = caches the fit, so rerunning the script does not refit five models
fit_het <- function(prior, label, family = gaussian(),
                    formula = formula_het_treat) {
  brm(formula,
      data = tw_rep,
      prior = prior,
      family = family,
      cores = 4,
      backend = "cmdstanr",
      control = list(adapt_delta = .95),
      file = paste0("data/tomz-weeks-rep/fits/tw-", label),
      refresh = 0)
}

dir.create("data/tomz-weeks-rep/fits", showWarnings = FALSE)

tw_treat_default <- fit_het(prior_default, "default")
tw_treat_hn1 <- fit_het(prior_hn1, "hn1")
tw_treat_hn5 <- fit_het(prior_hn5, "hn5")


### Bernoulli likelihood
# The outcome is binary and the paper concedes a logit fits better. The
# question is whether the regularization result depends on the linear
# probability specification.
tw_treat_bern <- fit_het(prior_default, "bern", family = bernoulli())


### coarser grouping
# Grouping is a researcher degree of freedom, so drop region and refit with 8
# cells rather than 32 to see how far the estimates move.
formula_coarse <- bf(
  force ~ lambda*alliance + controls,

  lambda ~ regime + stakes + costs + (1|treat_group_coarse),

  controls ~ white + male + hawk + intl +
    dem + rep + age + ed4,

  nl = TRUE
)

tw_treat_coarse <- fit_het(prior_default, "coarse", formula = formula_coarse)


### compile estimates and compare
est_sat <- group_effects_ols(lm_sat) %>% rename(ols_sat = est_med)
est_add <- group_effects_ols(lm_add) %>% rename(ols_add = est_med)
est_hier <- group_effects(tw_treat_default) %>%
  select(group, hier = est_med, hier_10 = est_10, hier_90 = est_90)
est_hn1 <- group_effects(tw_treat_hn1) %>% select(group, hier_hn1 = est_med)
est_hn5 <- group_effects(tw_treat_hn5) %>% select(group, hier_hn5 = est_med)
est_bern <- group_effects(tw_treat_bern) %>% select(group, hier_bern = est_med)

tw_compare <- est_sat %>%
  left_join(est_add, by = "group") %>%
  left_join(est_hier, by = "group") %>%
  left_join(est_hn1, by = "group") %>%
  left_join(est_hn5, by = "group") %>%
  left_join(est_bern, by = "group")


cat("\n### spread of group estimates by model\n")
print(as.data.frame(
  tw_compare %>%
    summarize(across(c(ols_sat, ols_add, hier, hier_hn1, hier_hn5, hier_bern),
                     list(sd = sd, min = min, max = max))) %>%
    pivot_longer(everything(),
                 names_to = c("model", "stat"), names_sep = "_(?=[a-z]+$)") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
), row.names = FALSE)


cat("\n### specification or pooling?\n")
cat("The additive model has the same covariates as the hierarchical model but\n")
cat("no pooling, so the first step is specification and the second is pooling.\n")
cat("Spread is the right measure here: displacement does not decompose,\n")
cat("because the two steps move individual groups in different directions.\n")
print(as.data.frame(
  tibble(
    step = c("fully crossed OLS", "additive OLS (drops the crossing)",
             "hierarchical (adds pooling)"),
    sd = round(c(sd(tw_compare$ols_sat), sd(tw_compare$ols_add),
                 sd(tw_compare$hier)), 3),
    range = round(c(diff(range(tw_compare$ols_sat)),
                    diff(range(tw_compare$ols_add)),
                    diff(range(tw_compare$hier))), 3)
  )
), row.names = FALSE)

# the specific claim in the manuscript: OLS finds one scenario where alliances
# reduce support by almost 20 points, and the hierarchical model treats it as
# noise. Where does the additive model, which has no pooling, put that group?
extreme <- tw_compare %>% slice_min(ols_sat, n = 3)
cat("\nthe three most negative groups under fully crossed OLS:\n")
print(as.data.frame(extreme %>%
  transmute(group, ols_sat = round(ols_sat, 3), ols_add = round(ols_add, 3),
            hier = round(hier, 3), hier_10 = round(hier_10, 3),
            hier_90 = round(hier_90, 3))), row.names = FALSE)


cat("\n### sign disagreement between models\n")
cat("No ground truth exists here, so these are disagreement rates, not error\n")
cat("rates. Error rates come from the simulation.\n")
print(as.data.frame(
  tw_compare %>%
    summarize(
      n_groups = n(),
      sat_vs_hier = sum(sign(ols_sat) != sign(hier)),
      add_vs_hier = sum(sign(ols_add) != sign(hier)),
      # groups where the saturated model is confident and the hierarchical
      # model puts the opposite sign inside its central interval
      sat_neg_hier_pos = sum(ols_sat < 0 & hier > 0),
      sat_pos_hier_neg = sum(ols_sat > 0 & hier < 0)
    )
), row.names = FALSE)


cat("\n### prior and likelihood sensitivity\n")
print(as.data.frame(
  tw_compare %>%
    summarize(
      hn1_vs_default = max(abs(hier_hn1 - hier)),
      hn5_vs_default = max(abs(hier_hn5 - hier)),
      bern_vs_default = max(abs(hier_bern - hier)),
      cor_hn5 = cor(hier_hn5, hier),
      cor_bern = cor(hier_bern, hier)
    ) %>%
    mutate(across(everything(), ~ round(.x, 4)))
), row.names = FALSE)


### sampler diagnostics
# R2 and R3 both note these are absent from the application.
diag_row <- function(model, label, sd_par = "sd_treat_group__lambda_Intercept") {
  np <- nuts_params(model)
  draws <- as_draws_df(model)
  sm <- posterior::summarise_draws(draws, rhat = posterior::rhat,
                                   ess_bulk = posterior::ess_bulk)
  sm <- sm[is.finite(sm$rhat), ]
  tibble(
    model = label,
    chains = posterior::nchains(draws),
    iter = posterior::niterations(draws),
    max_rhat = round(max(sm$rhat), 4),
    min_ess = round(min(sm$ess_bulk)),
    divergent = sum(np$Value[np$Parameter == "divergent__"]),
    sigma_theta = round(median(posterior::as_draws_matrix(draws)[, sd_par]), 3)
  )
}

cat("\n### sampler diagnostics\n")
print(as.data.frame(bind_rows(
  diag_row(tw_treat_default, "default prior (half-t)"),
  diag_row(tw_treat_hn1, "half-normal(0, 1)"),
  diag_row(tw_treat_hn5, "half-normal(0, .5)"),
  diag_row(tw_treat_bern, "Bernoulli likelihood"),
  diag_row(tw_treat_coarse, "coarse grouping",
           "sd_treat_group_coarse__lambda_Intercept")
)), row.names = FALSE)


### grouping sensitivity
# the coarse model has different cells, so compare the respondent-level
# predictions rather than the cell estimates
cat("\n### coarser grouping (8 cells rather than 32)\n")
coarse_est <- group_effects(tw_treat_coarse, "treat_group_coarse")
print(as.data.frame(coarse_est %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
cat(sprintf("SD of coarse estimates %.3f against %.3f for the 32-cell model\n",
            sd(coarse_est$est_med), sd(tw_compare$hier)))


### model comparison
cat("\n### leave-one-out comparison\n")
tw_treat_default <- add_criterion(tw_treat_default, "loo")
tw_treat_bern <- add_criterion(tw_treat_bern, "loo")
tw_treat_coarse <- add_criterion(tw_treat_coarse, "loo")
print(loo_compare(tw_treat_default, tw_treat_coarse))


### plot: where the models disagree
tw_compare_long <- tw_compare %>%
  select(group, ols_sat, ols_add, hier) %>%
  pivot_longer(cols = -group, names_to = "model", values_to = "estimate") %>%
  mutate(
    model = case_when(
      model == "ols_sat" ~ "OLS, fully crossed",
      model == "ols_add" ~ "OLS, additive interaction",
      model == "hier" ~ "Hierarchical"
    ),
    model = factor(model, levels = c("OLS, fully crossed",
                                     "OLS, additive interaction",
                                     "Hierarchical")),
    
    # take the group variables
    # paste(regime, stakes, costs, region_txt, sep = "_")
    exp_num = str_extract_all(group, "\\d+"),
    
    regime = as.numeric(map_chr(exp_num, 1, .default = NA)),
    stakes = as.numeric(map_chr(exp_num, 2, .default = NA)),
    costs = as.numeric(map_chr(exp_num, 3, .default = NA)),
    
    regime = case_when(
      regime == 0 ~ "Autocracy",
      regime == 1 ~ "Democracy"
    ),
    stakes = case_when(
      stakes == 0 ~ "Low Stakes",
      stakes == 1 ~ "High Stakes"
    ),
    costs = case_when(
      costs == 0 ~ "Low Costs",
      costs == 1 ~ "High Costs"
    ),
    
    region_txt = str_extract(group, "[^_]+$")
  )

ggplot(tw_compare_long, aes(x = estimate, y = reorder(group, estimate),
                            shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(alpha = .7) +
  scale_shape_manual(values = c(1, 2, 16)) +
  labs(
    x = "Estimated Alliance Treatment Effect",
    y = "",
    shape = "",
    title = "Pooling and Specification in Tomz and Weeks"
  ) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 6))


ggplot(tw_compare_long, aes(x = estimate, y = interaction(regime, costs,
                                                     sep = "\n"),
                         shape = model)) +
  facet_grid(stakes ~ region_txt) +
  geom_vline(xintercept = 0) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(1, 2, 16)) +
  labs(
    x = "Alliance Treatment Estimate",
    y = "",
    shape = "",
    title = "Alliance Treatment Estimates\nby Experimental Condition"
  ) +
  theme(legend.position = "bottom")
ggsave("figures/tw-het-treat-source.png", height = 6, width = 8)
