# Joshua Alley
# reanalyze Tomz and Weeks 2021
# Interaction model specification



# load data: appendix data with all controls
tw_rep <- read_dta("data/tomz-weeks-rep/2017-04-YouGov-extracted.dta")
glimpse(tw_rep)
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
              natl_sup = nat1 + nat2,
              ed4 = case_when(educ ==  1 | educ == 2 ~ 1,
                              educ ==  3 | educ == 4 ~ 2,
                              educ ==  5 ~ 3,
                              educ ==  6 ~ 4),
              age = (2016-birthyr) / 10,
              alliance = recode(alliance, `1` = 0, `2` = 1),
              regime = recode(regime, `1` = 0, `2` = 1),
              stakes = recode(stakes, `1` = 1, `2` = 0),
              costs = recode(costs, `1` = 0, `2` = 1),
              costs_num = as.numeric(costs),
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

              het_group = paste(dem, rep, high_newsint, white, male,
                                  sep = "_"),
              treat_group = paste(regime, stakes, costs, region_txt,
                            sep = "_"),
              treat_group_all = paste(alliance, regime, stakes, costs, region_txt,
                                  sep = "_")
            )
tw_rep$obs_id <- 1:nrow(tw_rep)
# number by group
sort(table(tw_rep$treat_group_all), decreasing = FALSE)
length(unique(tw_rep$treat_group_all))

sort(table(tw_rep$treat_group), decreasing = FALSE)
sort(table(tw_rep$het_group), decreasing = FALSE)
length(unique(tw_rep$het_group))



### model with respondent heterogeneity
formula_inter <- bf(
  force ~ (1 + alliance || het_group) +
    alliance*(dem + rep + high_newsint + white + male + hawk + intl) +
    regime + stakes + costs +
              africa + europe + asia +
              age + ed4)

inter_prior <- c(
  prior(normal(0, 1), class = "sigma"),
  prior(normal(0, .5), class = "b")
  )

tw_het <- brm(formula_inter,
                    data = tw_rep,
                    prior = inter_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 500
)
summary(tw_het)

# grab and plot the coefficients
coef_tw_het <- as.data.frame(fixef(tw_het, summary = TRUE))
coef_tw_het$variable <- rownames(coef_tw_het)
coef_tw_het <- coef_tw_het %>%
  mutate(
    equation = case_when(
      grepl("alliance", variable) ~ "Alliance Impact",
      variable %in% c("regime", "stakes", "costs",
                       "africa", "europe", "asia",
                       "age", "ed4") ~ "Controls",
      .default = NA_character_
    )
  ) %>%
  filter(!is.na(equation)) %>%
  mutate(
    variable = case_when(
          variable == "alliance" ~ "Baseline",
          variable == "alliance:dem" ~ "Democrat",
          variable == "alliance:rep" ~ "Republican",
          variable == "alliance:high_newsint" ~ "High News\nInterest",
          variable == "alliance:white" ~ "White",
          variable == "alliance:male" ~ "Male",
          variable == "alliance:hawk" ~ "Militant\nAssertiveness",
          variable == "alliance:intl" ~ "Internationalism",
          variable == "regime" ~ "Democracy",
          variable == "stakes" ~ "High Stakes",
          variable == "costs" ~ "High Costs",
          variable == "africa" ~ "Africa",
          variable == "europe" ~ "Eastern Europe",
          variable == "asia" ~ "Asia",
          variable == "age" ~ "Age",
          variable == "ed4" ~ "Education",
    ),
    variable = factor(variable, levels = c(
      "Baseline", "Democrat", "Republican", "High News\nInterest",
      "White", "Male", "Militant\nAssertiveness", "Internationalism",
      "Democracy",  "High Stakes", "High Costs",
      "Africa", "Eastern Europe", "Asia", "Age", "Education"), ordered = TRUE)
  )

ggplot(coef_tw_het, aes(x = Estimate, y = variable)) +
  facet_wrap(~ equation, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5)) +
  labs(y = "", x = "Estimate and 95% Interval",
       title = "Determinants of Support and Alliance Impact",
       subtitle = "Hierarchical Model of Respondent Heterogeneity") +
  theme_classic(base_size = 14)


# predicted outcomes
pred_het <- predictions(tw_het,
                              newdata = tw_rep) %>%
  posterior_draws() %>%
  group_by(obs_id) %>%
  select(obs_id, estimate, alliance) %>%
  summarise(across(everything(), list(pred_median = median)))

ggplot(pred_het, aes(x = estimate_pred_median, y = obs_id,
                           color = factor(alliance_pred_median))) +
  geom_point() +
  labs(x = "Predicted Support for Force",
       y = "",
       color = "Alliance")


# estimated slopes of alliance
slopes_het <- slopes(tw_het, variables = "alliance",
                     conf_level = 0.8)

tw_est <- slopes_het %>%
           filter(alliance == 1)

ggplot(tw_est, aes(x = estimate, y = rowid)) +
  geom_point() +
  labs(x = "Estimated Treatment Effect",
       y = "")

ggplot(tw_est, aes(x = estimate)) +
  geom_histogram()



### look at splits by variable
slopes_het_long <- tw_est %>%
  select(
    estimate, conf.low, conf.high,
    white, male, intl, hawk,
    dem, rep, high_newsint
  ) %>%
  pivot_longer(cols = -c(conf.low, conf.high,
                         estimate),
               names_to = "variable") %>%
  mutate(
    variable = case_when(
      variable == "intl" ~ "Internationalism",
      variable == "hawk" ~ "Militant Assertiveness",
      variable == "white" ~ "White",
      variable == "male" ~ "Male",
      variable == "dem" ~ "Democrat",
      variable == "rep" ~ "Republican",
      variable == "high_newsint" ~ "High News Interest"
    )
  )


ggplot(slopes_het_long, aes(x = factor(value), y = estimate)) +
  facet_wrap(~ variable, scales = "free_x") +
  geom_hline(yintercept = 0) +
    geom_boxplot(outlier.shape = NA) +
  labs(
    y = "Alliance Treatment Estimate",
    x = "Modifier Value",
    title = "Variation in Alliance Impact by Grouping Variable"
  )

# joint impacts
tw_est_sum <- tw_est %>%
               group_by(white, male, intl, hawk, high_newsint) %>%
               mutate(
                  nonzero_interval = ifelse(conf.low > 0, "Yes", "No")
               ) %>%
               summarise(
                 estimate = median(estimate),
                 nonzero_interval = first(nonzero_interval),
                 .groups = "drop")

ggplot(tw_est_sum, aes(x = factor(intl), y = factor(hawk),
                       z = estimate
                      )) +
  facet_grid(white + male ~ high_newsint, labeller = labeller(
                high_newsint = c(`0` = "Low News Interest", `1` = "High News Interest"),
                white = c(`0` = "Non-White", `1` = "White"),
                male = c(`0` = "Female", `1` = "Male")
              )
              ) +
  geom_tile(aes(fill = factor(nonzero_interval)), color = "white") +
  geom_text(aes(label = round(estimate, 2)), color = "white", size = 5) +
  scale_fill_grey(start = .6, end = .1,
  name = "Clear Positive Impact?") +
  theme_classic(base_size = 14) +
  labs(title = "Alliance Impact:",
       subtitle = "Foreign Policy Disposition, Gender, Race, and News Interest",
       x = "Internationalism",
       y = "Hawkishness") +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))



# comparison with ols
# model: OLS with interactions
lm_het <- lm(force ~
                     regime + stakes + costs +
                     asia + europe + africa +
                     alliance*(white + male + intl + hawk +
                                dem + rep + high_newsint),
                   data = tw_rep
)
summary(lm_het)

slopes_lm <- slopes(model = lm_het,
                    variables = "alliance")
hist(slopes_lm$estimate)

# comparison
slopes_het_comp <- tibble(
  hierarchical = slopes_het$estimate,
  ols = slopes_lm$estimate,
  alliance = slopes_het$alliance
) %>%
  mutate(
    diff = hierarchical - ols
  )
length(unique(slopes_het_comp$hierarchical))
length(unique(slopes_het_comp$ols))


slopes_het_comp_long <- slopes_het_comp %>%
  pivot_longer(
    names_to = "model",
    values_to = "te_est",
    cols = c(hierarchical, ols),
  ) %>%
  mutate(
    te_est = te_est * alliance,
    model = str_to_upper(model)
  )


ggplot(slopes_het_comp, aes(y = hierarchical, x = ols)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point()

slopes_diff_char_sum <- slopes_het_comp %>%
  group_by(ols) %>%
  summarise(
    mean_diff = mean(diff),
    median_diff = median(diff),
    sd_diff = sd(diff)
  ) %>%
  ungroup()

diff_char <- ggplot(slopes_diff_char_sum, aes(y = median_diff, x = ols)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point() +
  labs(
    x = "OLS Treatment Estimate",
    y = "(Hierarchical - OLS)"
  )
diff_char

disp_char <- ggplot(filter(slopes_het_comp_long, alliance == 1),
   aes(x = te_est)) +
  facet_wrap(~ model) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_histogram(color = "white") +
    labs(x = "Treatment Estimate",
       y = "Estimates",
  )

disp_char / diff_char +
  plot_annotation(
    title = "Comparison of Heterogeneous Treatment Estimates",
    subtitle = "Heterogeneity from Respondent Characteristics"
  )





### model with heterogeneous treatments
treat_inter_prior <- c(
  prior(normal(0, 1), class = "sigma"),
  prior(normal(0, 1), class = "b")
)

formula_treat_inter <- bf(
  force ~ (1 + alliance || treat_group) +
    alliance*(regime*stakes*costs*region_txt) +
    white + male + hawk + intl +
    dem + rep + age + ed4)

tw_het_treat <- brm(formula_treat_inter,
                    data = tw_rep,
                    prior = treat_inter_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 500
)
summary(tw_het_treat)


slopes_exp_data <- tw_rep %>%
                    distinct(treat_group,
                    regime, stakes, costs, region_txt) %>%
                    mutate(
                     dem = 0,
                     rep = 0,
                     male = 0,
                     white = 1,
                     intl = 3,
                     hawk = 2,
                     age = 4.5,
                     ed4 = 3
                    )


# estimated slopes
tw_est_exp <- slopes(tw_het_treat, variables = "alliance",
                       conf_level = 0.9,
                      newdata = slopes_exp_data)

ggplot(tw_est_exp, aes(x = estimate, y = treat_group)) +
  geom_point() +
  labs(x = "Estimated Treatment Effect",
       y = "")


### look at splits by variable
slopes_het_treat_long <- tw_est_exp %>%
  select(
    estimate, conf.low, conf.high,
    regime, stakes, costs, region_txt
  ) %>%
  mutate(
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
    )
  ) %>%
  pivot_longer(cols = -c(conf.low, conf.high,
                         estimate),
               names_to = "variable") %>%
  mutate(
    variable = case_when(
      variable == "regime" ~ "Regime",
      variable == "stakes" ~ "Stakes",
      variable == "costs" ~ "Costs",
      variable == "region_txt" ~ "Region"
    )
  )

ggplot(slopes_het_treat_long, aes(x = factor(value), y = estimate)) +
  facet_wrap(~ variable, scales = "free_x") +
  geom_hline(yintercept = 0) +
    geom_point(position = position_jitter(width = .25), alpha = .5) +
  geom_boxplot(outlier.shape = NA) +
  labs(
    y = "Alliance Treatment Estimate",
    x = "Modifier Value",
    title = "Variation in Alliance Impact Across Experimental Conditions"
  )



# comparison with ols
lm_het_treat <- lm(force ~
                     alliance*(regime*stakes*costs*region_txt) +
                         white + male + hawk + intl +
                        dem + rep + age + ed4,
                   data = tw_rep)
summary(lm_het_treat)

slopes_lm_treat <- slopes(model = lm_het_treat,
                    variables = "alliance",
                  newdata = slopes_exp_data)

# comparison
slopes_het_treat_comp <- tw_est_exp %>%
      rename("hierarchical" = estimate) %>%
      select(treat_group, hierarchical, regime, stakes, costs, region_txt) %>%
  left_join(select(slopes_lm_treat, estimate, treat_group)) %>%
    rename("ols" = estimate) %>%
  mutate(
    diff = hierarchical - ols
  )
length(unique(slopes_het_treat_comp$hierarchical))
length(unique(slopes_het_treat_comp$ols))


slopes_het_treat_comp_long <- slopes_het_treat_comp %>%
  pivot_longer(
    names_to = "model",
    values_to = "te_est",
    cols = c(hierarchical, ols),
  ) %>%
  mutate(
    model = str_to_upper(model),
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
    )
  )


tw_fig3_comp <- slopes_het_treat_comp_long %>%
                   group_by(regime, stakes, costs, model, region_txt) %>%
                   filter(te_est != 0) %>%
                   summarize(
                    te_est = mean(te_est),
                    .groups = "drop"
                   )

# rough approx of TW figure 3
ggplot(tw_fig3_comp, aes(x = te_est, y = interaction(regime, costs,
sep = "\n"),
             color = model)) +
  facet_grid(stakes ~ region_txt) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  scale_color_grey(start = .1, end = .6) +
  labs(
    x = "Alliance Treatment Estimate",
    y = "Experimental Condition",
    color = "Model Type",
    title = "Average Alliance Treatment Estimates\nby Experimental Condition"
  ) +
  theme(legend.position = "bottom")


slopes_diff_treat_sum <- slopes_het_treat_comp %>%
  group_by(ols) %>%
  summarise(
    mean_diff = mean(diff),
    median_diff = median(diff),
    sd_diff = sd(diff)
  ) %>%
  ungroup()



diff_treat <- ggplot(slopes_diff_treat_sum, aes(y = median_diff, x = ols)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point() +
  labs(
    x = "OLS Treatment Estimate",
    y = "(Hierarchical - OLS)"
  )
diff_treat

disp_treat <- ggplot(slopes_het_treat_comp_long,
                                  aes(x = te_est)) +
  facet_wrap(~ model) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_histogram() +
  labs(x = "Treatment Estimate",
       y = "Estimates"
  )


disp_treat / diff_treat +
  plot_annotation(
    title = "Comparison of Heterogeneous Treatment Estimates",
    subtitle = "Heterogeneity from Experimental Conditions"
  )

