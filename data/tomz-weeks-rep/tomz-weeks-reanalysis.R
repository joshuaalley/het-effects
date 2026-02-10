# Joshua Alley
# Tomz and Weeks reanalysis - nonlinear lambda specification
# see tw-analysis-interaction.R for interaction-based analysis


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
              dem = ifelse(pid7 >= 2, 1, 0),
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



### respondent heterogeneity- 
formula_pred <- bf(
  force ~ lambda*alliance + controls,

  lambda ~ (dem + rep + high_newsint + white + male) + (1| het_group),

  controls ~ regime + stakes + costs +
              africa + europe + asia +
              age + ed4,

  nl = TRUE
)

pred_prior <- c(
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, .5), nlpar = "controls")
  )

tw_het_nl <- brm(formula_pred,
                    data = tw_rep,
                    prior = pred_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 500
)
summary(tw_het_nl)

# coefficients
coef_tw_het_nl <- as.data.frame(fixef(tw_het_nl, summary = TRUE))
coef_tw_het_nl$variable <- rownames(coef_tw_het_nl)
coef_tw_het_nl <- coef_tw_het_nl %>%
  mutate(
    equation = case_when(
      grepl("lambda", variable) ~ "Alliance Impact",
      grepl("controls", variable) ~ "Controls"
    ),
    variable = case_when(
          variable == "lambda_Intercept" ~ "Intercept",
          variable == "lambda_dem" ~ "Democrat",
          variable == "lambda_rep" ~ "Republican",
          variable == "lambda_pid7" ~ "Party Position",
          variable == "lambda_high_newsint" ~ "High News\nInterest",
          variable == "lambda_white" ~ "White",
          variable == "lambda_male" ~ "Male",
          variable == "lambda_hawk" ~ "Militant\nAssertiveness",
          variable == "lambda_intl" ~ "Internationalism",
          variable == "controls_regime" ~ "Democracy",
          variable == "controls_stakes" ~ "High Stakes",
          variable == "controls_costs" ~ "High Costs",
          variable == "controls_africa" ~ "Africa",
          variable == "controls_europe" ~ "Eastern Europe",
          variable == "controls_asia" ~ "Asia",
          variable == "controls_age" ~ "Age",
          variable == "controls_ed4" ~ "Education",
    variable == "controls_Intercept" ~ "Intercept",
    ),
    variable = factor(variable, levels = c(
      "Intercept", "Democrat", "Republican", "High News\nInterest",
      "White", "Male", "Militant\nAssertiveness", "Internationalism",
      "Democracy",  "High Stakes", "High Costs",
      "Africa", "Eastern Europe", "Asia", "Age", "Education"), ordered = TRUE)
  )

ggplot(coef_tw_het_nl, aes(x = Estimate, y = variable)) +
  facet_wrap(~ equation, scales = "free") +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5)) +
  labs(y = "", x = "Estimate and 95% Interval",
       title = "Determinants of Support for Force and Alliance Impact",
       subtitle = "Hierarchical Model") +
  theme_classic(base_size = 14)
ggsave("figures/tw-treat-het-source.png", height = 6, width = 8)

# lambda estimates
lambda_est <- posterior_epred(tw_het_nl, nlpar = "lambda")
lambda_quant <- apply(lambda_est, 2,
                      function(x) quantile(x, probs = c(.1, .5, .9)))
rownames(lambda_quant) <- c("treat_10", "treat_med", "treat_90")

tw_est <- bind_cols(tw_rep, t(lambda_quant)) %>%
           distinct(
            het_group, treat_med, treat_10, treat_90,
             dem, rep, high_newsint,
             white, male, hawk, intl
           )
length(unique(tw_est$treat_med))

ggplot(tw_est, aes(x = treat_med, y = het_group)) +
  geom_point() +
  labs(x = "Estimated Treatment Effect",
       y = "")

ggplot(tw_est, aes(x = treat_med)) +
  geom_histogram()


### splits by variable
slopes_het_long <- tw_est %>%
  select(
    treat_med, treat_10, treat_90,
    white, male, #intl,
    dem, rep, high_newsint
  ) %>%
  pivot_longer(cols = -c(treat_10, treat_90,
                         treat_med),
               names_to = "variable") %>%
  mutate(
    variable = case_when(
      variable == "intl" ~ "Internationalism",
      variable == "hawk" ~ "Militant Assertiveness",
      variable == "white" ~ "White",
      variable == "male" ~ "Male",
      variable == "dem" ~ "Democrat",
      variable == "rep" ~ "Republican",
      variable == "pid7" ~ "Party ID",
      variable == "high_newsint" ~ "High News Interest"
    )
  )


ggplot(slopes_het_long, aes(x = factor(value), y = treat_med)) +
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
               group_by(white, male, dem, rep, 
                high_newsint) %>%
               mutate(
                  nonzero_interval = ifelse(treat_10 > 0, "Yes", "No")
               ) %>%
               summarise(
                 treat_med = median(treat_med),
                 nonzero_interval = first(nonzero_interval),
                 .groups = "drop") %>%
               mutate(
                party = case_when(
                  rep == 1 ~ "Republican",
                  dem == 1 ~ "Democrat", 
                  .default = "Independent"
                ),
                white = case_when(
                  white == 0 ~ "Non-White",
                  white == 1 ~ "White"),
                male = case_when(
                  male == 0 ~ "Female", 
                  male == 1 ~ "Male")
               )

ggplot(tw_est_sum, aes(y = white, x = male,
                       z = treat_med
                      )) +
  facet_grid(party ~ high_newsint, labeller = labeller(
                high_newsint = c(`0` = "Low News Interest", `1` = "High News Interest"),
                white = c(`0` = "Non-White", `1` = "White"),
                male = c(`0` = "Female", `1` = "Male")
              )
              ) +
  geom_tile(fill = "darkgrey", color = "white") +
  geom_text(aes(label = round(treat_med, 2)), color = "black", size = 5) +
  #scale_fill_grey(start = .6, end = .1,
  #name = "Clear Positive Impact?") +
  theme_classic(base_size = 14) +
  labs(title = "Alliance Impact:",
       subtitle = "Party, Gender, Race, and News Interest",
       x = "Gender",
       y = "Race") +
  theme(legend.position = "bottom",
        legend.key.width = unit(2, "cm"))
ggsave("figures/tw-treat-het-joint.png", height = 6, width = 8)

# comparison with ols
lm_het <- lm(force ~
                     regime + stakes + costs +
                     asia + europe + africa +
                     alliance*(dem + rep + high_newsint + white + male),
                   data = tw_rep
)
summary(lm_het)

slopes_data <- tw_rep %>%
             distinct(
             het_group,
             dem, rep, high_newsint,
             white, male, 
           )   %>%
           mutate(
            regime = 1,
            stakes = 1,
            costs = 0,
            asia = 0,
            africa = 0,
            europe = 1
           )

slopes_lm <- slopes(model = lm_het,
                    variables = "alliance",
                  newdata = slopes_data)

length(unique(slopes_lm$estimate))
ggplot(slopes_lm, aes(x = estimate)) +
  geom_histogram()

# comparison
slopes_het_comp <- tw_est %>%
  left_join(slopes_lm) %>%
  rename("hierarchical" = treat_med,
         "ols" = estimate) %>%
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

disp_char <- ggplot(slopes_het_comp_long, aes(x = te_est)) +
                facet_wrap(~ model) +
                geom_vline(xintercept = 0) +
                geom_hline(yintercept = 0) +
                geom_histogram(color = "white") +
                labs(x = "Treatment Estimate",
                     y = "Estimates")

disp_char / diff_char +
         plot_annotation(
           title = "Comparison of Heterogeneous Treatment Estimates",
           subtitle = "Heterogeneity from Respondent Characteristics")
ggsave("figures/tw-treat-het-comp.png", height = 6, width = 8)


### treatment heterogeneity- nonlinear lambda model
het_treat_prior <- c(
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, 1), nlpar = "controls")
)

formula_het_treat <- bf(
  force ~ lambda*alliance + controls,

  lambda ~ regime*stakes*costs*region_txt + (1|treat_group),

  controls ~ white + male + hawk + intl +
    dem + rep + age + ed4,

  nl = TRUE
)

tw_het_treat <- brm(formula_het_treat,
                    data = tw_rep,
                    prior = het_treat_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    control = list(adapt_delta = .95),
                    refresh = 500)
summary(tw_het_treat)


# lambda estimates
lambda_est_exp <- posterior_epred(tw_het_treat, nlpar = "lambda")
lambda_quant_exp <- apply(lambda_est_exp, 2,
                      function(x) quantile(x, probs = c(.1, .5, .9)))
rownames(lambda_quant_exp) <- c("treat_10", "treat_med", "treat_90")

tw_est_exp <- bind_cols(tw_rep, t(lambda_quant_exp)) %>%
           distinct(treat_group,
          treat_med, treat_10, treat_90,
          regime, stakes, costs, region_txt
          )

ggplot(tw_est_exp, aes(x = treat_med, y = treat_group)) +
  geom_point() +
  labs(x = "Estimated Treatment Effect",
       y = "")


### splits by variable
slopes_het_treat_long <- tw_est_exp %>%
  select(
    treat_med, treat_10, treat_90,
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
  pivot_longer(cols = -c(treat_10, treat_90,
                         treat_med),
               names_to = "variable") %>%
  mutate(
    variable = case_when(
      variable == "regime" ~ "Regime",
      variable == "stakes" ~ "Stakes",
      variable == "costs" ~ "Costs",
      variable == "region_txt" ~ "Region"
    )
  )

ggplot(slopes_het_treat_long, aes(x = factor(value), y = treat_med)) +
  facet_wrap(~ variable, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitter(width = .25), alpha = .5) +
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
                     age = 3.5,
                     ed4 = 2
                    )

slopes_lm_exp <- slopes(model = lm_het_treat,
                    variables = "alliance",
                  newdata = slopes_exp_data)

# comparison
slopes_het_treat_comp <- tw_est_exp %>%
  left_join(slopes_lm_exp) %>%
  rename("hierarchical" = treat_med,
         "ols" = estimate) %>%
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
    title = "Alliance Treatment Estimates\nby Experimental Condition"
  ) +
  theme(legend.position = "bottom")
ggsave("figures/tw-het-treat-source.png", height = 6, width = 8)

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
ggsave("figures/tw-het-treat-comp.png", height = 6, width = 8)
