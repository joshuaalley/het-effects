# Joshua Alley
# Bush and Prather reanalysis - nonlinear lambda specification
# see tomz-weeks-reanalysis.R for parallel analysis


# load data: appendix data with all controls
bp.us <- read_dta("data/bush-prather-rep/bush-prather-us-data-all.dta")
glimpse(bp.us)


# modify some variables to match filters in do-file
bp.us <- bp.us %>%
          mutate(
            # actively employed
            employ_dum = ifelse(employ == 1, 1, 0),
            # investment condition
            invest_cond = ifelse(treat_ipetype == 1, 1, 0),
            # 4-year college or more dummy
            college_educ = ifelse(educ >= 5, 1, 0),
            democrat = ifelse(party == 1, 1, 0),
            republican = ifelse(party == 2, 1, 0)
          )


# pull data and complete cases
bp.us.key <- bp.us %>%
               select(ipe_support,
                      treat_germrus, w2_vote_hill, treat_ipesidetaking,
                      woman, employ_dum, invest_cond,
                      college_educ, democrat, republican,
                      polknowledge, polint) %>%
              drop_na() %>%
              mutate(
                pol_engage = polknowledge + polint,
                high_pol_engage = ifelse(pol_engage >= median(pol_engage, na.rm = TRUE),
                                          1, 0),
                # group for heterogeneity equation
                het_group = paste(w2_vote_hill, woman, high_pol_engage, invest_cond,
                                  sep = "_")
              )
bp.us.key <- sjlabelled::remove_all_labels(bp.us.key)



### respondent heterogeneity - nonlinear lambda specification
formula_bp <- bf(
  ipe_support ~ lambda * treat_germrus + controls,

  lambda ~ (w2_vote_hill + woman + high_pol_engage + invest_cond) + (1 | het_group),

  controls ~ college_educ + employ_dum + woman + high_pol_engage,

  nl = TRUE
)

bp_prior <- c(
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, .5), nlpar = "controls")
)

bp_het_nl <- brm(formula_bp,
                 data = bp.us.key,
                 prior = bp_prior,
                 family = gaussian(),
                 cores = 4,
                 backend = "cmdstanr",
                 refresh = 500
)
summary(bp_het_nl)


# coefficients
coef_bp_het_nl <- as.data.frame(fixef(bp_het_nl, summary = TRUE))
coef_bp_het_nl$variable <- rownames(coef_bp_het_nl)
coef_bp_het_nl <- coef_bp_het_nl %>%
  mutate(
    equation = case_when(
      grepl("lambda", variable) ~ "Side-Taking Impact",
      grepl("controls", variable) ~ "Controls"
    ),
    variable = case_when(
      variable == "lambda_Intercept" ~ "Intercept",
      variable == "lambda_w2_vote_hill" ~ "Clinton Voter",
      variable == "lambda_woman" ~ "Woman",
      variable == "lambda_high_pol_engage" ~ "High Political\nEngagement",
      variable == "lambda_invest_cond" ~ "Investment\nCondition",
      variable == "controls_Intercept" ~ "Intercept",
      variable == "controls_college_educ" ~ "College\nEducation",
      variable == "controls_employ_dum" ~ "Employed",
      TRUE ~ variable
    )
  ) %>%
  filter(equation == "Side-Taking Impact",
         variable != "Intercept")

ggplot(coef_bp_het_nl, aes(x = Estimate, y = variable)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5)) +
  labs(
    x = "Coefficient Estimate",
    y = "",
    title = "Sources of Heterogeneity in Side-Taking Effects"
  )
ggsave("appendix/bp-het-source.png", height = 6, width = 8)


# group-level treatment effects
grid_bp <- bp.us.key %>%
  select(w2_vote_hill, woman, high_pol_engage, invest_cond, het_group) %>%
  distinct() %>%
  mutate(
    treat_germrus = 1,
    college_educ = 0,
    employ_dum = 0
  )

lambda_draws_bp <- posterior_epred(bp_het_nl, nlpar = "lambda",
                                   newdata = grid_bp)
lambda_summary_bp <- data.frame(
  het_group = grid_bp$het_group,
  w2_vote_hill = grid_bp$w2_vote_hill,
  woman = grid_bp$woman,
  high_pol_engage = grid_bp$high_pol_engage,
  invest_cond = grid_bp$invest_cond,
  median = apply(lambda_draws_bp, 2, median),
  lower = apply(lambda_draws_bp, 2, quantile, probs = 0.025),
  upper = apply(lambda_draws_bp, 2, quantile, probs = 0.975)
) %>%
  mutate(
    voter = ifelse(w2_vote_hill == 1, "Clinton Voter", "Trump Voter"),
    gender = ifelse(woman == 1, "Female", "Male"),
    engagement = ifelse(high_pol_engage == 1, "High Engagement", "Low Engagement"),
    econ_tie = ifelse(invest_cond == 1, "Investment", "Trade")
  )


ggplot(lambda_summary_bp, aes(x = median, y = reorder(het_group, median))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  facet_wrap(~ voter, scales = "free_y") +
  labs(
    x = "Effect of German Side-Taking",
    y = "",
    title = "Heterogeneous Effects of Side-Taking on Economic Engagement"
  )


# joint plot by group characteristics
ggplot(lambda_summary_bp, aes(x = median, y = engagement,
                              color = econ_tie, shape = gender)) +
  facet_wrap(~ voter) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper),
                  position = position_dodge(width = 0.5)) +
  scale_color_grey(start = 0.2, end = 0.6, name = "Economic Tie") +
  scale_shape_discrete(name = "Gender") +
  labs(
    x = "Effect of German Side-Taking",
    y = "Political Engagement",
    title = "Heterogeneous Effects by Voter Type"
  ) +
  theme(legend.position = "bottom")
ggsave("appendix/bp-het-joint.png", height = 6, width = 8)
