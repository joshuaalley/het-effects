# Joshua Alley
# reanalyze Tomz and Weeks 2021




# load data: appendix data with all controls
tw.rep <- read_dta("data/tomz-weeks-rep/2017-04-YouGov-extracted.dta")
glimpse(tw.rep)
tw.rep <- sjlabelled::remove_all_labels(tw.rep)

# clean data- 0/1 for treatments
tw.rep <- tw.rep %>%
            mutate(
              force = ifelse(pref >= 4, 1, 0),
              white = ifelse(race == 1, 1, 0),
              male = ifelse(gender == 1, 1, 0),
              hawk = abs(6 - hawk),
              natl.sup = nat1 + nat2,
              ed4 = case_when(educ ==  1 | educ == 2 ~ 1,
                              educ ==  3 | educ == 4 ~ 2,
                              educ ==  5 ~ 3,
                              educ ==  6 ~ 4),
              age = (2016-birthyr) / 10,
              alliance = factor(recode(alliance, `1` = 0, `2` = 1)),
              regime = factor(recode(regime, `1` = 0, `2` = 1)),
              stakes = factor(recode(stakes, `1` = 1, `2` = 0)),
              costs = factor(recode(costs, `1` = 0, `2` = 1)),
              region = as.integer(region),
              region.txt = case_when(
                region == 1 ~ "Africa",
                region == 2 ~ "Asia",
                region == 3 ~ "Eastern Europe",
                region == 4 ~ "South America"
              ),
              treat.group = paste(regime, stakes, costs, region.txt,
                            sep = "_")
            ) 




### model with treatment heterogeneity
tw.rep$het.group <- paste(tw.rep$white, tw.rep$male,
                          tw.rep$intl, tw.rep$hawk,
                            sep = "_")
treat.het.prior <- c(
  prior(normal(0, .75), class = "Intercept"),
  prior(normal(0, .25), class = "b"),
  prior(normal(0, .75), class = "sd"),
  prior(normal(0, 1), class = "sigma")
)
tw.treat.het <- brm(bf(force ~ 1 +
                         regime + stakes + costs + region.txt +
                         alliance +
                         (1 + alliance | white*male*intl*hawk)
                         
                         #white + male + intl + hawk +
                         #alliance*(white + male + intl + hawk) +
                         #  (1 + alliance | het.group) 
                         
                       #   (1 + alliance | white:male:intl:hawk) +
                       #   (1 + alliance | white:male) +
                       # (1 + alliance | white) +
                       # (1 + alliance | male) +
                       #   (1 + alliance | intl:hawk) +
                       #   (1 + alliance | intl) +
                       # (1 + alliance | hawk) 
                         

                       ),
                    data = tw.rep,
                    prior = treat.het.prior,
                    family = gaussian(),
                    cores = 4,
                    control = list(adapt_delta = .99,
                                   max_treedepth = 20),
                    backend = "cmdstanr",
                    refresh = 500
)
  summary(tw.treat.het)


# parameters
# modelplot(tw.treat.het,
#           coef_map =
#             c("b_alliance:hawk" = "Hawkishness and\nAlliance Impact",
#               "b_alliance:intl" = "Internationalism and\nAlliance Impact",
#               "b_alliance:male" = "Male and\nAlliance Impact",
#               "b_alliance:white" = "White and\nAlliance Impact",
#               "b_alliance" = "Alliance"),
#           size = 1, linewidth = 2 # to geom_pointrange
#           ) +
#   geom_vline(xintercept = 0) +
#   theme_bw(base_size = 14) +
#   labs(title = "Demographic Sources of Heterogeneous Alliance Effects",
#       x = "Estimate and 95% Credible Intervals")



# predictions
# get data 
grid.treat.het <- tw.rep %>%
                      select(alliance, white, male,
                           intl, hawk, het.group) %>%
                      group_by(het.group) %>%
                      mutate(
                        n = n()
                      ) %>%
                      distinct() %>%
                      mutate(
                        regime = 1,
                        stakes = 1,
                        region.txt = "Africa",
                        costs =  1
                      )

grid.treat.het$regime <- factor(grid.treat.het$regime)
grid.treat.het$stakes <- factor(grid.treat.het$stakes)
grid.treat.het$costs <- factor(grid.treat.het$costs)

# number of respondents
treat.het.num <- tw.rep %>%
                  group_by(het.group) %>%
                  summarize(
                    n = n()
                  )

ggplot(treat.het.num, aes(x = n)) +
  geom_bar() +
  ylim(0, 10) +
  labs(y = "Number of Groups",
       x = "Size of Group")

pred.treat.het <- predictions(tw.treat.het,
                              newdata = grid.treat.het) %>%
  posterior_draws()

ggplot(pred.treat.het, aes(x = draw, y = het.group, 
                          color = factor(alliance))) +
  stat_pointinterval() +
  labs(x = "Predicted Support for Force",
       y = "",
       color = "Alliance")

# slopes- create groups
slopes.treat.het <- slopes(model = tw.treat.het,
                           variables = "alliance",
                           newdata = grid.treat.het) 
slopes.treat.het 
slopes.treat.het <- left_join(slopes.treat.het, treat.het.num)

summary(slopes.treat.het$estimate)
sd(slopes.treat.het$estimate)

ggplot(slopes.treat.het, aes(y = estimate, x = factor(male),
                          color = factor(white))) +
  facet_grid(hawk ~ intl,
             labeller = label_both ) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 1)) +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c(`0` = "Female",
                              `1` = "Male")) +
  scale_color_grey(name = "Race",
                   labels = c(`0` = "Non-White",
                              `1` = "White")) +
  scale_shape_discrete(name = "Gender",
                       labels = c(`0` = "Female",
                                  `1` = "Male")) +
  theme(legend.position = "bottom") +
  labs(title = "Alliance Treatment Heterogeneity",
       subtitle = "Internationalism, Hawkishness, Race and Gender",
       x = "Gender", 
       y = "Marginal Effect of Alliance")
ggsave("figures/tw-treat-het.png", height = 8, width = 10)

# number and median
ggplot(slopes.treat.het, aes(x = n, y = estimate,
                             group = het.group)) +
  geom_hline(yintercept = 0) +
    geom_hline(yintercept = median(slopes.treat.het$estimate),
             linetype = "dashed") +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 1)) 

# alternative presentation
treat.het.all <- posterior_draws(slopes.treat.het)

# relevant summary info
summary(slopes.treat.het$estimate)
summary(treat.het.all$draw)
sd(treat.het.all$draw)

ggplot(treat.het.all, aes(x = draw, 
                          fill = het.group,
                          group = het.group)) +
  geom_density(alpha = .25) +
  scale_fill_grey() +
  annotate("text", label = "Mimimum Effect: -0.08\n(-.27, .21)",
           x = -.05, y = 4) +
  annotate("text", label = "Maximum Effect: .55\n(.36, .67)",
           x = .8, y = 4.25) +
  annotate("text", label = "Median Effect: .31\n(.11, .49)",
           x = .31, y = 6) +
  annotate("text", label = "SD of All Draws: .19",
           x = -.4, y = 6) +
  # annotate("text", label = "Median of All Draws: .31",
  #          x = -.4, y = 7.25) +
  # annotate("text", label = "Unexplained Variation: .05 (.00, .13)",
  #          x = -.35, y = 7.25) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  labs(x = "Marginal Effect of Alliance", y = "",
       title = "Posterior Distributions of Alliance Impact")
ggsave("figures/tw-treat-het-sum.png", height = 6, width = 9)


### look at splits by variable
slopes.treat.het.long <- slopes.treat.het %>% 
                          select(
                            estimate, conf.low, conf.high,
                            white, male, intl, hawk
                          ) %>%
                          pivot_longer(cols = -c(conf.low, conf.high,
                                                 estimate),
                                       names_to = "variable") %>%
                          mutate(
                            variable = case_when(
                              variable == "intl" ~ "Internationalism",
                              variable == "hawk" ~ "Militant Assertiveness",
                              variable == "white" ~ "White",
                              variable == "male" ~ "Male"
                            )
                          )

ggplot(slopes.treat.het.long, aes(x = factor(value), y = estimate)) +
  facet_wrap(~ variable, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_jitter(width = .33)) +
  labs(
    y = "Estimate and 95% Credible Interval",
    x = "Modifier Value"
  )


ggplot(slopes.treat.het.long, aes(x = factor(value), y = estimate)) +
  facet_wrap(~ variable, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = .25)) +
  labs(
    y = "Alliance Treatment Estimate",
    x = "Modifier Value",
    title = "Variation in Alliance Impact by Grouping Variable"
  )
ggsave("figures/tw-het-source.png", height = 6, width = 8)


# model: OLS with interactions 
lm.treat.het <- lm(force ~ 
    regime + stakes + costs + region.txt +
    alliance*(white*male*intl*hawk),
  data = tw.rep
)
summary(lm.treat.het)

slopes.lm <- slopes(model = lm.treat.het,
       variables = "alliance",
       newdata = grid.treat.het)

# comparison
slopes.treat.het.comp <- bind_rows(
  "Hierarchical"= slopes.treat.het,
  "OLS Interactions" = slopes.lm,
  .id = "model"
) %>%
  mutate(
    estimate = round(estimate, digits = 2),
    size_n = case_when(
      n <= fivenum(n)[2] ~ "1st Quartile",
      n > fivenum(n)[2] & n <= fivenum(n)[3] ~ "2nd Quartile",
      n > fivenum(n)[3] & n <= fivenum(n)[4] ~ "3rd Quartile",
      n > fivenum(n)[4] ~ "4th Quartile",
    )
  )

ggplot(slopes.treat.het.comp, aes(y = het.group,
                                  color = model,
                                  x = estimate)) +
  facet_wrap(~ model) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  size = .35, linewidth = 1) +
  #scale_color_manual(values = wesanderson::wes_palette("Royal1")) +
  scale_color_grey(start = .2, end = .4) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.y = element_blank()) +
  labs(title = "Alliance Treatment Heterogeneity",
       subtitle = c("Divided By Experimental Group"),
       x = "Estimate and 95% Credible Interval", 
       y = "Respondent Group")
ggsave("figures/tw-treat-het-comp.png", height = 8, width = 8)


# look at regularization by group size
ggplot(slopes.treat.het.comp, aes(y = het.group,
                                  color = model,
                                  x = estimate)) +
  facet_wrap(~ n, scales = "free_y") +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  size = .4, linewidth = 1,
                  position = position_dodge(width = 1)) +
  scale_color_manual(values = wesanderson::wes_palette("Royal1")) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank()) +
  labs(title = "Heterogeneous Alliance Treatments",
       subtitle = c("Divided By Experimental Group"),
       x = "Estimate and 95% Credible Interval", 
       y = "Group")

ggplot(slopes.treat.het.comp, aes(x = n,
                                  group = interaction(model, estimate),
                                  color = model,
                                  y = estimate)) +
  facet_wrap(~ size_n, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  size = .75, linewidth = 1.5,
                  position = position_dodge(width = .5)) +
  scale_color_manual(values = wesanderson::wes_palette("Royal1")) +
  theme(legend.position = "bottom") +
  labs(title = "Heterogeneous Alliance Treatments",
       subtitle = c("Divided By Heterogeneity Group"),
       y = "Estimate and 95% Credible Interval", 
       x = "Group Size")




# check regularization 


# RE with pred
tw.treat.het.pred <- brm(bf(force ~ 1 +
                         regime + stakes + costs + region.txt +
                         alliance +
                           alliance*(white + male + intl + hawk) +
                         (1 + alliance | white*male*intl*hawk)),
data = tw.rep,
prior = treat.het.prior,
family = gaussian(),
cores = 4,
control = list(adapt_delta = .99,
               max_treedepth = 20),
backend = "cmdstanr",
refresh = 500
)
summary(tw.treat.het.re)

slopes.re.pred <- slopes(model = tw.treat.het.pred,
                    variables = "alliance",
                    newdata = grid.treat.het)

# model with predictors, simple RE
tw.treat.het.pred <- brm(bf(force ~ 1 +
                            regime + stakes + costs + region.txt +
                            alliance*(white + male + intl + hawk) +
                            (1 + alliance | white:male:intl:hawk) 
),
data = tw.rep,
prior = treat.het.prior,
family = gaussian(),
cores = 4,
control = list(adapt_delta = .99,
               max_treedepth = 20),
backend = "cmdstanr",
refresh = 500
)
summary(tw.treat.het.pred)

slopes.pred <- slopes(model = tw.treat.het.pred,
                    variables = "alliance",
                    newdata = grid.treat.het)

# comparison- use RE median as baseline 
slopes.comp <- bind_rows("More RE with Group Pred" = slopes.re.pred,
                         "Group Pred and\n(white:male:hawk:intl) RE" = slopes.pred,
                         "More RE: (alliance | white), etc" = slopes.treat.het,
                        "OLS with Inter" = slopes.lm,
                         .id = "model")

ggplot(slopes.comp, aes(y = estimate,
                        x = het.group,
                        color = model,
                        shape = model)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = .25)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ggplot(slopes.comp, aes(y = estimate,
                        x = het.group)) +
  facet_wrap(~ model) +
  geom_hline(yintercept = median(slopes.re$estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = .25)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave("figures/RE-spec.png", height = 6, width = 8)






### model with heterogeneous treatments
het.treat.prior <- c(
  prior(normal(0, .5), class = "b"),
  prior(normal(0, 1), class = "sd")
)
tw.het.treat <- brm(bf(force ~ 1 + white + male + hawk + intl + 
                         pid7 + age + ed4 +
                         (1 + alliance | regime*stakes*costs*region.txt) 
                       
                         # (1 + alliance | regime) +
                         # (1 + alliance | stakes) +
                         # (1 + alliance | costs) +
                         # (1 + alliance | regime:stakes) +
                         # (1 + alliance | regime:costs) +
                         # (1 + alliance | regime:stakes:costs) +
                         # (1 + alliance | costs:stakes) +
                         # (1 + alliance | region.txt) 
                       ),
                    data = tw.rep,
                    prior = het.treat.prior,
                    family = gaussian(),
                    cores = 4,
                    control = list(adapt_delta = .99,
                                   max_treedepth = 15),
                    backend = "cmdstanr",
                    refresh = 500
)
summary(tw.het.treat)

# Using marginaleffects out of the box 
grid.het.treat <- tw.rep %>%
  select(alliance, regime, stakes,
         costs, region.txt, treat.group,
         region) %>%
  distinct() %>%
  mutate(
    white = median(tw.rep$white),
    male = median(tw.rep$male),
    hawk = median(tw.rep$hawk),
    intl =  median(tw.rep$intl),
    pid7 = median(tw.rep$pid7),
    age = median(tw.rep$age),
    ed4 = median(tw.rep$ed4)
  )

# predictions 
pred.het.treat <- predictions(tw.het.treat,
                              newdata = grid.het.treat) %>%
  posterior_draws()

ggplot(pred.het.treat, aes(x = draw, y = treat.group, 
                           fill = factor(alliance))) +
  stat_halfeye(slab_alpha = .5) +
  labs(x = "Predicted Support for Force",
       y = "",
       fill = "Alliance")

# slopes- create groups
slopes.het.treat <- slopes(model = tw.het.treat,
                           variables = "alliance",
                           newdata = grid.het.treat)
slopes.het.treat

ggplot(slopes.het.treat, aes(y = estimate, x = factor(regime),
                             shape = factor(costs),
                             color = factor(stakes))) +
  facet_wrap(~ region.txt, ncol = 5,
             labeller = labeller(regime = c(`0` = "Autocracy",
                                            `1` = "Democracy"))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  size = .75, linewidth = 1.5,
                  position = position_dodge(width = 1)) +
  scale_x_discrete(labels = c(`0` = "Autocracy",
                              `1` = "Democracy")) +
  scale_color_grey(name = "Stakes",
                   labels = c(`0` = "Low",
                              `1` = "High")) +
  scale_shape_discrete(name = "Costs",
                       labels = c(`0` = "High",
                                  `1` = "Low")) +
  labs(title = "Heterogeneous Treatments",
       subtitle = c("Region, Regime, Stakes, Cost"),
       x = "Regime", 
       y = "Marginal Effect of Alliance")
ggsave("appendix/tw-het-treat.png", height = 8, width = 10)

# give posterior mass
het.treat.all <- posterior_draws(slopes.het.treat)

ggplot(het.treat.all, aes(x = draw, y = treat.group)) +
  stat_halfeye() +
  labs(x = "Marginal effect", y = "")



# set up OLS model to compare
tw.het.treat.ols <- lm(force ~ white + male + hawk + intl + 
                         pid7 + age + ed4 +
                         alliance*regime*stakes*costs*region.txt,
                       data = tw.rep)
summary(tw.het.treat.ols)

# slopes- create groups
grid.het.treat.context <- tw.rep %>%
  select(alliance, regime, stakes,
         costs, region.txt, treat.group) %>%
  distinct() %>%
  mutate(
    white = median(tw.rep$white),
    male = median(tw.rep$male),
    hawk = median(tw.rep$hawk),
    intl =  median(tw.rep$intl),
    pid7 = median(tw.rep$pid7),
    age = median(tw.rep$age),
    ed4 = median(tw.rep$ed4)
  )

slopes.het.treat.lm <- slopes(model = tw.het.treat.ols,
                           variables = "alliance",
                           newdata = grid.het.treat.context)

slopes.het.treat.vs <- slopes(model = tw.het.treat,
                              variables = "alliance",
                              newdata = grid.het.treat.context) 
# comparison
slopes.het.treat.comp <- bind_rows(
  "Hierarchical"= slopes.het.treat.vs,
  "OLS Interactions" = slopes.het.treat.lm,
  .id = "model"
) %>%
  mutate(
    estimate = round(estimate, digits = 2)
  )


ggplot(slopes.het.treat.comp, aes(y = estimate, x = model,
                                  label = estimate)) +
  facet_grid(region.txt ~ stakes + regime + costs, 
             labeller = labeller(regime = c(`0` = "Autocracy",
                                     `1` = "Democracy"),
                          stakes = c(`0` = "Low Stakes",
                                     `1` = "High Stakes"),
                          costs = c(`0` = "High Costs",
                                    `1` = "Low Costs"))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  size = .75, linewidth = 1.5) +
  geom_label() +
  labs(title = "Heterogeneous Treatments",
       subtitle = c("Region, Regime, Stakes, Cost"),
       x = "Regime", 
       y = "Marginal Effect of Alliance")
ggsave("figures/tw-het-treat-comp1.png", height = 8, width = 10)


ggplot(slopes.het.treat.comp, aes(y = as.numeric(factor(treat.group)),
                                  x = estimate)) +
  facet_grid(region.txt ~ model) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  size = .75, linewidth = 1.5,
                  position = position_dodge(width = 1)) +
  #scale_color_manual(values = wesanderson::wes_palette("Royal1")) +
  theme(legend.position = "bottom") +
  labs(title = "Heterogeneous Alliance Treatments",
       subtitle = c("Divided By Experimental Group"),
       x = "Estimate and 95% Credible Interval", 
       y = "Treatment Group")
ggsave("figures/tw-het-treat-comp2.png", height = 8, width = 10)
