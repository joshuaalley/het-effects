# Joshua Alley
# reanalyze Tomz and Weeks 2021


# key packages
library(tidyverse)
library(haven)
library(brms)
library(bayesplot)
library(marginaleffects)
library(ggdist)

# set ggplot theme
theme_set(theme_bw(base_size = 14))


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
              alliance = as.integer(recode(alliance, `1` = 0, `2` = 1)),
              regime = as.integer(recode(regime, `1` = 0, `2` = 1)),
              stakes = as.integer(recode(stakes, `1` = 1, `2` = 0)),
              costs = as.integer(recode(costs, `1` = 0, `2` = 1)),
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


# model with heterogeneous treatments
het.treat.prior <- c(
  prior(normal(0, .5), class = "b"),
  prior(normal(0, 1), class = "sd")
)
tw.het.treat <- brm(bf(force ~ 1 + white + male + hawk + intl + natl.sup +
                         alliance*(regime + stakes + costs + region.txt) +
                      (1 + alliance | treat.group) ),
                    data = tw.rep,
                    prior = het.treat.prior,
                    family = gaussian(),
                    cores = 4,
                    control = list(adapt_delta = .95),
                    backend = "cmdstanr",
                    refresh = 500
                    )
summary(tw.het.treat)

# Using marginaleffects out of the box 

grid.het.treat <- tw.rep %>%
  select(alliance, regime, stakes,
         costs, region.txt, treat.group) %>%
  distinct() %>%
  mutate(
    white = median(tw.rep$white),
    male = median(tw.rep$male),
    hawk = median(tw.rep$hawk),
    intl =  median(tw.rep$intl),
    natl.sup = median(tw.rep$natl.sup)
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
       y = "Marginal effect of Alliance")
ggsave("figures/tw-het-treat.png", height = 6, width = 8)

# give posterior mass
het.treat.all <- posterior_draws(slopes.het.treat)

ggplot(het.treat.all, aes(x = draw, y = treat.group)) +
  stat_halfeye() +
  labs(x = "Marginal effect", y = "")



### model with treatment heterogeneity
tw.rep$het.group <- paste(tw.rep$white, tw.rep$male,
                          tw.rep$intl, tw.rep$hawk,
                            sep = "_")
treat.het.prior <- c(
  prior(normal(0, .5), class = "b"),
  prior(normal(0, 1), class = "sd")
)
tw.treat.het <- brm(bf(force ~ 1 +
                         regime + stakes + costs + region.txt +
                         alliance*(white + male + intl + hawk) +
                         (1 + alliance | het.group) ),
                    data = tw.rep,
                    prior = treat.het.prior,
                    family = gaussian(),
                    cores = 4,
                    control = list(adapt_delta = .95),
                    backend = "cmdstanr",
                    refresh = 500
)
summary(tw.treat.het)


# predictions
# get data 
grid.treat.het <- tw.rep %>%
                      select(alliance, white, male,
                           intl, hawk, het.group) %>%
                      distinct() %>%
                      mutate(
                        regime = median(tw.rep$regime),
                        stakes = median(tw.rep$stakes),
                        region.txt = "Africa",
                        costs =  median(tw.rep$costs)
                      )

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
  labs(title = "Treatment Heterogeneity",
       subtitle = "Internationalism, Hawkishness, Race and Gender",
       x = "Gender", 
       y = "Marginal effect of Alliance")
ggsave("figures/tw-treat-het.png", height = 6, width = 8)

# alternative presenation
treat.het.all <- posterior_draws(slopes.treat.het)

ggplot(treat.het.all, aes(x = draw, y = het.group)) +
  stat_pointinterval() +
  labs(x = "Marginal effect of Alliance", y = "")

