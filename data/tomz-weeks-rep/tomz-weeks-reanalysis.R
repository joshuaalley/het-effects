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
              treat.group = paste(regime, stakes, costs, region,
                            sep = "_")
            ) 


# model with heterogeneous treatments
het.treat.prior <- c(
  prior(normal(0, .5), class = "b"),
  prior(normal(0, 1), class = "sd")
)
tw.het.treat <- brm(bf(force ~ 1 + white + male + hawk + intl + natl.sup +
                         alliance*(regime + stakes + costs + region) +
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
# predictions 
pred.het.treat <- predictions(tw.het.treat,
            newdata = datagrid(model = tw.het.treat,
                       alliance = c(0, 1),
                       treat.group = unique(tw.rep$treat.group))) %>%
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
                           newdata = datagrid(
                             treat.group = unique(tw.rep$treat.group))
)
slopes.het.treat

het.treat.all <- posterior_draws(slopes.het.treat)

ggplot(het.treat.all, aes(x = draw, y = treat.group)) +
  stat_halfeye() +
  labs(x = "Marginal effect", y = "")



### model with treatment heterogeneity
tw.rep$het.group <- paste(tw.rep$white, tw.rep$male, tw.rep$hawk,
                            sep = "_")
treat.het.prior <- c(
  prior(normal(0, .5), class = "b"),
  prior(normal(0, 1), class = "sd")
)
tw.treat.het <- brm(bf(force ~ 1 +
                         regime + stakes + costs + region +
                         alliance*(white + male + hawk) +
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
pred.treat.het <- predictions(tw.treat.het,
                              newdata = datagrid(model = tw.treat.het,
                                                 alliance = c(0, 1),
                                                 het.group = unique(tw.rep$het.group))) %>%
  posterior_draws()

ggplot(pred.treat.het, aes(x = draw, y = het.group, 
                           fill = factor(alliance))) +
  stat_halfeye(slab_alpha = .5) +
  labs(x = "Predicted Support for Force",
       y = "",
       fill = "Alliance")

# slopes- create groups
slopes.treat.het <- slopes(model = tw.treat.het,
                           variables = "alliance",
                           newdata = datagrid(model = tw.treat.het,
                                      het.group = unique(tw.rep$het.group))
)
slopes.treat.het

treat.het.all <- posterior_draws(slopes.treat.het)

ggplot(treat.het.all, aes(x = draw, y = het.group)) +
  stat_halfeye() +
  labs(x = "Marginal effect of Alliance", y = "")
