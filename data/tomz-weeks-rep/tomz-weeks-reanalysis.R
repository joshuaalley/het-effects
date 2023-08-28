# Joshua Alley
# reanalyze Tomz and Weeks 2021


# key packages
library(tidyverse)
library(haven)
library(brms)
library(bayesplot)
library(marginaleffects)

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


# results
coef.het.treat <- coef(tw.het.treat)

# all draws
draws.het.treat <- prepare_predictions(tw.het.treat)

# fixed/systematic params
fixed.het.treat <- draws.het.treat[["dpars"]][["mu"]][["fe"]][["b"]]
var.het.treat <- draws.het.treat[["dpars"]][["mu"]][["re"]][["r"]][["regime:stakes:costs:region"]]


# slopes- create groups
slopes.het.treat <- slopes(model = tw.het.treat,
                variables = "alliance",
                 newdata = datagrid(
                   costs = c(0, 1),
                   regime = c(0, 1),
                   stakes = c(0, 1),
                   region = c(1, 2, 3, 4),
                   treat.group = unique(tw.rep$treat.group))
                 )
slopes.het.treat

draws.het.treat <- posterior_draws(slopes.het.treat)

ggplot(draws.het.treat, aes(x = draw, y = treat.group)) +
  ggdist::stat_halfeye() +
  labs(x = "Marginal effect", y = "")


### model with treatment heterogeneity 
treat.het.prior <- c(
  prior(normal(0, .5), class = "b"),
  prior(normal(0, 1), class = "sd")
)
tw.treat.het <- brm(bf(force ~ 1 +
                         regime + stakes + costs + region +
                         alliance*(white + male + hawk + intl) +
                         (1 + alliance | white:male:hawk:intl) ),
                    data = tw.rep,
                    prior = treat.het.prior,
                    family = gaussian(),
                    cores = 4,
                    control = list(adapt_delta = .95),
                    backend = "cmdstanr",
                    refresh = 500
)
summary(tw.treat.het)
