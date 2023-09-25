# Joshua Alley
# reanalyze Bush and Prather 2020


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
bp.us <- read_dta("data/bush-prather-rep/bush-prather-us-data-all.dta")
glimpse(bp.us)


# modify some variables to match filters in do-file
table(bp.us$treat_germrus)
table(bp.us$treat_ipesidetaking)
table(bp.us$treat_ipesidetaking, bp.us$treat_germrus)
table(bp.us$w2_vote_hill)
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
table(bp.us$w2_vote_hill, bp.us$democrat)
table(bp.us$w2_vote_hill, bp.us$treat_germrus)

# for het effects equation: gender, employment, pol knowlege, pol interest,
# investment/trade, interact Clinton voter dummy (w2_vote_hill) w/ treat_germrus


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
                                          1, 0)
              )

bp.us.key$treat_gr <- bp.us.key %>%
  group_by(
    w2_vote_hill,
    woman, pol_engage, invest_cond
    ) %>%
  group_indices()
class(bp.us.key)
bp.us.key <- sjlabelled::remove_all_labels(bp.us.key)

# simple brms model
formula.bp <- bf(ipe_support ~ 1 + treat_germrus +
                   treat_germrus*(
                   w2_vote_hill + pol_engage +
                   invest_cond +
                   woman) +
                   (1 + treat_germrus | 
                      treat_gr))
bp.mod.vars <- brm(formula.bp, 
                   data = bp.us.key,
                   family = gaussian(link = "identity"),
                   backend = "cmdstanr",
                   cores = 4,
                   control = list(adapt_delta = .9),
                   refresh = 500)
summary(bp.mod.vars)


# predictions 
# new data

grid.het.treat <- bp.us.key %>%
  select(treat_germrus, invest_cond,
         w2_vote_hill, woman, pol_engage,
         treat_gr) %>%
  distinct() 

pred.bp <- predictions(bp.mod.vars, newdata = grid.het.treat) %>%
  posterior_draws()

ggplot(pred.bp, aes(x = draw, y = treat_gr, 
                           fill = factor(treat_germrus))) +
  stat_halfeye(slab_alpha = .5) +
  labs(x = "Predicted Support for Engagement",
       y = "",
       fill = "Interference")

# slopes- create groups
slopes.bp <- slopes(model = bp.mod.vars,
                           variables = "treat_germrus",
                           newdata = grid.het.treat)
slopes.bp

ggplot(slopes.bp, aes(y = estimate, x = pol_engage,
                             color = factor(invest_cond),
                      shape = factor(woman))) +
  facet_wrap(~ w2_vote_hill, ncol = 5,
             labeller = labeller(w2_vote_hill = c(`0` = "Trump Voter",
                                            `1` = "Clinton Voter"))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(width = 1)) +
  scale_color_grey(name = "Econ. Tie",
                   labels = c(`0` = "Trade",
                              `1` = "Investment")) +
  scale_shape_discrete(name = "Gender",
                       labels = c(`0` = "Male",
                                  `1` = "Female")) +
  theme(legend.position = "bottom") +
  labs(title = "Heterogeneous Treatments",
       subtitle = c("Political Affiliation, Engagement, Economic Tie, Gender"),
       x = "Political Engagement", 
       y = "Marginal Effect of Electoral Endorsement")
ggsave("appendix/bp-het-est.png", height = 6, width = 8)

bp.me <- posterior_draws(slopes.bp)

ggplot(bp.me, aes(x = draw, y = treat_gr)) +
  stat_halfeye() +
  labs(x = "Marginal effect of Endorsement", y = "")

