# Joshua Alley
# reanalyze Bush and Prather 2020


# key packages
library(tidyverse)
library(haven)
library(cmdstanr)
library(bayesplot)

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
                      treat_germrus, treat_ipesidetaking,
                      woman, employ_dum, invest_cond,
                      college_educ, democrat, republican,
                      w2_vote_hill) %>%
              drop_na()
# down to 468 obs w/ focus on trump/Clinton voters 

# create a group indicator
# groups with particular combat experiences and demographics 
bp.us.split <- bp.us.key %>%
  group_split(treat_ipesidetaking) 
bp.us.sides <- bp.us.split[[2]]

bp.us.sides$treat_gr <- bp.us.sides %>%
  group_by(
    treat_germrus,
    woman, employ_dum, invest_cond,
    college_educ, 
    w2_vote_hill) %>%
  group_indices()

# return these
bp.us.clean <- bind_rows(bp.us.sides, bp.us.split[[1]])

# number of combat to zero for non-combat
bp.us.clean$treat_gr[bp.us.clean$treat_ipesidetaking == 0] <- 0
bp.us.clean$treat_gr <- bp.us.clean$treat_gr + 1
table(bp.us.clean$treat_gr)

# set up data: het effects matrix
het.bp.mat <- bp.us.clean %>%
  filter(treat_ipesidetaking == 1) %>%
  mutate(
    intercept = 1,
    clinton_germany = treat_germrus * w2_vote_hill,
  ) %>%
  select(intercept, treat_gr,
         treat_germrus, w2_vote_hill, clinton_germany,
         woman, employ_dum, invest_cond,
         college_educ, 
         ) %>%
  arrange(treat_gr) %>%
  distinct() %>%
  select(-treat_gr)

# add zeros for IVs in control group 
het.bp.mat <- rbind(rep(0, times = ncol(het.bp.mat)),
                     het.bp.mat)

# set up data
# data list
data.bp.het <- list(
  N = nrow(bp.us.clean),
  y = bp.us.clean$ipe_support,
  T = max(bp.us.clean$treat_gr),
  treat = bp.us.clean$treat_gr,
  L = ncol(het.bp.mat),
  M = het.bp.mat
)

# compile model
het.mod.bp <- cmdstan_model(stan_file = "data/bush-prather-rep/het-effects-bp.stan",
                         cpp_options = list(stan_threads = TRUE))


# fit model 
fit.het.bp <- het.mod.bp$sample(
  data = data.bp.het,
  chains = 4, 
  parallel_chains = 4,
  threads_per_chain = 2,
  seed = 12,
  #max_treedepth = 20, 
  adapt_delta = .99,
  refresh = 200
)

#  diagnose 
fit.het.bp$cmdstan_diagnose()

diagnostics <- fit.het.bp$diagnostic_summary()
print(diagnostics)

draws.bp <- fit.het.bp$draws(format = "df")

# heterogeneous effects model parameters
mcmc_intervals(draws.bp, regex_pars = "lambda") +
  scale_y_discrete(labels = colnames(data.bp.het$M))

color_scheme_set("gray")
mcmc_areas(draws.bp, regex_pars = "lambda",
               prob = .9) +
  scale_y_discrete(labels = c("Intercept",
                              "German Side-Taking\n(Against Trump)",
                              "Clinton Voter",
                              "Clinton Voter &\nGerman Side-Taking",
                              "Woman",
                              "Employed",
                              "Investment",
                              "College Education")) +
  labs(title = "Predictors of how Side-Taking in US Elections\nImpacts Support for Economic Engagement",
       x = "Estimated Shift in Treatment Effect",
       y = "")
ggsave("figures/bp-lambda.png", height = 6, width = 8)



# treatment group parameters
mcmc_intervals(draws.bp, regex_pars = "theta") 


# summary with predictors
theta.sum.bp <- mcmc_intervals(draws.bp, regex_pars = "theta")$data %>%
  filter(str_detect(parameter, "std|mu|sigma", negate = TRUE)) %>%
  bind_cols(het.bp.mat)


# plot it 
ggplot(theta.sum, aes(x = m, y = hh,
                      color = factor(treat_germrus),
                      #size = n,
                      shape = factor(w2_vote_hill))) +
  facet_wrap(~ invest_cond,
             labeller = labeller(invest_cond = c(`0` = "Trade",
                                                 `1` = "Investment"))) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(y = m, ymin = ll,
                      ymax = hh),
                  alpha = .75,
                  size = 1,
                  linewidth = 1) +
  scale_color_grey(start = .1, end = .6,
                   name = "Side-Taker",
                   labels = c("0" = "Russia: Support Trump", "1" = "Germany: Oppose Trump")) +
  scale_shape_discrete(name = "Voting Intention",
                       labels=c("0" = "Trump", "1" = "Clinton")) +
  labs(title = "Political Preference and Impact of Side-Taking on
       Support for Foreign Economic Engagement",
       x = "Median Estimated Impact of Side-Taking",
       y = "95% Credible Interval") +
  theme(legend.position = "bottom")
ggsave("figures/bp-theta-est.png", height = 6, width = 8)
