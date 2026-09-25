# Todd Sechser and Joshua Alley
# Analysis of YouGov Experiment


# check balance
balance.check <- glm(casualties_dummy ~ educ + age + republican + democrat +
                       hawkishness + isolationism + nationalism,
                     data = twdata)
summary(balance.check)


# check balance
balance.check.pass <- glm(casualties_dummy ~ educ + age + republican + democrat +
                       hawkishness + isolationism + nationalism,
                     data = twdata.pass)
summary(balance.check.pass)

check.models <- list(balance.check, balance.check.pass)
names(check.models) <- c("Full Sample", "Passed Manip.\nChecks")

modelplot(check.models,
          coef_rename = c("(Intercept)" = "Intercept",
                          "educ" = "Education",
                          "age" = "Age",
                          "republican" = "Republican",
                          "democrat" = "Democrat",
                          "hawkishness" = "Mil. Assertiveness",
                          "isolationism" = "Isolationism",
                          "nationalism" = "Natl. Chauvinism"),
          coef_omit = "Intercept") +
         labs(title = "Randomization Check",
              color = "Sample")
ggsave("appendix/randomization-check.png", height = 6, width = 8)



# results
# tabulate support for force
table(twdata$use.force, twdata$casualties)
table(twdata$use.force, twdata$casualties_dummy)
table(twdata$use.force, twdata$ally)
table(twdata$use.force, twdata$nuclear)
table(twdata$use.force, twdata$democracy)


# raw by treatment 
ggplot(twdata,
       aes(x = treat.text, y = use.force,
           group = factor(casualties),
           fill = factor(casualties))) +
  geom_col() +
  labs(y = "Number of Military Intervention Supporters",
       x = "",
       title = "Military Intervention Support by Condition",
       fill = "Casualties") +
  scale_fill_manual(values = wes_palette("GrandBudapest1")) +
  theme(legend.position = "bottom") +
  coord_flip() 
ggsave("appendix/all-treat-sum.png", height = 6, width = 8)


### Logit models of support for using force
# all treatments
interv.model <- glm(use.force ~ casualties_dummy +
                      nuclear + ally + democracy,
                    family = binomial(link = "logit"),
                    data = twdata.pass)
summary(interv.model)
modelplot(interv.model,
          coef_rename = c("(Intercept)" = "Intercept",
                          "casualties_dummy" = "US Casualties", 
                          "nuclear" = "Nuclear Adversary",
                          "ally" = "Ally",
                          "democracy" = "Democracy")) +
  geom_vline(xintercept = 0) +
  labs(
    title = "Estimated Treatment Effects: YouGov Experiment"
  )


# treatments with casualties magnitude
interv.model.num <- glm(use.force ~ casualties_9 + casualties_50 + casualties_250 +
                      nuclear + ally + democracy,
                    family = binomial(link = "logit"),
                    data = twdata.pass)
summary(interv.model.num)
modelplot(interv.model.num,
          coef_rename = c("(Intercept)" = "Intercept",
                          "casualties_9" = "9 US Casualties", 
                          "casualties_50" = "50 US Casualties", 
                          "casualties_250" = "250 US Casualties", 
                          "nuclear" = "Nuclear Adversary",
                          "ally" = "Ally",
                          "democracy" = "Democracy")) +
  geom_vline(xintercept = 0) +
  labs(
    title = "Estimated Treatment Effects with Number of Casualties: YouGov Experiment"
  )


# casulaties dummy w/ controls 
interv.model.control <- glm(use.force ~ casualties_dummy +
                      nuclear + ally + democracy + female +
                      educ + age + republican + democrat +
                      hawkishness + isolationism + nationalism,
                    family = binomial(link = "logit"),
                    data = twdata.pass)
summary(interv.model.control)




### policy choice: for appendix


# simple descriptive plot
table(twdata.pass$full.res, twdata.pass$casualties_dummy)
cas.policy <- data.frame(table(twdata.pass$full.res, twdata.pass$casualties_dummy))
colnames(cas.policy) <- c("policy.res", "Casualties", "number")

# create frequency variable 
cas.policy$perc <- NA
cas.policy$perc[cas.policy$Casualties == 1] <- 
                         cas.policy$number[cas.policy$Casualties == 1] / 
                               sum(twdata.pass$casualties_dummy)
cas.policy$perc[cas.policy$Casualties == 0] <- 
  cas.policy$number[cas.policy$Casualties == 0] / 
  (nrow(twdata.pass) - sum(twdata.pass$casualties_dummy))

# casualties text
cas.policy$cas.text <- ifelse(cas.policy$Casualties == 1, 
                              "Yes", "No")

# recode policy res
cas.policy$policy.res <- recode(cas.policy$policy.res,
                      "Do nothing" = "Nothing",
                      "Lodge a diplomatic protest" = "Diplomatic Protest",
                      "Impose economic sanctions on the attacker" = "Sanctions",
                      "Launch airstrikes against Country A's forces" = "Airstrikes",
                      "Send ground troops to attack Country A's forces" = "Ground Troops",
                      "Use nuclear weapons against Country A's forces" = "Nuclear")

# military vs not for facet
cas.policy$mil <- factor(ifelse(cas.policy$policy.res == "Airstrikes" |
                                      cas.policy$policy.res == "Ground Troops" |
                                      cas.policy$policy.res == "Nuclear",
                                    "Military", "Non-Military"),
                             ordered = TRUE,
                             levels = c("Non-Military", "Military"))


cas.policy <- cas.policy %>%
  group_by(mil) %>%
  mutate(
    n.mil = sum(number, na.rm = TRUE),
    prop.mil = number / n.mil
  )


# raw count 
ggplot(cas.policy, aes(y = number,
                       x = policy.res,
                      fill = Casualties)) + 
  geom_bar(position="dodge", stat = "identity") +
  scale_fill_grey(limits = rev) +
  labs(x = "Preffered Response",
       y = "Number of Respondents")

ggplot(cas.policy, aes(y = prop.mil,
                       x = policy.res,
                       fill = cas.text)) + 
  facet_wrap(~ mil, scales = "free_x") +
  geom_bar(position="dodge", stat = "identity") +
  scale_fill_grey(limits = rev) +
  labs(x = "Response",
       fill = "Casualties",
       y = "Percentage of Respondents in Treatment Group")
ggsave("appendix/policy-res-cas.png", height = 6, width = 8)


# Plot with number of casualties
table(twdata.pass$full.res, twdata.pass$casualties)
cas.policy.num <- data.frame(table(twdata$full.res, twdata$casualties))
colnames(cas.policy.num) <- c("policy.res", "Casualties", "number")
cas.policy.num$prop.cas <- cas.policy.num$number / 800 # 800 per casualties treatment

# recode policy res
cas.policy.num$policy.res <- recode(cas.policy.num$policy.res,
                                "Do nothing" = "Nothing",
                                "Lodge a diplomatic protest" = "Diplomatic Protest",
                                "Impose economic sanctions on the attacker" = "Sanctions",
                                "Launch airstrikes against Country A's forces" = "Airstrikes",
                                "Send ground troops to attack Country A's forces" = "Ground Troops",
                                "Use nuclear weapons against Country A's forces" = "Nuclear")

# military vs not for facet
cas.policy.num$mil <- factor(ifelse(cas.policy.num$policy.res == "Airstrikes" |
                            cas.policy.num$policy.res == "Ground Troops" |
                            cas.policy.num$policy.res == "Nuclear",
                            "Military", "Non-Military"),
                            ordered = TRUE,
                            levels = c("Non-Military", "Military"))
cas.policy.num$casualties_dummy <- ifelse(cas.policy.num$Casualties != 0, "Yes", "No")

cas.policy.num <- cas.policy.num %>%
                   group_by(mil, Casualties) %>%
                    mutate(
                      n.mil = sum(number, na.rm = TRUE),
                      prop.mil = number / n.mil
                    )

# raw count 
ggplot(cas.policy.num, aes(y = number,
                       x = policy.res,
                       fill = Casualties)) +
  facet_wrap(~ mil, scales = "free_x") +
  geom_bar(position="dodge", stat = "identity") +
  scale_fill_grey(limits = rev) +
  labs(title = "Casualties and Preferred Invasion Response",
       x = "Policy",
       y = "Number of Respondents")
ggsave("appendix/policy-res-ncas.png", height = 6, width = 8)




# proportion of each treatment group
ggplot(cas.policy.num, aes(y = prop.cas,
                           x = policy.res,
                           fill = Casualties)) +
  facet_wrap(~ mil, scales = "free_x") +
  geom_bar(position="dodge", stat = "identity") +
  scale_fill_grey(limits = rev) +
  labs(title = "Casualties and Preferred Invasion Response",
       x = "Policy",
       y = "Share of Respondents in Casulaties Treatment Group")
ggsave("appendix/policy-res-cas.png", height = 6, width = 8)



# proportion within each policy
ggplot(cas.policy.num, aes(y = prop.mil,
                           x = policy.res,
                           fill = Casualties)) +
  facet_wrap(~ mil, scales = "free_x") +
  geom_bar(position="dodge", stat = "identity") +
  scale_fill_grey(limits = rev) +
  labs(title = "Casualties and Preferred Invasion Response",
       x = "Policy",
       y = "Share of Respondents in Military Choice Group")



# sum proportions
cas.policy.num %>% 
  group_by(policy.res) %>%
   summarize(
     sum.prop = sum(number) / 3200,
     .groups = "keep"
   )




# military policy response
table(twdata.pass$mil.res, twdata.pass$casualties)
table(twdata.pass$mil.res, twdata.pass$casualties_dummy)

# all these analyses focus on respondents who selected mil/non-mil
# due to list-wise deletion on the outcome variable 

# assume equal
mil.action.reg <- lm(as.numeric(mil.res) ~ casualties_dummy + ally + 
                       nuclear + democracy,
                      data = twdata.pass)
summary(mil.action.reg)

# assume equal: different dummies
mil.action.regn <- lm(as.numeric(mil.res) ~ casualties_9 + casualties_50 + casualties_250 + 
                       ally + 
                       nuclear + democracy,
                     data = twdata.pass)
summary(mil.action.regn)

# ordinal regression model
mil.action.ord <- polr(mil.res ~ casualties_dummy + 
                         ally + 
                         nuclear + democracy,
                       data = twdata.pass)
summary(mil.action.ord)
brant(mil.action.ord)




# generalized ordinal model
gen.milact <- clm(mil.res ~ democracy + ally + nuclear,
                      nominal = ~ casualties_dummy,
                      data = twdata.pass)
summary(gen.milact)
gen.mil.choice <- tidy(gen.milact)
class(gen.mil.choice) <- "data.frame"
# relative risk
exp(coef(gen.milact))


# generalized ordinal model: split casualties
gen.milact.ncas <- clm(mil.res ~ democracy + 
                         ally + nuclear,
                      nominal = ~ casualties_9 + casualties_50 + casualties_250 ,
                      data = twdata.pass)
summary(gen.milact.ncas)
mil.choice.ncas <- broom::tidy(gen.milact.ncas)
class(mil.choice.ncas) <- "data.frame"
# relative risk
exp(coef(gen.milact.ncas))


# military action models plot
mil.action.models <- list(mil.action.reg, mil.action.ord)
names(mil.action.models) <- c("OLS", "Ordinal\nLogit")
modelplot(mil.action.models,
          coef_map = c(
                          "casualties_dummy" = "US Casualties",
                          # "casualties_9" = "9 US Casualties", 
                          # "casualties_50" = "50 US Casualties", 
                          # "casualties_250" = "250 US Casualties", 
                          "nuclear" = "Nuclear Adversary",
                          "ally" = "Alliance",
                          "democracy" = "Democracy")) +
  labs(
    color = "Model",
    title = "Models of Military Policy Choice"
  )
ggsave("appendix/mil-choice-models.png", height = 6, width = 8)



# nonmilitary policy response
table(twdata$nonmil.res, twdata$casualties)
table(twdata$nonmil.res, twdata$casualties_dummy)

# assume equal
nonmil.action.reg <- lm(as.numeric(nonmil.res) ~ casualties_dummy + ally + 
                       nuclear + democracy,
                     data = twdata.pass)
summary(nonmil.action.reg)

# assume equal: different dummies
nonmil.action.regn <- lm(as.numeric(nonmil.res) ~ casualties_9 + casualties_50 + casualties_250 + 
                        ally + 
                        nuclear + democracy,
                      data = twdata.pass)
summary(nonmil.action.regn)

# ordinal regression model
nonmil.action.ord <- polr(nonmil.res ~ casualties_dummy + ally + 
                         nuclear + democracy,
                       data = twdata.pass)
summary(nonmil.action.ord)
brant(nonmil.action.ord)

# ordinal regression model
nonmil.action.ordn <- polr(nonmil.res ~ casualties_9 + casualties_50 + casualties_250 + 
                             ally + 
                            nuclear + democracy,
                          data = twdata.pass)
summary(nonmil.action.ordn)
brant(nonmil.action.ordn)

# generalized ordinal model
gen.nonmilact <- clm(nonmil.res ~ democracy + nuclear + ally,
                  nominal = ~ casualties_dummy,
                  data = twdata.pass)
summary(gen.nonmilact)
gen.nonmil.choice <- broom::tidy(gen.nonmilact)
class(gen.nonmil.choice) <- "data.frame"
# relative risk
exp(coef(gen.nonmilact))


# generalized ordinal model: split casualties
gen.nonmilact.ncas <- clm(nonmil.res ~ democracy +
                            ally + nuclear,
                       nominal = ~ casualties_9 + casualties_50 + casualties_250,
                       data = twdata.pass)
summary(gen.nonmilact.ncas)
nonmil.choice.ncas <- broom::tidy(gen.nonmilact.ncas)
class(nonmil.choice.ncas) <- "data.frame"
# relative risk
exp(coef(gen.nonmilact.ncas))


# military action models plot
nonmil.action.models <- list(nonmil.action.reg, nonmil.action.ord)
names(nonmil.action.models) <- c("OLS", "Ordinal\nLogit")
modelplot(nonmil.action.models,
          coef_map = c(
                          "casualties_dummy" = "US Casualties",
                          # "casualties_9" = "9 US Casualties", 
                          # "casualties_50" = "50 US Casualties", 
                          # "casualties_250" = "250 US Casualties", 
                          "nuclear" = "Nuclear Adversary",
                          "ally" = "Alliance",
                          "democracy" = "Democracy")) +
  labs(
    color = "Model",
    title = "Models of Non-Military Policy Choice"
  )
ggsave("appendix/nonmil-choice-models.png", height = 6, width = 8)


### treatment interactions: 

# subgroup means, please
twdata %>%
   group_by(treat.text) %>%
   summarize(
     use_force_perc = sum(use.force, na.rm = TRUE) / n()
   ) %>%
  ggplot(aes(x = use_force_perc, y = treat.text)) +
  geom_point() +
  labs(x = "Percentage Supporting Force",
       y = "")


interv.model.inter <- brm(use.force ~ 
                            #casualties_dummy*(nuclear + democracy + ally) +
                            casualties_9*(nuclear + democracy + ally) +
                            casualties_50*(nuclear + democracy + ally) +
                            casualties_250*(nuclear + democracy +ally) +
                        female +
                        educ + age + republican + democrat +
                        hawkishness + isolationism + nationalism +
                         (1 | treat.text.het),
                    family = bernoulli(link = "logit"),
                    backend = "cmdstanr",
                    chains = 4, cores = 4,
                    data = twdata.pass)
summary(interv.model.inter)


# grid of modifier values for the marginal effects below.
# treat.text.het has to be rebuilt from ally/nuclear/democracy: datagrid() holds
# unspecified variables at their mode, which pins every row to "Alliance | Nuclear |"
# and leaves 7 of 8 profiles internally inconsistent. Because slope = "eydx" on a
# logit depends on the baseline probability, that mismatch distorts the estimates.
het.datagrid <- function(model){
  datagrid(model = model,
           nuclear = c(0, 1),
           democracy = c(0, 1),
           ally = c(0, 1)) %>%
    mutate(
      treat.text.het = paste(ifelse(ally == 1, "Alliance |", ""),
                             ifelse(nuclear == 1, "Nuclear |", ""),
                             ifelse(democracy == 1, "Democracy", ""),
                             sep = " ")
    )
}


slopes.inter.pass <- slopes(interv.model.inter,
                            variables = c("casualties_9",
                                          "casualties_50",
                                          "casualties_250"),
                            slope = "eydx",
                      newdata = het.datagrid(interv.model.inter))


# full sample
interv.model.inter.full <- brm(use.force ~ 
                                 casualties_9*(nuclear + democracy + ally) +
                                 casualties_50*(nuclear + democracy + ally) +
                                 casualties_250*(nuclear + democracy +ally) +
                                 female +
                                 educ + age + republican + democrat +
                                 hawkishness + isolationism + nationalism +
                                 (1 | treat.text.het),
                               family = bernoulli(link = "logit"),
                               backend = "cmdstanr",
                               chains = 4, cores = 4,
                               data = twdata)
summary(interv.model.inter.full)


slopes.inter.full <- slopes(interv.model.inter.full, 
                            variables = c("casualties_9",
                                          "casualties_50",
                                          "casualties_250"),
                            slope = "eydx",
                            newdata = het.datagrid(interv.model.inter.full))
                             

slopes.inter.all <- bind_rows("All Respondents" = slopes.inter.full,
                              "All Manipulation\nChecks Passed" = slopes.inter.pass,
                              .id = "sample") %>%
  mutate(
    nuclear = case_when(
      nuclear == 1 ~ "Nuclear Adversary",
      nuclear == 0 ~ "Non-Nuclear"
    ),
    ally = case_when(
      ally == 1 ~ "Allied Partner",
      ally == 0 ~ "Non-Allied Partner"
    ),
    democracy = case_when(
      democracy == 1 ~ "Democratic Partner",
      democracy == 0 ~ "Non-Democratic Partner"
    ),
    casualties = case_when(
      term == "casualties_9" ~ "9 Casualties",
      term == "casualties_50" ~ "50 Casualties",
      term == "casualties_250" ~ "250 Casualties"
    )
  ) 

ggplot(slopes.inter.all, aes(x = estimate, 
                             y = interaction(ally, democracy, sep = "\n"),
                             color = nuclear)) +
  facet_grid(sample ~ fct_rev(casualties)) +
  geom_vline(xintercept = 0) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = .8),
                  size = .75, linewidth = 1.5) +
  scale_color_grey(start = .6, end = 0) +
  scale_x_continuous(labels = scales::percent) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    y = "",
    color = "",
    x = "Estimated Effect and 95% Confidence Interval",
    title = "Heterogeneous Effects of Casualties",
    subtitle = "Partner Democracy, Alliance, and Adversary Nuclear Capabilities"
  )
ggsave("figures/het-effects.png", height = 6, width = 8)



# quick model with any casualties dummy, fit into explicit het effects framework
formula_pred <- bf(
  use.force ~ lambda*casualties_9 + eta*casualties_50 + tau*casualties_250 + controls,
  
  lambda ~ (nuclear + democracy + ally) + (1 |g| treat.text.het),
  eta ~ (nuclear + democracy + ally) + (1 |g| treat.text.het),
  tau ~ (nuclear + democracy + ally) + (1 |g| treat.text.het),

  # modifiers enter the outcome equation too; otherwise lambda absorbs
  # differences in baseline support across groups.
  # this is alpha_g in the outcome equation: with the modifiers and the group
  # intercept omitted here, baseline support is forced to one constant across
  # all eight cells (0.30 to 0.75 in the raw data), and lambda/eta/tau end up
  # estimating each cell's level relative to the pooled baseline rather than
  # the within-cell effect of casualties. that reverses the ranking of the
  # groups: the effect estimate then correlates about +0.9 with baseline
  # support instead of about -0.8.
  controls ~ nuclear + democracy + ally + (1 |g| treat.text.het) +
    female +
    educ + age + republican + democrat +
    hawkishness + isolationism + nationalism,

  nl = TRUE
)

# every nlpar needs its own prior. previously only lambda and controls were named,
# so eta and tau got improper flat priors on their population-level effects (no
# lprior statement at all in the generated Stan), while lambda was regularized
# N(0, 1) -- the three were not comparable. the sd() terms did fall back to brms'
# default half-student_t(3, 0, 2.5); exponential(2) is set here deliberately,
# since 8 groups carry little information about the variance
pred_prior <- c(
  prior(normal(0, 1), nlpar = "lambda"),
  prior(normal(0, 1), nlpar = "eta"),
  prior(normal(0, 1), nlpar = "tau"),
  prior(normal(0, .5), nlpar = "controls"),
  prior(exponential(2), class = "sd", nlpar = "lambda"),
  prior(exponential(2), class = "sd", nlpar = "eta"),
  prior(exponential(2), class = "sd", nlpar = "tau"),
  prior(exponential(2), class = "sd", nlpar = "controls")
)

# note: this fits the full sample while interv.model.inter above uses
# twdata.pass, so the two are not directly comparable until one is switched
tw_het_nl <- brm(formula_pred,
                 data = twdata,
                 prior = pred_prior,
                 family = bernoulli(link = "logit"),
                 cores = 4,
                 backend = "cmdstanr",
                 refresh = 500,
                 # 11 divergent transitions at the default; only 8 groups inform
                 # each sd, so the funnel needs smaller steps
                 control = list(adapt_delta = 0.99)
)
summary(tw_het_nl)


# just take the lambdas 
# lambda estimates
lambda_est <- posterior_epred(tw_het_nl, nlpar = "lambda")
lambda_quant <- apply(lambda_est, 2,
                      function(x) quantile(x, probs = c(.1, .5, .9)))
rownames(lambda_quant) <- c("treat_10", "treat_med", "treat_90")

lambda_res <- bind_cols(tw_het_nl$data,
                        t(lambda_quant)) %>%
  distinct(
    treat.text.het, treat_med, treat_10, treat_90,
    nuclear, democracy, ally
  )
length(unique(lambda_res$treat_med))

ggplot(lambda_res, aes(x = treat_med, y = treat.text.het)) +
  geom_pointrange(aes(xmin = treat_10, xmax = treat_90)) +
  labs(x = "Estimated Treatment Effect",
       y = "")


# nuclear only
# negative w/ nukes if full sample
interv.model.nuke <- glm(use.force ~ 
                                 casualties_9*(nuclear) +
                                 casualties_50*(nuclear) +
                                 casualties_250*(nuclear) +
                                 nuclear + ally + democracy + female +
                                 educ + age + republican + democrat +
                                 hawkishness + isolationism + nationalism,
                               family = binomial(link = "logit"),
                               data = twdata)
summary(interv.model.nuke)



### FP disposition interactions

# raw data 
dispo.data <- twdata.pass %>%
  select(hawkishness, isolationism,
         nationalism, casualties,
         use.force) %>% 
  pivot_longer( # longer to bring in the groups
    cols = c("hawkishness", "isolationism",
             "nationalism"),
    names_to = "dispo",
    values_to = "value"
  ) %>%
  group_by(casualties, dispo, value) %>% 
  summarize(
    number = sum(use.force, na.rm = TRUE),
    n = n(),
    perc = number / n,
    .groups = "keep"
  )

# plot dispositions by treatment
dispo.labs <- c("hawkishness" = "Mil. Assertiveness",
                         "isolationism" = "Isolationism",
                         "nationalism" = "Natl. Chauvinism")
ggplot(dispo.data, aes(y =  n, x = value,
                           fill = factor(casualties))) +
  facet_wrap(~ dispo, 
             labeller = labeller(dispo = dispo.labs)) +
  geom_bar(stat = "identity",
         position = position_dodge(width = .5)) 

# plot force use 
ggplot(dispo.data, aes(x = value, y = perc,
                       color = factor(casualties))) +
  facet_wrap(~ dispo, 
             labeller = labeller(dispo = dispo.labs)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = wes_palette("Zissou1", n = 4)) +
  labs(
    color = "Casualties",
    y = "Percentage Supporting Force",
    x = "Disposition"
  )


# grouped summaries for plotting 
dispo.data.sum <- twdata.pass %>% 
  select(hawkishness, isolationism,
         nationalism, casualties_dummy,
         use.force) %>% 
  pivot_longer( # longer to bring in the groups
    cols = c("hawkishness", "isolationism",
             "nationalism"),
    names_to = "dispo",
    values_to = "value"
  ) %>% 
     group_by(dispo, value,
              casualties_dummy) %>%
               summarize(
                 n = n(),
                 force.perc = (sum(use.force, na.rm = T)) / n,
                 .groups = "keep"
                 )

# plot- percentages by disposition 
ggplot(dispo.data.sum, aes(x = value, y = force.perc,
                           color = factor(casualties_dummy),
                           group = factor(casualties_dummy))) +
  facet_wrap(~ dispo,
             labeller = labeller(dispo = dispo.labs)) +
  geom_point(size = 3) +
  geom_line(linewidth = 2) +
  scale_color_grey(start = .6, end = .4,
                   labels = c("0" = "No",
                              "1" = "Yes")) +
  labs(
    color = "Casualties",
    y = "Percentage Supporting Force",
    x = "Disposition Score (Rescaled)"
  )
ggsave("appendix/dispo-interactions-perc.png", 
       height = 6, width = 8)

# differences at each disposition level- less informative
dispo.data.diff <- dispo.data.sum %>% 
                    select(dispo, value, casualties_dummy, force.perc) %>%
                    pivot_wider(names_from = c("casualties_dummy"),
                                values_from = "force.perc") %>%
                    mutate(
                      diff.supp = `1` - `0`
                    )
ggplot(dispo.data.diff, aes(x = value, y = diff.supp)) +
  facet_wrap(~ dispo,
             labeller = labeller(dispo = dispo.labs)) +
    geom_bar(stat = "identity")


# interact w/ foreign policy dispositions- liner prob model 
interv.ols.disp <- lm(use.force ~ casualties_dummy +
                           casualties_dummy:hawkishness +
                           casualties_dummy:isolationism +
                           casualties_dummy:nationalism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         data = twdata.pass)
summary(interv.ols.disp)


# interact w/ foreign policy dispositions- ols, mil. assertiveness
interv.ols.hawk <- lm(use.force ~ casualties_dummy +
                           casualties_dummy:hawkishness +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         data = twdata.pass)
summary(interv.ols.hawk)

# interact w/ foreign policy dispositions- logit, isolationism
interv.ols.isol <- lm(use.force ~ casualties_dummy +
                           casualties_dummy:isolationism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         data = twdata.pass)
summary(interv.ols.isol)

# interact w/ foreign policy dispositions- logit, natl chauvinism
interv.ols.natl <- lm(use.force ~ casualties_dummy +
                           casualties_dummy:nationalism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         data = twdata.pass)
summary(interv.ols.natl)


# table for the appendix
ols.tab <- modelsummary(list(interv.ols.hawk, 
                               interv.ols.isol,
                               interv.ols.natl),
                          output = 'latex',
                          statistic = "({conf.low}, {conf.high})",
                          coef_map = c(
                            'casualties_dummy' = 'Casualties',
                            'casualties_9' = '9 Casualties',
                            'casualties_50' = '50 Casualties',
                            'casualties_250' = '250 Casualties',
                            'hawkishness' = 'Mil. Assert',
                            'isolationism' = 'Isolationism',
                            'nationalism' = 'National Chauv.',
                            'casualties_dummy:hawkishness' =
                              'Casualties x Mil. Assert',
                            'casualties_dummy:isolationism' =
                              'Casualties x Isolationism',
                            'casualties_dummy:nationalism' =
                              'Casualties x National Chauv.',
                            'nuclear' = 'Nuclear Adversary',
                            'ally' = 'Alliance',
                            'democracy' = 'Democracy',
                            'educ' = 'Education',
                            'age' = 'Age',
                            'female' = 'Female',
                            'republican' = 'Republican',
                            'democrat' = 'Democrat'
                          ),
                          #output = 
                          gof_omit = 'AIC|Log.Lik|F|RMSE|R2',
                          latex_options = "scale_down",
                          notes = paste0("95\\\\% Confidence Intervals in Parentheses"),
                          caption = "Coefficient estimates from linear probability models of support for military intervention in the attentive sample.\\label{tab:ols-inter-coefs}") %>%
  kableExtra::kable_styling(font_size = 12, full_width = FALSE)
kableExtra::save_kable(ols.tab, "appendix/ols-inter-coefs.tex")

# plot ME 
# plot each: first mil assert
mil.assert.ols <- margins.dispo.func(interv.ols.hawk, "hawkishness", "Militant Assertiveness")
# then isolationism
isol.ols <- margins.dispo.func(interv.ols.isol, "isolationism", "Isolationism")
# last, national chauvinism
natl.ols <- margins.dispo.func(interv.ols.natl, "nationalism", "National Chauvinism")


# plot ME: all together
inter.est.ols <- bind_rows("Militant Assertiveness" = mil.assert.ols[[1]],
                       "Isolationism" = isol.ols[[1]],
                       "National Chauvinism" = natl.ols[[1]],
                       .id = "Disposition")


ggplot(inter.est.ols, aes(x = xvals,
                      y = yvals)) +
  facet_wrap(~ Disposition) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, 
                  ymax = upper),
              alpha = .5) +
  labs(x = "Modifying Variable Range (Rescaled)",
       y = "Estimated Marginal Effect of Casualties")
ggsave("appendix/disp-interactions-ols.png",
       height = 6, width = 8)



### Logit for the interactions

# interact w/ foreign policy dispositions- logit
interv.model.disp <- glm(use.force ~ casualties_dummy +
                           casualties_dummy:hawkishness +
                           casualties_dummy:isolationism +
                           casualties_dummy:nationalism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         family = binomial(link = "logit"),
                         data = twdata.pass)
summary(interv.model.disp)


# interact w/ foreign policy dispositions- logit, mil. assertiveness
interv.model.hawk <- glm(use.force ~ casualties_dummy +
                           casualties_dummy:hawkishness +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         family = binomial(link = "logit"),
                         data = twdata.pass)
summary(interv.model.hawk)

# interact w/ foreign policy dispositions- logit, isolationism
interv.model.isol <- glm(use.force ~ casualties_dummy +
                           casualties_dummy:isolationism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         family = binomial(link = "logit"),
                         data = twdata.pass)
summary(interv.model.isol)

# interact w/ foreign policy dispositions- logit, natl chauvinism
interv.model.natl <- glm(use.force ~ casualties_dummy +
                           casualties_dummy:nationalism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         family = binomial(link = "logit"),
                         data = twdata.pass)
summary(interv.model.natl)


# table for the appendix
logit.tab <- modelsummary(list(interv.model, 
                  interv.model.control,
                  interv.model.num,
                  interv.model.hawk, 
                  interv.model.isol,
                  interv.model.natl),
                  output = 'latex',
             statistic = "({conf.low}, {conf.high})",
             coef_map = c(
               'casualties_dummy' = 'Casualties',
               'casualties_9' = '9 Casualties',
               'casualties_50' = '50 Casualties',
               'casualties_250' = '250 Casualties',
               'hawkishness' = 'Mil. Assert',
               'isolationism' = 'Isolationism',
               'nationalism' = 'National Chauv.',
               'casualties_dummy:hawkishness' =
                 'Casualties x Mil. Assert',
               'casualties_dummy:isolationism' =
                 'Casualties x Isolationism',
               'casualties_dummy:nationalism' =
                 'Casualties x National Chauv.',
               'nuclear' = 'Nuclear Adversary',
               'ally' = 'Alliance',
               'democracy' = 'Democracy',
               'educ' = 'Education',
               'age' = 'Age',
               'female' = 'Female',
               'republican' = 'Republican',
               'democrat' = 'Democrat'
               ),
             exponentiate = TRUE,
             #output = 
             gof_omit = 'AIC|Log.Lik|F|RMSE',
             latex_options = "scale_down",
             notes = paste0("95\\\\% Confidence Intervals in Parentheses"),
             caption = "Exponentiated coefficient estimates from logistic regression models of support for military intervention in the attentive sample.\\label{tab:logit-coefs}") %>%
             kableExtra::kable_styling(font_size = 9, full_width = FALSE)
kableExtra::save_kable(logit.tab, "appendix/logit-coefs.tex")





# plot each: first mil assert
mil.assert.marg <- margins.dispo.func(interv.model.hawk, "hawkishness", "Militant Assertiveness")
# then isolationism
isol.marg <- margins.dispo.func(interv.model.isol, "isolationism", "Isolationism")
# last, national chauvinism
natl.marg <- margins.dispo.func(interv.model.natl, "nationalism", "National Chauvinism")


# combine all three:
grid.arrange(mil.assert.marg[[2]], natl.marg[[2]], isol.marg[[2]], nrow = 1)
margins.disp <- arrangeGrob(mil.assert.marg[[2]], isol.marg[[2]], natl.marg[[2]],
                            nrow = 1)

# combine data 
inter.est <- bind_rows("Militant Assertiveness" = mil.assert.marg[[1]],
                       "Isolationism" = isol.marg[[1]],
                       "National Chauvinism" = natl.marg[[1]],
                       .id = "Disposition")


ggplot(inter.est, aes(x = xvals,
                       y = yvals)) +
  facet_wrap(~ Disposition) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, 
                      ymax = upper),
              alpha = .5) +
  labs(x = "Modifying Variable Range (Rescaled)",
    y = "Estimated Marginal Effect of Casualties")
ggsave("appendix/disp-interactions.png",
       height = 6, width = 8)


# combine all three: dist of modifying vars
grid.arrange(mil.assert.marg[[3]], natl.marg[[3]], isol.marg[[3]], nrow = 1,
             top = "Distribution of Modifying Variables")
margins.dist <- arrangeGrob(mil.assert.marg[[3]], isol.marg[[3]], natl.marg[[3]], nrow = 1,
                            top = "Distribution of Modifying Variables")
# ggsave("appendix/disp-interactions.png", margins.dist,
#        height = 6, width = 8)



# interflex militant assertiveness
bin.mil.assert <- interflex(estimator = "binning",
          Y = "use.force", 
          D = "casualties_dummy", X = "hawkishness", 
          data = twdata.pass,
          method = "logit",
          weights = NULL, 
          Ylabel = "Support for Force", 
          Dlabel = "Casualties", 
          Xlabel = "Militant Assertiveness", 
          bin.labs = FALSE,
          na.rm = TRUE)
plot(bin.mil.assert)

kernel.mil.assert <- interflex(estimator = "kernel",
                            Y = "use.force", 
                            D = "casualties_dummy", X = "hawkishness", 
                            data = twdata.pass,
                            method = "logit",
                            weights = NULL, 
                            Ylabel = "Support for Force", 
                            Dlabel = "Casualties", 
                            Xlabel = "Militant Assertiveness", 
                            main = "Raw Plot",
                            na.rm = TRUE)
plot(kernel.mil.assert) 


# interflex isolationism
bin.isolation <- interflex(estimator = "binning",
                            Y = "use.force", 
                            D = "casualties_dummy", X = "isolationism", 
                            data = twdata.pass,
                            nbins = 2,
                            method = "logit",
                            weights = NULL, 
                            Ylabel = "Support for Force", 
                            Dlabel = "Casualties", 
                            Xlabel = "Isolationism", 
                            bin.labs = FALSE,
                            na.rm = TRUE)
plot(bin.isolation)

kernel.isolation <- interflex(estimator = "kernel",
                               Y = "use.force", 
                               D = "casualties_dummy", X = "isolationism", 
                               data = twdata.pass,
                               method = "logit",
                               weights = NULL, 
                               Ylabel = "Support for Force", 
                               Dlabel = "Casualties", 
                               Xlabel = "Isolationism", 
                               main = "Raw Plot",
                               na.rm = TRUE)
plot(kernel.isolation) 

# interflex national chauvinism
bin.natl.chauv <- interflex(estimator = "binning",
                           Y = "use.force", 
                           D = "casualties_dummy", X = "nationalism", 
                           data = twdata.pass,
                           nbins = 2,
                           method = "logit",
                           weights = NULL, 
                           Ylabel = "Support for Force", 
                           Dlabel = "Casualties", 
                           Xlabel = "National Chauvinism", 
                           bin.labs = FALSE,
                           na.rm = TRUE)
plot(bin.natl.chauv)

kernel.natl.chauv <- interflex(estimator = "kernel",
                              Y = "use.force", 
                              D = "casualties_dummy", X = "nationalism", 
                              data = twdata.pass,
                              method = "logit",
                              weights = NULL, 
                              Ylabel = "Support for Force", 
                              Dlabel = "Casualties", 
                              Xlabel = "National Chauvinism", 
                              main = "Raw Plot",
                              na.rm = TRUE)
plot(kernel.natl.chauv) 



# combine plots: binning
grid.arrange(plot(bin.mil.assert),
             plot(bin.isolation),
             plot(bin.natl.chauv),
             nrow = 1)

bin.inter.res <- arrangeGrob(plot(bin.mil.assert),
                              plot(bin.isolation),
                              plot(bin.natl.chauv),
                              nrow = 1)
ggsave("appendix/bin-inter-res.png", bin.inter.res,
       height = 8, width = 10)

# combine plots: kernel
grid.arrange(plot(kernel.mil.assert),
             plot(kernel.isolation),
             plot(kernel.natl.chauv),
             nrow = 1)

kernel.inter.res <- arrangeGrob(plot(kernel.mil.assert),
                             plot(kernel.isolation),
                             plot(kernel.natl.chauv),
                             nrow = 1)
ggsave("appendix/kernel-inter-res.png", kernel.inter.res,
       height = 8, width = 10)



# interact w/ foreign policy dispositions- ols
interv.disp.lm <- lm(use.force ~ casualties_dummy +
                           casualties_dummy:hawkishness +
                           casualties_dummy:isolationism +
                           casualties_dummy:nationalism +
                           nuclear + ally + democracy + female +
                           educ + age + republican + democrat +
                           hawkishness + isolationism + nationalism,
                         data = twdata.pass)
summary(interv.disp.lm)


### gender interaction
### treatment interactions

# estimate model: stronger impact on men
interv.model.gender <- glm(use.force ~ casualties_dummy*female +
                            nuclear + ally + democracy + 
                            educ + age + republican + democrat +
                            hawkishness + isolationism + nationalism,
                          family = binomial(link = "logit"),
                          data = twdata.pass)
summary(interv.model.gender)



### July 2023 responses
t.test(force ~ casualties, data = july.23.data)
table(july.23.data$response.txt)
ggplot(july.23.data, aes(x = response.txt)) +
  geom_bar()


july.23.data.pct <- july.23.data %>%
  group_by(scenario, response.txt) %>%
  summarise(country = first(country),
            casualties = median(casualties),
            count = n(), .groups = 'drop') %>%
  group_by(scenario) %>%
  mutate(percentage = count / sum(count) * 100)

ggplot(july.23.data.pct, aes(x = response.txt,
                             y = percentage,
                             group = factor(casualties),
                             fill = factor(casualties))) +
  facet_wrap(~ country, ncol = 3) +
  ylim(0, 45) +
  coord_flip() +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(percentage, 0)),
            position = position_dodge2(width = .9),
            hjust = -0.2) +
  scale_fill_grey(start = .5, end = .2, 
                  name = "Casualties",
                  labels = c(`0` = "None",
                             `1` = "50")) +
  labs(x = "",
       y = "Percentage of Respondents",
       title = "Policy Response Preferences in Three Scenarios") +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave("figures/july-23-res.png", height = 6, width = 8)


# randomization check
# check balance
balance.check.23 <- glm(casualties ~ male + college_grad + age + 
                            conservative +
                            mil_assert + isolationism + natl_chauv,
                          data = july.23.data)
summary(balance.check.23)


modelplot(balance.check.23,
          coef_rename = c("(Intercept)" = "Intercept",
                          "male" = "Male",
                          "college_grad" = "College Graduate",
                          "age" = "Age",
                          "conservative" = "Republican",
                          "mil_assert" = "Mil. Assertiveness",
                          "isolationism" = "Isolationism",
                          "natl_chauv" = "Natl. Chauvinism"),
          coef_omit = "Intercept") +
  labs(title = "Randomization Check: July 2023 Experiment")
ggsave("appendix/randomization-check-23.png", height = 6, width = 8)
