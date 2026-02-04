# Joshua Alley
# compare model specs 


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
              
              het_group = paste(dem, rep, high_newsint, white, male, hawk, intl,
                                  sep = "_"),

              het_group_small = paste(intl, white, male,
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


# 1- brms reg spec
formula_pred <- bf(
  force ~ lambda*alliance + controls,
  
  lambda ~ (dem + rep + high_newsint + white + male + hawk + intl) + (1|het_group), 
  
  controls ~ regime + stakes + costs + 
              africa + europe + asia +
              age + ed4,
  
  nl = TRUE
)
formula_pred

pred_prior <- c(
  prior(normal(0, .5), nlpar = "lambda"),
  prior(normal(0, .5), nlpar = "controls")
  )


het_reg <- brm(formula_pred,
                    data = tw_rep,
                    prior = pred_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 500
)
summary(het_reg)
sum(subset(nuts_params(het_reg), Parameter == "divergent__")$Value)

# varying intercepts
formula_vi <- bf(
  force ~ (1 + alliance || het_group) +
    dem + rep + high_newsint + white + male + hawk + intl + 
    regime + stakes + costs + 
              africa + europe + asia +
              age + ed4)
formula_vi

vi_prior <- c(
  prior(normal(0, 1), class = "sigma"),
  prior(normal(0, .5), class = "b")
  )


het_vi <- brm(formula_vi,
                    data = tw_rep,
                    prior = vi_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 500
)
summary(het_vi)
sum(subset(nuts_params(het_vi), Parameter == "divergent__")$Value)


# interaction
formula_inter <- bf(
  force ~ (1 + alliance || het_group) +
    alliance*(dem + rep + high_newsint + white + male + hawk + intl) +
    regime + stakes + costs + 
              africa + europe + asia +
              age + ed4)
formula_inter

inter_prior <- c(
  prior(normal(0, 1), class = "sigma"),
  prior(normal(0, .5), class = "b")
  )


het_inter <- brm(formula_inter,
                    data = tw_rep,
                    prior = inter_prior,
                    family = gaussian(),
                    cores = 4,
                    backend = "cmdstanr",
                    refresh = 500
)
summary(het_inter)
table(tw_rep$het_group_small, tw_rep$treat_group)
sum(subset(nuts_params(het_inter), Parameter == "divergent__")$Value)


# ols with interactions 
het_inter_ols <- lm(force ~ 
     alliance*(dem + rep + high_newsint + white + male + hawk + intl) +
    regime + stakes + costs + 
              africa + europe + asia +
              age + ed4,
          data = tw_rep)
summary(het_inter_ols)


# pull models and compare 
coef_reg <- as.data.frame(fixef(het_reg, summary = TRUE))
coef_reg$variable <- rownames(coef_reg)


# Look at slopes
lambda_est <- posterior_epred(het_reg, nlpar = "lambda")
lambda_quant <- apply(lambda_est, 2, 
                      function(x) quantile(x, probs = c(.1, .5, .9)))
rownames(lambda_quant) <- c("treat_10", "treat_med", "treat_90")

reg_est <- bind_cols(tw_rep, t(lambda_quant)) %>%
            #distinct(het_group_small, intl, white, male, treat_med) %>%
            distinct(het_group, treat_med, dem, rep, high_newsint, white, male, hawk, intl) %>%
            rename(
              estimate_reg = treat_med
            )


# create newdata
pred_data <- tw_rep %>%
              #distinct(het_group_small, intl, white, male) %>%
  distinct(het_group, dem, rep, high_newsint, white, male, hawk, intl) %>%
              mutate(
                regime = 1,
                stakes = 1, 
                costs = 0,
                africa = 0,
                europe = 1,
                asia = 1,
                ed4 = 2,
                age = 4.5
              )



# slopes by the VI model
vi_est <- slopes(het_vi, variables = "alliance",
                  newdata = pred_data)
vi_est <- vi_est %>%
            rename(
              estimate_vi = estimate
            ) %>%
            select(
              het_group, estimate_vi
            )
table(vi_est$estimate_vi)
sum(length(unique(vi_est$estimate_vi)))

# slopes by the interaction model
inter_est <- slopes(het_inter, variables = "alliance",
                newdata = pred_data) %>%
            rename(
              estimate_inter = estimate
            ) %>%
            select(
              het_group, estimate_inter
            )
table(inter_est$estimate_inter)
sum(length(unique(inter_est$estimate_inter)))


# ols interaction estimates 
ols_est <- slopes(het_inter_ols, variables = "alliance",
                  newdata = pred_data) %>%
            rename(
              estimate_ols = estimate
            ) %>%
            select(
              het_group, estimate_ols
            )

# combine these three things
est_all <- left_join(reg_est, inter_est) %>%
             left_join(vi_est) %>%
             left_join(ols_est)
summary(est_all$estimate_reg) 
sd(est_all$estimate_reg)
summary(est_all$estimate_vi)
sd(est_all$estimate_vi)
summary(est_all$estimate_inter)
sd(est_all$estimate_inter)
summary(est_all$estimate_ols)
sd(est_all$estimate_ols)


est_all_long <- est_all %>%
             pivot_longer(
              cols = c(estimate_reg, estimate_vi,
                 estimate_inter, estimate_ols),
              names_to = "model",
              values_to = "estimate"
             )


ggplot(est_all_long, aes(y = het_group, x = estimate)) +
  geom_point(aes(color = model)) +
  theme(legend.position = "bottom")

ggplot(est_all_long, aes(y = het_group, x = estimate)) +
  facet_wrap(~ model) +
  geom_point() +
  theme(legend.position = "bottom")


ggplot(filter(est_all_long, model == "estimate_inter" | model == "estimate_ols"),
 aes(y = het_group, x = estimate)) +
  facet_wrap(~ model, ncol = 1) +
  geom_point(aes(color = model)) +
  theme(legend.position = "bottom")
