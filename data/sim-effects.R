# Joshua Alley
# analysis with simulated data- recover true values


# vary sample size
n.sim <- c(1000, 2000, 3000)
# vary SD of coefficients
sd.coef <- c(.05, .25, .75)


# create coefficient vector outside of estimation function
# 3 vectors- one for each SD
beta.list <- vector(mode = "list", length = 3)

for(i in 1:length(beta.list)){
beta.list[[i]] <- rnorm(32, mean = 0, sd = sd.coef[i])
}

# simulation data
sim.data <- vector(mode = "list", 
                   length = length(beta.list))
for(i in 1:length(beta.list)){
  
sample.size <- n.sim[i]  
# predictors
predictors.sim <- data.frame( 
  group1 = rbinom(sample.size, size = 1, prob = .3),
  group2 = rbinom(sample.size, size = 1, prob = .4),
  group3 = rbinom(sample.size, size = 1, prob = .5),
  group4 = rbinom(sample.size, size = 1, prob = .6),
  
  treat = rbinom(sample.size, size = 1, prob = .5)
)


# Create interaction terms
predictors.sim.inter <- as.data.frame(
  model.matrix( ~ treat*group1*group2*group3*group4,
                data = predictors.sim)) 

# outcome- normally distributed
mu.y <- as.numeric(as.matrix(predictors.sim.inter) %*% beta.list[[i]])

y <- rnorm(sample.size, mean = mu.y, sd = .75)


# simulated data
sim.data[[i]] <- bind_cols(y = y, mu.y = mu.y, predictors.sim)

}


# function for simulation
simulation.inter.comp <- function(list){
  
sim.data <- as.data.frame(list$data)
beta <- unlist(list$beta)

# differences in means by groups
grid.calc <- sim.data %>%
  group_by(treat,
           group1, group2, group3, group4) %>%
  summarize(
    n = n(),
    mu.y = mean(mu.y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    id_cols = c(group1, group2, group3, group4),
    names_from = treat,
    values_from = c(n, mu.y)
  ) %>%
  mutate(# treatment effect for each group
    true = mu.y_1 - mu.y_0 
  )

# hypothetical data grid for slopes
grid.sim <- sim.data %>%
  ungroup() %>% 
  select(group1, group2, group3, group4) %>%
  distinct() 


### OLS regression
# fit model
ols.inter <- lm(y ~ treat*group1*group2*group3*group4,
                data = sim.data)
summary(ols.inter)
 
# get slopes
slopes.ols <- slopes(model = ols.inter,
                           variables = "treat",
                           newdata = grid.sim)  %>%
               left_join(grid.calc) %>%
               mutate(
                bias = estimate - true,
                in.interval = ifelse(true > conf.low & true < conf.high, 1, 0)
               )
table(slopes.ols$bias)
sum(slopes.ols$in.interval)

# look at constituent terms
bias.ols.coef <- coef(ols.inter) - beta


### hierarchical regression
# VS model
sim.het.prior <- c(
  prior(normal(0, .75), class = "Intercept"),
  prior(normal(0, 1), class = "sd"),
  prior(normal(0, 1), class = "sigma")
)

start.time <- Sys.time()
vs.mod <- brm(bf(y ~ 1 +
                         (1 + treat | group1*group2*group3*group4)),
            data = sim.data,
            prior = sim.het.prior,
            family = gaussian(),
            cores = 4,
            control = list(adapt_delta = .99,
               max_treedepth = 20),
            backend = "cmdstanr",
           refresh = 500
)
summary(vs.mod)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# slopes again
slopes.vs <- slopes(model = vs.mod,
                    conf_level = .9,
                     variables = "treat",
                     newdata = grid.sim,
                    allow_new_levels = TRUE) %>%
  left_join(grid.calc) %>%
  mutate(
    bias = estimate - true,
    in.interval = ifelse(true > conf.low & true < conf.high, 1, 0)
  )

# output everything in a list
sim.input <- list("outcome" = sim.data$y, "data" = sim.data,
                  "beta" = beta, "group.effects" = grid.calc$true)
ols.res <- list("ols.mod" = ols.inter, "slopes.ols" = slopes.ols)
vs.res <- list("vs.mod" = vs.mod, "slopes.vs" = slopes.vs,
               "time.vs" = time.taken)

output  <- list(sim.input, ols.res, vs.res)
output 

# end comparison function
}


# set combinations
combos <- expand.grid(sim.data, beta.list)
colnames(combos) <- c("data", "beta")
combos$sample.size <- rep(n.sim, times = 3)
combos$sd.coef <- rep(sd.coef, each = 3)
combos$combo <- rownames(combos)
combos$pair <- paste0(combos$sample.size, "_", combos$sd.coef)  

combos.list <- split(combos, f = combos$combo)

# run it
simulation.all <- lapply(combos.list,
                         simulation.inter.comp)

names(simulation.all) <- combos$pair
# save results
saveRDS(simulation.all, file = "data/simulation-res.rds")
# load 
simulation.all <- readRDS(file = "data/simulation-res.rds")


# return data
combos.ret <- bind_rows(combos.list)

# calculate RMSE and group coef bias 
sim.res <- vector(mode = "list", length = length(simulation.all))
sim.slopes.res <- vector(mode = "list", length = length(simulation.all))

# for loop given nested lists
for(i in 1:length(simulation.all)){
  
  est <- unlist(simulation.all[[i]], recursive = FALSE)
  
  data <- est$data
  data$group.var <- data %>%
                      select(-c(y, mu.y)) %>%
                      group_by(treat, group1, group2, 
                               group3, group4) %>%
                      group_indices()
  
  group.sum <- data %>%
                  group_by(group.var) %>%
                  summarize(
                    n = n()
                  )
  
  distinct(est$slopes.vs, group1, group2, group3,
            .keep_all = TRUE) %>%
      glimpse()

  # predictions
  pred.vs <- posterior_epred(object = est$vs.mod)
  pred.vs.med <- apply(pred.vs, 2, median)

  
  # RMSE OLS:
  # outcome
  rmse.ols <- sqrt(mean(est$ols.mod$residuals^2))
  # group coefs
  ols.est <- est$slopes.ols$estimate
  bias.ols <- est$group.effects - ols.est
  rmse.coef.ols <- sqrt(mean((est$group.effects - ols.est)^2))
  mean(bias.ols)
  
  # RMSE VS:
  # outcome
  resid.vs <- est$outcome - pred.vs.med
  rmse.vs <- sqrt(mean((est$outcome - pred.vs.med)^2))

  # group coefs
  bias.vs <- est$group.effects - est$slopes.vs$estimate
  rmse.coef.vs <- sqrt(mean((est$group.effects - est$slopes.vs$estimate)^2))
  mean(bias.vs)

  # results
  model <- c("OLS", "Hierarchical")
  rmse.out <- c(rmse.ols, rmse.vs)
  bias <- c(bias.ols, bias.vs)
  rmse.coef <- c(rmse.coef.ols, rmse.coef.vs)

  res <- data.frame(model = model,
                    rmse.out = rmse.out,
                    bias = bias,
                    rmse.coef = rmse.coef,
                    sample.size = combos.ret$sample.size[i],
                    sd.coef = combos.ret$sd.coef[i]
                    )

  sim.res[[i]] <- res
  sim.slopes.res[[i]] <- bind_rows("OLS" = est$slopes.ols,
                             "Hierarchical" = est$slopes.vs,
                             .id = "model") %>%
                    mutate(
                      sample.size = combos.ret$sample.size[i],
                      sd.coef = combos.ret$sd.coef[i],
                      n = group.sum$n
                      )
} # finish comparison


# pull together comparison data
sim.res.data <- bind_rows(sim.res,
                          .id = "sim") %>%
                mutate(
                  sd.coef = paste0("Coef SD=", sd.coef),
                  scen = paste0("N=", sample.size, ",\n",
                                sd.coef)
                ) %>%
                select(-bias) %>%
                distinct() %>%
                group_by(scen) %>%
                mutate(
                  change.rmse.coef = rmse.coef - lag(rmse.coef),
                  change.rmse.out = rmse.out - lag(rmse.out)
                )

ggplot(sim.res.data, aes(x = factor(sample.size), y = rmse.out,
                         color = model)) +
  facet_wrap(~ sd.coef, scales = "free_x") +
  geom_point(size = 3) +
  labs(title = "RMSE of Varying Slopes and OLS",
       x = "Sample Size",
       y = "RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")

ggplot(sim.res.data, aes(x = factor(sample.size), 
                         y = change.rmse.out)) +
  facet_wrap(~ sd.coef) +
  geom_point(size = 3) +
  labs(title = "Coefficient RMSE",
       subtitle = "Improvement with Hierarchical Model",
       x = "Sample Size",
       y = "Change in Coefficient RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")

# present this in a different way- change in R
ggplot(sim.res.data, aes(x = factor(sample.size), 
                         y = change.rmse.coef)) +
  facet_wrap(~ sd.coef) +
  geom_point(size = 3) +
  labs(title = "Coefficient RMSE",
  subtitle = "Improvement with Hierarchical Model",
       x = "Sample Size",
       y = "Change in Coefficient RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")
ggsave("figures/sim-rmse-coef-nsd.png", height = 6, width = 8)



# futher treatment data- slopes 
sim.slopes.data <- bind_rows(sim.slopes.res,
                          .id = "sim") %>%
  mutate(
    sd.impact = paste0("SD Coef=", sd.coef),
    scen = paste0("N=", sample.size, ",\n",
                  sd.impact)
  )

ggplot(sim.slopes.data, aes(x = rowid, y = bias,
                            color = model)) +
  geom_hline(yintercept = 0) +
  facet_grid(sd.impact ~ sample.size,
             scales = "free_y") +
  geom_point() 

ggplot(sim.slopes.data, aes(x = bias,
                            fill = model)) +
  geom_hline(yintercept = 0) +
  facet_grid(sd.impact ~ sample.size,
             scales = "free_y") +
  geom_histogram(#position = position_dodge(width = .5),
                 bins = 16) 



ggplot(sim.slopes.data, aes(x = factor(in.interval),
                            fill = model)) +
  facet_grid(sd.impact ~ sample.size) +
  geom_bar(position = position_dodge(width = 1)) 




