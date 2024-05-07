# Joshua Alley
# analysis with simulated data- recover true values


# vary sample size
n.sim <- c(1000, 2500, 5000)
# vary SD of coefficients
sd.coef <- c(.05, .25, .75)



simulation.inter.comp <- function(data){
  
sample.size = data$sample.size
sd.impact = data$sd.impact

# predictors
predictors.sim <- data.frame( 
# group1 = factor(
#             sample(size = sample.size, 
#                x = c(0, 1, 2, 3, 4),
#                replace = TRUE,
#                prob = c(.2, .2, .2, .2, .2))
#             ),
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

# coefficients will want to vary SD 
beta <- rnorm(ncol(predictors.sim.inter), mean = 0, sd = sd.impact)

# group indices
# predictors.sim.inter$group.ind <- predictors.sim.inter %>%
#                             group_by(group1, group2, group3) %>%
#                             group_indices()
# table(predictors.sim.inter$group.ind)
# 
# n.group <- max(predictors.sim.inter$group.ind)
# print(n.group)

# outcome- normally distributed
mu.y <- as.numeric(as.matrix(predictors.sim.inter) %*% beta)

y <- rnorm(sample.size, mean = mu.y, sd = .25)


# simulated data
sim.data <- bind_cols(y = y, mu.y = mu.y, predictors.sim)


# fit model
ols.inter <- lm(y ~ treat*group1*group2*group3*group4,
                data = sim.data)
summary(ols.inter)

# hypothetical data
grid.sim <- sim.data %>%
   ungroup() %>% 
   select(group1, group2, group3, group4) %>%
   distinct() 

#top-end betas
 beta.group <- rowSums(as.matrix(grid.sim) * beta)
 
#  list(grid.sim, predictors.sim.inter, beta, beta.group)
# } 

# get slopes
slopes.ols <- slopes(model = ols.inter,
                           variables = "treat",
                           newdata = grid.sim)  %>%
               mutate(
                true = beta.group,
                bias = estimate - true,
                in.interval = ifelse(true > conf.low & true < conf.high, 1, 0)
               )
table(slopes.ols$bias)
sum(slopes.ols$in.interval)

# look at constituent terms
bias.ols.coef <- coef(ols.inter) - beta
print(mean(bias.ols.coef))

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
  mutate(
    true = beta.group,
    bias = estimate - true,
    in.interval = ifelse(true > conf.low & true < conf.high, 1, 0)
  )

# output everything in a list
sim.input <- list("outcome" = y, "data" = sim.data,
                  "beta" = beta, "group.effects" = beta.group)
ols.res <- list("ols.mod" = ols.inter, "slopes.ols" = slopes.ols)
vs.res <- list("vs.mod" = vs.mod, "slopes.vs" = slopes.vs,
               "time.vs" = time.taken)

output  <- list(sim.input, ols.res, vs.res)
output 

# end comparison function
}

combos <- expand.grid(n.sim, sd.coef)
colnames(combos) <- c("sample.size", "sd.impact")
combos$combo <- rownames(combos)
combos$pair <- paste0(combos$sample.size, "_", combos$sd.impact)  

combos.list <- split(combos, f = combos$combo)

simulation.all <- lapply(combos.list,
                         simulation.inter.comp)

names(simulation.all) <- combos$pair
# save results
saveRDS(simulation.all, file = "data/simulation-res.rds")


# calculate RMSE and group coef bias 
sim.res <- vector(mode = "list", length = length(simulation.all))
sim.slopes.res <- vector(mode = "list", length = length(simulation.all))

# for loop given nested lists
for(i in 1:length(simulation.all)){
  
  est <- unlist(simulation.all[[i]], recursive = FALSE)
  
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
  bias.ols <- sqrt(mean((est$group.effects - ols.est)^2))
  mean(bias.ols)
  
  # RMSE VS: 
  # outcome
  resid.vs <- est$outcome - pred.vs.med
  rmse.vs <- sqrt(mean((est$outcome - pred.vs.med)^2))
  
  # group coefs
  bias.vs <- sqrt(mean((est$group.effects - est$slopes.vs$estimate)^2))
  mean(bias.vs)
  
  # results
  model <- c("OLS", "VS")
  rmse <- c(rmse.ols, rmse.vs)
  bias <- c(bias.ols, bias.vs)

  res <- data.frame(model = model,
                    rmse = rmse,
                    bias = bias,
                    sample.size = combos$sample.size[i],
                    sd.impact = str_remove(as.character(
                      combos$sd.impact[i]),
                      "^0+"))
  
  sim.res[[i]] <- res
  sim.slopes.res[[i]] <- bind_rows("OLS" = est$slopes.ols,
                             "VS" = est$slopes.vs,
                             .id = "model") %>%
                    mutate(
                      sample.size = combos$sample.size[i],
                      sd.impact = str_remove(as.character(
                      combos$sd.impact[i]),
                      "^0+"))
}

sim.res.data <- bind_rows(sim.res,
                          .id = "sim") %>%
                mutate(
                  sd.impact = paste0("SD Coef=", sd.impact),
                  scen = paste0("N=", sample.size, ",\n",
                                sd.impact)
                )

ggplot(sim.res.data, aes(x = factor(sample.size), y = rmse,
                         color = model)) +
  facet_wrap(~ sd.impact, scales = "free_x") +
  geom_point(size = 3) +
  labs(title = "RMSE of Varying Slopes and OLS",
       x = "Sample Size",
       y = "RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")
ggsave("figures/sim-rmse-out.png", height = 6, width = 8)

ggplot(sim.res.data, aes(x = factor(sample.size), y = bias,
                         color = model)) +
  facet_wrap(~ sd.impact) +
  geom_point(size = 3) +
  labs(title = "RMSE Treatment Estimate: Varying Slopes and OLS",
       x = "Sample Size",
       y = "Estimate RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")
ggsave("figures/sim-rmse-coef.png", height = 6, width = 8)



# futher treatment data- slopes 
sim.slopes.data <- bind_rows(sim.slopes.res,
                          .id = "sim") %>%
  mutate(
    sd.impact = paste0("SD Coef=", sd.impact),
    scen = paste0("N=", sample.size, ",\n",
                  sd.impact)
  )

ggplot(sim.slopes.data, aes(x = rowid, y = bias,
                            color = model)) +
  geom_hline(yintercept = 0) +
  facet_wrap(sd.impact ~ sample.size,
             scales = "free_y") +
  geom_point() 


ggplot(sim.slopes.data, aes(x = factor(in.interval),
                            fill = model)) +
  facet_grid(sd.impact ~ sample.size) +
  geom_bar(position = position_dodge(width = 1)) 

