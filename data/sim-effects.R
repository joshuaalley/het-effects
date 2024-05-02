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
group1 = rmultinom(n = sample.size, size = 4, 
                    prob = c(.25, .25, .25, .25))[1 ,],
group2 = rbinom(sample.size, size = 1, prob = .3),
group3 = rbinom(sample.size, size = 1, prob = .7),

treat = rbinom(sample.size, size = 1, prob = .5)
)


# Create interaction terms
predictors.sim.inter <- as.data.frame(model.matrix( ~ treat*group1*group2*group3,
                                      data = predictors.sim))

# coefficients will want to vary SD 
beta <- rnorm(ncol(predictors.sim.inter), mean = .25, sd = sd.impact)

# group indices
predictors.sim.inter$group.ind <- predictors.sim.inter %>%
                            group_by(group1, group2, group3) %>%
                            group_indices()
table(predictors.sim.inter$group.ind)

n.group <- max(predictors.sim.inter$group.ind)

# outcome- normally distributed
mu.y <- as.numeric(as.matrix(select(predictors.sim.inter,
                                    -group.ind)) %*% beta)

y <- rnorm(sample.size, mean = mu.y, sd = .25)


# simulated data
sim.data <- bind_cols(y = y, mu.y = mu.y, predictors.sim.inter)




# fit model
ols.inter <- lm(y ~ treat*group1*group2*group3,
                data = sim.data)
summary(ols.inter)

# hypothetical data
grid.sim <- predictors.sim.inter %>%
  # select(treat, group1, group2, group3,
  #        group.ind) %>%
  group_by(group.ind) %>%
  mutate(
    n = n()
  ) %>%
  distinct() %>%
  ungroup() %>%
  filter(treat == 1)

# top-end betas  
beta.group <- rowSums(as.matrix(select(grid.sim, -c(n, group.ind))) * beta)

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
mean(bias.ols.coef)


# VS model
sim.het.prior <- c(
  prior(normal(0, .75), class = "Intercept"),
  prior(normal(0, 1), class = "sd"),
  prior(normal(0, 1), class = "sigma")
)

start.time <- Sys.time()
vs.mod <- brm(bf(y ~ 1 +
                         (1 + treat | group1*group2*group3)),
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
table(slopes.vs$bias)
mean(slopes.vs$bias)
sum(slopes.vs$in.interval)

# predictions
pred.vs <- posterior_epred(object = vs.mod) 
pred.vs.med <- apply(pred.vs, 2, median)


# RMSE OLS:
# outcome
sqrt(mean(ols.inter$residuals^2))
# group coefs
sqrt(mean((beta.group - slopes.ols$estimate)^2))
mean(slopes.ols$bias)

# RMSE VS: 
# outcome
sqrt(mean((sim.data$y - pred.vs.med)^2))
# group coefs
sqrt(mean((beta.group - slopes.vs$estimate)^2))
mean(slopes.vs$bias)


# 90% interval coverage
ggplot(slopes.vs, aes(x = true, y = group.ind)) +
  geom_vline(xintercept = .25,
             linetype = "dashed") +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high))


# 90% interval coverage
ggplot(slopes.ols, aes(x = true, y = group.ind)) +
  geom_vline(xintercept = .25,
             linetype = "dashed") +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high))

# output everything in a list
sim.data <- list("outcome" = y, "beta"= beta, "group.effects" = beta.group)
ols.res <- list(ols.inter, slopes.ols)
vs.res <- list(vs.mod, slopes.vs, time.taken)

output  <- list(sim.data, ols.res, vs.res)
output 

# end comparison function
}

results.list <- vector(mode = "list",
                       length = length(n.sim)*length(sd.coef))

combos <- expand.grid(n.sim, sd.coef)
colnames(combos) <- c("sample.size", "sd.impact")
combos$combo <- rownames(combos)
  

combos.list <- split(combos, f = combos$combo)

simulation.all <- lapply(combos.list,
                         simulation.inter.comp)
# save results
saveRDS(simulation.all, file = "data/simulation-res.rds")