# Joshua Alley
# vary number of groups

### Simulation that varies number of groups, sets group effects
# express as treat | group1 + group2 ... 

# key parameter here
num.group <- c(2, 3, 4, 5, 6)
# vary SD of coefficients
sd.coef <- c(.05, .25, .75)


group.data.sim <- function(list){
  
  num.group <- list$num.group
  sd.coef <- list$sd.coef
  
  data <- data.frame(id = 1:1000, value = rnorm(1000))
  
  # add group numbers after shuffling
  data <- as.data.frame(data[sample(nrow(data)), ])
  
  data$group <- sjmisc::split_var(data$id, n = num.group)
  
  data$treat <- rbinom(n = nrow(data), size = 1,
                       prob = .5)
  
  group_dums <- as.data.frame(
    model.matrix( ~ group - 1,
                  data = data)) 
  
  data <- bind_cols(data, group_dums)
  
  # combine group dummies in eqn
  group.preds <- str_flatten(colnames(group_dums)[1:(length(colnames(group_dums)) - 1)], 
                             collapse = "*")
  
  formula.out <- as.formula(paste("~ treat *(", group.preds, ")",
                                  collapse = " "))
  pred.mat <- model.matrix(formula.out,
                           data = data)
  
  beta <- rnorm(n = ncol(pred.mat), sd = .25)
  beta[2] <- .5
  
  # outcome- normally distributed
  mu.y <- as.numeric(pred.mat %*% beta)
  
  y <- rnorm(nrow(data), mean = mu.y, sd = sd.coef)
  
  # add outcome to data
  data <- bind_cols(data, "mu.y" = mu.y,
                    "y" = y) %>%
    select(-value)
  
  
  grid.calc <- data %>%
    group_by(treat, group) %>%
    summarize(
      n = n(),
      mu.y = mean(mu.y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      id_cols = c(group),
      names_from = treat,
      values_from = c(n, mu.y)
    ) %>%
    mutate( # treatment effect for each group
      true = mu.y_1 - mu.y_0
    ) %>%
    select(
      group, true
    )
  
  #print(head(data))
  formula.ols <- as.formula(paste("y ~ treat *(", group.preds, ")",
                                  collapse = " "))
  
  formula.vs <- as.formula(paste("y ~ (1 + treat |", group.preds, ")",
                                 collapse = " "))
  
  out = list("data" = data, 
             "coef" = beta,
             "treat.effect" = grid.calc,
             "formula.ols" = formula.ols,
             "formula.vs" = formula.vs)
  
}

# set combinations
combos.gr <- expand.grid(num.group, sd.coef)
colnames(combos.gr) <- c("num.group", "sd.coef")
combos.gr$combo <- rownames(combos.gr)
combos.gr$pair <- paste0(combos.gr$num.group, "_", combos.gr$sd.coef)  

combos.gr.list <- split(combos.gr, f = combos.gr$pair)

# run it
sim.data.group <- lapply(combos.gr.list,
                         group.data.sim)


### model grouped data
simulation.gr.model <- function(list){
  
  data <- list$data
  coef <- list$coef
  treat.effect <- list$treat.effect
  formula.ols <- list$formula.ols
  formula.vs <- list$formula.vs
  
  
  # hypothetical data grid for slopes
  grid.sim <- data %>%
    select(starts_with("group")) %>%
    distinct() 
  
  
  ### OLS regression
  # fit model
  ols.inter <- lm(formula.ols,
                  data = data)
  summary(ols.inter)
  
  # get slopes
  slopes.ols <- slopes(model = ols.inter,
                       variables = "treat",
                       newdata = grid.sim)  %>%
    left_join(treat.effect) %>%
    mutate(
      bias = estimate - true,
      in.interval = ifelse(true > conf.low & true < conf.high, 1, 0)
    )
  
  
  ### hierarchical regression
  # VS model
  sim.het.prior <- c(
    prior(normal(0, .75), class = "Intercept"),
    prior(normal(0, 1), class = "sd"),
    prior(normal(0, 1), class = "sigma")
  )

  start.time <- Sys.time()
  vs.mod <- brm(bf(formula.vs),
                data = data,
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
    left_join(treat.effect) %>%
    mutate(
      bias = estimate - true,
      in.interval = ifelse(true > conf.low & true < conf.high, 1, 0)
    )

  # output everything in a list
  sim.input <- list("data" = data,
                    "coef" = coef, "treat.effect" = treat.effect)
  ols.res <- list("ols.mod" = ols.inter, "slopes.ols" = slopes.ols)
  vs.res <- list("vs.mod" = vs.mod, "slopes.vs" = slopes.vs,
                 "time.vs" = time.taken)
  
  output  <- list(sim.input, ols.res, vs.res)
  output 
  
  
}


# run it
res.sim.group <- lapply(sim.data.group,
                        simulation.gr.model)
# save results
saveRDS(res.sim.group, file = "data/simulation-res-group.rds")
# load 
res.sim.group <- readRDS(file = "data/simulation-res-group.rds")


# return data
combos.gr.ret <- bind_rows(combos.gr.list)

# for loop given nested lists
sim.gr.res <- vector(mode = "list", length = length(res.sim.group))
slopes.gr.res <- vector(mode = "list", length = length(res.sim.group)) 

for(i in 1:length(res.sim.group)){
  
  est <- unlist(res.sim.group[[i]], recursive = FALSE)
  
  # predictions
  pred.vs <- posterior_epred(object = est$vs.mod)
  pred.vs.med <- apply(pred.vs, 2, median)
  
  
  # RMSE OLS:
  # outcome
  rmse.ols <- sqrt(mean(est$ols.mod$residuals^2))
  # group coefs
  ols.est <- est$slopes.ols$estimate
  bias.ols <- est$slopes.ols$bias
  rmse.coef.ols <- sqrt(mean((est$slopes.ols$true - ols.est)^2))
  mean(bias.ols)
  
  # RMSE VS:
  # outcome
  resid.vs <- est$outcome - pred.vs.med
  rmse.vs <- sqrt(mean((est$outcome - pred.vs.med)^2))
  
  # group coefs
  bias.vs <-  est$slopes.vs$true - est$slopes.vs$estimate
  rmse.coef.vs <- sqrt(mean((est$slopes.vs$true - est$slopes.vs$estimate)^2))
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
                    num.group = combos.gr.ret$num.group[i],
                    sd.coef = combos.gr.ret$sd.coef[i]
  )
  
  sim.gr.res[[i]] <- res
  slopes.gr.res[[i]] <- bind_rows("OLS" = est$slopes.ols,
                                  "Hierarchical" = est$slopes.vs,
                                  .id = "model") %>%
    mutate(
      num.group = combos.gr.ret$num.group[i],
      sd.coef = combos.gr.ret$sd.coef[i]
    )
}


sim.gr.res <- bind_rows(sim.gr.res,
                        .id = "sim") %>%
  mutate(
    sd.coef = paste0("Coef SD=", sd.coef),
    scen = paste0("Groups=", num.group, ",\n",
                  sd.coef)
  ) %>%
  select(-bias) %>%
  distinct() %>%
  group_by(scen) %>%
  mutate(
    change.rmse.coef = rmse.coef - lag(rmse.coef),
    change.rmse.out = rmse.out - lag(rmse.out)
  )

ggplot(sim.gr.res, aes(x = factor(num.group), y = change.rmse.coef)) +
  facet_wrap(~ sd.coef, scales = "free_x") +
  geom_hline(yintercept = 0) +
  geom_point(size = 3) +
  labs(title = "Change in RMSE of Varying Slopes Compared to OLS",
       x = "Number of Groups",
       y = "RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")
ggsave("figures/sim-rmse-coef.png", height = 6, width = 8)


ggplot(sim.gr.res, aes(x = factor(num.group), 
                       color = model,
                       y = rmse.coef)) +
  facet_wrap(~ sd.coef) +
  geom_point(size = 3) +
  labs(title = "RMSE Treatment Estimate: Varying Slopes and OLS",
       x = "Number of Groups",
       y = "Estimate RMSE",
       color = "Model") +
  scale_color_grey(start = .6, end = .1) +
  theme(legend.position = "bottom")
ggsave("figures/sim-rmse-coef-level.png", height = 6, width = 8)


# futher treatment data- slopes 
sim.slopes.gr <- bind_rows(slopes.gr.res,
                           .id = "sim") %>%
  mutate(
    sd.impact = paste0("SD Coef=", sd.coef),
    scen = paste0("Groups=", num.group, ",\n",
                  sd.impact)
  )  %>%
  group_by(scen, rowid) %>%
  mutate(
    change.bias = bias - lag(bias)
  )

ggplot(sim.slopes.gr, aes(x = group, y = bias,
                          fill = model)) +
  geom_hline(yintercept = 0) +
  facet_grid(sd.impact ~ num.group,
             scales = "free_y") +
  geom_bar(stat = "identity",
           width = .5,
           position = position_dodge(width = .5)) 

ggplot(sim.slopes.gr, aes(x = group, y = change.bias)) +
  geom_hline(yintercept = 0) +
  facet_grid(sd.impact ~ num.group,
             scales = "free_y") +
  geom_point()


ggplot(sim.slopes.gr, aes(x = factor(in.interval),
                          fill = model
)) +
  facet_grid(sd.impact ~ num.group) +
  geom_bar(position = position_dodge(width = 1)) +
  theme(legend.position = "bottom")


ggplot(sim.slopes.gr, aes(y = factor(group),
                          x = estimate,
                          shape = sd.impact,
                          color = model)) +
  facet_wrap(~ sd.impact + num.group, scales = "free",
             ncol = 9) +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
                  position = position_dodge(width = 1)) +
  geom_text(aes(x = true), label = "x", color = "black") +
  theme(legend.position = "bottom")
ggsave("figures/sim-rmse-coef-coverage.png", height = 10, width = 12)
