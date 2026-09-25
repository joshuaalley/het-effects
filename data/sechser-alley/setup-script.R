# Joshua Alley and Todd Sechser
# Testing Trip Wires
# Set up script


# load packages
library(conflicted)
#library(rIP)
library(tidyverse)
library(systemfit)
library(nnet)
library(MASS)
library(gridExtra)
library(modelsummary)
library(marginaleffects)
library(brant)
library(ordinal)
library(haven)
library(interplot)
library(margins)
library(grid)
library(interflex)
library(wesanderson)
#library(ggcarly)
library(modelsummary)
library(brms)


# set seed
set.seed(12)

# manage conflicts 
conflict_scout()
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("Position", "ggplot2")
conflict_prefer("recode", "dplyr")
conflict_prefer("expand", "tidyr")
conflict_prefer("pack", "tidyr")
conflict_prefer("unpack", "tidyr")
conflict_prefer("combine", "dplyr")

# set default ggplot theme
theme_set(theme_bw(base_size = 14))
options("modelsummary_format_numeric_latex" = "plain")

# load information on treatments
short.vignette <- read.csv("data/vignette-content-short.csv")
long.vignette <- read.csv("data/vignette-content-long.csv")


# margins plot for each disposition
# use a function that inputs modifying variable and label 
margins.dispo.func <- function(model, factor, label){
  
  inter.plot <- cplot(model, 
                      x = factor, dx = "casualties_dummy", 
                      what = "effect",
                      rug = TRUE,
                      draw = FALSE)
  
  
  # ggplot nicely
  plot <- ggplot(inter.plot, aes(x = xvals,
                                 y = yvals)) +
    geom_blank() +
    geom_hline(yintercept = 0) +
    geom_pointrange(aes(ymin = lower, 
                        ymax = upper)) +
    xlim(c(0, 4)) + 
    ggtitle(label) +
    labs(x = label,
         y = "Estimated Marginal Effect of Casualties")
  
  # # add bar plot via: https://stackoverflow.com/questions/40136528/interaction-marginal-effects-plot-with-overlay-histogram-using-ggplot
  
  plot2 <- ggplot(data = model$model,
                  aes_string(x = factor)) +
    geom_blank() +
    geom_bar(fill="gray30", alpha = 0.75,
             width = .3) +
    labs(x = paste("Distribution of", label),
         y = "Number of Respondents")
  
  grid.arrange(plot, plot2, nrow = 2)
  plots <- list(inter.plot, plot, plot2)
  
}