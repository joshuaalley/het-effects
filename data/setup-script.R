# Joshua Alley
# set up script for het effects paper


# key packages
library(tidyverse)
library(haven)
library(brms)
library(bayesplot)
library(marginaleffects)
library(modelsummary)
library(ggdist)
library(conflicted)
library(Matrix)

# set ggplot theme
theme_set(theme_bw(base_size = 14))

# set seed
set.seed(12)

# manage conflicts
conflict_scout()
conflict_prefer("lag", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("expand", "tidyr")
