# Joshua Alley
# Alternative applications: separate pairwise interactions versus one model
#
# Published experiments that estimate several modifiers of one treatment in
# separate interaction models. None of them fully crosses its modifiers, so the
# comparisons are the authors' separate models, one additive model with every
# modifier, and the hierarchical model. For each paper:
#   1. reproduce the authors' separate models
#   2. fit one additive interaction model with every modifier
#   3. fit the hierarchical model in lme4, and again in brms, since REML
#      reports a group SD of exactly zero whenever the likelihood peaks there
#   4. split-half check: estimate the authors' subgroup effects, and effects in
#      the hierarchical model's groups, on half the sample and score them
#      against raw differences in the other half
# Assumes data/setup-script.R has been sourced.
#
# data, from Harvard Dataverse on 2026-09-15
#   ../barber-pope-rep/Ideology_Trump.dta               Barber & Pope, APSR 2019, doi:10.7910/DVN/38BFML
#   buzas/studyI_data.csv                          Buzas & Bassan-Nygate, AJPS 2024, doi:10.7910/DVN/NWDZJB
#   webster/anger_social.RData                     Webster, Connors & Sinclair, JOP 2021, doi:10.7910/DVN/0XS5NX
#   jacob/recoded_data.csv                         Jacob, AJPS 2023, doi:10.7910/DVN/XERHRF
#   reeves-rogowski/RR-AJPS-taps-processed.RData   Reeves & Rogowski, AJPS 2017, doi:10.7910/DVN/XRJ43D

ALT <- "data/alt-examples"
N_SPLITS <- 200
N_SPLITS_WEBSTER <- 100 # sixteen outcomes
RUN <- c("barber", "buzas", "webster", "jacob", "reeves")
WORKERS <- max(1, min(8, parallelly::availableCores() - 1))
dir.create(file.path(ALT, "fits"), showWarnings = FALSE)


### tools

# interaction terms of a fit, whichever side of the colon the treatment is on
int_rows <- function(fit, treat, vc = NULL, label = "") {
  b <- if (inherits(fit, "merMod")) lme4::fixef(fit) else coef(fit)
  if (is.null(vc)) vc <- as.matrix(vcov(fit))
  pat <- paste0("(^", treat, ":)|(:", treat, "$)")
  nm <- names(b)[grepl(pat, names(b)) & !is.na(b)]
  se <- sqrt(diag(vc)[nm])
  tibble(model = label, term = gsub(pat, "", nm), est = unname(b[nm]),
         se = unname(se), p = 2 * pnorm(-abs(unname(b[nm]) / unname(se))))
}

# posterior mean and 95% interval for the same terms in a brms fit
brms_int_rows <- function(fit, treat) {
  fe <- brms::fixef(fit)
  pat <- paste0("(^", treat, ":)|(:", treat, "$)")
  nm <- rownames(fe)[grepl(pat, rownames(fe))]
  tibble(term = gsub(pat, "", nm), brms_est = round(fe[nm, "Estimate"], 4),
         brms_lo = round(fe[nm, "Q2.5"], 4), brms_hi = round(fe[nm, "Q97.5"], 4))
}

# effect of the treatment for every row, from any fit with a predict method
row_effects <- function(fit, d, treat) {
  d1 <- d
  d1[[treat]] <- 1
  d0 <- d
  d0[[treat]] <- 0
  suppressWarnings(predict(fit, newdata = d1) - predict(fit, newdata = d0))
}

# raw treated-minus-control difference within each set of rows
raw_diffs <- function(d, y, treat, idx) {
  vapply(idx, function(i) {
    t1 <- d[[treat]][i] == 1
    if (!any(t1) || all(t1)) return(NA_real_)
    mean(d[[y]][i][t1]) - mean(d[[y]][i][!t1])
  }, numeric(1))
}

subgroup_rows <- function(d, defs) lapply(defs, function(e) which(eval(e, d)))

fit_formulas <- function(spec) {
  mods <- paste(spec$mods, collapse = " + ")
  list(
    joint = as.formula(paste(spec$y, "~", spec$treat, "* (", mods, ")", spec$controls)),
    hier = as.formula(paste(spec$y, "~", spec$treat, "* (", mods, ")", spec$controls,
                            spec$re_extra, "+ (1 +", spec$treat, "|| grp)"))
  )
}

fit_hier <- function(f, d) {
  suppressMessages(suppressWarnings(lme4::lmer(
    f, data = d, REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE))))
}

sigma_theta <- function(fit, treat) {
  vc <- as.data.frame(lme4::VarCorr(fit))
  vc$sdcor[vc$var1 %in% treat & is.na(vc$var2)][1]
}

# the same hierarchical model in brms, default priors, fits cached
fit_brms <- function(f, d, label) {
  brm(f, data = d, family = gaussian(),
      cores = 4, backend = "cmdstanr",
      control = list(adapt_delta = .95),
      seed = 12, refresh = 0,
      file = file.path(ALT, "fits", label), file_refit = "on_change")
}

# group SD for the treatment: REML estimate beside the posterior
sigma_compare <- function(brms_fit, lme4_fit, treat, label) {
  s <- posterior::as_draws_df(brms_fit)[[paste0("sd_grp__", treat)]]
  np <- brms::nuts_params(brms_fit)
  tibble(model = label, n = nrow(brms_fit$data),
         groups = length(unique(brms_fit$data$grp)),
         lme4 = round(sigma_theta(lme4_fit, treat), 4),
         brms_median = round(median(s), 4),
         brms_05 = round(unname(quantile(s, .05)), 4),
         brms_95 = round(unname(quantile(s, .95)), 4),
         divergent = sum(np$Value[np$Parameter == "divergent__"]),
         max_rhat = round(max(brms::rhat(brms_fit), na.rm = TRUE), 3))
}

# subgroup effects in the full sample: separate models, joint model, lme4
subgroup_estimates <- function(d, spec) {
  f <- fit_formulas(spec)
  sg <- subgroup_rows(d, spec$subgroups)
  avg <- function(e) vapply(sg, function(i) mean(e[i]), numeric(1))
  fit_h <- fit_hier(f$hier, d)
  sep <- setNames(rep(NA_real_, length(sg)), names(sg))
  for (s in spec$separate) {
    fs <- as.formula(paste(spec$y, "~", spec$treat, "*", s$term, spec$controls, s$extra))
    e <- row_effects(stats::lm(fs, data = d), d, spec$treat)
    for (nm in s$subgroups) sep[nm] <- mean(e[sg[[nm]]])
  }
  list(
    fit = fit_h,
    table = tibble(subgroup = names(sg), n = lengths(sg),
                   raw = raw_diffs(d, spec$y, spec$treat, sg), separate = sep,
                   joint = avg(row_effects(stats::lm(f$joint, data = d), d, spec$treat)),
                   hier = avg(row_effects(fit_h, d, spec$treat)))
  )
}

# one split: halves drawn within treatment-by-group strata, at the level of
# the respondent when respondents answer several items
split_once <- function(split_id, d, spec) {
  f <- fit_formulas(spec)
  cl <- if (is.null(spec$cluster)) seq_len(nrow(d)) else d[[spec$cluster]]
  first <- !duplicated(cl)
  strata <- paste(d[[spec$treat]], d$grp)[first]
  half_of <- integer(sum(first))
  for (s in unique(strata)) {
    k <- which(strata == s)
    half_of[k] <- sample(rep(1:2, length.out = length(k)))
  }
  h <- half_of[match(cl, cl[first])]
  A <- d[h == 1, ]
  B <- d[h == 2, ]

  sg_a <- subgroup_rows(A, spec$subgroups)
  sg_b <- subgroup_rows(B, spec$subgroups)
  g_a <- split(seq_len(nrow(A)), as.character(A$grp))
  g_b <- split(seq_len(nrow(B)), as.character(B$grp))
  g_names <- intersect(names(g_a), names(g_b))
  tgt_sg <- raw_diffs(B, spec$y, spec$treat, sg_b)
  tgt_g <- raw_diffs(B, spec$y, spec$treat, g_b[g_names])
  ok_sg <- !is.na(tgt_sg) & !is.na(raw_diffs(A, spec$y, spec$treat, sg_a))
  ok_g <- !is.na(tgt_g) & !is.na(raw_diffs(A, spec$y, spec$treat, g_a[g_names]))

  e_joint <- row_effects(stats::lm(f$joint, data = A), A, spec$treat)
  fit_h <- fit_hier(f$hier, A)
  e_hier <- row_effects(fit_h, A, spec$treat)
  sep <- setNames(rep(NA_real_, length(sg_a)), names(sg_a))
  for (s in spec$separate) {
    fs <- as.formula(paste(spec$y, "~", spec$treat, "*", s$term, spec$controls, s$extra))
    e <- row_effects(stats::lm(fs, data = A), A, spec$treat)
    for (nm in s$subgroups) sep[nm] <- mean(e[sg_a[[nm]]])
  }

  avg <- function(e, idx) vapply(idx, function(i) mean(e[i]), numeric(1))
  err <- function(est, tgt, ok) (est[ok] - tgt[ok])^2
  list(
    sg = data.frame(split = split_id, subgroup = names(sg_a)[ok_sg],
                    separate = err(sep, tgt_sg, ok_sg),
                    joint = err(avg(e_joint, sg_a), tgt_sg, ok_sg),
                    hier = err(avg(e_hier, sg_a), tgt_sg, ok_sg)),
    g = data.frame(split = split_id, group = g_names[ok_g],
                   joint = err(avg(e_joint, g_a[g_names]), tgt_g, ok_g),
                   hier = err(avg(e_hier, g_a[g_names]), tgt_g, ok_g)),
    sigma_theta = sigma_theta(fit_h, spec$treat)
  )
}

run_splits <- function(d, spec, n) {
  future.apply::future_lapply(seq_len(n), split_once, d = d, spec = spec,
                              future.seed = 12L)
}

mc_se <- function(x) sd(x) / sqrt(length(x))

# held-out error by approach. Every score includes the noise in the held-out
# half, which is the same for every approach, so differences are differences
# in mean squared error.
split_summary <- function(res) {
  sg <- bind_rows(lapply(res, `[[`, "sg"))
  g <- bind_rows(lapply(res, `[[`, "g"))
  sgs <- sg %>% group_by(split) %>%
    summarize(across(c(separate, joint, hier), mean), .groups = "drop")
  gs <- g %>% group_by(split) %>%
    summarize(n_groups = n(), across(c(joint, hier), mean), .groups = "drop")
  sig <- vapply(res, `[[`, numeric(1), "sigma_theta")
  list(
    overall = tibble(
      sg_separate = mean(sgs$separate), sg_joint = mean(sgs$joint), sg_hier = mean(sgs$hier),
      sg_hier_vs_sep = mean(sgs$hier - sgs$separate), sg_hier_vs_sep_se = mc_se(sgs$hier - sgs$separate),
      sg_hier_beats_sep = mean(sgs$hier < sgs$separate),
      sg_joint_beats_sep = mean(sgs$joint < sgs$separate),
      sg_hier_beats_joint = mean(sgs$hier < sgs$joint),
      g_n = mean(gs$n_groups), g_joint = mean(gs$joint), g_hier = mean(gs$hier),
      g_hier_vs_joint = mean(gs$hier - gs$joint), g_hier_vs_joint_se = mc_se(gs$hier - gs$joint),
      g_hier_beats_joint = mean(gs$hier < gs$joint),
      sigma_theta_median = median(sig), sigma_theta_zero = mean(sig < 1e-6)
    ),
    by_subgroup = sg %>% group_by(subgroup) %>%
      summarize(splits = n(), separate = mean(separate), joint = mean(joint),
                hier = mean(hier), hier_vs_sep = mean(hier - separate),
                hier_beats_sep = mean(hier < separate), .groups = "drop") %>%
      arrange(hier_vs_sep)
  )
}

print_split_summary <- function(s, label) {
  o <- s$overall
  cat(sprintf("\n--- %s: split-half ---\n", label))
  cat(sprintf("authors' subgroup effects: separate %.5f  joint %.5f  hierarchical %.5f\n",
              o$sg_separate, o$sg_joint, o$sg_hier))
  cat(sprintf("  hierarchical - separate %+.5f (MC SE %.5f); hierarchical better in %.0f%% of splits, joint better than separate in %.0f%%\n",
              o$sg_hier_vs_sep, o$sg_hier_vs_sep_se, 100 * o$sg_hier_beats_sep, 100 * o$sg_joint_beats_sep))
  cat(sprintf("  hierarchical better than joint in %.0f%% of splits\n", 100 * o$sg_hier_beats_joint))
  cat(sprintf("hierarchical model's groups (%.0f scored per split): joint %.5f  hierarchical %.5f\n",
              o$g_n, o$g_joint, o$g_hier))
  cat(sprintf("  hierarchical - joint %+.5f (MC SE %.5f); better in %.0f%% of splits\n",
              o$g_hier_vs_joint, o$g_hier_vs_joint_se, 100 * o$g_hier_beats_joint))
  cat(sprintf("half-sample lme4 sigma_theta: median %.4f, zero in %.0f%% of splits\n",
              o$sigma_theta_median, 100 * o$sigma_theta_zero))
  cat("by subgroup (mean squared held-out error):\n")
  print(as.data.frame(s$by_subgroup %>%
                        mutate(across(c(separate, joint, hier, hier_vs_sep), ~ signif(.x, 3)),
                               hier_beats_sep = round(hier_beats_sep, 2))), row.names = FALSE)
}

# one row per analysis, for papers with many outcomes
split_rows <- function(s, label) {
  s$overall %>%
    transmute(analysis = label,
              sg_separate = signif(sg_separate, 3), sg_joint = signif(sg_joint, 3),
              sg_hier = signif(sg_hier, 3),
              hier_beats_sep = round(sg_hier_beats_sep, 2),
              joint_beats_sep = round(sg_joint_beats_sep, 2),
              groups_hier_beats_joint = round(g_hier_beats_joint, 2),
              lme4_sigma_zero = round(sigma_theta_zero, 2))
}

# separate, joint and hierarchical interaction estimates side by side
side_by_side <- function(sep, joint, hier) {
  sep %>% select(term, sep_est = est, sep_p = p) %>%
    full_join(joint %>% select(term, joint_est = est, joint_p = p), by = "term") %>%
    full_join(hier %>% select(term, hier_est = est, hier_p = p), by = "term") %>%
    mutate(change = case_when(
      sep_p < .05 & joint_p >= .05 ~ "significant alone only",
      sep_p >= .05 & joint_p < .05 ~ "significant jointly only",
      sign(sep_est) != sign(joint_est) & pmin(sep_p, joint_p) < .10 ~ "sign flips",
      TRUE ~ ""),
      across(ends_with("_est"), ~ round(.x, 4)),
      across(ends_with("_p"), ~ round(.x, 3)))
}

future::plan(future::multisession, workers = WORKERS)


### Barber & Pope (APSR 2019): Trump cues and policy support
if ("barber" %in% RUN) {
  cat("\n\n############ Barber & Pope (APSR 2019) ############\n")
  bp <- readstata13::read.dta13("data/barber-pope-rep/Ideology_Trump.dta")
  bp_arms <- bp %>% filter(contrump == 1 | self == 1 | libtrump == 1)

  # Table 1: one interacted modifier per model, the others as controls.
  # Respondents answer ten items, so SEs are clustered by respondent here;
  # the published table treats the 9,000 rows as independent.
  bp_tab1 <- list(
    knowledge = lm(Support ~ libtrump*knowledge + contrump*knowledge + trump_approve +
                     ideo5b + republican + party_strength + race_white, data = bp_arms),
    party_strength = lm(Support ~ libtrump*party_strength + contrump*party_strength +
                          knowledge + ideo5b + trump_approve + republican + race_white,
                        data = bp_arms %>% filter(pid7 %in% 4:7)),
    trump_approve = lm(Support ~ libtrump*trump_approve + contrump*trump_approve +
                         knowledge + ideo5b + republican + party_strength + race_white,
                       data = bp_arms),
    ideo5b = lm(Support ~ libtrump*ideo5b + contrump*ideo5b + knowledge + trump_approve +
                  republican + party_strength + race_white, data = bp_arms)
  )
  bp_sep <- function(treat) bind_rows(lapply(names(bp_tab1), function(m) {
    fit <- bp_tab1[[m]]
    int_rows(fit, treat, vc = sandwich::vcovCL(fit, cluster = ~caseid), label = m) %>%
      mutate(p_published = 2 * pnorm(-abs(est / sqrt(diag(vcov(fit))[
        grepl(paste0("(^", treat, ":)|(:", treat, "$)"), names(coef(fit)))]))))
  }))

  # The joint model uses the directional seven-point party scale. Party dummies
  # plus the folded strength measure would force strong Democrats and strong
  # Republicans to respond to a Trump cue the same way; the authors avoid that
  # by fitting party strength only among independents and Republicans.
  bp_mods <- c("knowledge", "pid7", "trump_approve", "ideo5b")
  bp_cc <- bp_arms %>%
    filter(pid7 %in% 1:7) %>%
    mutate(party3 = case_when(pid7 >= 5 ~ "R", pid7 <= 3 ~ "D", TRUE ~ "I"),
           approve3 = cut(trump_approve, c(0, 2, 3, 5), labels = c("dis", "neu", "app")),
           know2 = as.numeric(knowledge >= 5)) %>%
    select(caseid, Support, libtrump, contrump, race_white, republican, democrat,
           party_strength, party3, approve3, know2, all_of(bp_mods)) %>%
    na.omit() %>%
    mutate(grp = interaction(party3, approve3, know2, drop = TRUE))
  cat(sprintf("party3 from pid7 agrees with the authors' republican and democrat dummies: %.3f\n",
              mean((bp_cc$party3 == "R") == (bp_cc$republican == 1) &
                     (bp_cc$party3 == "D") == (bp_cc$democrat == 1))))
  bp_f_joint <- as.formula(paste("Support ~ (libtrump + contrump) * (",
                                 paste(bp_mods, collapse = " + "), ") + race_white"))
  bp_f_hier <- as.formula(paste("Support ~ (libtrump + contrump) * (",
                                paste(bp_mods, collapse = " + "),
                                ") + race_white + (1 | caseid) + (1 + libtrump + contrump || grp)"))
  bp_joint_fit <- lm(bp_f_joint, data = bp_cc)
  bp_joint_vc <- sandwich::vcovCL(bp_joint_fit, cluster = ~caseid)
  bp_hier_fit <- fit_hier(bp_f_hier, bp_cc)
  bp_brms_fit <- fit_brms(bp_f_hier, bp_cc, "bp-both-cues")
  for (treat in c("libtrump", "contrump")) {
    cat(sprintf("\n%s interactions: separate (Table 1; clustered p, published p) vs joint vs lme4 vs brms (95%% interval)\n", treat))
    cat("party_strength is fit among independents and Republicans; pid7 is its full-sample counterpart\n")
    sep <- bp_sep(treat)
    print(as.data.frame(
      side_by_side(sep, int_rows(bp_joint_fit, treat, vc = bp_joint_vc), int_rows(bp_hier_fit, treat)) %>%
        left_join(sep %>% select(term, p_published) %>% mutate(p_published = round(p_published, 3)),
                  by = "term") %>%
        left_join(brms_int_rows(bp_brms_fit, treat), by = "term")), row.names = FALSE)
  }

  # do the correlated modifiers moderate each cue jointly, even when none does alone?
  bp_wald <- function(terms) {
    b <- coef(bp_joint_fit)[terms]
    stat <- as.numeric(t(b) %*% solve(bp_joint_vc[terms, terms]) %*% b)
    sprintf("chi2(%d) = %.2f, p = %.4f", length(terms), stat,
            pchisq(stat, length(terms), lower.tail = FALSE))
  }
  bp_nm <- names(coef(bp_joint_fit))
  for (treat in c("libtrump", "contrump")) {
    pick <- function(m) bp_nm[grepl(paste0("(^", treat, ":", m, "$)|(^", m, ":", treat, "$)"), bp_nm)]
    cat(sprintf("%s, joint test of approval + ideology + party: %s; all four modifiers: %s\n", treat,
                bp_wald(c(pick("trump_approve"), pick("ideo5b"), pick("pid7"))),
                bp_wald(c(pick("knowledge"), pick("trump_approve"), pick("ideo5b"), pick("pid7")))))
  }

  cat("\ngroup SD for each cue: lme4 and brms\n")
  print(as.data.frame(bind_rows(
    sigma_compare(bp_brms_fit, bp_hier_fit, "libtrump", "liberal cue"),
    sigma_compare(bp_brms_fit, bp_hier_fit, "contrump", "conservative cue"))), row.names = FALSE)

  # split-half, one cue at a time against the no-cue condition
  bp_spec <- list(
    y = "Support", treat = "trt", mods = bp_mods, controls = "+ race_white",
    re_extra = "+ (1 | caseid)", cluster = "caseid",
    subgroups = list(
      republicans = quote(party3 == "R"), democrats = quote(party3 == "D"),
      independents = quote(party3 == "I"),
      disapprove = quote(trump_approve <= 2), neutral = quote(trump_approve == 3),
      approve = quote(trump_approve >= 4),
      low_knowledge = quote(knowledge <= 3), mid_knowledge = quote(knowledge %in% 4:5),
      high_knowledge = quote(knowledge >= 6),
      liberal = quote(ideo5b <= 2), moderate = quote(ideo5b == 3), conservative = quote(ideo5b >= 4)),
    # party strength is left out: the authors fit it only among independents
    # and Republicans, so there is no full-sample separate estimate to compare
    separate = list(
      list(term = "factor(party3)", extra = "", subgroups = c("republicans", "democrats", "independents")),
      list(term = "trump_approve", extra = "+ knowledge + ideo5b + republican + party_strength",
           subgroups = c("disapprove", "neutral", "approve")),
      list(term = "knowledge", extra = "+ trump_approve + ideo5b + republican + party_strength",
           subgroups = c("low_knowledge", "mid_knowledge", "high_knowledge")),
      list(term = "ideo5b", extra = "+ knowledge + trump_approve + republican + party_strength",
           subgroups = c("liberal", "moderate", "conservative")))
  )
  for (arm in c("libtrump", "contrump")) {
    d_arm <- bp_cc %>% filter(.data[[arm]] == 1 | (libtrump == 0 & contrump == 0)) %>%
      mutate(trt = .data[[arm]])
    full <- subgroup_estimates(d_arm, bp_spec)
    cat(sprintf("\n%s vs no cue: subgroup effects in the full sample (n = %d rows)\n", arm, nrow(d_arm)))
    print(as.data.frame(full$table %>% mutate(across(raw:hier, ~ round(.x, 3)))), row.names = FALSE)
    print_split_summary(split_summary(run_splits(d_arm, bp_spec, N_SPLITS)),
                        paste("Barber & Pope,", arm))
  }
}


### Buzas & Bassan-Nygate (AJPS 2024): shaming Israel, Study I
if ("buzas" %in% RUN) {
  cat("\n\n############ Buzas & Bassan-Nygate (AJPS 2024), Study I ############\n")
  agree5 <- function(x) case_when(
    x == "Strongly agree" ~ 100, x == "Somewhat agree" ~ 75,
    x == "Neither agree nor disagree" ~ 50, x == "Somewhat disagree" ~ 25,
    x == "Strongly disagree" ~ 0)
  religions <- c("Protestant", "Other", "Atheist/agnostic", "Nothing in particular", "Buddhist",
                 "Roman Catholic", "Mormon", "Orthodox (Greek or Russian)", "Jewish", "Muslim", "Hindu")
  parties <- c("Strong Republican", "Republican", "Independent, but Lean Republican", "Independent",
               "Independent, but Lean Democrat", "Democrat", "Strong Democrat")
  # recodes follow analysis_I.R; leaners count as partisans
  bb <- read.csv(file.path(ALT, "buzas", "studyI_data.csv"))
  bb <- bb[3:nrow(bb), ] %>%
    mutate(
      shaming = case_when(condition == "1" ~ 0, condition == "2" ~ 1),
      support_Israel = (agree5(israel_sentiment_1) + (100 - agree5(israel_sentiment_2)) +
                          agree5(israel_sentiment_3) + (100 - agree5(israel_sentiment_4)) +
                          (100 - agree5(israel_sentiment_5)) + (100 - agree5(israel_sentiment_6))) / 6,
      antisemitism = (agree5(antisemitism_1) + agree5(antisemitism_2) + agree5(antisemitism_3) +
                        agree5(antisemitism_4) + agree5(antisemitism_5)) / 5,
      Republican = ifelse(partisanship %in% parties, as.numeric(partisanship %in% parties[1:3]), NA),
      Democrat = ifelse(partisanship %in% parties, as.numeric(partisanship %in% parties[5:7]), NA),
      Independent = ifelse(partisanship %in% parties, as.numeric(partisanship == "Independent"), NA),
      Christian = ifelse(religion %in% religions, as.numeric(religion %in% c("Protestant", "Roman Catholic")), NA),
      Jewish = ifelse(religion %in% religions, as.numeric(religion == "Jewish"), NA),
      Muslim = ifelse(religion %in% religions, as.numeric(religion == "Muslim"), NA),
      party3 = case_when(Republican == 1 ~ "R", Democrat == 1 ~ "D", Independent == 1 ~ "I"),
      religion4 = case_when(Christian == 1 ~ "Christian", Jewish == 1 ~ "Jewish",
                            Muslim == 1 ~ "Muslim", !is.na(Christian) ~ "other")
    )
  bb_mods6 <- c("Democrat", "Republican", "Independent", "Christian", "Jewish", "Muslim")
  bb_sigma <- list()
  for (y in c("support_Israel", "antisemitism")) {
    # Tables D1 and D3: each modifier in its own model
    sep <- bind_rows(lapply(bb_mods6, function(m) {
      int_rows(lm(as.formula(paste(y, "~ shaming *", m)), data = bb), "shaming", label = m)
    }))
    d_y <- bb %>%
      select(all_of(y), shaming, party3, religion4, all_of(bb_mods6)) %>%
      na.omit() %>%
      mutate(grp = interaction(party3, religion4, drop = TRUE))
    spec_bb <- list(
      y = y, treat = "shaming", mods = c("Democrat", "Republican", "Christian", "Jewish", "Muslim"),
      controls = "", re_extra = "", cluster = NULL,
      subgroups = setNames(lapply(bb_mods6, function(m) bquote(.(as.name(m)) == 1)), bb_mods6),
      separate = lapply(bb_mods6, function(m) list(term = m, extra = "", subgroups = m))
    )
    f <- fit_formulas(spec_bb)
    full <- subgroup_estimates(d_y, spec_bb)
    brms_fit <- fit_brms(f$hier, d_y, paste0("bb-", y))
    bb_sigma[[y]] <- sigma_compare(brms_fit, full$fit, "shaming", y)
    cat(sprintf("\n%s: separate (Tables D1, D3) vs joint vs lme4 vs brms; n = %d, %d groups\n",
                y, nrow(d_y), nlevels(d_y$grp)))
    print(as.data.frame(side_by_side(sep, int_rows(lm(f$joint, data = d_y), "shaming"),
                                     int_rows(full$fit, "shaming")) %>%
                          left_join(brms_int_rows(brms_fit, "shaming"), by = "term")), row.names = FALSE)
    print(as.data.frame(full$table %>% mutate(across(raw:hier, ~ round(.x, 2)))), row.names = FALSE)
    print_split_summary(split_summary(run_splits(d_y, spec_bb, N_SPLITS)),
                        paste("Buzas & Bassan-Nygate,", y))
  }
  cat("\ngroup SD for shaming: lme4 and brms\n")
  print(as.data.frame(bind_rows(bb_sigma)), row.names = FALSE)
}


### Webster, Connors & Sinclair (JOP 2021): anger and social distance
if ("webster" %in% RUN) {
  cat("\n\n############ Webster, Connors & Sinclair (JOP 2021) ############\n")
  load(file.path(ALT, "webster", "anger_social.RData")) # dat
  freq5 <- c("Never", "Once in a while", "Some of the time", "Most of the time", "Always")
  qual5 <- c("Very poor", "Poor", "Fair", "Good", "Excellent")
  # named so it cannot be masked by the ideo7 column inside mutate()
  ideo7_levels <- c("Very liberal", "Liberal", "Slightly liberal", "Moderate; middle of the road",
                    "Slightly conservative", "Conservative", "Very conservative")
  wcs_levels <- list(
    fourpack_1 = c("Always", "Most of the time", "About half the time", "Sometimes", "Never"),
    talking = c("Continue talking to them, including about politics",
                "Continue talking to them, but not about politics",
                "Try to find a polite way out of the conversation",
                "Leave the conversation without worrying about being polite",
                "Attack their political views"),
    meal = c("Certainly go", "Go if I had nothing better to do", "Try to find a polite way to say no",
             "Say no", "Say no and attack their political views"),
    club = c("Certainly go", "Go if I had nothing better to do", "Try to find a polite way to say no",
             "Say no", "Say no and talk badly about those who are attending to your friends"),
    date = c("Certainly go", "Go if I had nothing better to do", "Try to find a polite way to say no",
             "Say no", "Say no and talk badly about the person to your friends"),
    family_date = c("Certainly", "Yes, but only if the person was otherwise a good person",
                    "Probably not", "Certainly not"),
    friendship = c("I would change nothing about my friendship",
                   "I would remain friends with them but would not discuss politics",
                   "I would remain friends with them but I would attack their political beliefs",
                   "I would end the friendship"),
    family = c("I would treat them the exact same",
               "I would treat them the same but not talk about politics",
               "I would distance myself a bit",
               "I would cut them out of my life as much as possible",
               "I would cut them out of my life as much as possible, and I would attack their political beliefs"),
    divorce_gop = c("Very sad", "Sad", "Neither sad nor happy", "Happy", "Very happy"),
    friends_comfort = c("Extremely comfortable", "Moderately comfortable", "Slightly comfortable",
                        "Neither comfortable nor uncomfortable", "Slightly uncomfortable",
                        "Moderately uncomfortable", "Extremely uncomfortable"),
    marry_upset = c("Not at all upset", "Not too upset", "Somewhat upset", "Upset", "Extremely upset"),
    news_share = c("Extremely unlikely", "Moderately unlikely", "Slightly unlikely",
                   "Neither unlikely nor likely", "Slightly likely", "Moderately likely",
                   "Extremely likely")
  )
  wcs_levels$fourpack_2 <- wcs_levels$fourpack_3 <- wcs_levels$fourpack_4 <- wcs_levels$fourpack_1
  wcs_levels$neighbors_comfort <- wcs_levels$friends_comfort

  # recodes follow anger_social_replication.R
  wcs <- dat %>%
    filter(!(pid7 %in% c("Completely Independent", "")), Finished == "True") %>%
    mutate(
      pid7 = as.character(pid7),
      democrat = as.numeric(pid7 %in% c("Independent but lean Democrat", "Weak Democrat", "Strong Democrat")),
      female = case_when(gender == "Female" ~ 1, gender == "Male" ~ 0),
      nonwhite = ifelse(as.character(race_eth) == "", NA,
                        as.numeric(as.character(race_eth) != "White, non-Hispanic")),
      ideology = match(as.character(ideo7), ideo7_levels),
      selfmonitoring = (match(as.character(sm1), freq5) - 1) +
        (match(as.character(sm2), freq5) - 1) + (match(as.character(sm3), qual5) - 1),
      treated = as.numeric(experimentfordemocrats_DO == "treatment_dems" |
                             experimentforrepublicans_DO == "treatment_gop"),
      strong_partisan = as.numeric(pid7 %in% c("Strong Democrat", "Strong Republican")),
      ideological_extremity = abs(ideology - 4)
    )
  for (v in names(wcs_levels)) {
    x <- match(as.character(wcs[[v]]), wcs_levels[[v]]) - 1
    wcs[[paste0(v, "_01")]] <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  wcs_outcomes <- paste0(c("fourpack_1", "fourpack_2", "fourpack_3", "fourpack_4", "talking", "meal",
                           "club", "date", "family_date", "friendship", "family", "divorce_gop",
                           "friends_comfort", "neighbors_comfort", "marry_upset", "news_share"), "_01")
  wcs_mods <- c("selfmonitoring", "strong_partisan", "ideological_extremity", "female", "nonwhite")
  sm_median <- median(wcs$selfmonitoring, na.rm = TRUE)
  spec_wcs <- function(y) list(
    y = y, treat = "treated", mods = wcs_mods, controls = "+ democrat", re_extra = "", cluster = NULL,
    subgroups = list(
      strong = quote(strong_partisan == 1), weak = quote(strong_partisan == 0),
      female = quote(female == 1), male = quote(female == 0),
      nonwhite = quote(nonwhite == 1), white = quote(nonwhite == 0),
      extreme = quote(ideological_extremity >= 2), not_extreme = quote(ideological_extremity <= 1),
      high_self_monitoring = bquote(selfmonitoring > .(sm_median)),
      low_self_monitoring = bquote(selfmonitoring <= .(sm_median))),
    separate = list(
      list(term = "strong_partisan", extra = "", subgroups = c("strong", "weak")),
      list(term = "female", extra = "", subgroups = c("female", "male")),
      list(term = "nonwhite", extra = "", subgroups = c("nonwhite", "white")),
      list(term = "ideological_extremity", extra = "", subgroups = c("extreme", "not_extreme")),
      list(term = "selfmonitoring", extra = "", subgroups = c("high_self_monitoring", "low_self_monitoring")))
  )

  wcs_counts <- list()
  wcs_splits <- list()
  wcs_by_sg <- list()
  wcs_sigma <- list()
  for (y in wcs_outcomes) {
    sep <- bind_rows(lapply(wcs_mods, function(m) {
      int_rows(lm(as.formula(paste(y, "~ treated *", m, "+ democrat")), data = wcs), "treated", label = m)
    }))
    d_y <- wcs %>%
      select(all_of(y), treated, democrat, all_of(wcs_mods)) %>%
      na.omit() %>%
      mutate(grp = interaction(strong_partisan, female, nonwhite, ideological_extremity >= 2,
                               selfmonitoring > sm_median, democrat, drop = TRUE))
    spec <- spec_wcs(y)
    f <- fit_formulas(spec)
    fit_h <- fit_hier(f$hier, d_y)
    brms_fit <- fit_brms(f$hier, d_y, paste0("wcs-", y))
    wcs_sigma[[y]] <- sigma_compare(brms_fit, fit_h, "treated", y)
    wcs_counts[[y]] <- side_by_side(sep, int_rows(lm(f$joint, data = d_y), "treated"),
                                    int_rows(fit_h, "treated")) %>%
      left_join(brms_int_rows(brms_fit, "treated"), by = "term") %>%
      mutate(outcome = y)
    s <- split_summary(run_splits(d_y, spec, N_SPLITS_WEBSTER))
    wcs_splits[[y]] <- split_rows(s, y)
    wcs_by_sg[[y]] <- s$by_subgroup %>% mutate(outcome = y)
  }
  wcs_counts <- bind_rows(wcs_counts)
  cat(sprintf("\n%d separate interaction models (5 modifiers x 16 outcomes)\n", nrow(wcs_counts)))
  cat(sprintf("p < .05: separate %d, joint %d, lme4 %d; brms 95%% interval excludes zero %d\n",
              sum(wcs_counts$sep_p < .05), sum(wcs_counts$joint_p < .05, na.rm = TRUE),
              sum(wcs_counts$hier_p < .05, na.rm = TRUE),
              sum(wcs_counts$brms_lo > 0 | wcs_counts$brms_hi < 0, na.rm = TRUE)))
  print(as.data.frame(wcs_counts %>%
                        filter(sep_p < .05 | joint_p < .05 | hier_p < .05 | brms_lo > 0 | brms_hi < 0) %>%
                        select(outcome, term, sep_est, sep_p, joint_est, joint_p, hier_p, brms_lo, brms_hi, change)),
        row.names = FALSE)
  cat("\ngroup SD for the anger treatment: lme4 and brms\n")
  print(as.data.frame(bind_rows(wcs_sigma)), row.names = FALSE)
  cat("\n--- Webster et al.: split-half by outcome ---\n")
  print(as.data.frame(bind_rows(wcs_splits)), row.names = FALSE)
  cat("averaged over outcomes, by subgroup:\n")
  print(as.data.frame(bind_rows(wcs_by_sg) %>% group_by(subgroup) %>%
                        summarize(separate = signif(mean(separate), 3), hier = signif(mean(hier), 3),
                                  hier_beats_sep = round(mean(hier_beats_sep), 2), .groups = "drop") %>%
                        arrange(hier - separate)), row.names = FALSE)
}


### Jacob (AJPS 2023): sanctioning attacks on democracy in Poland
if ("jacob" %in% RUN) {
  cat("\n\n############ Jacob (AJPS 2023) ############\n")
  jac <- read.csv(file.path(ALT, "jacob", "recoded_data.csv"))
  # The paper scores liberal and majoritarian notions of democracy with an
  # ordinal CFA in lavaan (03_results_appendix.R), which is not installed here.
  # Standardized means of the same three items stand in for the factor scores.
  jac <- jac %>%
    mutate(liberal = as.numeric(scale((L1 + L2 + L3) / 3)),
           majoritarian = as.numeric(scale((M1 + M2 + M3) / 3)),
           # the undemocratic response: attacking after winning, refusing to
           # concede after losing
           attack = as.numeric(behavior == "undemocratic"),
           pis = as.numeric(choice_pre == "Law and Justice (PiS)"),
           knowledgeable = as.numeric(knowledge == "knowledgeable"))
  jac_mods <- c("pis", "knowledgeable", "liberal", "majoritarian", "inparty_strength",
                "outparty_strength", "monthly_income")
  lib_med <- median(jac$liberal)
  maj_med <- median(jac$majoritarian)
  spec_jac <- function(y) list(
    y = y, treat = "attack", mods = jac_mods, controls = "", re_extra = "", cluster = NULL,
    subgroups = list(
      pis = quote(pis == 1), ko = quote(pis == 0),
      knowledgeable = quote(knowledgeable == 1), unknowledgeable = quote(knowledgeable == 0),
      liberal_high = bquote(liberal > .(lib_med)), liberal_low = bquote(liberal <= .(lib_med)),
      majoritarian_high = bquote(majoritarian > .(maj_med)), majoritarian_low = bquote(majoritarian <= .(maj_med)),
      strong_inparty = quote(inparty_strength >= 6), weak_inparty = quote(inparty_strength <= 5),
      outparty_attached = quote(outparty_strength >= 3), outparty_detached = quote(outparty_strength <= 2),
      income_high = quote(monthly_income >= 4), income_low = quote(monthly_income <= 3)),
    separate = list(
      list(term = "pis", extra = "", subgroups = c("pis", "ko")),
      list(term = "knowledgeable", extra = "", subgroups = c("knowledgeable", "unknowledgeable")),
      list(term = "liberal", extra = "", subgroups = c("liberal_high", "liberal_low")),
      list(term = "majoritarian", extra = "", subgroups = c("majoritarian_high", "majoritarian_low")),
      list(term = "inparty_strength", extra = "", subgroups = c("strong_inparty", "weak_inparty")),
      list(term = "outparty_strength", extra = "", subgroups = c("outparty_attached", "outparty_detached")),
      list(term = "monthly_income", extra = "", subgroups = c("income_high", "income_low")))
  )

  jac_counts <- list()
  jac_sigma <- list()
  jac_splits <- list()
  jac_by_sg <- list()
  for (sc in c("win", "lose")) {
    for (y in c("approval", "democratic", "copartisan_shift")) {
      label <- paste(sc, y)
      d_sc <- jac %>% filter(scenario == sc)
      # appendix models: each modifier in its own model, within each scenario
      sep <- bind_rows(lapply(jac_mods, function(m) {
        int_rows(lm(as.formula(paste(y, "~ attack *", m)), data = d_sc), "attack", label = m)
      }))
      d_y <- d_sc %>%
        select(all_of(y), attack, all_of(jac_mods)) %>%
        na.omit() %>%
        mutate(grp = interaction(pis, knowledgeable, inparty_strength >= 6, liberal > lib_med,
                                 majoritarian > maj_med, drop = TRUE))
      spec <- spec_jac(y)
      f <- fit_formulas(spec)
      fit_h <- fit_hier(f$hier, d_y)
      brms_fit <- fit_brms(f$hier, d_y, paste0("jacob-", sc, "-", y))
      jac_sigma[[label]] <- sigma_compare(brms_fit, fit_h, "attack", label)
      jac_counts[[label]] <- side_by_side(sep, int_rows(lm(f$joint, data = d_y), "attack"),
                                          int_rows(fit_h, "attack")) %>%
        left_join(brms_int_rows(brms_fit, "attack"), by = "term") %>%
        mutate(analysis = label)
      s <- split_summary(run_splits(d_y, spec, N_SPLITS))
      jac_splits[[label]] <- split_rows(s, label)
      jac_by_sg[[label]] <- s$by_subgroup %>% mutate(analysis = label)
    }
  }
  jac_counts <- bind_rows(jac_counts)
  cat(sprintf("\n%d separate interaction models (7 modifiers x 3 outcomes x 2 scenarios)\n", nrow(jac_counts)))
  cat(sprintf("p < .05: separate %d, joint %d, lme4 %d; brms 95%% interval excludes zero %d\n",
              sum(jac_counts$sep_p < .05), sum(jac_counts$joint_p < .05, na.rm = TRUE),
              sum(jac_counts$hier_p < .05, na.rm = TRUE),
              sum(jac_counts$brms_lo > 0 | jac_counts$brms_hi < 0, na.rm = TRUE)))
  print(as.data.frame(jac_counts %>%
                        filter(sep_p < .05 | joint_p < .05 | hier_p < .05 | brms_lo > 0 | brms_hi < 0) %>%
                        select(analysis, term, sep_est, sep_p, joint_est, joint_p, hier_p, brms_lo, brms_hi, change)),
        row.names = FALSE)
  cat("\ngroup SD for the undemocratic response: lme4 and brms\n")
  print(as.data.frame(bind_rows(jac_sigma)), row.names = FALSE)
  cat("\n--- Jacob: split-half by analysis ---\n")
  print(as.data.frame(bind_rows(jac_splits)), row.names = FALSE)
  cat("averaged over analyses, by subgroup:\n")
  print(as.data.frame(bind_rows(jac_by_sg) %>% group_by(subgroup) %>%
                        summarize(separate = signif(mean(separate), 3), hier = signif(mean(hier), 3),
                                  hier_beats_sep = round(mean(hier_beats_sep), 2), .groups = "drop") %>%
                        arrange(hier - separate)), row.names = FALSE)
}


### Reeves & Rogowski (AJPS 2017): the public cost of unilateral action
if ("reeves" %in% RUN) {
  cat("\n\n############ Reeves & Rogowski (AJPS 2017) ############\n")
  rr_env <- new.env()
  load(file.path(ALT, "reeves-rogowski", "RR-AJPS-taps-processed.RData"), envir = rr_env)
  # unilateral action against legislation; the control arm is not part of the contrast
  rr <- rr_env$tapsData %>%
    filter(!is.na(treatment2)) %>%
    mutate(unilateral = treatment2,
           age = as.character(age),
           age_30_44 = as.numeric(age == "30-44"),
           age_45_59 = as.numeric(age == "45-59"),
           age_60plus = as.numeric(age == "60+"),
           know = sum.correct)
  rr_mods <- c("att", "college", "age_30_44", "age_45_59", "age_60plus", "female", "know")
  spec_rr <- list(
    y = "y", treat = "unilateral", mods = rr_mods, controls = "", re_extra = "", cluster = NULL,
    subgroups = list(
      supporters = quote(att > 0), neutral = quote(att == 0), opponents = quote(att < 0),
      college = quote(college == 1), no_college = quote(college == 0),
      age_18_29 = quote(age == "18-29"), age_30_44 = quote(age == "30-44"),
      age_45_59 = quote(age == "45-59"), age_60plus = quote(age == "60+"),
      female = quote(female == 1), male = quote(female == 0),
      high_knowledge = quote(know > 0), low_knowledge = quote(know <= 0)),
    separate = list(
      list(term = "att", extra = "", subgroups = c("supporters", "neutral", "opponents")),
      list(term = "college", extra = "", subgroups = c("college", "no_college")),
      list(term = "(age_30_44 + age_45_59 + age_60plus)", extra = "",
           subgroups = c("age_18_29", "age_30_44", "age_45_59", "age_60plus")),
      list(term = "female", extra = "", subgroups = c("female", "male")),
      list(term = "know", extra = "", subgroups = c("high_knowledge", "low_knowledge")))
  )
  rr_sep_terms <- c("att", "college", "(age_30_44 + age_45_59 + age_60plus)", "female", "know")

  rr_counts <- list()
  rr_sigma <- list()
  rr_splits <- list()
  rr_by_sg <- list()
  for (dom in c("pot", "tax", "defense")) {
    for (type in c("candidate", "handling")) {
      label <- paste(dom, type)
      d_dom <- rr %>% mutate(y = .data[[paste0(dom, ".", type, ".binary")]],
                             att = .data[[paste0(dom, ".attitudes")]])
      # Table 2 and Tables A14-A17: weighted logits, one modifier per model
      sep <- bind_rows(lapply(rr_sep_terms, function(m) {
        int_rows(glm(as.formula(paste("y ~ unilateral *", m)), data = d_dom, weights = oct2015wt1,
                     family = binomial(link = "logit")), "unilateral", label = m)
      }))
      d_y <- d_dom %>%
        select(y, unilateral, age, oct2015wt1, all_of(rr_mods)) %>%
        na.omit() %>%
        mutate(grp = interaction(college, female, age, sign(att), drop = TRUE))
      f <- fit_formulas(spec_rr)
      # the joint model in the authors' own family, so a change is not a change of model
      joint_logit <- glm(f$joint, data = d_y, weights = oct2015wt1, family = binomial(link = "logit"))
      fit_h <- fit_hier(f$hier, d_y)
      brms_fit <- fit_brms(f$hier, d_y, paste0("rr-", dom, "-", type))
      rr_sigma[[label]] <- sigma_compare(brms_fit, fit_h, "unilateral", label)
      rr_counts[[label]] <- side_by_side(sep, int_rows(joint_logit, "unilateral"),
                                         int_rows(fit_h, "unilateral")) %>%
        left_join(brms_int_rows(brms_fit, "unilateral"), by = "term") %>%
        mutate(analysis = label)
      s <- split_summary(run_splits(d_y, spec_rr, N_SPLITS))
      rr_splits[[label]] <- split_rows(s, label)
      rr_by_sg[[label]] <- s$by_subgroup %>% mutate(analysis = label)
    }
  }
  rr_counts <- bind_rows(rr_counts)
  cat("\nseparate and joint models are weighted logits, as published; lme4 and brms are linear probability models\n")
  cat(sprintf("%d interaction terms across 6 outcomes\n", nrow(rr_counts)))
  cat(sprintf("p < .05: separate %d, joint %d, lme4 %d; brms 95%% interval excludes zero %d\n",
              sum(rr_counts$sep_p < .05, na.rm = TRUE), sum(rr_counts$joint_p < .05, na.rm = TRUE),
              sum(rr_counts$hier_p < .05, na.rm = TRUE),
              sum(rr_counts$brms_lo > 0 | rr_counts$brms_hi < 0, na.rm = TRUE)))
  print(as.data.frame(rr_counts %>%
                        filter(sep_p < .05 | joint_p < .05 | hier_p < .05 | brms_lo > 0 | brms_hi < 0) %>%
                        select(analysis, term, sep_est, sep_p, joint_est, joint_p, hier_p, brms_lo, brms_hi, change)),
        row.names = FALSE)
  cat("\ngroup SD for unilateral action: lme4 and brms\n")
  print(as.data.frame(bind_rows(rr_sigma)), row.names = FALSE)
  cat("\n--- Reeves & Rogowski: split-half by outcome ---\n")
  print(as.data.frame(bind_rows(rr_splits)), row.names = FALSE)
  cat("averaged over outcomes, by subgroup:\n")
  print(as.data.frame(bind_rows(rr_by_sg) %>% group_by(subgroup) %>%
                        summarize(separate = signif(mean(separate), 3), hier = signif(mean(hier), 3),
                                  hier_beats_sep = round(mean(hier_beats_sep), 2), .groups = "drop") %>%
                        arrange(hier - separate)), row.names = FALSE)
}

future::plan(future::sequential)
