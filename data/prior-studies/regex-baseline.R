# Joshua Alley
# Corrected regular-expression baseline for Figure 1.
#
# This re-runs the published pattern-matching pipeline so that the LLM
# classification has an honest "before" to be compared against. Three bugs in
# the original tally are fixed here; the patterns themselves are left alone,
# because changing them would mean comparing against a method the paper never
# used.
#
# Fixed:
#   1. Year. The original took the year from `published_at`, the Dataverse
#      DEPOSIT date. Journals adopted deposit mandates around 2015, so older
#      articles were archived retroactively: by deposit date 2012-2014 hold
#      1, 4 and 6 papers and 2015 holds 274, against true article counts of
#      29, 66, 71 and 133. The article year is recoverable exactly --- it
#      parses out of the Dataverse citation string for all 2,080 records and
#      agrees with `search_year` 100% of the time --- so `search_year` is used.
#   2. Aggregation. A paper is counted once if any of its files matches, rather
#      than summing matches across files and patterns.
#   3. Corpus. Only author-written analysis code is read. `.ado` files are
#      downloaded third-party Stata packages whose help text contains the very
#      syntax being searched for; including them inflates the count with
#      matches from library documentation. This mirrors the file filter used by
#      the LLM pipeline so the two run on exactly the same corpus.
#
# Two variants are reported:
#   as_published  every pattern from the original script
#   syntax_only   only patterns that identify a model actually being FIT,
#                 dropping (a) prose terminology, which matches comments and
#                 cited titles rather than code, and (b) the fixest pattern
#                 `feols(...|`, which matches the fixed-effects separator and
#                 is not an interaction at all --- Reviewer 1's example.
#
# The gap between the two is the share of the published count that rests on
# matching English rather than statistics.

library(dplyr)
library(tidyr)
library(readr)
library(stringr)

ROOT <- "data/prior-studies"
REPFILES <- file.path(ROOT, "replication-files")
CODE_EXT <- c("r", "do", "py", "stan", "jags")


# ---------------------------------------------------------------------------
# patterns, exactly as published
# ---------------------------------------------------------------------------

INTERACTION_PATTERNS <- list(
  r_lm_interaction = "(lm|glm|felm|feols|plm|ivreg)\\s*\\([^)]*~[^)]*[*:]",
  r_fixest_interaction = "feols?\\s*\\([^)]*\\|",
  stata_interaction_hash = "[ci]\\.\\w+#",
  stata_interaction_double = "\\w+##\\w+",
  margins_dydx = "dydx\\s*\\(",
  margins_at = "margins\\s*,.*at\\s*\\(",
  marginaleffects_r = "marginaleffects|slopes\\s*\\(|comparisons\\s*\\(",
  heterogeneous_effect = "heterogen(eous|eity)\\s*(effect|treatment)?",
  conditional_effect = "conditional\\s*(average)?\\s*(treatment)?\\s*effect",
  subgroup_analysis = "subgroup\\s*analysis|stratified\\s*analysis",
  moderation = "moderat(ion|ing|or)\\s*(effect|analysis|variable)?",
  interaction_effect = "interaction\\s*(effect|term|model)",
  interaction_plot = "interplot|interact_plot|marginsplot|coefplot.*#"
)

ML_PATTERNS <- list(
  causal_forest = "causal_forest\\s*\\(",
  grf_package = "library\\s*\\(\\s*grf\\s*\\)|require\\s*\\(\\s*grf\\s*\\)",
  grf_cate = "average_treatment_effect|predict.*causal_forest",
  bartcause = "bartc\\s*\\(|library\\s*\\(\\s*bartCause\\s*\\)",
  bcf = "bcf\\s*\\(|library\\s*\\(\\s*bcf\\s*\\)",
  cate_explicit = "\\bCATE\\b|conditional\\s+average\\s+treatment\\s+effect",
  ite_explicit = "\\bITE\\b|individual(ized)?\\s+treatment\\s+effect",
  heterogeneous_ml = "(causal|treatment)\\s*(forest|tree|learning)",
  double_ml = "DoubleML|double.*machine.*learning|dml_",
  meta_learner = "[STXR][-_]?learner|metalearner"
)

HIER_PATTERNS <- list(
  lme4 = "lmer\\s*\\(|glmer\\s*\\(",
  brms = "brm\\s*\\(|brms::",
  rstanarm = "stan_lmer|stan_glmer|stan_glm\\s*\\(",
  varying_slopes_r = "\\([^|]+\\+[^|]+\\|",
  varying_slopes_explicit = "\\|\\|",
  stan_model = "\\.stan$|\\.jags$|stan_code|stan_model|rstan::",
  stata_mixed = "mixed\\s+|xtmixed|meglm|melogit|xtmelogit",
  stata_random_slope = "\\|\\|\\s*\\w+:",
  partial_pooling = "partial\\s*pool|shrinkage|borrow.*strength",
  varying_effects = "varying\\s*(slope|effect|coefficient|intercept)",
  random_effects = "random\\s*(slope|effect|coefficient)",
  multilevel = "multilevel|multi-level|hierarchical\\s*(model|linear|regression)",
  bayesian_het = "posterior.*effect|effect.*posterior|credible.*interval"
)

# Patterns that match prose rather than a fitted model, plus the fixest
# separator. These are what `syntax_only` drops.
PROSE_OR_WRONG <- c(
  "r_fixest_interaction",
  "heterogeneous_effect", "conditional_effect", "subgroup_analysis",
  "moderation", "interaction_effect",
  "cate_explicit", "ite_explicit", "heterogeneous_ml", "double_ml",
  "partial_pooling", "varying_effects", "random_effects", "multilevel",
  "bayesian_het"
)


# ---------------------------------------------------------------------------
# scan
# ---------------------------------------------------------------------------

files <- read_csv(file.path(ROOT, "code_files_list.csv"),
                  show_col_types = FALSE) %>%
  mutate(
    safe_doi = str_replace_all(dataset_doi, "[:/]", "_"),
    ext = tolower(str_extract(filename, "(?<=\\.)[A-Za-z0-9]+$")),
    path = file.path(REPFILES, safe_doi, filename)
  ) %>%
  filter(ext %in% CODE_EXT, file.exists(path))

message("scanning ", nrow(files), " files in ",
        n_distinct(files$dataset_doi), " papers")

hit <- function(pat, txt) {
  m <- gregexpr(pat, txt, ignore.case = TRUE, perl = TRUE)[[1]]
  length(m) > 0 && m[1] != -1
}

all_pats <- c(INTERACTION_PATTERNS, ML_PATTERNS, HIER_PATTERNS)
pat_group <- c(rep("interaction", length(INTERACTION_PATTERNS)),
               rep("ml", length(ML_PATTERNS)),
               rep("hier", length(HIER_PATTERNS)))
names(pat_group) <- names(all_pats)

t0 <- Sys.time()
res <- vapply(seq_len(nrow(files)), function(i) {
  txt <- tryCatch(
    paste(readLines(files$path[i], warn = FALSE, encoding = "UTF-8"),
          collapse = "\n"),
    error = function(e) ""
  )
  vapply(all_pats, hit, logical(1), txt = txt)
}, logical(length(all_pats)))
message("scan elapsed: ", format(Sys.time() - t0))

file_hits <- as.data.frame(t(res))
names(file_hits) <- names(all_pats)
file_hits$dataset_doi <- files$dataset_doi


# ---------------------------------------------------------------------------
# aggregate to paper: any file matching counts the paper once
# ---------------------------------------------------------------------------

summarise_variant <- function(keep_pats, label) {
  ip <- intersect(names(INTERACTION_PATTERNS), keep_pats)
  mp <- intersect(names(ML_PATTERNS), keep_pats)
  hp <- intersect(names(HIER_PATTERNS), keep_pats)

  file_hits %>%
    group_by(dataset_doi) %>%
    summarize(
      interactions = any(across(all_of(ip)) %>% rowSums() > 0),
      ml_het = any(across(all_of(mp)) %>% rowSums() > 0),
      hierarchical = any(across(all_of(hp)) %>% rowSums() > 0),
      .groups = "drop"
    ) %>%
    mutate(variant = label,
           # any of the three, counted once --- the original summed these,
           # so a paper using both an interaction and a mixed model was
           # counted twice
           any_het = interactions | ml_het | hierarchical)
}

keep_all <- names(all_pats)
keep_syntax <- setdiff(names(all_pats), PROSE_OR_WRONG)

baseline <- bind_rows(
  summarise_variant(keep_all, "as_published"),
  summarise_variant(keep_syntax, "syntax_only")
)

# article year, not deposit date
meta <- read_csv(file.path(ROOT, "dataverse_datasets.csv"),
                 show_col_types = FALSE) %>%
  transmute(dataset_doi = global_id,
            year = search_year,
            deposit_year = as.numeric(substr(published_at, 1, 4)))

baseline <- baseline %>% left_join(meta, by = "dataset_doi")
write_csv(baseline, file.path(ROOT, "regex-baseline.csv"))


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

cat("\n=========== Papers flagged, by variant ===========\n")
print(as.data.frame(
  baseline %>%
    group_by(variant) %>%
    summarize(papers = n(),
              interactions = sum(interactions), ml = sum(ml_het),
              hierarchical = sum(hierarchical), any = sum(any_het),
              share_any = sprintf("%.1f%%", 100 * mean(any_het)))
), row.names = FALSE)

cat("\n=== How much of the published count rests on prose, not code? ===\n")
w <- baseline %>%
  select(dataset_doi, variant, any_het) %>%
  pivot_wider(names_from = variant, values_from = any_het)
cat(sprintf("  flagged by published patterns : %d\n", sum(w$as_published)))
cat(sprintf("  flagged by syntax patterns    : %d\n", sum(w$syntax_only)))
cat(sprintf("  flagged ONLY by prose/fixest  : %d (%.1f%% of the published count)\n",
            sum(w$as_published & !w$syntax_only),
            100 * sum(w$as_published & !w$syntax_only) / sum(w$as_published)))

cat("\n=== Year variable: deposit date vs article year ===\n")
yr <- baseline %>%
  filter(variant == "as_published") %>%
  select(year, deposit_year, any_het)
cmp <- full_join(
  yr %>% count(year, name = "papers_by_article_year"),
  yr %>% count(deposit_year, name = "papers_by_deposit") %>%
    rename(year = deposit_year),
  by = "year"
) %>% arrange(year)
print(as.data.frame(cmp), row.names = FALSE)

cat("\n=== Share of papers estimating heterogeneous effects, by year ===\n")
print(as.data.frame(
  baseline %>%
    group_by(variant, year) %>%
    summarize(papers = n(), flagged = sum(any_het),
              share = sprintf("%.1f%%", 100 * mean(any_het)), .groups = "drop") %>%
    pivot_wider(names_from = variant, values_from = c(papers, flagged, share)) %>%
    select(year, papers = papers_as_published,
           flagged_as_published, share_as_published,
           flagged_syntax_only, share_syntax_only) %>%
    arrange(year)
), row.names = FALSE)

cat("\nWrote", file.path(ROOT, "regex-baseline.csv"), "\n")
