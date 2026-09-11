# Joshua Alley
# Figure 1: how often political scientists estimate heterogeneous effects, and
# how many they estimate at once.
#
# Replaces llm-aggregate.R and llm-interaction-counts.R, which did these two
# halves separately. Assumes data/setup-script.R has been sourced.
#
# Inputs, both written by llm-classify.py and archived in full:
#   llm-output/gpt-oss.classify-v2.jsonl   one judgment per file-chunk
#   llm-output/gpt-oss.extract-v1.jsonl    interaction inventory, run only over
#                                          the chunks the classifier flagged
#
# Panel A counts papers: does this paper estimate an effect that varies?
# Panel B counts interactions within those papers: how many distinct ones?
# The first is a claim about breadth, the second about proliferation, and they
# answer different halves of the motivating argument.

ROOT <- "data/prior-studies"
OUT <- "figures/het-effects-prev.png"

# A run is a model paired with a prompt version; archives are keyed on both, so
# revisions sit beside their predecessors rather than replacing them. The
# earlier runs are read only to report reliability, not to produce the figure.
HEADLINE <- c(tag = "gpt-oss", prompt = "classify-v2")
COMPARISONS <- list(
  gpt_oss_v1 = c(tag = "gpt-oss", prompt = "classify-v1"),
  llama_v1   = c(tag = "llama",   prompt = "classify-v1")
)
YEARS <- 2012:2024


# ---------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------

# Parse line by line. A run killed mid-write leaves a truncated line, and one
# bad line should cost one record rather than the whole archive.
read_jsonl <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  ok <- vapply(lines, function(l)
    !inherits(try(fromJSON(l, simplifyVector = FALSE), silent = TRUE),
              "try-error"), logical(1), USE.NAMES = FALSE)
  if (any(!ok)) message("  ", basename(path), ": skipped ", sum(!ok),
                        " unparseable line(s)")
  raw <- stream_in(textConnection(lines[ok]), verbose = FALSE)
  raw[!sapply(raw$result, is.null), ]
}

archive <- function(run) {
  file.path(ROOT, "llm-output",
            paste0(run[["tag"]], ".", run[["prompt"]], ".jsonl"))
}

# the last successful attempt for a chunk wins; failures were retried and both
# attempts are kept in the archive for audit
latest <- function(d) {
  d %>% arrange(file_id, chunk, ts) %>%
    group_by(file_id, chunk) %>% slice_tail(n = 1) %>% ungroup()
}

read_classify <- function(run) {
  raw <- read_jsonl(archive(run))
  if (is.null(raw)) return(NULL)
  res <- raw$result
  has <- function(tag) vapply(res$methods, function(m) tag %in% unlist(m),
                              logical(1))
  tibble(
    dataset_doi = raw$dataset_doi, file_id = raw$file_id,
    filename = raw$filename, chunk = raw$chunk, ts = raw$ts,
    het = vapply(res$estimates_heterogeneous_effects, isTRUE, logical(1)),
    interaction = has("interaction_term"),
    subgroup = has("subgroup_split"),
    ml = has("ml_heterogeneity")
  ) %>% latest()
}

to_paper <- function(d) {
  d %>%
    group_by(dataset_doi) %>%
    summarize(across(c(het, interaction, subgroup, ml), any), .groups = "drop") %>%
    # the reported measure: an interaction or a machine learning method.
    # Subgroup splitting is excluded --- see the appendix; it is the least
    # reliably measured category, kappa 0.49 against 0.92 for this one.
    mutate(estimates_het = interaction | ml)
}

meta <- read_csv(file.path(ROOT, "dataverse_datasets.csv"),
                 show_col_types = FALSE) %>%
  transmute(dataset_doi = global_id, year = search_year)

message("reading classifications")
paper <- read_classify(HEADLINE) %>% to_paper() %>%
  left_join(meta, by = "dataset_doi") %>%
  filter(!is.na(year), year %in% YEARS)


# ---------------------------------------------------------------------------
# reliability: agreement with the runs this one supersedes
# ---------------------------------------------------------------------------

kappa <- function(a, b) {
  po <- mean(a == b)
  pe <- mean(a) * mean(b) + (1 - mean(a)) * (1 - mean(b))
  (po - pe) / (1 - pe)
}

cat("=== agreement with superseded runs (reliability, not used in the figure) ===\n")
for (nm in names(COMPARISONS)) {
  other <- read_classify(COMPARISONS[[nm]])
  if (is.null(other)) next
  o <- to_paper(other) %>% select(dataset_doi, other = estimates_het)
  j <- inner_join(paper, o, by = "dataset_doi")
  cat(sprintf("  %-12s this %.1f%%  that %.1f%%  agree %.1f%%  kappa %.3f  (n=%d)\n",
              nm, 100 * mean(j$estimates_het), 100 * mean(j$other),
              100 * mean(j$estimates_het == j$other),
              kappa(j$estimates_het, j$other), nrow(j)))
}


# ---------------------------------------------------------------------------
# interaction inventory, with the names verified against the source
# ---------------------------------------------------------------------------

message("reading interaction inventory")
ex <- read_jsonl(file.path(ROOT, "llm-output", "gpt-oss.extract-v1.jsonl"))

inv <- tibble(
  dataset_doi = ex$dataset_doi, file_id = ex$file_id,
  filename = ex$filename, chunk = ex$chunk, ts = ex$ts,
  pairs = lapply(ex$result$interactions,
                 function(x) if (is.null(x)) character(0) else unlist(x)),
  complete = vapply(ex$result$list_complete,
                    function(x) if (is.null(x)) TRUE else isTRUE(x), logical(1))
) %>% latest()

# Strip a transformation wrapper so log(x) and i.state verify against x, state.
core <- function(x) {
  x <- str_replace(x, "^[a-z_.]*\\(", "")
  x <- str_replace(x, "\\)+$", "")
  x <- str_replace(x, "^(c|i|ib[0-9]*|l|f|d)\\.", "")
  str_replace(x, "[^\\w.].*$", "")
}

tok_cache <- new.env(hash = TRUE)
tokens_of <- function(path) {
  hit <- tok_cache[[path]]
  if (!is.null(hit)) return(hit)
  txt <- tryCatch(paste(readLines(path, warn = FALSE), collapse = "\n"),
                  error = function(e) "")
  tk <- unique(tolower(unlist(str_extract_all(txt, "[A-Za-z_][\\w.]*"))))
  assign(path, tk, envir = tok_cache)
  tk
}

message("verifying variable names against source files")
n_raw <- n_kept <- 0L
inv$kept <- vector("list", nrow(inv))
for (i in seq_len(nrow(inv))) {
  tk <- tokens_of(file.path(ROOT, "replication-files",
                            str_replace_all(inv$dataset_doi[i], "[:/]", "_"),
                            inv$filename[i]))
  keep <- character(0)
  for (pr in inv$pairs[[i]]) {
    n_raw <- n_raw + 1L
    v <- tolower(trimws(str_split(pr, ":", n = 2)[[1]]))
    if (length(v) != 2 || v[1] == v[2]) next
    present <- vapply(v, function(x)
      x %in% tk || (nzchar(core(x)) && core(x) %in% tk), logical(1))
    if (all(present)) keep <- c(keep, paste(sort(v), collapse = ":"))
  }
  inv$kept[[i]] <- unique(keep)
  n_kept <- n_kept + length(inv$kept[[i]])
}
cat(sprintf("\ninteraction pairs returned %d; both names found in source %d (%.1f%%)\n",
            n_raw, n_kept, 100 * n_kept / max(n_raw, 1)))
cat("Pairs naming a variable absent from the file are dropped, so the counts\n")
cat("below are a floor.\n")

counts <- inv %>%
  group_by(dataset_doi) %>%
  summarize(n_pairs = length(unique(unlist(kept))),
            n_vars = length(unique(unlist(str_split(unique(unlist(kept)), ":")))),
            censored = any(!complete), .groups = "drop") %>%
  left_join(meta, by = "dataset_doi") %>%
  filter(!is.na(year), year %in% YEARS, n_pairs > 0)

write_csv(counts, file.path(ROOT, "llm-interaction-counts.csv"))


# ---------------------------------------------------------------------------
# numbers
# ---------------------------------------------------------------------------

prev <- paper %>%
  group_by(year) %>%
  summarize(papers = n(), any_het = mean(estimates_het), ml = mean(ml),
            .groups = "drop")

cat("\n=== Panel A: share of papers estimating heterogeneous effects ===\n")
print(as.data.frame(prev %>%
  transmute(year, papers,
            interaction_or_ml = sprintf("%.0f%%", 100 * any_het),
            ml_only = sprintf("%.1f%%", 100 * ml))), row.names = FALSE)

cat("\n=== Panel B: distinct interactions per paper, among papers with any ===\n")
cat("papers:", nrow(counts), " right-censored at the 30-pair cap:",
    sum(counts$censored), sprintf("(%.1f%%)\n", 100 * mean(counts$censored)))
print(round(quantile(counts$n_pairs, c(.25, .5, .75, .9, .95, 1)), 1))
cat("distinct variables entering interactions:\n")
print(round(quantile(counts$n_vars, c(.25, .5, .75, .9, .95, 1)), 1))

fit <- lm(n_pairs ~ year, data = counts)
cat(sprintf("\ntrend: %+.3f interactions per paper per year (SE %.3f)\n",
            coef(fit)[2], sqrt(diag(vcov(fit)))[2]))
cat(sprintf("median crossed: a paper with %d modifiers implies 2^%d = %d cells\n",
            median(counts$n_vars), median(counts$n_vars),
            2^median(counts$n_vars)))


# ---------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------

pa <- prev %>%
  pivot_longer(c(any_het, ml), names_to = "measure", values_to = "share") %>%
  mutate(measure = factor(measure, levels = c("any_het", "ml"),
                          labels = c("Interaction or machine learning",
                                     "Machine learning only"))) %>%
  ggplot(aes(x = year, y = share, group = measure)) +
  geom_line(aes(linetype = measure)) +
  geom_point(aes(shape = measure), size = 2.2, fill = "white") +
  scale_shape_manual(values = c(21, 24)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_x_continuous(breaks = seq(min(YEARS), max(YEARS), 2),
                     limits = range(YEARS) + c(-0.6, 0.6)) +
  labs(x = NULL, y = "Share of papers", linetype = NULL, shape = NULL) +
  theme(legend.position = "bottom")
pa

# ML het is usually
prev %>%
  pivot_longer(c(any_het, ml), names_to = "measure", values_to = "share") %>% 
  select(measure, share) %>%
  print(share)

pb <- counts %>%
  ggplot(aes(x = year, y = n_pairs, group = year)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, quantile(counts$n_pairs, .95))) +
  scale_x_continuous(breaks = seq(min(YEARS), max(YEARS), 2),
                     limits = range(YEARS) + c(-0.6, 0.6)) +
  labs(x = "Year of Publication",
       y = "Distinct interactions\nper paper")

p_inter <- pa / pb + patchwork::plot_layout(heights = c(1.15, 1)) +
  patchwork::plot_annotation(
    title = "Heterogeneous Effect Estimation in Leading Journals",
    subtitle = "2012-2024",
    caption = paste0(
      "Lower panel covers papers estimating at least one interaction and drops outliers from view"),
    theme = theme(plot.caption = element_text(hjust = 0, size = 9)))
p_inter

ggsave(OUT, p_inter, height = 8, width = 8)
cat("\nwrote", OUT, "\n")
