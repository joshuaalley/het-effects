# Joshua Alley
# Analyze replication files to detect heterogeneous effects methods
#
# Approach:
# 1. Use dataverse package to find replication archives for APSR/AJPS/JOP
# 2. Download R/Stata/Python code files
# 3. Parse code to detect: interactions, ML heterogeneity, hierarchical models
# 4. Future: use LLM API for more sophisticated code parsing


# =============================================================================
# 1. CONFIGURATION
# =============================================================================

# Harvard Dataverse (primary repository for polisci)
DATAVERSE_SERVER <- "dataverse.harvard.edu"

# API Key: Set DATAVERSE_KEY environment variable before running
# Options:
#   1. Sys.setenv(DATAVERSE_KEY = "your-key-here") in R
#   2. Add to .Renviron file: DATAVERSE_KEY=your-key-here
#   3. Set in system environment variables
if (Sys.getenv("DATAVERSE_KEY") == "") {
  warning("DATAVERSE_KEY environment variable not set. Some files may be inaccessible.")
} else {
  message("Using Dataverse API key: ", substr(Sys.getenv("DATAVERSE_KEY"), 1, 8), "...")
}

# Target journals and their Dataverse collections (if they have dedicated ones)
# Many journals require deposit in Harvard Dataverse
JOURNAL_SEARCH_TERMS <- c(
  "American Political Science Review",
  "American Journal of Political Science",
  "Journal of Politics"
)

# Code file extensions to download
CODE_EXTENSIONS <- c(".R", ".r", ".do", ".ado", ".py", ".stan", ".jags")

# Output directories
OUTPUT_DIR <- "data/prior-studies"  # For CSV results
CODE_DIR <- "data/prior-studies/replication-files"  # For downloaded code

# Create directories if needed
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(CODE_DIR)) dir.create(CODE_DIR, recursive = TRUE)


# =============================================================================
# 2. SEARCH DATAVERSE FOR REPLICATION ARCHIVES
# =============================================================================

search_journal_replications <- function(journal_term, start_year = 2010,
                                         end_year = 2024, per_page = 100) {
  #' Search Dataverse for replication files from a journal
  #'
  #' @param journal_term Search term for journal
 #' @param start_year Start year
  #' @param end_year End year
  #' @param per_page Results per page
  #' @return Data frame of datasets

  message("Searching Dataverse for: ", journal_term)

  all_results <- list()

  # Search by year to avoid hitting limits
 for (year in start_year:end_year) {

    query <- paste0(
      '"', journal_term, '"',
      ' AND publicationDate:', year, '*'
    )

    tryCatch({
      results <- dataverse_search(
        query,
        type = "dataset",
        per_page = per_page,
        server = DATAVERSE_SERVER
      )

      if (length(results) > 0 && nrow(results) > 0) {
        results$search_year <- year
        results$search_term <- journal_term
        all_results[[length(all_results) + 1]] <- results
        message("  ", year, ": ", nrow(results), " datasets")
      }

      Sys.sleep(5)  # Rate limiting - 5 sec gap to avoid 403 blocking from Dataverse ISP

    }, error = function(e) {
      message("  ", year, ": Error - ", e$message)
    })
  }

  if (length(all_results) == 0) {
    return(tibble())
  }

  combined <- bind_rows(all_results)

  # Remove duplicates
  combined <- combined %>%
    distinct(global_id, .keep_all = TRUE)

  message("  Total unique datasets: ", nrow(combined))

  return(combined)
}


search_all_journals <- function(start_year = 2010, end_year = 2024) {
  #' Search for replications from all target journals
  #'
  #' @return Combined data frame

  all_datasets <- list()

  for (term in JOURNAL_SEARCH_TERMS) {
    results <- search_journal_replications(term, start_year, end_year)
    if (nrow(results) > 0) {
      all_datasets[[length(all_datasets) + 1]] <- results
    }
  }

  combined <- bind_rows(all_datasets) %>%
    distinct(global_id, .keep_all = TRUE)

  message("\nTotal unique datasets across all journals: ", nrow(combined))

  return(combined)
}


# =============================================================================
# 3. GET DATASET FILE LISTS
# =============================================================================

# Debug function - test a single DOI to see response structure
debug_dataset <- function(doi) {
  #' Debug function to inspect Dataverse response
  #'
  #' @param doi Dataset DOI
  #' @return List with raw responses

  message("Testing DOI: ", doi)

  results <- list()

  # Try dataset_files()
  message("\n1. Trying dataset_files()...")
  tryCatch({
    files <- dataset_files(doi, server = DATAVERSE_SERVER)
    message("   Success! Type: ", class(files))
    message("   Length: ", length(files))
    if (length(files) > 0) {
      message("   First element names: ", paste(names(files[[1]]), collapse = ", "))
    }
    results$dataset_files <- files
  }, error = function(e) {
    message("   Error: ", e$message)
  })

  # Try get_dataset()
  message("\n2. Trying get_dataset()...")
  tryCatch({
    dataset <- get_dataset(doi, server = DATAVERSE_SERVER)
    message("   Success! Type: ", class(dataset))
    message("   Names: ", paste(names(dataset), collapse = ", "))
    if (!is.null(dataset$files)) {
      message("   Has $files: ", length(dataset$files), " items")
      message("   $files type: ", class(dataset$files))
      if (length(dataset$files) > 0) {
        message("   First file names: ", paste(names(dataset$files[[1]]), collapse = ", "))
        if (!is.null(dataset$files[[1]]$dataFile)) {
          message("   dataFile names: ", paste(names(dataset$files[[1]]$dataFile), collapse = ", "))
        }
      }
    }
    results$get_dataset <- dataset
  }, error = function(e) {
    message("   Error: ", e$message)
  })

  return(results)
}


get_dataset_files <- function(doi) {
  #' Get list of files in a Dataverse dataset
  #'
  #' @param doi Dataset DOI (global_id)
  #' @return Data frame of files

  tryCatch({
    # Use get_dataset() which returns structured data
    dataset <- get_dataset(doi, server = DATAVERSE_SERVER)

    # Extract files - this is already a data.frame
    files <- dataset$files

    if (is.null(files) || nrow(files) == 0) {
      return(tibble())
    }

    # Select and rename columns we need
    file_df <- tibble(
      file_id = files$id,
      filename = files$filename,
      content_type = files$contentType,
      file_size = files$filesize,
      dataset_doi = doi
    )

    # Remove rows where we couldn't get filename
    file_df <- file_df %>% filter(!is.na(filename))

    return(file_df)

  }, error = function(e) {
    message("  Error getting files for ", doi, ": ", e$message)
    return(tibble())
  })
}


get_code_files <- function(datasets, max_datasets = NULL) {
  #' Get code file information for multiple datasets
  #'
  #' @param datasets Data frame with global_id column
  #' @param max_datasets Maximum datasets to process (NULL for all)
  #' @return Data frame of code files

  dois <- datasets$global_id

  if (!is.null(max_datasets)) {
    dois <- head(dois, max_datasets)
  }

  message("Getting file lists for ", length(dois), " datasets...")

  all_files <- list()

  for (i in seq_along(dois)) {
    files <- get_dataset_files(dois[i])

    if (nrow(files) > 0) {
      # Filter to code files
      code_files <- files %>%
        filter(
          str_detect(tolower(filename),
                     paste0("(", paste(CODE_EXTENSIONS, collapse = "|"), ")$"))
        )

      if (nrow(code_files) > 0) {
        all_files[[length(all_files) + 1]] <- code_files
      }
    }

    if (i %% 50 == 0) {
      message("  Processed ", i, " of ", length(dois), " datasets")
    }

    Sys.sleep(5)  # Rate limiting - 5 sec gap to avoid 403 blocking from Dataverse ISP
  }

  combined <- bind_rows(all_files)

  message("Found ", nrow(combined), " code files across ",
          n_distinct(combined$dataset_doi), " datasets")

  return(combined)
}


# =============================================================================
# 4. DOWNLOAD CODE FILES
# =============================================================================

download_code_file <- function(file_id, filename, dataset_doi, output_dir = CODE_DIR) {
  #' Download a single code file
  #'
  #' @param file_id Dataverse file ID
  #' @param filename Original filename
  #' @param dataset_doi Dataset DOI for organizing
  #' @param output_dir Output directory for code files
  #' @return Path to downloaded file or NA

  # Create safe directory name from DOI
  safe_doi <- gsub("[:/]", "_", dataset_doi)
  dataset_dir <- file.path(output_dir, safe_doi)

  if (!dir.exists(dataset_dir)) {
    dir.create(dataset_dir, recursive = TRUE)
  }

  output_path <- file.path(dataset_dir, filename)

  # Skip if already downloaded
  if (file.exists(output_path)) {
    return(output_path)
  }

  tryCatch({
    content <- get_file(file_id, server = DATAVERSE_SERVER)
    writeBin(content, output_path)
    return(output_path)

  }, error = function(e) {
    message("  Error downloading ", filename, ": ", e$message)
    return(NA_character_)
  })
}


download_code_files <- function(code_files, max_files = NULL) {
  #' Download multiple code files
  #'
  #' @param code_files Data frame from get_code_files
  #' @param max_files Maximum files to download (NULL for all)
  #' @return Data frame with download paths

  if (!dir.exists(CODE_DIR)) {
    dir.create(CODE_DIR, recursive = TRUE)
  }

  if (!is.null(max_files)) {
    code_files <- head(code_files, max_files)
  }

  message("Downloading ", nrow(code_files), " code files...")

  code_files$local_path <- NA_character_

  for (i in seq_len(nrow(code_files))) {
    code_files$local_path[i] <- download_code_file(
      code_files$file_id[i],
      code_files$filename[i],
      code_files$dataset_doi[i]
    )

    if (i %% 20 == 0) {
      message("  Downloaded ", i, " of ", nrow(code_files), " files")
    }

    Sys.sleep(5)  # Rate limiting - 5 sec gap to avoid 403 blocking from Dataverse ISP
  }

  n_success <- sum(!is.na(code_files$local_path))
  message("Successfully downloaded ", n_success, " files")

  return(code_files)
}


# =============================================================================
# 5. PARSE CODE FILES FOR METHODS
# =============================================================================

# Pattern definitions for detecting methods
# NOTE: These are designed to be specific to heterogeneous effects analysis,
# not just any use of the method. This means some false negatives but fewer
# false positives, which is better for estimating method prevalence.

INTERACTION_PATTERNS <- list(
  # R model formulas with interactions (in regression context)
  r_lm_interaction = "(lm|glm|felm|feols|plm|ivreg)\\s*\\([^)]*~[^)]*[*:]",
  r_fixest_interaction = "feols?\\s*\\([^)]*\\|",  # fixest with FE often has interactions

  # Stata interactions (factor variable notation)
  stata_interaction_hash = "[ci]\\.\\w+#",  # c.x#c.y or i.x#i.y notation
  stata_interaction_double = "\\w+##\\w+",  # Full factorial notation x##y

  # Marginal effects (strong indicator of interest in conditional effects)
  margins_dydx = "dydx\\s*\\(",
  margins_at = "margins\\s*,.*at\\s*\\(",
  marginaleffects_r = "marginaleffects|slopes\\s*\\(|comparisons\\s*\\(",

  # Explicit heterogeneity/moderation terminology
  heterogeneous_effect = "heterogen(eous|eity)\\s*(effect|treatment)?",
  conditional_effect = "conditional\\s*(average)?\\s*(treatment)?\\s*effect",
  subgroup_analysis = "subgroup\\s*analysis|stratified\\s*analysis",
  moderation = "moderat(ion|ing|or)\\s*(effect|analysis|variable)?",
  interaction_effect = "interaction\\s*(effect|term|model)",

  # Interaction plots (indicates substantive interest in interactions)
  interaction_plot = "interplot|interact_plot|marginsplot|coefplot.*#"
)

ML_HETEROGENEITY_PATTERNS <- list(
  # Causal forests - very specific to CATE
  causal_forest = "causal_forest\\s*\\(",
  grf_package = "library\\s*\\(\\s*grf\\s*\\)|require\\s*\\(\\s*grf\\s*\\)",
  grf_cate = "average_treatment_effect|predict.*causal_forest",

  # BART for causal inference specifically
  bartcause = "bartc\\s*\\(|library\\s*\\(\\s*bartCause\\s*\\)",
  bcf = "bcf\\s*\\(|library\\s*\\(\\s*bcf\\s*\\)",  # Bayesian Causal Forest

  # CATE-specific terminology
  cate_explicit = "\\bCATE\\b|conditional\\s+average\\s+treatment\\s+effect",
  ite_explicit = "\\bITE\\b|individual(ized)?\\s+treatment\\s+effect",
  heterogeneous_ml = "(causal|treatment)\\s*(forest|tree|learning)",

  # Double ML for heterogeneity
  double_ml = "DoubleML|double.*machine.*learning|dml_",

  # Meta-learners
  meta_learner = "[STXR][-_]?learner|metalearner"
)

HIERARCHICAL_PATTERNS <- list(
  # R mixed effects packages
  lme4 = "lmer\\s*\\(|glmer\\s*\\(",
  brms = "brm\\s*\\(|brms::",
  rstanarm = "stan_lmer|stan_glmer|stan_glm\\s*\\(",


  # Random/varying slopes syntax in R (the key indicator beyond random intercepts)
  varying_slopes_r = "\\([^|]+\\+[^|]+\\|",  # (1 + x | group) pattern
  varying_slopes_explicit = "\\|\\|",  # brms || notation for uncorrelated random effects

  # Stan/JAGS (Bayesian modeling)
  stan_model = "\\.stan$|\\.jags$|stan_code|stan_model|rstan::",

  # Stata mixed models
  stata_mixed = "mixed\\s+|xtmixed|meglm|melogit|xtmelogit",
  stata_random_slope = "\\|\\|\\s*\\w+:",  # Random slope syntax in Stata mixed

  # Explicit terminology
  partial_pooling = "partial\\s*pool|shrinkage|borrow.*strength",
  varying_effects = "varying\\s*(slope|effect|coefficient|intercept)",
  random_effects = "random\\s*(slope|effect|coefficient)",
  multilevel = "multilevel|multi-level|hierarchical\\s*(model|linear|regression)",

  # Bayesian heterogeneity
  bayesian_het = "posterior.*effect|effect.*posterior|credible.*interval"
)


parse_code_file <- function(file_path, verbose = FALSE) {
  #' Parse a code file for method patterns
  #'
  #' @param file_path Path to code file
  #' @param verbose If TRUE, return details of which patterns matched
  #' @return List of detected methods

  if (is.na(file_path) || !file.exists(file_path)) {
    return(list(interactions = 0, ml_het = 0, hierarchical = 0, error = TRUE))
  }

  tryCatch({
    # Read file content
    content <- paste(readLines(file_path, warn = FALSE), collapse = "\n")
    content_lower <- tolower(content)

    # Count pattern matches and optionally return which matched
    count_patterns <- function(patterns, text, return_names = FALSE) {
      matches <- sapply(patterns, function(p) {
        length(gregexpr(p, text, ignore.case = TRUE, perl = TRUE)[[1]]) > 0 &&
          gregexpr(p, text, ignore.case = TRUE, perl = TRUE)[[1]][1] != -1
      })

      if (return_names) {
        return(names(patterns)[matches])
      }
      sum(matches)  # Count unique patterns matched
    }

    result <- list(
      interactions = count_patterns(INTERACTION_PATTERNS, content),
      ml_het = count_patterns(ML_HETEROGENEITY_PATTERNS, content_lower),
      hierarchical = count_patterns(HIERARCHICAL_PATTERNS, content_lower),
      n_lines = length(readLines(file_path, warn = FALSE)),
      error = FALSE
    )

    if (verbose) {
      result$interaction_patterns <- count_patterns(INTERACTION_PATTERNS, content, TRUE)
      result$ml_patterns <- count_patterns(ML_HETEROGENEITY_PATTERNS, content_lower, TRUE)
      result$hierarchical_patterns <- count_patterns(HIERARCHICAL_PATTERNS, content_lower, TRUE)
    }

    return(result)

  }, error = function(e) {
    list(interactions = 0, ml_het = 0, hierarchical = 0, error = TRUE)
  })
}


parse_all_code_files <- function(code_files) {
  #' Parse all downloaded code files
  #'
  #' @param code_files Data frame with local_path column
  #' @return Data frame with method counts

  message("Parsing ", sum(!is.na(code_files$local_path)), " code files...")

  results <- map_df(seq_len(nrow(code_files)), function(i) {
    parsed <- parse_code_file(code_files$local_path[i])

    tibble(
      file_id = code_files$file_id[i],
      filename = code_files$filename[i],
      dataset_doi = code_files$dataset_doi[i],
      interactions = parsed$interactions,
      ml_het = parsed$ml_het,
      hierarchical = parsed$hierarchical,
      n_lines = parsed$n_lines %||% NA_integer_,
      parse_error = parsed$error
    )
  })

  # Summary
  message("\nMethod detection summary:")
  message("  Files with interactions: ",
          sum(results$interactions > 0, na.rm = TRUE))
  message("  Files with ML heterogeneity: ",
          sum(results$ml_het > 0, na.rm = TRUE))
  message("  Files with hierarchical models: ",
          sum(results$hierarchical > 0, na.rm = TRUE))

  return(results)
}


# =============================================================================
# 6. AGGREGATE TO DATASET/ARTICLE LEVEL
# =============================================================================

aggregate_to_dataset <- function(parsed_files, datasets) {
  #' Aggregate file-level results to dataset level
  #'
  #' @param parsed_files Results from parse_all_code_files
  #' @param datasets Original dataset metadata
  #' @return Data frame with one row per dataset

  dataset_summary <- parsed_files %>%
    group_by(dataset_doi) %>%
    summarise(
      n_code_files = n(),
      total_lines = sum(n_lines, na.rm = TRUE),
      has_interactions = any(interactions > 0, na.rm = TRUE),
      has_ml_het = any(ml_het > 0, na.rm = TRUE),
      has_hierarchical = any(hierarchical > 0, na.rm = TRUE),
      max_interactions = max(interactions, na.rm = TRUE),
      max_ml_het = max(ml_het, na.rm = TRUE),
      max_hierarchical = max(hierarchical, na.rm = TRUE),
      .groups = "drop"
    )

  # Join with dataset metadata
  dataset_summary <- dataset_summary %>%
    left_join(
      datasets %>% select(global_id, name, description, published_at),
      by = c("dataset_doi" = "global_id")
    )

  return(dataset_summary)
}


# =============================================================================
# 7. MAIN WORKFLOW
# =============================================================================

run_replication_analysis <- function(start_year = 2015, end_year = 2024,
                                      max_datasets = 100, max_files = 500) {
  #' Complete workflow for replication file analysis
  #'
  #' @param start_year Start year for search
  #' @param end_year End year
  #' @param max_datasets Max datasets to process
  #' @param max_files Max files to download
  #' @return List with all results

  message("=== Replication File Analysis ===\n")

  # 1. Search for datasets
  datasets <- search_all_journals(start_year, end_year)

  if (nrow(datasets) == 0) {
    stop("No datasets found")
  }

  # Save dataset list
  write_csv(datasets, file.path(OUTPUT_DIR, "dataverse_datasets.csv"))

  # 2. Get code file information
  code_files <- get_code_files(datasets, max_datasets = max_datasets)

  if (nrow(code_files) == 0) {
    stop("No code files found")
  }

  # Save code file list
  write_csv(code_files, file.path(OUTPUT_DIR, "code_files_list.csv"))

  # 3. Download code files
  code_files <- download_code_files(code_files, max_files = max_files)

  # 4. Parse code files
  parsed <- parse_all_code_files(code_files)

  # Save parsed results
  write_csv(parsed, file.path(OUTPUT_DIR, "parsed_code_files.csv"))

  # 5. Aggregate to dataset level
  dataset_summary <- aggregate_to_dataset(parsed, datasets)

  # Save dataset summary
  write_csv(dataset_summary, file.path(OUTPUT_DIR, "dataset_method_summary.csv"))

  # 6. Overall summary
  message("\n=== Overall Summary ===")
  message("Datasets searched: ", nrow(datasets))
  message("Datasets with code: ", n_distinct(code_files$dataset_doi))
  message("Code files analyzed: ", nrow(parsed))
  message("\nMethod prevalence (dataset level):")
  message("  Interactions: ", sum(dataset_summary$has_interactions, na.rm = TRUE),
          " (", round(100 * mean(dataset_summary$has_interactions, na.rm = TRUE), 1), "%)")
  message("  ML heterogeneity: ", sum(dataset_summary$has_ml_het, na.rm = TRUE),
          " (", round(100 * mean(dataset_summary$has_ml_het, na.rm = TRUE), 1), "%)")
  message("  Hierarchical: ", sum(dataset_summary$has_hierarchical, na.rm = TRUE),
          " (", round(100 * mean(dataset_summary$has_hierarchical, na.rm = TRUE), 1), "%)")

  return(list(
    datasets = datasets,
    code_files = code_files,
    parsed = parsed,
    dataset_summary = dataset_summary
  ))
}


### ------------------------------------- ### 
# Full analysis
# results <- run_replication_analysis(
#   start_year = 2010,
#   end_year = 2024,
#   max_datasets = NULL,
#   max_files = NULL
# )

# get results with all files downloaded

  # Load the saved code file list          
 code_files <- read_csv(file.path("data/prior-studies", "code_files_list.csv"))                                                                                                                    
  # Reconstruct local_path from dataset_doi and filename
  code_files <- code_files %>%
    mutate(
      safe_doi = gsub("[:/]", "_", dataset_doi),
      local_path = file.path("data/prior-studies/replication-files", safe_doi, filename)
    ) %>%
    # Only keep files that actually exist
    filter(file.exists(local_path))

  # Now parse
  parsed <- parse_all_code_files(code_files)

  # Load datasets too
  datasets <- read_csv(file.path(OUTPUT_DIR, "dataverse_datasets.csv"))

  # Aggregate
  dataset_summary <- aggregate_to_dataset(parsed, datasets)

  # Build results list
  results <- list(
    datasets = datasets,
    code_files = code_files,
    parsed = parsed,
    dataset_summary = dataset_summary
  )

# results
rep_results <- results[["dataset_summary"]] %>%
                     mutate(
                      year = as.numeric(substr(published_at, 1, 4))
                     )
glimpse(rep_results)

# summarize
rep_results_sum <- rep_results %>% 
                     group_by(year) %>%
                     filter(year <= 2024) %>%
                     summarize(
                        num_papers = n(),
                        used_interactions = sum(has_interactions),
                        used_ml = sum(has_ml_het),
                        total_interactions = sum(max_interactions),
                        share_interactions = used_interactions / num_papers,
                        inter_per_paper = total_interactions / used_interactions,
                        total_ml = sum(max_ml_het),
                        share_ml = used_ml / num_papers,
                        share_het = share_interactions + share_ml
                     ) 


# plot results
ggplot(rep_results_sum, aes(x = year, y = used_interactions)) +
  geom_line()

ggplot(rep_results_sum, aes(x = year, y = share_interactions)) +
  geom_line()

ggplot(rep_results_sum, aes(x = year, y = share_het)) +
  geom_line()

ggplot(rep_results_sum, aes(x = year, y = inter_per_paper)) +
  geom_line()

# paper prevalence
results_prev_long <- rep_results_sum %>%
                      select(year, num_papers, used_interactions, used_ml) %>%
                      pivot_longer(
                        cols = !year,
                        names_to = "papers"
                      ) %>%
                      mutate(
                        papers = case_when(
                          str_detect(papers, "num") ~ "Total",
                          str_detect(papers, "inter") ~ "Interactions",
                          str_detect(papers, "ml") ~ "ML",
                        )
                      )


ggplot(results_prev_long, aes(x = year, y = value,
             color = papers)) +
  geom_line() +
  scale_color_grey(
    start = .1, end = .6
  ) +
  labs(
    title = "Prevalence of Heterogeneity Techniques",
    subtitle = "APSR, AJPS and JOP- 2010 to 2024",
    x = "Year",
    y = "Number of Papers",
    color = ""
  ) +
  theme(legend.position = "bottom")
ggsave("figures/het-effects-prev.png", height = 6, width = 8)


# paper share
results_share_long <- rep_results_sum %>%
                      select(year, share_interactions, share_ml) %>%
                      pivot_longer(
                        cols = !year,
                        names_to = "papers"
                      ) %>%
                      mutate(
                        papers = case_when(
                          str_detect(papers, "inter") ~ "Interactions",
                          str_detect(papers, "ml") ~ "ML",
                        )
                      )

ggplot(results_share_long, aes(x = year, y = value,
             color = papers)) +
  geom_line() +
  scale_color_grey(
    start = .1, end = .6
  ) +
  labs(
    title = "Share of Papers with Different Heterogeneity Techniques",
    subtitle = "APSR, AJPS and JOP: 2010 to 2024",
    x = "Year",
    y = "Share of Papers",
    color = ""
  ) +
  theme(legend.position = "bottom")
