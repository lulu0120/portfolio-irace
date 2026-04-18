#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

resolve_qap_hydra_dir <- function() {
  find_qap_hydra_dir <- function(start) {
    current <- normalizePath(start, winslash = "/", mustWork = TRUE)
    if (!dir.exists(current)) current <- dirname(current)
    repeat {
      candidate <- file.path(current, "QAP-hydra")
      if (dir.exists(candidate)) {
        return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
      }
      if (basename(current) == "QAP-hydra") {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    NA_character_
  }

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    path <- if (!is.null(ctx)) ctx$path else ""
    if (is.character(path) && nzchar(path) && file.exists(path)) {
      resolved <- find_qap_hydra_dir(dirname(path))
      if (is.character(resolved) && !is.na(resolved)) return(resolved)
    }
  }
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    resolved <- find_qap_hydra_dir(dirname(sub("^--file=", "", file_flag[[1L]])))
    if (is.character(resolved) && !is.na(resolved)) return(resolved)
  }
  resolved <- find_qap_hydra_dir(getwd())
  if (is.character(resolved) && !is.na(resolved)) return(resolved)
  stop("Cannot resolve QAP-hydra directory. In RStudio, set working directory inside the portfolio-irace repo.")
}

script_dir <- resolve_qap_hydra_dir()
source(file.path(script_dir, "qap_hydra_interface.R"), local = .GlobalEnv)

get_arg <- function(args, flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

parse_int_vec <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(as.integer(default))
  as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
}

run_qap_hydra_from_args <- function(args = commandArgs(trailingOnly = TRUE), save_outputs = TRUE) {
  split_seeds <- parse_int_vec(get_arg(args, "--split-seeds", "11"), 11L)
  budgets <- parse_int_vec(get_arg(args, "--budgets", "1000"), 1000L)
  target_portfolio_sizes <- parse_int_vec(get_arg(args, "--target-portfolio-sizes", "2"), 2L)
  algo_seeds <- parse_int_vec(get_arg(args, "--algo-seeds", "1"), 1L)
  best_k <- as.integer(get_arg(args, "--best-k", "1"))
  race_seed_offset <- as.integer(get_arg(args, "--race-seed-offset", "700"))
  time_limit <- as.numeric(get_arg(args, "--time-limit", "5"))
  T_first <- as.integer(get_arg(args, "--T-first", "5"))
  T_each <- as.integer(get_arg(args, "--T-each", "1"))
  N_min <- as.integer(get_arg(args, "--N-min", "3"))
  name <- get_arg(args, "--name", NULL)
  out_prefix <- get_arg(args, "--out-prefix", "qap_hydra_run")
  quiet <- "--quiet" %in% args

  run_qap_hydra_overall_rstudio(
    split_seeds = split_seeds,
    budgets = budgets,
    target_portfolio_sizes = target_portfolio_sizes,
    algo_seeds = algo_seeds,
    best_k = best_k,
    race_seed_offset = race_seed_offset,
    time_limit = time_limit,
    T_first = T_first,
    T_each = T_each,
    N_min = N_min,
    verbose = !isTRUE(quiet),
    name = name,
    out_prefix = out_prefix,
    save_outputs = isTRUE(save_outputs)
  )
}

main <- function() {
  run_qap_hydra_from_args()
}

if (sys.nframe() == 0L) {
  res <- main()
} else {
  cat("Loaded run_qap_hydra_overall.R\n")
  cat("Use either:\n")
  cat("  run_qap_hydra_overall_rstudio(...)\n")
  cat("or:\n")
  cat("  run_qap_hydra_from_args(c('--split-seeds','11','--budgets','100,200,300','--target-portfolio-sizes','2','--algo-seeds','2,3,4,5,6','--name','qap_hydra_irace'))\n")
}
