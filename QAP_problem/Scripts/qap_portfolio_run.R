#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

script_dir <- local({
  find_scripts_dir <- function(start) {
    if (!is.character(start) || length(start) == 0L || is.na(start) || !nzchar(start)) {
      return(NA_character_)
    }
    current <- normalizePath(start[[1L]], winslash = "/", mustWork = TRUE)
    if (!dir.exists(current)) current <- dirname(current)
    repeat {
      candidate <- file.path(current, "QAP_problem", "Scripts")
      if (dir.exists(candidate)) {
        return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
      }
      if (basename(current) == "Scripts" && basename(dirname(current)) == "QAP_problem") {
        return(normalizePath(current, winslash = "/", mustWork = TRUE))
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    NA_character_
  }

  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    resolved <- find_scripts_dir(sub("^--file=", "", file_flag[[1L]]))
    if (is.character(resolved) && !is.na(resolved)) {
      resolved
    } else {
      dirname(normalizePath(sub("^--file=", "", file_flag[[1L]]), winslash = "/", mustWork = TRUE))
    }
  } else {
    frames <- sys.frames()
    ofiles <- vapply(frames, function(env) {
      if (exists("ofile", envir = env, inherits = FALSE)) as.character(get("ofile", envir = env, inherits = FALSE)) else NA_character_
    }, character(1))
    ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
    if (length(ofiles) > 0L) {
      dirname(normalizePath(tail(ofiles, 1L), winslash = "/", mustWork = TRUE))
    } else {
      resolved <- find_scripts_dir(getwd())
      if (is.character(resolved) && !is.na(resolved)) {
        resolved
      } else {
        stop("Cannot resolve QAP_problem/Scripts from current working directory.")
      }
    }
  }
})
source(file.path(script_dir, "qap_portfolio_interface.R"), local = .GlobalEnv)
script_dir <- qap_scripts_dir()

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

train_split <- get_arg("--train-split", file.path(script_dir, "train_instances.txt"))
budget <- as.integer(get_arg("--budget", "1000"))
port_size <- as.integer(get_arg("--port-size", "2"))
best_k <- as.integer(get_arg("--best-k", "1"))
seed0 <- as.integer(get_arg("--seed", "1"))
race_seed0 <- as.integer(get_arg("--race-seed", "700"))
time_limit <- as.numeric(get_arg("--time-limit", "5"))
enable_second_level_test <- !tolower(get_arg("--enable-second-level-test", "true")) %in% c("false", "f", "0", "no")
global_free_rider_elimination <- !tolower(get_arg("--global-free-rider-elimination", "true")) %in% c("false", "f", "0", "no")
enable_duplicate <- !tolower(get_arg("--enable-duplicate", "true")) %in% c("false", "f", "0", "no")

res <- run_qap_portfolio_once(
  train_split = train_split,
  budget = budget,
  port_size = port_size,
  best_k = best_k,
  enable_second_level_test = enable_second_level_test,
  global_free_rider_elimination = global_free_rider_elimination,
  enable_duplicate = enable_duplicate,
  seed0 = seed0,
  race_seed0 = race_seed0,
  time_limit = time_limit,
  verbose = TRUE
)

cat("best_name:", res$result$best_name, "\n")
cat("best_name_decoded:", res$best_name_decoded, "\n")
cat("total_eval_used:", res$result$total_eval_used, "\n")
cat("total_instances_used:", res$result$total_instances_used, "\n")
