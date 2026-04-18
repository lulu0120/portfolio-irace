#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

resolve_qap_scripts_dir <- function() {
  if (exists("script_dir", envir = .GlobalEnv, inherits = FALSE)) {
    existing <- get("script_dir", envir = .GlobalEnv, inherits = FALSE)
    if (is.character(existing) && length(existing) > 0L && nzchar(existing[[1L]]) && dir.exists(existing[[1L]])) {
      candidate <- normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
      if (basename(candidate) == "Scripts" && basename(dirname(candidate)) == "QAP_problem") {
        return(candidate)
      }
    }
  }

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    editor_path <- if (!is.null(ctx)) ctx$path else ""
    if (is.character(editor_path) && length(editor_path) > 0L && nzchar(editor_path) && file.exists(editor_path)) {
      current <- dirname(normalizePath(editor_path[[1L]], winslash = "/", mustWork = TRUE))
      repeat {
        candidate <- file.path(current, "QAP_problem", "Scripts")
        if (dir.exists(candidate)) {
          return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
        }
        if (basename(current) == "Scripts" && basename(dirname(current)) == "QAP_problem") {
          return(current)
        }
        parent <- dirname(current)
        if (identical(parent, current)) break
        current <- parent
      }
    }
  }

  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_flag[[1L]]), winslash = "/", mustWork = TRUE)))
  }

  frames <- sys.frames()
  ofiles <- vapply(frames, function(env) {
    if (exists("ofile", envir = env, inherits = FALSE)) {
      as.character(get("ofile", envir = env, inherits = FALSE))
    } else {
      NA_character_
    }
  }, character(1))
  ofiles <- ofiles[!is.na(ofiles) & nzchar(ofiles)]
  if (length(ofiles) > 0L) {
    return(dirname(normalizePath(tail(ofiles, 1L), winslash = "/", mustWork = TRUE)))
  }

  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    candidate <- file.path(current, "QAP_problem", "Scripts")
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
    if (basename(current) == "Scripts" && basename(dirname(current)) == "QAP_problem") {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }

  stop("Cannot resolve QAP_problem/Scripts. In RStudio, set working directory inside the portfolio-irace repo.")
}

script_dir <- resolve_qap_scripts_dir()
source(file.path(resolve_qap_scripts_dir(), "qap_portfolio_interface.R"), local = .GlobalEnv)

get_arg <- function(args, flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

parse_int_vec <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(as.integer(default))
  as.integer(strsplit(x, ",", fixed = TRUE)[[1L]])
}

parse_flag <- function(args, flag) {
  flag %in% args
}

parse_toggle_flag <- function(args, enable_flag, disable_flag, default = TRUE) {
  if (enable_flag %in% args) return(TRUE)
  if (disable_flag %in% args) return(FALSE)
  isTRUE(default)
}

run_qap_runoverall_from_args <- function(args = commandArgs(trailingOnly = TRUE),
                                         save_outputs = TRUE) {
  split_seeds <- parse_int_vec(get_arg(args, "--split-seeds", "11"), 11L)
  budgets <- parse_int_vec(get_arg(args, "--budgets", "1000"), 1000L)
  port_sizes <- parse_int_vec(get_arg(args, "--port-sizes", "1,2,4"), c(1L, 2L, 4L))
  algo_seeds <- parse_int_vec(get_arg(args, "--algo-seeds", "1"), 1L)
  best_k <- as.integer(get_arg(args, "--best-k", "1"))
  enable_second_level_test <- parse_toggle_flag(args, "--second-level-test", "--no-second-level-test", TRUE)
  global_free_rider_elimination <- parse_toggle_flag(args, "--global-free-rider-elimination", "--no-global-free-rider-elimination", TRUE)
  enable_duplicate <- parse_toggle_flag(args, "--duplicate", "--no-duplicate", TRUE)
  time_limit <- as.numeric(get_arg(args, "--time-limit", "5"))
  race_seed_offset <- as.integer(get_arg(args, "--race-seed-offset", "700"))
  max_T <- as.integer(get_arg(args, "--max-T", "120"))
  T_first <- as.integer(get_arg(args, "--T-first", "5"))
  T_each <- as.integer(get_arg(args, "--T-each", "1"))
  N_iter <- as.integer(get_arg(args, "--N-iter", "3"))
  N_min <- as.integer(get_arg(args, "--N-min", as.character(max(3L, floor(2 + log2(3))))))
  name <- get_arg(args, "--name", NULL)
  out_prefix <- get_arg(args, "--out-prefix", "qap_run_overall")
  debug <- parse_flag(args, "--debug")
  debug_eval_hit <- !("--no-debug-eval-hit" %in% args)
  debug_softrestart <- parse_flag(args, "--debug-softrestart")
  debug_duplicate <- parse_flag(args, "--debug-duplicate")
  quiet <- parse_flag(args, "--quiet")

  run_qap_overall_rstudio(
    split_seeds = split_seeds,
    budgets = budgets,
    port_sizes = port_sizes,
    algo_seeds = algo_seeds,
    best_k = best_k,
    enable_second_level_test = enable_second_level_test,
    global_free_rider_elimination = global_free_rider_elimination,
    enable_duplicate = enable_duplicate,
    race_seed_offset = race_seed_offset,
    time_limit = time_limit,
    max_T = max_T,
    T_first = T_first,
    T_each = T_each,
    N_iter = N_iter,
    N_min = N_min,
    debug = debug,
    debug_eval_hit = debug_eval_hit,
    debug_softrestart = debug_softrestart,
    debug_duplicate = debug_duplicate,
    verbose = !isTRUE(quiet),
    name = name,
    out_prefix = out_prefix,
    save_outputs = isTRUE(save_outputs)
  )
}

main <- function() {
  run_qap_runoverall_from_args()
}

if (sys.nframe() == 0L) {
  res <- main()
} else {
  cat("Loaded qap_run_overall.R\n")
  cat("Use either:\n")
  cat("  run_qap_overall_rstudio(...)\n")
  cat("or:\n")
  cat("  run_qap_runoverall_from_args(c('--split-seeds','11','--budgets','500','--port-sizes','2','--algo-seeds','21,22','--best-k','1','--name','portsize2_20seeds_budget500'))\n")
}
