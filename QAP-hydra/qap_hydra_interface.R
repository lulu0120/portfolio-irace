#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
  library(truncnorm)
}))

qap_hydra_find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  if (!dir.exists(current)) current <- dirname(current)
  repeat {
    if (dir.exists(file.path(current, "QAP_problem", "Scripts")) &&
        dir.exists(file.path(current, "benchmark", "hydra-irace"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Cannot resolve portfolio-irace repo root.")
}

repo_root <- local({
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    path <- if (!is.null(ctx)) ctx$path else ""
    if (is.character(path) && nzchar(path) && file.exists(path)) {
      return(qap_hydra_find_repo_root(dirname(path)))
    }
  }
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    return(qap_hydra_find_repo_root(dirname(sub("^--file=", "", file_flag[[1L]]))))
  }
  qap_hydra_find_repo_root(getwd())
})

source(file.path(repo_root, "QAP_problem", "Scripts", "qap_common.R"), local = .GlobalEnv)

load_qap_hydra_modules <- function(repo_root = NULL) {
  if (is.null(repo_root)) repo_root <- get("repo_root", envir = .GlobalEnv)
  multidim_dir <- file.path(repo_root, "Multiple_multidim")
  multidim_files <- c(
    "initial_and_env.R",
    "statistical_test.R",
    "elites_functions.R",
    "generate_discrete_input.R",
    "within_race.R",
    "wrapper_new.R"
  )
  hydra_dir <- file.path(repo_root, "benchmark", "hydra-irace")
  hydra_files <- c("hydra_wrapper.R")
  invisible(lapply(file.path(multidim_dir, multidim_files), function(f) source(f, local = .GlobalEnv)))
  invisible(lapply(file.path(hydra_dir, hydra_files), function(f) source(f, local = .GlobalEnv)))
}

rtnorm_trunc <- function(n, mean, sd, lower = LB, upper = UB, do_round = TRUE) {
  n <- as.integer(n)
  if (n <= 0L) return(numeric(0))
  if (length(mean) == 1L) mean <- rep(as.numeric(mean), n)
  if (length(sd) == 1L) sd <- rep(as.numeric(sd), n)
  if (length(mean) != n) stop("rtnorm_trunc: length(mean) must be 1 or n")
  if (length(sd) != n) stop("rtnorm_trunc: length(sd) must be 1 or n")
  mean <- as.numeric(mean)
  sd <- as.numeric(sd)
  x <- numeric(n)
  bad_sd <- !is.finite(sd) | sd <= 0
  if (any(bad_sd)) x[bad_sd] <- pmin(pmax(mean[bad_sd], lower), upper)
  ok <- !bad_sd
  if (any(ok)) x[ok] <- truncnorm::rtruncnorm(sum(ok), a = lower, b = upper, mean = mean[ok], sd = sd[ok])
  if (isTRUE(do_round)) x <- round(x)
  x
}

make_qap_param_config <- function(a_values = 0:9, p_values = 0:11, l_values = 0:5) {
  list(
    type = "multi",
    params = list(
      list(name = "a", type = "categorical", values = as.character(a_values)),
      list(name = "p", type = "categorical", values = as.character(p_values)),
      list(name = "l", type = "categorical", values = as.character(l_values))
    )
  )
}

read_qap_instance_table <- function(split_file) {
  paths <- readLines(split_file, warn = FALSE)
  paths <- paths[nzchar(trimws(paths))]
  data.frame(
    instance_id = seq_along(paths),
    instance_path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

get_qap_split_paths <- function(split_seed, repo_root = repo_root) {
  split_dir <- file.path(repo_root, "QAP_problem", "Scripts", "splits", paste0("seed_", as.integer(split_seed)))
  train_split <- file.path(split_dir, "train_instances.txt")
  test_split <- file.path(split_dir, "test_instances.txt")
  if (!file.exists(train_split) || !file.exists(test_split)) {
    stop("Split files not found for seed ", as.integer(split_seed), ": ", split_dir)
  }
  list(
    split_dir = split_dir,
    train_split = train_split,
    test_split = test_split
  )
}

qap_decode_cfg_vector <- function(c, state = NULL, param_config) {
  x <- if (!is.null(state) && !is.null(state$cfg_vec_by_id_env)) {
    key <- as.character(as.integer(c))
    if (exists(key, envir = state$cfg_vec_by_id_env, inherits = FALSE)) {
      as.numeric(get(key, envir = state$cfg_vec_by_id_env, inherits = FALSE))
    } else {
      as.numeric(c)
    }
  } else {
    as.numeric(c)
  }
  if (length(x) < length(param_config$params)) {
    x <- c(x, rep(1, length(param_config$params) - length(x)))
  }
  vals <- mapply(function(param, raw) {
    idx <- max(1L, min(length(param$values), as.integer(round(raw))))
    as.integer(param$values[[idx]])
  }, param_config$params, x, SIMPLIFY = TRUE, USE.NAMES = TRUE)
  out <- as.list(vals)
  names(out) <- vapply(param_config$params, function(param) as.character(param$name), character(1))
  out
}

decode_qap_cfg_label <- function(state, id, param_config) {
  vals <- tryCatch(
    qap_decode_cfg_vector(c = id, state = state, param_config = param_config),
    error = function(e) list(a = NA_integer_, p = NA_integer_, l = NA_integer_)
  )
  if (length(vals) == 0L || any(!c("a", "p", "l") %in% names(vals))) {
    return(sprintf("%d:[a=?,p=?,l=?]", as.integer(id)))
  }
  sprintf("%d:[a=%d,p=%d,l=%d]", as.integer(id), vals$a, vals$p, vals$l)
}

make_qap_score_functions <- function(instance_table, param_config, time_limit = 5, solver = NULL) {
  instance_lookup <- setNames(instance_table$instance_path, instance_table$instance_id)

  cost_config <- function(c, i,
                          w = NULL, p = NULL,
                          sigma2 = NULL,
                          instance_random = 0,
                          eps = NULL,
                          state = NULL) {
    decoded <- qap_decode_cfg_vector(c = c, state = state, param_config = param_config)
    instance_path <- instance_lookup[[as.character(as.integer(i))]]
    if (is.null(instance_path)) return(PENALTY)
    seed_val <- if (is.null(eps)) NA_integer_ else as.integer(abs(round(eps)))
    res <- run_qap_solver(
      instance = instance_path,
      a = decoded$a,
      p = decoded$p,
      l = decoded$l,
      t = time_limit,
      seed = seed_val,
      solver = solver
    )
    as.numeric(res$cost)
  }

  scores_configs_on_instance <- function(configs, i, base_seed,
                                         w = NULL, sigma1 = NULL, sigma2 = NULL,
                                         p = NULL,
                                         state = NULL) {
    configs <- as.numeric(configs)
    cfg_u <- sort(unique(configs))
    scores_u <- vapply(seq_along(cfg_u), function(k) {
      seed_val <- as.integer(base_seed + 1000L * as.integer(i) + k)
      decoded <- qap_decode_cfg_vector(c = cfg_u[[k]], state = state, param_config = param_config)
      instance_path <- instance_lookup[[as.character(as.integer(i))]]
      res <- run_qap_solver(
        instance = instance_path,
        a = decoded$a,
        p = decoded$p,
        l = decoded$l,
        t = time_limit,
        seed = seed_val,
        solver = solver
      )
      as.numeric(res$cost)
    }, numeric(1))
    scores_u[match(configs, cfg_u)]
  }

  list(
    cost_config = cost_config,
    scores_configs_on_instance = scores_configs_on_instance
  )
}

run_qap_hydra_once <- function(train_split,
                               budget,
                               target_portfolio_size = 2L,
                               best_k = 1L,
                               seed0 = 1L,
                               race_seed0 = 700L,
                               time_limit = 5,
                               a_values = 0:9,
                               p_values = 0:11,
                               l_values = 0:5,
                               T_first = 5L,
                               T_each = 1L,
                               N_min = 3L,
                               verbose = TRUE) {
  load_qap_hydra_modules(repo_root)

  param_config <- make_qap_param_config(a_values = a_values, p_values = p_values, l_values = l_values)
  instance_table <- read_qap_instance_table(train_split)
  scorer <- make_qap_score_functions(instance_table = instance_table, param_config = param_config, time_limit = time_limit)
  assign("cost_config", scorer$cost_config, envir = .GlobalEnv)
  assign("scores_configs_on_instance", scorer$scores_configs_on_instance, envir = .GlobalEnv)

  N_param <<- length(param_config$params)
  N_min_default <<- as.integer(N_min)
  T_first <<- as.integer(T_first)
  LB <<- 1
  UB <<- max(length(a_values), length(p_values), length(l_values))

  if (isTRUE(verbose)) {
    slot_budgets <- compute_hydra_budget_split(
      B_total = as.integer(budget),
      n_slots = as.integer(target_portfolio_size)
    )
    cat(sprintf(
      "[QAP-HYDRA] total_budget=%d | target_portfolio_size=%d | slot_budgets=%s | time_limit=%s\n",
      as.integer(budget),
      as.integer(target_portfolio_size),
      paste(as.integer(slot_budgets), collapse = ","),
      format(as.numeric(time_limit), scientific = FALSE)
    ))
  }

  res <- run_hydra_wrapper_irace(
    param_config = param_config,
    target_portfolio_size = as.integer(target_portfolio_size),
    B_total = as.integer(budget),
    training_instances = instance_table$instance_id,
    min_T_for_test = as.integer(T_first),
    T_each = as.integer(T_each),
    N_min = as.integer(N_min),
    seed0 = as.integer(seed0),
    race_seed0 = as.integer(race_seed0),
    verbose = isTRUE(verbose)
  )

  portfolio_ids <- as.integer(res$portfolio)
  decoded <- if (length(portfolio_ids) > 0L) {
    paste(vapply(portfolio_ids, function(id) decode_qap_cfg_label(res$state, id, param_config), character(1)), collapse = "; ")
  } else {
    "<none>"
  }
  best_name <- paste0("(", paste(portfolio_ids, collapse = ","), ")")

  list(
    result = res,
    best_name = best_name,
    best_name_decoded = decoded,
    best_k = as.integer(best_k),
    best_k_names = best_name,
    best_k_union = best_name,
    best_k_union_decoded = decoded,
    instance_table = instance_table,
    param_config = param_config
  )
}

run_qap_hydra_grid <- function(split_seeds,
                               budgets,
                               target_portfolio_sizes,
                               algo_seeds,
                               best_k = 1L,
                               race_seed_offset = 700L,
                               time_limit = 5,
                               T_first = 5L,
                               T_each = 1L,
                               N_min = 3L,
                               verbose = TRUE) {
  grid <- expand.grid(
    split_seed = as.integer(split_seeds),
    budget = as.integer(budgets),
    target_portfolio_size = as.integer(target_portfolio_sizes),
    algo_seed = as.integer(algo_seeds),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(grid))
  details <- vector("list", nrow(grid))

  for (idx in seq_len(nrow(grid))) {
    split_seed <- grid$split_seed[[idx]]
    budget <- grid$budget[[idx]]
    target_portfolio_size <- grid$target_portfolio_size[[idx]]
    algo_seed <- grid$algo_seed[[idx]]
    race_seed <- as.integer(race_seed_offset + algo_seed)
    split_paths <- get_qap_split_paths(split_seed = split_seed, repo_root = repo_root)

    if (isTRUE(verbose)) {
      cat(sprintf(
        "\n[QAP-HYDRA] split_seed=%d | algo_seed=%d | budget=%d | target_portfolio_size=%d\n",
        split_seed, algo_seed, budget, target_portfolio_size
      ))
    }

    one <- run_qap_hydra_once(
      train_split = split_paths$train_split,
      budget = budget,
      target_portfolio_size = target_portfolio_size,
      best_k = best_k,
      seed0 = algo_seed,
      race_seed0 = race_seed,
      time_limit = time_limit,
      T_first = T_first,
      T_each = T_each,
      N_min = N_min,
      verbose = verbose
    )

    rr <- one$result
    rows[[idx]] <- data.frame(
      split_seed = split_seed,
      algo_seed = algo_seed,
      race_seed0 = race_seed,
      budget = budget,
      target_portfolio_size = target_portfolio_size,
      time_limit = as.numeric(time_limit),
      train_split = normalizePath(split_paths$train_split, winslash = "/", mustWork = TRUE),
      test_split = normalizePath(split_paths$test_split, winslash = "/", mustWork = TRUE),
      train_instances = length(readLines(split_paths$train_split, warn = FALSE)),
      test_instances = length(readLines(split_paths$test_split, warn = FALSE)),
      best_k = as.integer(best_k),
      best_name = as.character(one$best_name),
      best_name_decoded = as.character(one$best_name_decoded),
      best_k_names = as.character(one$best_k_names),
      best_k_union = as.character(one$best_k_union),
      best_k_union_decoded = as.character(one$best_k_union_decoded),
      total_eval_used = as.integer(sum(vapply(rr$slot_results, function(x) x$solver_result$total_eval_used, numeric(1)))),
      total_instances_used = as.integer(sum(vapply(rr$slot_results, function(x) x$solver_result$total_instances_used, numeric(1)))),
      n_slots = length(rr$slot_results),
      stringsAsFactors = FALSE
    )
    details[[idx]] <- one
  }

  list(
    summary = do.call(rbind, rows),
    details = details
  )
}

run_qap_hydra_overall_rstudio <- function(split_seeds = 11L,
                                          budgets = c(1000L, 2500L),
                                          target_portfolio_sizes = c(2L, 4L),
                                          algo_seeds = 1L,
                                          best_k = 1L,
                                          race_seed_offset = 700L,
                                          time_limit = 5,
                                          T_first = 5L,
                                          T_each = 1L,
                                          N_min = 3L,
                                          verbose = TRUE,
                                          name = NULL,
                                          out_prefix = "qap_hydra_run",
                                          save_outputs = TRUE) {
  if (!is.null(name) && nzchar(as.character(name))) {
    out_prefix <- as.character(name)[[1L]]
  }

  cat("=== QAP Hydra-irace Run Overall ===\n")
  cat("split_seeds:", paste(as.integer(split_seeds), collapse = ","), "\n")
  cat("budgets:", paste(as.integer(budgets), collapse = ","), "\n")
  cat("target_portfolio_sizes:", paste(as.integer(target_portfolio_sizes), collapse = ","), "\n")
  cat("algo_seeds:", paste(as.integer(algo_seeds), collapse = ","), "\n")
  cat("best_k:", as.integer(best_k), "\n")
  cat("time_limit:", as.numeric(time_limit), "\n")

  res <- run_qap_hydra_grid(
    split_seeds = as.integer(split_seeds),
    budgets = as.integer(budgets),
    target_portfolio_sizes = as.integer(target_portfolio_sizes),
    algo_seeds = as.integer(algo_seeds),
    best_k = as.integer(best_k),
    race_seed_offset = as.integer(race_seed_offset),
    time_limit = as.numeric(time_limit),
    T_first = as.integer(T_first),
    T_each = as.integer(T_each),
    N_min = as.integer(N_min),
    verbose = isTRUE(verbose)
  )

  summary_df <- res$summary
  print(summary_df)

  out_dir <- file.path(repo_root, "QAP-hydra")
  csv_path <- file.path(out_dir, paste0(out_prefix, ".csv"))
  rds_path <- file.path(out_dir, paste0(out_prefix, ".Rds"))

  if (isTRUE(save_outputs)) {
    write.csv(summary_df, csv_path, row.names = FALSE)
    saveRDS(res, rds_path)
    cat("\nSaved summary:", csv_path, "\n")
    cat("Saved details:", rds_path, "\n")
  }

  invisible(res)
}
