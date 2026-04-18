#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

script_dir <- local({
  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    dirname(normalizePath(sub("^--file=", "", file_flag[[1L]]), winslash = "/", mustWork = TRUE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  }
})
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "newinstances_generate_split.R"), local = .GlobalEnv)
source(file.path(repo_root, "QAP-hydra", "qap_hydra_interface.R"), local = .GlobalEnv)

parse_int_vec_newinst_hydra <- function(x) {
  vals <- strsplit(as.character(x), ",", fixed = TRUE)[[1L]]
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  as.integer(vals)
}

run_qap_newinstances_hydra_grid <- function(split_seed = 11L,
                                            train_ratio = 0.5,
                                            budgets = 1000L,
                                            k = 2L,
                                            algo_seeds = 21:60,
                                            time_limit = 5,
                                            T_first = 5L,
                                            T_each = 1L,
                                            N_min = 3L,
                                            verbose = TRUE,
                                            out_prefix = NULL,
                                            save_outputs = TRUE) {
  split_paths <- generate_newinstances_split(split_seed = split_seed, train_ratio = train_ratio)

  grid <- expand.grid(
    budget = as.integer(budgets),
    algo_seed = as.integer(algo_seeds),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(grid))
  details <- vector("list", nrow(grid))

  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    out_prefix <- sprintf("qap_newinstances_hydra%d_budget%s",
                          as.integer(k), paste(as.integer(unique(budgets)), collapse = "_"))
  }
  out_dir <- file.path(newinstances_root_default(), "results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(out_dir, paste0(out_prefix, ".csv"))
  rds_path <- file.path(out_dir, paste0(out_prefix, ".Rds"))

  for (idx in seq_len(nrow(grid))) {
    budget <- grid$budget[[idx]]
    algo_seed <- grid$algo_seed[[idx]]
    race_seed <- as.integer(700L + algo_seed)

    if (isTRUE(verbose)) {
      slot_budgets <- if (exists("compute_hydra_budget_split", mode = "function", inherits = TRUE)) {
        compute_hydra_budget_split(B_total = as.integer(budget), n_slots = as.integer(k))
      } else {
        integer(0)
      }
      slot_part <- if (length(slot_budgets) > 0L) {
        sprintf(" | slot_budgets=%s", paste(as.integer(slot_budgets), collapse = ","))
      } else {
        ""
      }
      cat(sprintf(
        "\n[QAP-newinstances-hydra] split_seed=%d | algo_seed=%d | total_budget=%d | target_portfolio_size=%d%s\n",
        as.integer(split_seed), algo_seed, budget, as.integer(k), slot_part
      ))
    }

    one <- run_qap_hydra_once(
      train_split = split_paths$train,
      budget = budget,
      target_portfolio_size = as.integer(k),
      best_k = 1L,
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
      dataset = "NewInstances",
      split_seed = as.integer(split_seed),
      algo_seed = algo_seed,
      race_seed0 = race_seed,
      budget = budget,
      port_size = 1L,
      target_portfolio_size = as.integer(k),
      time_limit = as.numeric(time_limit),
      train_split = normalizePath(split_paths$train, winslash = "/", mustWork = TRUE),
      test_split = normalizePath(split_paths$test, winslash = "/", mustWork = TRUE),
      train_instances = length(readLines(split_paths$train, warn = FALSE)),
      test_instances = length(readLines(split_paths$test, warn = FALSE)),
      best_k = 1L,
      best_name = as.character(one$best_name),
      best_name_decoded = as.character(one$best_name_decoded),
      best_k_names = as.character(one$best_k_names),
      best_k_names_decoded = as.character(one$best_name_decoded),
      best_k_union = as.character(one$best_k_union),
      best_k_union_decoded = as.character(one$best_k_union_decoded),
      total_eval_used = as.integer(sum(vapply(rr$slot_results, function(x) x$solver_result$total_eval_used, numeric(1)))),
      total_instances_used = as.integer(sum(vapply(rr$slot_results, function(x) x$solver_result$total_instances_used, numeric(1)))),
      n_races = length(rr$slot_results),
      stringsAsFactors = FALSE
    )
    details[[idx]] <- one

    if (isTRUE(verbose)) {
      cat(sprintf(
        "[QAP-newinstances-hydra-done] split_seed=%d | algo_seed=%d | total_budget=%d | best_name=%s | best_decoded=%s\n",
        as.integer(split_seed),
        algo_seed,
        budget,
        as.character(one$best_name),
        as.character(one$best_name_decoded)
      ))
      cat(sprintf(
        "  best_k_union=%s\n  best_k_union_decoded=%s\n",
        as.character(one$best_k_union),
        as.character(one$best_k_union_decoded)
      ))
    }

    if (isTRUE(save_outputs)) {
      summary_partial <- do.call(rbind, rows[seq_len(idx)])
      details_partial <- details[seq_len(idx)]
      saveRDS(one, file.path(out_dir, sprintf("%s_seed%d_budget%d.Rds", out_prefix, algo_seed, budget)))
      utils::write.csv(summary_partial, csv_path, row.names = FALSE)
      saveRDS(list(summary = summary_partial, details = details_partial, split_paths = split_paths), rds_path)
    }
  }

  res <- list(summary = do.call(rbind, rows), details = details, split_paths = split_paths)
  if (isTRUE(save_outputs)) {
    utils::write.csv(res$summary, csv_path, row.names = FALSE)
    saveRDS(res, rds_path)
  }
  invisible(res)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[[idx + 1L]]
  }
  split_seed <- as.integer(get_arg("--split-seed", "11"))
  train_ratio <- as.numeric(get_arg("--train-ratio", "0.5"))
  budgets <- parse_int_vec_newinst_hydra(get_arg("--budgets", "1000"))
  k <- as.integer(get_arg("--k", "2"))
  algo_seeds <- parse_int_vec_newinst_hydra(get_arg("--algo-seeds", "21,22,23"))
  time_limit <- as.numeric(get_arg("--time-limit", "5"))
  out_prefix <- get_arg("--name", "")

  run_qap_newinstances_hydra_grid(
    split_seed = split_seed,
    train_ratio = train_ratio,
    budgets = budgets,
    k = k,
    algo_seeds = algo_seeds,
    time_limit = time_limit,
    out_prefix = out_prefix,
    save_outputs = TRUE,
    verbose = TRUE
  )
}
