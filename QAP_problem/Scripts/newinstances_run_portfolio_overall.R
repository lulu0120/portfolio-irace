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
source(file.path(script_dir, "qap_portfolio_interface.R"), local = .GlobalEnv)
source(file.path(script_dir, "newinstances_generate_split.R"), local = .GlobalEnv)

parse_int_vec_newinst <- function(x) {
  vals <- strsplit(as.character(x), ",", fixed = TRUE)[[1L]]
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  as.integer(vals)
}

run_qap_newinstances_portfolio_grid <- function(method = c("portsize", "bestk"),
                                                split_seed = 11L,
                                                train_ratio = 0.5,
                                                budgets = 1000L,
                                                k = 2L,
                                                algo_seeds = 21:60,
                                                time_limit = 5,
                                                race_seed_offset = 700L,
                                                max_T = 120L,
                                                T_first = 5L,
                                                T_each = 1L,
                                                N_iter = 3L,
                                                N_min = 3L,
                                                enable_second_level_test = TRUE,
                                                global_free_rider_elimination = TRUE,
                                                enable_duplicate = TRUE,
                                                verbose = TRUE,
                                                out_prefix = NULL,
                                                save_outputs = TRUE) {
  method <- match.arg(method)
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
    out_prefix <- sprintf("qap_newinstances_%s%d_budget%s",
                          method, as.integer(k), paste(as.integer(unique(budgets)), collapse = "_"))
  }
  out_dir <- file.path(newinstances_root_default(), "results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  csv_path <- file.path(out_dir, paste0(out_prefix, ".csv"))
  rds_path <- file.path(out_dir, paste0(out_prefix, ".Rds"))

  port_size <- if (identical(method, "portsize")) as.integer(k) else 1L
  best_k <- if (identical(method, "bestk")) as.integer(k) else 1L

  for (idx in seq_len(nrow(grid))) {
    budget <- grid$budget[[idx]]
    algo_seed <- grid$algo_seed[[idx]]
    race_seed <- as.integer(race_seed_offset + algo_seed)

    if (isTRUE(verbose)) {
      cat(sprintf(
        "\n[QAP-newinstances-%s] split_seed=%d | algo_seed=%d | budget=%d | k=%d\n",
        method, as.integer(split_seed), algo_seed, budget, as.integer(k)
      ))
    }

    one <- run_qap_portfolio_once(
      train_split = split_paths$train,
      budget = budget,
      port_size = port_size,
      best_k = best_k,
      enable_second_level_test = isTRUE(enable_second_level_test),
      global_free_rider_elimination = isTRUE(global_free_rider_elimination),
      enable_duplicate = isTRUE(enable_duplicate),
      seed0 = algo_seed,
      race_seed0 = race_seed,
      time_limit = time_limit,
      max_T = max_T,
      T_first = T_first,
      T_each = T_each,
      N_iter = N_iter,
      N_min = N_min,
      verbose = verbose
    )

    rr <- one$result
    last_race <- if (length(rr$races) > 0L) rr$races[[length(rr$races)]]$result else NULL
    decode_state <- if (!is.null(rr$state)) rr$state else if (!is.null(last_race)) last_race$state else NULL
    top_k_names_raw <- get_top_k_qap_names_from_result(rr, k = best_k)
    best_name_chr <- as.character(rr$best_name)
    if (!is.na(best_name_chr) && nzchar(best_name_chr)) {
      rest <- top_k_names_raw[top_k_names_raw != best_name_chr]
      top_k_names <- c(best_name_chr, rest)
    } else {
      top_k_names <- top_k_names_raw
    }
    top_k_names <- top_k_names[seq_len(min(length(top_k_names), max(1L, as.integer(best_k))))]
    top_k_union <- sort(unique(unlist(lapply(top_k_names, parse_qap_port_name))))
    top_k_union <- top_k_union[is.finite(top_k_union)]
    top_k_union_name <- format_qap_port_name(top_k_union)
    top_k_names_decoded <- if (length(top_k_names) > 0L) {
      paste(vapply(top_k_names, function(nm) decode_qap_port_name(nm, decode_state, one$param_config), character(1)), collapse = " || ")
    } else {
      "<none>"
    }
    top_k_union_decoded <- if (length(top_k_union) > 0L) {
      decode_qap_port_name(top_k_union_name, decode_state, one$param_config)
    } else {
      "<none>"
    }

    rows[[idx]] <- data.frame(
      dataset = "NewInstances",
      split_seed = as.integer(split_seed),
      algo_seed = algo_seed,
      race_seed0 = race_seed,
      budget = budget,
      port_size = port_size,
      time_limit = as.numeric(time_limit),
      train_split = normalizePath(split_paths$train, winslash = "/", mustWork = TRUE),
      test_split = normalizePath(split_paths$test, winslash = "/", mustWork = TRUE),
      train_instances = length(readLines(split_paths$train, warn = FALSE)),
      test_instances = length(readLines(split_paths$test, warn = FALSE)),
      best_k = best_k,
      best_name = as.character(rr$best_name),
      best_name_decoded = as.character(one$best_name_decoded),
      best_k_names = paste(top_k_names, collapse = ";"),
      best_k_names_decoded = top_k_names_decoded,
      best_k_union = as.character(top_k_union_name),
      best_k_union_decoded = top_k_union_decoded,
      total_eval_used = as.integer(rr$total_eval_used),
      total_instances_used = as.integer(rr$total_instances_used),
      n_races = length(rr$races),
      stringsAsFactors = FALSE
    )
    details[[idx]] <- one

    if (isTRUE(verbose)) {
      cat(sprintf(
        "[QAP-newinstances-%s-done] split_seed=%d | algo_seed=%d | budget=%d | best_name=%s | best_decoded=%s\n",
        method,
        as.integer(split_seed),
        algo_seed,
        budget,
        as.character(rr$best_name),
        as.character(one$best_name_decoded)
      ))
      if (identical(method, "bestk")) {
        cat(sprintf(
          "  top_k_names=%s\n  top_k_union=%s\n  top_k_union_decoded=%s\n",
          paste(top_k_names, collapse = ";"),
          as.character(top_k_union_name),
          as.character(top_k_union_decoded)
        ))
      }
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
  method <- get_arg("--method", "portsize")
  split_seed <- as.integer(get_arg("--split-seed", "11"))
  train_ratio <- as.numeric(get_arg("--train-ratio", "0.5"))
  budgets <- parse_int_vec_newinst(get_arg("--budgets", "1000"))
  k <- as.integer(get_arg("--k", "2"))
  algo_seeds <- parse_int_vec_newinst(get_arg("--algo-seeds", "21,22,23"))
  time_limit <- as.numeric(get_arg("--time-limit", "5"))
  out_prefix <- get_arg("--name", "")

  run_qap_newinstances_portfolio_grid(
    method = method,
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
