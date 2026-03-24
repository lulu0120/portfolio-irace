#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

source("/Users/liuxiaolu/Desktop/portfolio-irace/QAP_problem/Scripts/qap_portfolio_interface.R", local = .GlobalEnv)

qap_verbose_mode <- function(verbose) {
  if (isFALSE(verbose)) return("none")
  if (isTRUE(verbose)) return("light")
  if (is.character(verbose) && length(verbose) == 1L) {
    x <- tolower(trimws(verbose))
    if (x %in% c("none", "false", "f", "0")) return("none")
    if (x %in% c("light", "true", "t", "1")) return("light")
    if (x %in% c("full", "verbose")) return("full")
  }
  stop("verbose must be FALSE, TRUE/'light', or 'full'.")
}

read_instances_list <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
  x <- x[nzchar(x)]
  normalizePath(x, winslash = "/", mustWork = TRUE)
}

normalize_qap_config <- function(config) {
  if (is.list(config)) {
    config <- unlist(config, use.names = TRUE)
  }
  if (is.data.frame(config)) {
    if (nrow(config) != 1L) stop("config data.frame must have exactly one row.")
    config <- unlist(config[1, , drop = TRUE], use.names = TRUE)
  }

  if (is.null(names(config)) || !all(c("a", "p", "l") %in% names(config))) {
    if (length(config) != 3L) {
      stop("config must be named c(a=..., p=..., l=...) or length-3 in order (a, p, l).")
    }
    names(config) <- c("a", "p", "l")
  }

  out <- c(
    a = as.integer(config[["a"]]),
    p = as.integer(config[["p"]]),
    l = as.integer(config[["l"]])
  )
  if (any(!is.finite(out))) {
    stop("config contains non-finite values.")
  }
  out
}

assign_qap_duplicate_ticks <- function(cfg_df) {
  cfg_df <- as.data.frame(cfg_df, stringsAsFactors = FALSE)
  keys <- vapply(seq_len(nrow(cfg_df)), function(i) {
    qap_config_key(cfg_df$a[[i]], cfg_df$p[[i]], cfg_df$l[[i]])
  }, character(1))
  seen <- new.env(parent = emptyenv())
  dup_tick <- integer(nrow(cfg_df))
  for (i in seq_len(nrow(cfg_df))) {
    key <- keys[[i]]
    current <- if (exists(key, envir = seen, inherits = FALSE)) get(key, envir = seen, inherits = FALSE) else 0L
    dup_tick[[i]] <- current
    assign(key, current + 1L, envir = seen)
  }
  cfg_df$duplicate_tick <- as.integer(dup_tick)
  cfg_df$member_index <- seq_len(nrow(cfg_df))
  cfg_df
}

normalize_qap_config_set <- function(config_set) {
  if (is.data.frame(config_set)) {
    need_cols <- c("a", "p", "l")
    if (!all(need_cols %in% names(config_set))) {
      stop("config_set data.frame must contain columns a, p, l.")
    }
    rows <- lapply(seq_len(nrow(config_set)), function(i) normalize_qap_config(config_set[i, , drop = FALSE]))
  } else if (is.list(config_set) && !is.null(names(config_set)) && all(c("a", "p", "l") %in% names(config_set))) {
    rows <- list(normalize_qap_config(config_set))
  } else if (is.atomic(config_set)) {
    rows <- list(normalize_qap_config(config_set))
  } else if (is.list(config_set)) {
    rows <- lapply(config_set, normalize_qap_config)
  } else {
    stop("config_set must be a single config, a list of configs, or a data.frame with columns a, p, l.")
  }

  if (length(rows) == 0L) {
    stop("config_set must contain at least one configuration.")
  }

  out <- do.call(rbind, lapply(rows, function(z) {
    data.frame(
      a = as.integer(z[["a"]]),
      p = as.integer(z[["p"]]),
      l = as.integer(z[["l"]]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  assign_qap_duplicate_ticks(out)
}

qap_config_key <- function(a, p, l, duplicate_tick = 0L) {
  sprintf(
    "a=%d|p=%d|l=%d|dup=%d",
    as.integer(a), as.integer(p), as.integer(l), as.integer(duplicate_tick)
  )
}

qap_eval_seed <- function(instance_index, duplicate_tick = 0L) {
  as.integer(100000L * as.integer(duplicate_tick) + as.integer(instance_index))
}

score_qap_config_on_instances <- function(config, instances, time_limit = 5, verbose = TRUE) {
  config <- normalize_qap_config(config)
  verbose_mode <- qap_verbose_mode(verbose)
  rows <- lapply(seq_along(instances), function(i) {
    inst <- instances[[i]]
    rr <- run_qap_solver(
      instance = inst,
        a = config[["a"]],
        p = config[["p"]],
        l = config[["l"]],
        t = time_limit,
        seed = qap_eval_seed(i, 0L)
      )
    if (verbose_mode != "none") {
      cat(sprintf(
        "[test-eval] %d/%d | a=%d p=%d l=%d | cost=%s\n",
        i, length(instances), config[["a"]], config[["p"]], config[["l"]], format(rr$cost, scientific = FALSE)
      ))
    }
    data.frame(
      instance = inst,
      cost = as.numeric(rr$cost),
      status = as.integer(rr$status),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

score_qap_config_grid_on_instances <- function(cfg_df, instances, time_limit = 5, verbose = TRUE) {
  cfg_df <- normalize_qap_config_set(cfg_df)
  verbose_mode <- qap_verbose_mode(verbose)
  score_mat <- matrix(NA_real_, nrow = length(instances), ncol = nrow(cfg_df))
  colnames(score_mat) <- vapply(seq_len(nrow(cfg_df)), function(k) {
    qap_config_key(cfg_df$a[[k]], cfg_df$p[[k]], cfg_df$l[[k]], cfg_df$duplicate_tick[[k]])
  }, character(1))

  for (i in seq_along(instances)) {
    inst <- instances[[i]]
    for (k in seq_len(nrow(cfg_df))) {
      rr <- run_qap_solver(
        instance = inst,
        a = cfg_df$a[[k]],
        p = cfg_df$p[[k]],
        l = cfg_df$l[[k]],
        t = time_limit,
        seed = qap_eval_seed(i, cfg_df$duplicate_tick[[k]])
      )
      score_mat[i, k] <- as.numeric(rr$cost)
      if (verbose_mode == "light" && k == 1L) {
        cat(sprintf(
          "[grid-test-eval] instance %d/%d\n",
          i, length(instances)
        ))
      }
      if (verbose_mode == "full") {
        if (k == 1L) {
          cat(sprintf(
            "[grid-test-eval] instance %d/%d\n",
            i, length(instances)
          ))
        }
        cat(sprintf(
          "  cfg %d/%d | a=%d p=%d l=%d dup=%d | cost=%s\n",
          k, nrow(cfg_df),
          cfg_df$a[[k]], cfg_df$p[[k]], cfg_df$l[[k]], cfg_df$duplicate_tick[[k]],
          format(score_mat[i, k], scientific = FALSE)
        ))
      }
    }
  }

  list(
    configs = cfg_df,
    instances = instances,
    score_matrix = score_mat
  )
}

score_qap_portfolio_on_instances <- function(config_set, instances, time_limit = 5, verbose = TRUE) {
  cfg_df <- normalize_qap_config_set(config_set)
  verbose_mode <- qap_verbose_mode(verbose)

  rows <- lapply(seq_along(instances), function(i) {
    inst <- instances[[i]]
    cfg_scores <- vapply(seq_len(nrow(cfg_df)), function(k) {
      rr <- run_qap_solver(
        instance = inst,
        a = cfg_df$a[[k]],
        p = cfg_df$p[[k]],
        l = cfg_df$l[[k]],
        t = time_limit,
        seed = qap_eval_seed(i, cfg_df$duplicate_tick[[k]])
      )
      as.numeric(rr$cost)
    }, numeric(1))

    best_idx <- which.min(cfg_scores)
    if (verbose_mode != "none") {
      cat(sprintf(
        "[portfolio-test-eval] %d/%d | size=%d | best_cfg=(a=%d,p=%d,l=%d) | cost=%s\n",
        i, length(instances), nrow(cfg_df),
        cfg_df$a[[best_idx]], cfg_df$p[[best_idx]], cfg_df$l[[best_idx]],
        format(min(cfg_scores), scientific = FALSE)
      ))
    }
    data.frame(
      instance = inst,
      portfolio_cost = min(cfg_scores),
      best_member_index = as.integer(best_idx),
      best_member_duplicate_tick = cfg_df$duplicate_tick[[best_idx]],
      best_member_a = cfg_df$a[[best_idx]],
      best_member_p = cfg_df$p[[best_idx]],
      best_member_l = cfg_df$l[[best_idx]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

score_qap_portfolio_from_grid <- function(config_set, pooled_scores, verbose = TRUE) {
  cfg_df <- normalize_qap_config_set(config_set)
  verbose_mode <- qap_verbose_mode(verbose)
  pooled_cfg <- pooled_scores$configs
  instances <- pooled_scores$instances
  score_mat <- pooled_scores$score_matrix
  keys_all <- vapply(seq_len(nrow(pooled_cfg)), function(k) {
    qap_config_key(pooled_cfg$a[[k]], pooled_cfg$p[[k]], pooled_cfg$l[[k]], pooled_cfg$duplicate_tick[[k]])
  }, character(1))
  keys_port <- vapply(seq_len(nrow(cfg_df)), function(k) {
    qap_config_key(cfg_df$a[[k]], cfg_df$p[[k]], cfg_df$l[[k]], cfg_df$duplicate_tick[[k]])
  }, character(1))
  idx <- match(keys_port, keys_all)
  if (anyNA(idx)) {
    stop("Portfolio contains config not found in pooled score grid.")
  }

  rows <- lapply(seq_along(instances), function(i) {
    cfg_scores <- score_mat[i, idx]
    best_idx <- which.min(cfg_scores)
    if (verbose_mode != "none") {
      cat(sprintf(
        "[portfolio-from-grid] %d/%d | size=%d | best_cfg=(a=%d,p=%d,l=%d,dup=%d) | cost=%s\n",
        i, length(instances), nrow(cfg_df),
        cfg_df$a[[best_idx]], cfg_df$p[[best_idx]], cfg_df$l[[best_idx]], cfg_df$duplicate_tick[[best_idx]],
        format(min(cfg_scores), scientific = FALSE)
      ))
    }
    data.frame(
      instance = instances[[i]],
      portfolio_cost = min(cfg_scores),
      best_member_index = as.integer(best_idx),
      best_member_duplicate_tick = cfg_df$duplicate_tick[[best_idx]],
      best_member_a = cfg_df$a[[best_idx]],
      best_member_p = cfg_df$p[[best_idx]],
      best_member_l = cfg_df$l[[best_idx]],
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

evaluate_qap_config_testset <- function(config,
                                        split_seed = 11L,
                                        test_split = NULL,
                                        time_limit = 5,
                                        out_prefix = NULL,
                                        save_outputs = TRUE,
                                        verbose = TRUE) {
  config <- normalize_qap_config(config)
  verbose_mode <- qap_verbose_mode(verbose)
  if (is.null(test_split) || !nzchar(test_split)) {
    split_paths <- get_qap_split_paths(
      split_seed = as.integer(split_seed),
      script_dir = qap_scripts_dir()
    )
    test_split <- split_paths$test_split
  }
  test_split <- normalizePath(test_split, winslash = "/", mustWork = TRUE)
  instances <- read_instances_list(test_split)
  instance_scores <- score_qap_config_on_instances(
    config = config,
    instances = instances,
    time_limit = time_limit,
    verbose = verbose_mode
  )

  summary_df <- data.frame(
    split_seed = as.integer(split_seed),
    test_split = test_split,
    a = as.integer(config[["a"]]),
    p = as.integer(config[["p"]]),
    l = as.integer(config[["l"]]),
    time_limit = as.numeric(time_limit),
    test_instances = nrow(instance_scores),
    test_mean_cost = mean(instance_scores$cost),
    test_median_cost = median(instance_scores$cost),
    test_sd_cost = if (nrow(instance_scores) > 1L) stats::sd(instance_scores$cost) else NA_real_,
    test_min_cost = min(instance_scores$cost),
    test_max_cost = max(instance_scores$cost),
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_df,
    instance_scores = instance_scores
  )

  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    out_prefix <- file.path(
      qap_scripts_dir(),
      sprintf("qap_test_eval_a%d_p%d_l%d_seed%d", config[["a"]], config[["p"]], config[["l"]], as.integer(split_seed))
    )
  }

  if (isTRUE(save_outputs)) {
    csv_path <- paste0(out_prefix, ".csv")
    rds_path <- paste0(out_prefix, ".Rds")
    write.csv(summary_df, csv_path, row.names = FALSE)
    saveRDS(out, rds_path)
    if (verbose_mode != "none") {
      cat("Saved summary:", csv_path, "\n")
      cat("Saved details:", rds_path, "\n")
    }
  }

  invisible(out)
}

evaluate_qap_portfolio_testset <- function(config_set,
                                           split_seed = 11L,
                                           test_split = NULL,
                                           time_limit = 5,
                                           out_prefix = NULL,
                                           save_outputs = TRUE,
                                           verbose = TRUE) {
  cfg_df <- normalize_qap_config_set(config_set)
  verbose_mode <- qap_verbose_mode(verbose)
  if (is.null(test_split) || !nzchar(test_split)) {
    split_paths <- get_qap_split_paths(
      split_seed = as.integer(split_seed),
      script_dir = qap_scripts_dir()
    )
    test_split <- split_paths$test_split
  }
  test_split <- normalizePath(test_split, winslash = "/", mustWork = TRUE)
  instances <- read_instances_list(test_split)
  instance_scores <- score_qap_portfolio_on_instances(
    config_set = cfg_df,
    instances = instances,
    time_limit = time_limit,
    verbose = verbose_mode
  )

  summary_df <- data.frame(
    split_seed = as.integer(split_seed),
    test_split = test_split,
    portfolio_size = nrow(cfg_df),
    portfolio_label = paste(apply(cfg_df, 1, function(z) sprintf("(a=%d,p=%d,l=%d,dup=%d)", z[["a"]], z[["p"]], z[["l"]], z[["duplicate_tick"]])), collapse = "; "),
    time_limit = as.numeric(time_limit),
    test_instances = nrow(instance_scores),
    test_mean_cost = mean(instance_scores$portfolio_cost),
    test_median_cost = median(instance_scores$portfolio_cost),
    test_sd_cost = if (nrow(instance_scores) > 1L) stats::sd(instance_scores$portfolio_cost) else NA_real_,
    test_min_cost = min(instance_scores$portfolio_cost),
    test_max_cost = max(instance_scores$portfolio_cost),
    stringsAsFactors = FALSE
  )

  out <- list(
    summary = summary_df,
    configs = cfg_df,
    instance_scores = instance_scores
  )

  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    out_prefix <- file.path(
      qap_scripts_dir(),
      sprintf("qap_portfolio_test_eval_size%d_seed%d", nrow(cfg_df), as.integer(split_seed))
    )
  }

  if (isTRUE(save_outputs)) {
    csv_path <- paste0(out_prefix, ".csv")
    rds_path <- paste0(out_prefix, ".Rds")
    write.csv(summary_df, csv_path, row.names = FALSE)
    saveRDS(out, rds_path)
    if (verbose_mode != "none") {
      cat("Saved summary:", csv_path, "\n")
      cat("Saved details:", rds_path, "\n")
    }
  }

  invisible(out)
}

normalize_qap_portfolio_set <- function(portfolios) {
  if (is.data.frame(portfolios)) {
    need_cols <- c("a", "p", "l")
    if (!all(need_cols %in% names(portfolios))) {
      stop("portfolio data.frame must contain columns a, p, l.")
    }
    return(list(portfolio_1 = normalize_qap_config_set(portfolios)))
  }

  if (!is.list(portfolios) || length(portfolios) == 0L) {
    stop("portfolios must be a non-empty list of portfolio config sets.")
  }

  is_single_named_config <- !is.null(names(portfolios)) && all(c("a", "p", "l") %in% names(portfolios))
  if (is_single_named_config) {
    return(list(portfolio_1 = normalize_qap_config_set(portfolios)))
  }

  out <- lapply(portfolios, normalize_qap_config_set)
  nm <- names(portfolios)
  if (is.null(nm) || any(!nzchar(nm))) {
    nm <- paste0("portfolio_", seq_along(out))
  }
  names(out) <- nm
  out
}

parse_qap_decoded_portfolio <- function(label) {
  if (is.null(label) || length(label) != 1L || !nzchar(trimws(label))) {
    stop("Decoded portfolio label must be a non-empty string.")
  }
  mm <- gregexpr("a=([0-9]+),p=([0-9]+),l=([0-9]+)", label, perl = TRUE)
  hits <- regmatches(label, mm)[[1L]]
  if (length(hits) == 0L) {
    stop(sprintf("Could not parse decoded portfolio label: %s", label))
  }
  rows <- lapply(hits, function(hit) {
    vals <- sub("^a=([0-9]+),p=([0-9]+),l=([0-9]+)$", "\\1,\\2,\\3", hit, perl = TRUE)
    vals <- as.integer(strsplit(vals, ",", fixed = TRUE)[[1L]])
    data.frame(
      a = vals[[1L]],
      p = vals[[2L]],
      l = vals[[3L]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

aggregate_budget_metric <- function(row_df, value_col, prefix) {
  mean_df <- stats::aggregate(
    row_df[[value_col]],
    by = list(method = row_df$method, budget = row_df$budget),
    FUN = mean
  )
  names(mean_df)[names(mean_df) == "x"] <- sprintf("%s_mean", prefix)

  median_df <- stats::aggregate(
    row_df[[value_col]],
    by = list(method = row_df$method, budget = row_df$budget),
    FUN = median
  )
  names(median_df)[names(median_df) == "x"] <- sprintf("%s_median", prefix)

  sd_df <- stats::aggregate(
    row_df[[value_col]],
    by = list(method = row_df$method, budget = row_df$budget),
    FUN = function(x) if (length(x) > 1L) stats::sd(x) else NA_real_
  )
  names(sd_df)[names(sd_df) == "x"] <- sprintf("%s_sd", prefix)

  n_df <- stats::aggregate(
    row_df[[value_col]],
    by = list(method = row_df$method, budget = row_df$budget),
    FUN = length
  )
  names(n_df)[names(n_df) == "x"] <- sprintf("%s_n", prefix)

  out <- merge(mean_df, median_df, by = c("method", "budget"), all = TRUE, sort = TRUE)
  out <- merge(out, sd_df, by = c("method", "budget"), all = TRUE, sort = TRUE)
  out <- merge(out, n_df, by = c("method", "budget"), all = TRUE, sort = TRUE)
  out[[sprintf("%s_se", prefix)]] <- out[[sprintf("%s_sd", prefix)]] / sqrt(out[[sprintf("%s_n", prefix)]])
  out
}

evaluate_qap_benchmark_df <- function(benchmark_df,
                                      portfolio_col = c("best_name_decoded", "best_k_union_decoded"),
                                      method_name = NULL,
                                      out_prefix = NULL,
                                      save_outputs = TRUE,
                                      verbose = TRUE) {
  portfolio_col <- match.arg(portfolio_col)
  verbose_mode <- qap_verbose_mode(verbose)
  df <- as.data.frame(benchmark_df, stringsAsFactors = FALSE)
  if (!portfolio_col %in% names(df)) {
    stop(sprintf("benchmark_df is missing column '%s'.", portfolio_col))
  }
  need_cols <- c("budget", "algo_seed", "split_seed", "test_split", "time_limit")
  miss <- setdiff(need_cols, names(df))
  if (length(miss) > 0L) {
    stop(sprintf("benchmark_df is missing required columns: %s", paste(miss, collapse = ", ")))
  }
  if (is.null(method_name) || !nzchar(method_name)) {
    method_name <- portfolio_col
  }

  df$.eval_name <- sprintf(
    "%s|b%s|s%s|row%s",
    method_name,
    df$budget,
    df$algo_seed,
    seq_len(nrow(df))
  )

  group_key <- paste(df$split_seed, df$test_split, df$time_limit, sep = "|")
  groups <- split(df, group_key, drop = TRUE)

  row_results <- vector("list", length(groups))
  detail_map <- list()
  pooled_map <- list()
  idx <- 1L

  for (group_name in names(groups)) {
    grp <- groups[[group_name]]
    portfolios <- setNames(
      lapply(grp[[portfolio_col]], parse_qap_decoded_portfolio),
      grp$.eval_name
    )
    if (verbose_mode != "none") {
      cat(sprintf(
        "[benchmark-eval] method=%s | group=%s | portfolios=%d | split_seed=%s\n",
        method_name, group_name, length(portfolios), grp$split_seed[[1L]]
      ))
    }
    eval_res <- evaluate_many_qap_portfolios_testset(
      portfolios = portfolios,
      split_seed = as.integer(grp$split_seed[[1L]]),
      test_split = grp$test_split[[1L]],
      time_limit = as.numeric(grp$time_limit[[1L]]),
      save_outputs = FALSE,
      verbose = verbose
    )
    summary_i <- eval_res$summary
    names(summary_i)[names(summary_i) == "portfolio_name"] <- ".eval_name"
    merged <- merge(grp, summary_i, by = ".eval_name", sort = FALSE, all.x = TRUE)
    merged$method <- method_name
    row_results[[idx]] <- merged
    detail_map[[group_name]] <- eval_res$details
    pooled_map[[group_name]] <- eval_res$pooled_configs
    idx <- idx + 1L
  }

  row_df <- do.call(rbind, row_results)
  row_df$.eval_name <- NULL
  rownames(row_df) <- NULL

  budget_mean_df <- aggregate_budget_metric(
    row_df = row_df,
    value_col = "test_mean_cost",
    prefix = "test_mean_cost"
  )
  budget_median_df <- aggregate_budget_metric(
    row_df = row_df,
    value_col = "test_median_cost",
    prefix = "test_median_cost"
  )
  budget_df <- merge(
    budget_mean_df,
    budget_median_df,
    by = c("method", "budget"),
    all = TRUE,
    sort = TRUE
  )
  keep_cols <- c(
    "method",
    "budget",
    "test_mean_cost_mean",
    "test_mean_cost_sd",
    "test_mean_cost_n",
    "test_mean_cost_se",
    "test_median_cost_mean",
    "test_median_cost_sd",
    "test_median_cost_n",
    "test_median_cost_se"
  )
  budget_df <- budget_df[, keep_cols, drop = FALSE]
  rownames(budget_df) <- NULL

  out <- list(
    row_results = row_df,
    budget_summary = budget_df,
    details = detail_map,
    pooled_configs = pooled_map
  )

  if (!is.null(out_prefix) && nzchar(out_prefix) && isTRUE(save_outputs)) {
    row_csv <- paste0(out_prefix, "_rows.csv")
    budget_csv <- paste0(out_prefix, "_by_budget.csv")
    rds_path <- paste0(out_prefix, ".Rds")
    write.csv(row_df, row_csv, row.names = FALSE)
    write.csv(budget_df, budget_csv, row.names = FALSE)
    saveRDS(out, rds_path)
    if (verbose_mode != "none") {
      cat("Saved row summary:", row_csv, "\n")
      cat("Saved budget summary:", budget_csv, "\n")
      cat("Saved details:", rds_path, "\n")
    }
  }

  invisible(out)
}

evaluate_qap_benchmark_csv <- function(csv_path,
                                       portfolio_col = c("best_name_decoded", "best_k_union_decoded"),
                                       method_name = NULL,
                                       out_prefix = NULL,
                                       save_outputs = TRUE,
                                       verbose = TRUE) {
  csv_path <- normalizePath(csv_path, winslash = "/", mustWork = TRUE)
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (is.null(method_name) || !nzchar(method_name)) {
    method_name <- tools::file_path_sans_ext(basename(csv_path))
  }
  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    out_prefix <- file.path(qap_scripts_dir(), method_name)
  }
  evaluate_qap_benchmark_df(
    benchmark_df = df,
    portfolio_col = portfolio_col,
    method_name = method_name,
    out_prefix = out_prefix,
    save_outputs = save_outputs,
    verbose = verbose
  )
}

plot_qap_budget_comparison <- function(budget_summary,
                                       out_prefix = NULL,
                                       main = NULL,
                                       show_median = TRUE) {
  df <- as.data.frame(budget_summary, stringsAsFactors = FALSE)
  methods <- unique(df$method)
  budgets <- sort(unique(df$budget))
  cols <- seq_along(methods)
  names(cols) <- methods
  if (is.null(main)) main <- "QAP Test Score By Budget"

  if (!is.null(out_prefix) && nzchar(out_prefix)) {
    png(filename = paste0(out_prefix, ".png"), width = 960, height = 640)
    on.exit(dev.off(), add = TRUE)
  }

  y_lo <- pmin(
    df$test_mean_cost_mean - ifelse(is.na(df$test_mean_cost_sd), 0, df$test_mean_cost_sd),
    df$test_median_cost_mean,
    na.rm = TRUE
  )
  y_hi <- pmax(
    df$test_mean_cost_mean + ifelse(is.na(df$test_mean_cost_sd), 0, df$test_mean_cost_sd),
    df$test_median_cost_mean,
    na.rm = TRUE
  )
  plot(
    range(budgets),
    range(c(y_lo, y_hi), na.rm = TRUE),
    type = "n",
    xlab = "Budget",
    ylab = "Average Test Score Across Seeds",
    main = main
  )
  for (m in methods) {
    sub <- df[df$method == m, , drop = FALSE]
    sub <- sub[order(sub$budget), , drop = FALSE]
    lines(sub$budget, sub$test_mean_cost_mean, type = "b", lwd = 2, col = cols[[m]], pch = 19, lty = 1)
    if (isTRUE(show_median)) {
      lines(sub$budget, sub$test_median_cost_mean, type = "b", lwd = 2, col = cols[[m]], pch = 1, lty = 2)
    }
    err <- ifelse(is.na(sub$test_mean_cost_sd), 0, sub$test_mean_cost_sd)
    ok <- is.finite(err) & err > 0
    if (any(ok)) {
      arrows(
        x0 = sub$budget[ok],
        y0 = sub$test_mean_cost_mean[ok] - err[ok],
        x1 = sub$budget[ok],
        y1 = sub$test_mean_cost_mean[ok] + err[ok],
        angle = 90,
        code = 3,
        length = 0.05,
        col = cols[[m]],
        lwd = 1.5
      )
    }
  }
  if (isTRUE(show_median)) {
    legend(
      "topright",
      legend = c(
        paste(methods, "mean"),
        paste(methods, "median")
      ),
      col = c(cols[methods], cols[methods]),
      lwd = 2,
      pch = c(rep(19, length(methods)), rep(1, length(methods))),
      lty = c(rep(1, length(methods)), rep(2, length(methods))),
      bty = "n"
    )
  } else {
    legend(
      "topright",
      legend = paste(methods, "mean"),
      col = cols[methods],
      lwd = 2,
      pch = 19,
      lty = 1,
      bty = "n"
    )
  }
  invisible(df)
}

compare_qap_benchmark_csvs <- function(csv_paths,
                                       portfolio_cols,
                                       method_names = NULL,
                                       out_prefix = NULL,
                                       save_outputs = TRUE,
                                       verbose = TRUE,
                                       show_median = TRUE) {
  if (length(csv_paths) != length(portfolio_cols)) {
    stop("csv_paths and portfolio_cols must have the same length.")
  }
  if (is.null(method_names)) {
    method_names <- tools::file_path_sans_ext(basename(csv_paths))
  }
  if (length(method_names) != length(csv_paths)) {
    stop("method_names must have the same length as csv_paths.")
  }

  verbose_mode <- qap_verbose_mode(verbose)
  benchmark_list <- lapply(seq_along(csv_paths), function(i) {
    csv_path <- normalizePath(csv_paths[[i]], winslash = "/", mustWork = TRUE)
    df <- read.csv(csv_path, stringsAsFactors = FALSE)
    if (!portfolio_cols[[i]] %in% names(df)) {
      stop(sprintf("benchmark_df is missing column '%s'.", portfolio_cols[[i]]))
    }
    need_cols <- c("budget", "algo_seed", "split_seed", "test_split", "time_limit")
    miss <- setdiff(need_cols, names(df))
    if (length(miss) > 0L) {
      stop(sprintf("benchmark_df is missing required columns: %s", paste(miss, collapse = ", ")))
    }
    df$.portfolio_col <- portfolio_cols[[i]]
    df$method <- method_names[[i]]
    df
  })

  all_df <- do.call(rbind, benchmark_list)
  rownames(all_df) <- NULL
  all_df$.eval_name <- sprintf(
    "%s|b%s|s%s|row%s",
    all_df$method,
    all_df$budget,
    all_df$algo_seed,
    seq_len(nrow(all_df))
  )

  group_key <- paste(all_df$split_seed, all_df$test_split, all_df$time_limit, sep = "|")
  groups <- split(all_df, group_key, drop = TRUE)
  row_results <- vector("list", length(groups))
  detail_map <- list()
  pooled_map <- list()
  idx <- 1L

  for (group_name in names(groups)) {
    grp <- groups[[group_name]]
    portfolios <- setNames(
      lapply(seq_len(nrow(grp)), function(i) {
        parse_qap_decoded_portfolio(grp[[grp$.portfolio_col[[i]]]][[i]])
      }),
      grp$.eval_name
    )
    if (verbose_mode != "none") {
      cat(sprintf(
        "[benchmark-compare] group=%s | methods=%s | portfolios=%d | split_seed=%s\n",
        group_name,
        paste(sort(unique(grp$method)), collapse = ","),
        length(portfolios),
        grp$split_seed[[1L]]
      ))
    }
    eval_res <- evaluate_many_qap_portfolios_testset(
      portfolios = portfolios,
      split_seed = as.integer(grp$split_seed[[1L]]),
      test_split = grp$test_split[[1L]],
      time_limit = as.numeric(grp$time_limit[[1L]]),
      save_outputs = FALSE,
      verbose = verbose
    )
    summary_i <- eval_res$summary
    names(summary_i)[names(summary_i) == "portfolio_name"] <- ".eval_name"
    merged <- merge(grp, summary_i, by = ".eval_name", sort = FALSE, all.x = TRUE)
    row_results[[idx]] <- merged
    detail_map[[group_name]] <- eval_res$details
    pooled_map[[group_name]] <- eval_res$pooled_configs
    idx <- idx + 1L
  }

  row_df <- do.call(rbind, row_results)
  row_df$.eval_name <- NULL
  row_df$.portfolio_col <- NULL
  rownames(row_df) <- NULL

  budget_mean_df <- aggregate_budget_metric(
    row_df = row_df,
    value_col = "test_mean_cost",
    prefix = "test_mean_cost"
  )
  budget_median_df <- aggregate_budget_metric(
    row_df = row_df,
    value_col = "test_median_cost",
    prefix = "test_median_cost"
  )
  budget_df <- merge(
    budget_mean_df,
    budget_median_df,
    by = c("method", "budget"),
    all = TRUE,
    sort = TRUE
  )
  keep_cols <- c(
    "method",
    "budget",
    "test_mean_cost_mean",
    "test_mean_cost_sd",
    "test_mean_cost_n",
    "test_mean_cost_se",
    "test_median_cost_mean",
    "test_median_cost_sd",
    "test_median_cost_n",
    "test_median_cost_se"
  )
  budget_df <- budget_df[, keep_cols, drop = FALSE]
  rownames(budget_df) <- NULL

  by_method <- split(row_df, row_df$method, drop = TRUE)
  out <- list(
    row_results = row_df,
    budget_summary = budget_df,
    by_method = by_method,
    details = detail_map,
    pooled_configs = pooled_map
  )

  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    out_prefix <- file.path(qap_scripts_dir(), "qap_benchmark_compare")
  }

  if (isTRUE(save_outputs)) {
    row_csv <- paste0(out_prefix, "_rows.csv")
    budget_csv <- paste0(out_prefix, "_by_budget.csv")
    rds_path <- paste0(out_prefix, ".Rds")
    write.csv(row_df, row_csv, row.names = FALSE)
    write.csv(budget_df, budget_csv, row.names = FALSE)
    saveRDS(out, rds_path)
    if (verbose_mode != "none") {
      cat("Saved row summary:", row_csv, "\n")
      cat("Saved budget summary:", budget_csv, "\n")
      cat("Saved details:", rds_path, "\n")
    }
  }

  plot_qap_budget_comparison(
    budget_summary = budget_df,
    out_prefix = out_prefix,
    show_median = isTRUE(show_median)
  )

  invisible(out)
}

evaluate_many_qap_portfolios_testset <- function(portfolios,
                                                 split_seed = 11L,
                                                 test_split = NULL,
                                                 time_limit = 5,
                                                 out_prefix = NULL,
                                                 save_outputs = TRUE,
                                                 verbose = TRUE) {
  verbose_mode <- qap_verbose_mode(verbose)
  portfolio_list <- normalize_qap_portfolio_set(portfolios)
  if (is.null(test_split) || !nzchar(test_split)) {
    split_paths <- get_qap_split_paths(
      split_seed = as.integer(split_seed),
      script_dir = qap_scripts_dir()
    )
    test_split <- split_paths$test_split
  }
  test_split <- normalizePath(test_split, winslash = "/", mustWork = TRUE)
  instances <- read_instances_list(test_split)
  pooled_cfg_df <- unique(do.call(rbind, lapply(portfolio_list, function(df) {
    df[, c("a", "p", "l", "duplicate_tick"), drop = FALSE]
  })))
  rownames(pooled_cfg_df) <- NULL
  if (verbose_mode != "none") {
    cat(sprintf(
      "[multi-portfolio-test-eval] pooled unique configs: %d across %d portfolios\n",
      nrow(pooled_cfg_df), length(portfolio_list)
    ))
  }
  pooled_scores <- score_qap_config_grid_on_instances(
    cfg_df = pooled_cfg_df,
    instances = instances,
    time_limit = time_limit,
    verbose = verbose_mode
  )

  result_list <- lapply(seq_along(portfolio_list), function(i) {
    name_i <- names(portfolio_list)[[i]]
    cfg_df <- portfolio_list[[i]]
    if (verbose_mode != "none") {
      cat(sprintf(
        "[multi-portfolio-test-eval] %d/%d | %s | size=%d\n",
        i, length(portfolio_list), name_i, nrow(cfg_df)
      ))
    }
    instance_scores <- score_qap_portfolio_from_grid(
      config_set = cfg_df,
      pooled_scores = pooled_scores,
      verbose = verbose_mode
    )
        summary_df <- data.frame(
      split_seed = as.integer(split_seed),
      test_split = test_split,
      portfolio_size = nrow(cfg_df),
      portfolio_label = paste(apply(cfg_df, 1, function(z) sprintf("(a=%d,p=%d,l=%d,dup=%d)", z[["a"]], z[["p"]], z[["l"]], z[["duplicate_tick"]])), collapse = "; "),
      time_limit = as.numeric(time_limit),
      test_instances = nrow(instance_scores),
      test_mean_cost = mean(instance_scores$portfolio_cost),
      test_median_cost = median(instance_scores$portfolio_cost),
      test_sd_cost = if (nrow(instance_scores) > 1L) stats::sd(instance_scores$portfolio_cost) else NA_real_,
      test_min_cost = min(instance_scores$portfolio_cost),
      test_max_cost = max(instance_scores$portfolio_cost),
      portfolio_name = name_i,
      stringsAsFactors = FALSE
    )
    list(summary = summary_df, configs = cfg_df, instance_scores = instance_scores)
  })

  summary_df <- do.call(rbind, lapply(result_list, function(x) x$summary))
  summary_df <- summary_df[, c(
    "portfolio_name", "split_seed", "test_split", "portfolio_size", "portfolio_label",
    "time_limit", "test_instances", "test_mean_cost", "test_median_cost",
    "test_sd_cost", "test_min_cost", "test_max_cost"
  )]
  rownames(summary_df) <- NULL

  out <- list(
    summary = summary_df,
    portfolios = lapply(result_list, function(x) x$configs),
    details = setNames(lapply(result_list, function(x) x$instance_scores), summary_df$portfolio_name),
    pooled_configs = pooled_cfg_df
  )

  if (is.null(out_prefix) || !nzchar(out_prefix)) {
    out_prefix <- file.path(
      qap_scripts_dir(),
      sprintf("qap_many_portfolios_test_eval_n%d_seed%d", length(portfolio_list), as.integer(split_seed))
    )
  }

  if (isTRUE(save_outputs)) {
    csv_path <- paste0(out_prefix, ".csv")
    rds_path <- paste0(out_prefix, ".Rds")
    write.csv(summary_df, csv_path, row.names = FALSE)
    saveRDS(out, rds_path)
    if (verbose_mode != "none") {
      cat("Saved summary:", csv_path, "\n")
      cat("Saved details:", rds_path, "\n")
    }
  }

  invisible(out)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  get_arg <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[[idx + 1L]]
  }

  a <- as.integer(get_arg("--a", NA_character_))
  p <- as.integer(get_arg("--p", NA_character_))
  l <- as.integer(get_arg("--l", NA_character_))
  split_seed <- as.integer(get_arg("--split-seed", "11"))
  test_split <- get_arg("--test-split", NULL)
  time_limit <- as.numeric(get_arg("--time-limit", "5"))
  out_prefix <- get_arg("--out-prefix", NULL)

  if (!is.finite(a) || !is.finite(p) || !is.finite(l)) {
    stop("Usage: Rscript evaluate_qap_config_testset.R --a <int> --p <int> --l <int> [--split-seed <int>] [--test-split <path>] [--time-limit <num>] [--out-prefix <prefix>]")
  }

  evaluate_qap_config_testset(
    config = c(a = a, p = p, l = l),
    split_seed = split_seed,
    test_split = test_split,
    time_limit = time_limit,
    out_prefix = out_prefix,
    save_outputs = TRUE,
    verbose = TRUE
  )
}

if (sys.nframe() == 0L) {
  main()
} else {
  cat("Loaded evaluate_qap_config_testset.R\n")
  cat("Use:\n")
  cat("  evaluate_qap_config_testset(c(a=2,p=5,l=1), split_seed=11, time_limit=5, verbose=TRUE)\n")
  cat("  evaluate_qap_portfolio_testset(data.frame(a=c(1,2),p=c(7,9),l=c(0,0)), split_seed=11, time_limit=5, verbose=\"full\")\n")
  cat("  evaluate_many_qap_portfolios_testset(list(p1=data.frame(a=c(1,2),p=c(7,9),l=c(0,0)), p2=data.frame(a=c(1,1),p=c(9,6),l=c(0,0))), split_seed=11, time_limit=5, verbose=FALSE)\n")
  cat("  evaluate_qap_benchmark_csv(\"/abs/path/qap_bestk2_budget200.csv\", portfolio_col=\"best_k_union_decoded\", method_name=\"bestk2\", out_prefix=\"/abs/path/qap_bestk2_eval\", verbose=TRUE)\n")
  cat("  compare_qap_benchmark_csvs(csv_paths=c(\"/abs/path/qap_bestk2_budget200.csv\",\"/abs/path/qap_port_size2_budget200.csv\"), portfolio_cols=c(\"best_k_union_decoded\",\"best_name_decoded\"), method_names=c(\"bestk2\",\"port2\"), out_prefix=\"/abs/path/qap_compare\", verbose=TRUE)\n")
}
