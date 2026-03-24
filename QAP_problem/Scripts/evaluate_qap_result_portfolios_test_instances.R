#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

source("/Users/liuxiaolu/Desktop/portfolio-irace/QAP_problem/Scripts/evaluate_qap_config_testset.R", local = .GlobalEnv)

resolve_group_compare_test_split <- function(family, split_seed, script_dir = qap_scripts_dir()) {
  if (is.null(family) || !nzchar(trimws(family))) {
    stop("family is required when test_split is not provided.")
  }
  split_dir <- file.path(
    script_dir,
    "group_compare",
    "splits",
    as.character(family),
    paste0("seed_", as.integer(split_seed))
  )
  test_split <- file.path(split_dir, "test_instances.txt")
  if (!file.exists(test_split)) {
    stop("Test split not found: ", test_split)
  }
  normalizePath(test_split, winslash = "/", mustWork = TRUE)
}

resolve_group_compare_test_splits <- function(families, split_seed, script_dir = qap_scripts_dir()) {
  families <- as.character(families)
  families <- trimws(families)
  families <- families[nzchar(families)]
  if (length(families) == 0L) {
    stop("families must contain at least one non-empty family name.")
  }
  unique(vapply(families, function(fam) {
    resolve_group_compare_test_split(
      family = fam,
      split_seed = split_seed,
      script_dir = script_dir
    )
  }, character(1)))
}

read_instances_from_splits <- function(test_splits) {
  test_splits <- as.character(test_splits)
  test_splits <- test_splits[nzchar(trimws(test_splits))]
  if (length(test_splits) == 0L) {
    stop("test_splits must contain at least one path.")
  }
  unique(unlist(lapply(test_splits, read_instances_list), use.names = FALSE))
}

compute_qap_average_rank_by_method <- function(instance_scores_wide, source_metadata, detail_names) {
  if (is.null(instance_scores_wide) || nrow(instance_scores_wide) == 0L) {
    return(list(by_method = data.frame(), by_instance_method = data.frame(), by_portfolio = data.frame()))
  }
  if (is.null(source_metadata) || nrow(source_metadata) == 0L) {
    stop("source_metadata is required to compute method-level ranks.")
  }
  detail_names <- as.character(detail_names)
  score_cols <- make.names(detail_names, unique = TRUE)
  missing_cols <- setdiff(score_cols, names(instance_scores_wide))
  if (length(missing_cols) > 0L) {
    stop("Missing portfolio score columns in instance_scores_wide: ", paste(missing_cols, collapse = ", "))
  }

  rank_mat <- t(apply(instance_scores_wide[, score_cols, drop = FALSE], 1, function(x) {
    rank(as.numeric(x), ties.method = "average")
  }))
  rank_df <- as.data.frame(rank_mat, stringsAsFactors = FALSE)
  names(rank_df) <- score_cols
  rank_df$instance <- instance_scores_wide$instance

  portfolio_map <- data.frame(
    portfolio_name = detail_names,
    score_col = score_cols,
    stringsAsFactors = FALSE
  )
  portfolio_map <- merge(
    portfolio_map,
    unique(source_metadata[, c("portfolio_name", "method")]),
    by = "portfolio_name",
    all.x = TRUE,
    sort = FALSE
  )

  long_rows <- lapply(seq_len(nrow(portfolio_map)), function(i) {
    data.frame(
      instance = rank_df$instance,
      portfolio_name = portfolio_map$portfolio_name[[i]],
      method = portfolio_map$method[[i]],
      rank = as.numeric(rank_df[[portfolio_map$score_col[[i]]]]),
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_rows)

  by_portfolio <- stats::aggregate(
    rank ~ method + portfolio_name,
    data = long_df,
    FUN = mean
  )
  names(by_portfolio)[names(by_portfolio) == "rank"] <- "average_rank"
  rownames(by_portfolio) <- NULL

  by_instance_method <- stats::aggregate(
    rank ~ instance + method,
    data = long_df,
    FUN = mean
  )
  names(by_instance_method)[names(by_instance_method) == "rank"] <- "average_rank"
  rownames(by_instance_method) <- NULL

  by_method <- stats::aggregate(
    average_rank ~ method,
    data = by_instance_method,
    FUN = mean
  )
  names(by_method)[names(by_method) == "average_rank"] <- "average_rank_over_instances"
  by_method$n_instances <- as.integer(vapply(split(by_instance_method$instance, by_instance_method$method), length, integer(1)))
  rownames(by_method) <- NULL
  by_method <- by_method[order(by_method$average_rank_over_instances), , drop = FALSE]

  list(
    by_method = by_method,
    by_instance_method = by_instance_method,
    by_portfolio = by_portfolio
  )
}

plot_qap_method_rank_distribution <- function(res,
                                              out_path = NULL,
                                              main = "Method Rank Distribution",
                                              ylab = "Per-instance Average Rank",
                                              point_alpha = 0.55) {
  if (is.null(res$rank_summary$by_instance_method) || nrow(res$rank_summary$by_instance_method) == 0L) {
    stop("res$rank_summary$by_instance_method is empty.")
  }
  df <- res$rank_summary$by_instance_method
  method_order <- res$rank_summary$by_method$method
  if (is.null(method_order) || length(method_order) == 0L) {
    method_order <- sort(unique(df$method))
  }
  df$method <- factor(df$method, levels = method_order)
  cols <- c(bestk2 = "#E69F00", portsize2 = "#0072B2", hydra2 = "#009E73")
  fill_cols <- cols[as.character(method_order)]
  fill_cols[is.na(fill_cols)] <- "grey70"

  if (!is.null(out_path) && nzchar(trimws(out_path))) {
    grDevices::pdf(out_path, width = 7, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  op <- graphics::par(mar = c(5, 5, 3, 1) + 0.1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::boxplot(
    average_rank ~ method,
    data = df,
    col = fill_cols,
    border = "grey25",
    outline = FALSE,
    main = main,
    xlab = "",
    ylab = ylab
  )
  for (i in seq_along(method_order)) {
    vals <- df$average_rank[df$method == method_order[[i]]]
    xj <- stats::runif(length(vals), min = i - 0.14, max = i + 0.14)
    grDevices::adjustcolor(fill_cols[[i]], alpha.f = point_alpha)
    graphics::points(
      xj,
      vals,
      pch = 16,
      cex = 0.8,
      col = grDevices::adjustcolor(fill_cols[[i]], alpha.f = point_alpha)
    )
  }
  invisible(df)
}

plot_qap_portfolio_rank_distribution <- function(res,
                                                 out_path = NULL,
                                                 main = "Portfolio Rank Distribution",
                                                 ylab = "Average Rank Across Instances",
                                                 las = 2) {
  if (is.null(res$rank_summary$by_portfolio) || nrow(res$rank_summary$by_portfolio) == 0L) {
    stop("res$rank_summary$by_portfolio is empty.")
  }
  df <- res$rank_summary$by_portfolio
  ord <- order(df$average_rank, df$method, df$portfolio_name)
  df <- df[ord, , drop = FALSE]
  cols <- c(bestk2 = "#E69F00", portsize2 = "#0072B2", hydra2 = "#009E73")
  point_cols <- cols[df$method]
  point_cols[is.na(point_cols)] <- "grey40"

  if (!is.null(out_path) && nzchar(trimws(out_path))) {
    grDevices::pdf(out_path, width = 10, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  op <- graphics::par(mar = c(9, 5, 3, 1) + 0.1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::plot(
    seq_len(nrow(df)),
    df$average_rank,
    pch = 16,
    col = point_cols,
    xaxt = "n",
    xlab = "",
    ylab = ylab,
    main = main
  )
  graphics::axis(
    side = 1,
    at = seq_len(nrow(df)),
    labels = df$portfolio_name,
    las = las,
    cex.axis = 0.7
  )
  graphics::legend(
    "topright",
    legend = unique(df$method),
    col = cols[unique(df$method)],
    pch = 16,
    bty = "n"
  )
  invisible(df)
}

compute_qap_method_instance_score <- function(res, aggregate_fun = c("mean", "median")) {
  aggregate_fun <- match.arg(aggregate_fun)
  if (is.null(res$instance_scores_wide) || nrow(res$instance_scores_wide) == 0L) {
    stop("res$instance_scores_wide is empty.")
  }
  if (is.null(res$source_metadata) || nrow(res$source_metadata) == 0L) {
    stop("res$source_metadata is empty.")
  }

  detail_names <- names(res$details)
  score_cols <- make.names(detail_names, unique = TRUE)
  score_df <- res$instance_scores_wide[, c("instance", score_cols), drop = FALSE]
  portfolio_map <- data.frame(
    portfolio_name = detail_names,
    score_col = score_cols,
    stringsAsFactors = FALSE
  )
  portfolio_map <- merge(
    portfolio_map,
    unique(res$source_metadata[, c("portfolio_name", "method")]),
    by = "portfolio_name",
    all.x = TRUE,
    sort = FALSE
  )

  long_rows <- lapply(seq_len(nrow(portfolio_map)), function(i) {
    data.frame(
      instance = score_df$instance,
      portfolio_name = portfolio_map$portfolio_name[[i]],
      method = portfolio_map$method[[i]],
      performance_value = as.numeric(score_df[[portfolio_map$score_col[[i]]]]),
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_rows)

  fun <- if (aggregate_fun == "mean") mean else median
  by_instance_method <- stats::aggregate(
    performance_value ~ instance + method,
    data = long_df,
    FUN = fun
  )
  rownames(by_instance_method) <- NULL

  by_method <- stats::aggregate(
    performance_value ~ method,
    data = by_instance_method,
    FUN = mean
  )
  names(by_method)[names(by_method) == "performance_value"] <- "average_performance_value"
  by_method$n_instances <- as.integer(vapply(split(by_instance_method$instance, by_instance_method$method), length, integer(1)))
  rownames(by_method) <- NULL
  by_method <- by_method[order(by_method$average_performance_value), , drop = FALSE]

  list(
    by_instance_method = by_instance_method,
    by_method = by_method,
    by_portfolio = long_df,
    aggregate_fun = aggregate_fun
  )
}

plot_qap_method_score_distribution <- function(res,
                                               aggregate_fun = c("mean", "median"),
                                               out_path = NULL,
                                               main = NULL,
                                               ylab = "Per-instance Aggregated Performance Value",
                                               point_alpha = 0.55) {
  aggregate_fun <- match.arg(aggregate_fun)
  score_summary <- compute_qap_method_instance_score(res, aggregate_fun = aggregate_fun)
  df <- score_summary$by_instance_method
  method_order <- score_summary$by_method$method
  df$method <- factor(df$method, levels = method_order)
  cols <- c(bestk2 = "#E69F00", portsize2 = "#0072B2", hydra2 = "#009E73")
  fill_cols <- cols[as.character(method_order)]
  fill_cols[is.na(fill_cols)] <- "grey70"

  if (is.null(main)) {
    main <- sprintf("Method Performance Distribution (%s)", aggregate_fun)
  }
  if (!is.null(out_path) && nzchar(trimws(out_path))) {
    grDevices::pdf(out_path, width = 7, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  op <- graphics::par(mar = c(5, 5, 3, 1) + 0.1)
  on.exit(graphics::par(op), add = TRUE)
  graphics::boxplot(
    performance_value ~ method,
    data = df,
    col = fill_cols,
    border = "grey25",
    outline = FALSE,
    main = main,
    xlab = "",
    ylab = ylab
  )
  for (i in seq_along(method_order)) {
    vals <- df$performance_value[df$method == method_order[[i]]]
    xj <- stats::runif(length(vals), min = i - 0.14, max = i + 0.14)
    graphics::points(
      xj,
      vals,
      pch = 16,
      cex = 0.8,
      col = grDevices::adjustcolor(fill_cols[[i]], alpha.f = point_alpha)
    )
  }
  invisible(score_summary)
}

plot_qap_method_score_ecdf <- function(res,
                                       aggregate_fun = c("mean", "median"),
                                       out_path = NULL,
                                       main = NULL,
                                       xlab = "Per-instance Aggregated Performance Value",
                                       ylab = "ECDF") {
  aggregate_fun <- match.arg(aggregate_fun)
  score_summary <- compute_qap_method_instance_score(res, aggregate_fun = aggregate_fun)
  df <- score_summary$by_instance_method
  method_order <- score_summary$by_method$method
  cols <- c(bestk2 = "#E69F00", portsize2 = "#0072B2", hydra2 = "#009E73")
  line_cols <- cols[as.character(method_order)]
  line_cols[is.na(line_cols)] <- "grey40"

  if (is.null(main)) {
    main <- sprintf("Method Performance ECDF (%s)", aggregate_fun)
  }
  if (!is.null(out_path) && nzchar(trimws(out_path))) {
    grDevices::pdf(out_path, width = 7, height = 5)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  x_range <- range(df$performance_value, finite = TRUE)
  graphics::plot(
    NA,
    xlim = x_range,
    ylim = c(0, 1),
    xlab = xlab,
    ylab = ylab,
    main = main
  )
  for (i in seq_along(method_order)) {
    method_i <- method_order[[i]]
    vals <- df$performance_value[df$method == method_i]
    stats::plot.ecdf(vals, add = TRUE, verticals = TRUE, do.points = FALSE, col = line_cols[[i]], lwd = 2)
  }
  graphics::legend(
    "bottomright",
    legend = method_order,
    col = line_cols,
    lwd = 2,
    bty = "n"
  )
  invisible(score_summary)
}

make_portfolio_eval_label <- function(result_path, row_df, row_index) {
  tag <- tools::file_path_sans_ext(basename(result_path))
  sprintf(
    "%s|seed=%s|budget=%s|port=%s|row=%d",
    tag,
    as.integer(row_df$algo_seed[[1L]]),
    as.integer(row_df$budget[[1L]]),
    as.integer(row_df$port_size[[1L]]),
    as.integer(row_index)
  )
}

extract_qap_portfolios_from_result <- function(result_path,
                                               portfolio_col = c("best_k_union_decoded", "best_name_decoded"),
                                               rows = NULL,
                                               drop_empty = TRUE) {
  portfolio_col <- match.arg(portfolio_col)
  result_path <- normalizePath(result_path, winslash = "/", mustWork = TRUE)
  obj <- readRDS(result_path)
  if (is.null(obj$summary) || !is.data.frame(obj$summary)) {
    stop("Result file does not contain a summary data.frame: ", result_path)
  }
  summary_df <- as.data.frame(obj$summary, stringsAsFactors = FALSE)
  if (!portfolio_col %in% names(summary_df)) {
    stop("Column missing in summary: ", portfolio_col)
  }
  if (is.null(rows)) {
    rows <- seq_len(nrow(summary_df))
  }
  rows <- as.integer(rows)
  rows <- rows[rows >= 1L & rows <= nrow(summary_df)]
  if (length(rows) == 0L) {
    stop("No valid rows selected for: ", result_path)
  }

  portfolios <- list()
  meta_rows <- vector("list", length(rows))
  out_idx <- 1L
  for (i in seq_along(rows)) {
    row_idx <- rows[[i]]
    row_df <- summary_df[row_idx, , drop = FALSE]
    decoded <- as.character(row_df[[portfolio_col]][[1L]])
    if (is.na(decoded) || !nzchar(trimws(decoded)) || identical(decoded, "<none>")) {
      if (isTRUE(drop_empty)) next
      stop("Selected row has empty decoded portfolio: ", result_path, " row ", row_idx)
    }
    cfg_df <- parse_qap_decoded_portfolio(decoded)
    label <- make_portfolio_eval_label(result_path, row_df, row_idx)
    portfolios[[label]] <- cfg_df
    row_df$result_path <- result_path
    row_df$result_row <- as.integer(row_idx)
    row_df$portfolio_name <- label
    meta_rows[[out_idx]] <- row_df
    out_idx <- out_idx + 1L
  }
  meta_rows <- Filter(Negate(is.null), meta_rows)
  meta_df <- if (length(meta_rows) > 0L) do.call(rbind, meta_rows) else data.frame()
  list(
    portfolios = portfolios,
    metadata = meta_df
  )
}

extract_qap_portfolios_from_results <- function(result_paths,
                                                portfolio_col = c("best_k_union_decoded", "best_name_decoded"),
                                                rows_by_result = NULL,
                                                drop_empty = TRUE) {
  portfolio_col <- match.arg(portfolio_col)
  result_paths <- as.character(result_paths)
  result_paths <- result_paths[nzchar(trimws(result_paths))]
  if (length(result_paths) == 0L) {
    stop("result_paths must contain at least one path.")
  }

  portfolios <- list()
  meta_rows <- list()
  for (i in seq_along(result_paths)) {
    rows_i <- NULL
    if (!is.null(rows_by_result)) {
      if (is.list(rows_by_result) && length(rows_by_result) >= i) {
        rows_i <- rows_by_result[[i]]
      } else if (!is.list(rows_by_result)) {
        rows_i <- rows_by_result
      }
    }
    one <- extract_qap_portfolios_from_result(
      result_path = result_paths[[i]],
      portfolio_col = portfolio_col,
      rows = rows_i,
      drop_empty = drop_empty
    )
    portfolios <- c(portfolios, one$portfolios)
    if (nrow(one$metadata) > 0L) {
      meta_rows[[length(meta_rows) + 1L]] <- one$metadata
    }
  }
  list(
    portfolios = portfolios,
    metadata = if (length(meta_rows) > 0L) do.call(rbind, meta_rows) else data.frame()
  )
}

evaluate_qap_portfolio_list_test_instances <- function(portfolios,
                                                       split_seed = 11L,
                                                       families = c("RandomRandom", "GridStructured"),
                                                       test_splits = NULL,
                                                       time_limit = 5,
                                                       out_prefix = NULL,
                                                       verbose = TRUE,
                                                       save_outputs = TRUE) {
  verbose_mode <- qap_verbose_mode(verbose)
  portfolio_list <- normalize_qap_portfolio_set(portfolios)
  if (is.null(test_splits) || length(test_splits) == 0L) {
    test_splits <- resolve_group_compare_test_splits(
      families = families,
      split_seed = split_seed,
      script_dir = qap_scripts_dir()
    )
  } else {
    test_splits <- normalizePath(as.character(test_splits), winslash = "/", mustWork = TRUE)
  }

  instances <- read_instances_from_splits(test_splits)
  pooled_cfg_df <- unique(do.call(rbind, lapply(portfolio_list, function(df) {
    df[, c("a", "p", "l", "duplicate_tick"), drop = FALSE]
  })))
  rownames(pooled_cfg_df) <- NULL

  if (verbose_mode != "none") {
    cat(sprintf(
      "[result-portfolio-test-eval] pooled unique configs: %d across %d portfolios on %d instances\n",
      nrow(pooled_cfg_df), length(portfolio_list), length(instances)
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
        "[result-portfolio-test-eval] %d/%d | %s | size=%d\n",
        i, length(portfolio_list), name_i, nrow(cfg_df)
      ))
    }
    instance_scores <- score_qap_portfolio_from_grid(
      config_set = cfg_df,
      pooled_scores = pooled_scores,
      verbose = "none"
    )
    summary_df <- data.frame(
      portfolio_name = name_i,
      split_seed = as.integer(split_seed),
      families = paste(families, collapse = ","),
      portfolio_size = nrow(cfg_df),
      portfolio_label = paste(apply(cfg_df, 1, function(z) sprintf(
        "(a=%d,p=%d,l=%d,dup=%d)", z[["a"]], z[["p"]], z[["l"]], z[["duplicate_tick"]]
      )), collapse = "; "),
      time_limit = as.numeric(time_limit),
      test_instances = nrow(instance_scores),
      test_mean_cost = mean(instance_scores$portfolio_cost),
      test_median_cost = median(instance_scores$portfolio_cost),
      test_sd_cost = if (nrow(instance_scores) > 1L) stats::sd(instance_scores$portfolio_cost) else NA_real_,
      test_min_cost = min(instance_scores$portfolio_cost),
      test_max_cost = max(instance_scores$portfolio_cost),
      stringsAsFactors = FALSE
    )
    list(summary = summary_df, configs = cfg_df, instance_scores = instance_scores)
  })

  summary_df <- do.call(rbind, lapply(result_list, function(x) x$summary))
  rownames(summary_df) <- NULL

  detail_list <- setNames(lapply(result_list, function(x) x$instance_scores), summary_df$portfolio_name)
  wide_cost_df <- data.frame(instance = instances, stringsAsFactors = FALSE)
  for (nm in names(detail_list)) {
    safe_nm <- make.names(nm, unique = TRUE)
    wide_cost_df[[safe_nm]] <- detail_list[[nm]]$portfolio_cost
  }
  cost_cols <- setdiff(names(wide_cost_df), "instance")
  cost_mat <- as.matrix(wide_cost_df[, cost_cols, drop = FALSE])
  best_idx <- max.col(-cost_mat, ties.method = "first")
  wide_cost_df$min_portfolio_cost <- vapply(seq_len(nrow(cost_mat)), function(i) min(cost_mat[i, ]), numeric(1))
  wide_cost_df$best_portfolio_name <- names(detail_list)[best_idx]

  out <- list(
    summary = summary_df,
    portfolios = lapply(result_list, function(x) x$configs),
    details = detail_list,
    pooled_configs = pooled_cfg_df,
    instance_scores_wide = wide_cost_df,
    test_splits = test_splits
  )

  if (isTRUE(save_outputs) && !is.null(out_prefix) && nzchar(trimws(out_prefix))) {
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(summary_df, paste0(out_prefix, "_summary.csv"), row.names = FALSE)
    utils::write.csv(wide_cost_df, paste0(out_prefix, "_instance_scores.csv"), row.names = FALSE)
    saveRDS(out, paste0(out_prefix, ".Rds"))
  }

  out
}

evaluate_qap_result_portfolios_on_combined_group_compare_testset <- function(result_paths,
                                                                              families = c("RandomRandom", "GridStructured"),
                                                                              split_seed = 11L,
                                                                              portfolio_col = c("best_k_union_decoded", "best_name_decoded"),
                                                                              rows_by_result = NULL,
                                                                              time_limit = 5,
                                                                              out_prefix = NULL,
                                                                              verbose = TRUE,
                                                                              save_outputs = TRUE,
                                                                              drop_empty = TRUE) {
  portfolio_col <- match.arg(portfolio_col)
  extracted <- extract_qap_portfolios_from_results(
    result_paths = result_paths,
    portfolio_col = portfolio_col,
    rows_by_result = rows_by_result,
    drop_empty = drop_empty
  )
  eval_res <- evaluate_qap_portfolio_list_test_instances(
    portfolios = extracted$portfolios,
    families = families,
    split_seed = split_seed,
    time_limit = time_limit,
    out_prefix = out_prefix,
    verbose = verbose,
    save_outputs = save_outputs
  )
  summary_df <- merge(
    extracted$metadata,
    eval_res$summary,
    by = "portfolio_name",
    all.y = TRUE,
    sort = FALSE
  )
  eval_res$source_metadata <- extracted$metadata
  eval_res$summary <- summary_df
  if (isTRUE(save_outputs) && !is.null(out_prefix) && nzchar(trimws(out_prefix))) {
    utils::write.csv(summary_df, paste0(out_prefix, "_summary.csv"), row.names = FALSE)
    saveRDS(eval_res, paste0(out_prefix, ".Rds"))
  }
  eval_res
}

portfolio_mode_to_col <- function(portfolio_mode = c("best_k_union", "best_k_names", "best_name")) {
  portfolio_mode <- match.arg(portfolio_mode)
  switch(
    portfolio_mode,
    best_name = "best_name_decoded",
    best_k_names = "best_k_names_decoded",
    best_k_union = "best_k_union_decoded"
  )
}

infer_method_name_from_csv <- function(csv_path) {
  nm <- tolower(tools::file_path_sans_ext(basename(csv_path)))
  if (grepl("bestk", nm)) return("bestk2")
  if (grepl("portsize2", nm)) return("portsize2")
  if (grepl("hydra2", nm)) return("hydra2")
  nm
}

normalize_portfolio_mode_by_method <- function(portfolio_mode_by_method = NULL,
                                               default_mode = "best_name") {
  if (is.null(portfolio_mode_by_method)) {
    portfolio_mode_by_method <- c(
      bestk2 = "best_k_union",
      portsize2 = default_mode,
      hydra2 = default_mode
    )
  }
  if (is.list(portfolio_mode_by_method)) {
    portfolio_mode_by_method <- unlist(portfolio_mode_by_method, use.names = TRUE)
  }
  if (is.null(names(portfolio_mode_by_method)) || any(!nzchar(names(portfolio_mode_by_method)))) {
    stop("portfolio_mode_by_method must be a named vector/list keyed by method name.")
  }
  portfolio_mode_by_method <- as.character(portfolio_mode_by_method)
  for (i in seq_along(portfolio_mode_by_method)) {
    portfolio_mode_by_method[[i]] <- match.arg(
      portfolio_mode_by_method[[i]],
      c("best_k_union", "best_k_names", "best_name")
    )
  }
  portfolio_mode_by_method
}

extract_qap_portfolios_from_csv <- function(csv_path,
                                            portfolio_col,
                                            method_name = NULL,
                                            rows = NULL,
                                            drop_empty = TRUE) {
  csv_path <- normalizePath(csv_path, winslash = "/", mustWork = TRUE)
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (!portfolio_col %in% names(df)) {
    stop("Column missing in csv: ", portfolio_col, " | ", csv_path)
  }
  if (is.null(method_name) || !nzchar(method_name)) {
    method_name <- infer_method_name_from_csv(csv_path)
  }
  if (is.null(rows)) {
    rows <- seq_len(nrow(df))
  }
  rows <- as.integer(rows)
  rows <- rows[rows >= 1L & rows <= nrow(df)]
  if (length(rows) == 0L) {
    stop("No valid rows selected for: ", csv_path)
  }

  portfolios <- list()
  meta_rows <- list()
  out_idx <- 1L
  for (row_idx in rows) {
    row_df <- df[row_idx, , drop = FALSE]
    decoded <- as.character(row_df[[portfolio_col]][[1L]])
    if (is.na(decoded) || !nzchar(trimws(decoded)) || identical(decoded, "<none>")) {
      if (isTRUE(drop_empty)) next
      stop("Selected row has empty decoded portfolio: ", csv_path, " row ", row_idx)
    }
    cfg_df <- parse_qap_decoded_portfolio(decoded)
    label <- sprintf(
      "%s|seed=%s|budget=%s|port=%s|row=%d",
      method_name,
      as.integer(row_df$algo_seed[[1L]]),
      as.integer(row_df$budget[[1L]]),
      as.integer(row_df$port_size[[1L]]),
      as.integer(row_idx)
    )
    portfolios[[label]] <- cfg_df
    row_df$source_csv <- csv_path
    row_df$result_row <- as.integer(row_idx)
    row_df$method <- method_name
    row_df$portfolio_name <- label
    meta_rows[[out_idx]] <- row_df
    out_idx <- out_idx + 1L
  }

  list(
    portfolios = portfolios,
    metadata = if (length(meta_rows) > 0L) do.call(rbind, meta_rows) else data.frame()
  )
}

evaluate_qap_result_csv_folder_on_combined_group_compare_testset <- function(results_dir,
                                                                              families = c("GridRandom", "GridStructured", "RandomStructured"),
                                                                              split_seed = 11L,
                                                                              portfolio_mode = c("best_k_union", "best_k_names", "best_name"),
                                                                              portfolio_mode_by_method = NULL,
                                                                              rows_by_csv = NULL,
                                                                              time_limit = 5,
                                                                              out_prefix = NULL,
                                                                              verbose = TRUE,
                                                                              save_outputs = TRUE,
                                                                              drop_empty = TRUE) {
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)
  portfolio_mode <- match.arg(portfolio_mode)
  portfolio_col <- portfolio_mode_to_col(portfolio_mode)
  portfolio_mode_by_method <- normalize_portfolio_mode_by_method(
    portfolio_mode_by_method = portfolio_mode_by_method,
    default_mode = portfolio_mode
  )

  csv_paths <- list.files(results_dir, pattern = "\\.csv$", full.names = TRUE)
  csv_paths <- csv_paths[!grepl("_rows\\.csv$|_by_budget\\.csv$|_summary\\.csv$|_instance_scores\\.csv$|_average_rank_by_method\\.csv$|_average_rank_by_instance_method\\.csv$|_average_rank_by_portfolio\\.csv$", csv_paths)]
  csv_paths <- sort(csv_paths)
  if (length(csv_paths) == 0L) {
    stop("No CSV result files found in: ", results_dir)
  }

  portfolios <- list()
  meta_rows <- list()
  for (i in seq_along(csv_paths)) {
    rows_i <- NULL
    if (!is.null(rows_by_csv)) {
      if (is.list(rows_by_csv) && length(rows_by_csv) >= i) {
        rows_i <- rows_by_csv[[i]]
      } else if (!is.list(rows_by_csv)) {
        rows_i <- rows_by_csv
      }
    }
    one <- extract_qap_portfolios_from_csv(
      csv_path = csv_paths[[i]],
      portfolio_col = {
        method_i <- infer_method_name_from_csv(csv_paths[[i]])
        mode_i <- if (method_i %in% names(portfolio_mode_by_method)) {
          portfolio_mode_by_method[[method_i]]
        } else {
          portfolio_mode
        }
        portfolio_mode_to_col(mode_i)
      },
      method_name = infer_method_name_from_csv(csv_paths[[i]]),
      rows = rows_i,
      drop_empty = drop_empty
    )
    portfolios <- c(portfolios, one$portfolios)
    if (nrow(one$metadata) > 0L) {
      meta_rows[[length(meta_rows) + 1L]] <- one$metadata
    }
  }

  metadata <- if (length(meta_rows) > 0L) do.call(rbind, meta_rows) else data.frame()
  eval_res <- evaluate_qap_portfolio_list_test_instances(
    portfolios = portfolios,
    families = families,
    split_seed = split_seed,
    time_limit = time_limit,
    out_prefix = out_prefix,
    verbose = verbose,
    save_outputs = save_outputs
  )
  row_df <- merge(
    metadata,
    eval_res$summary,
    by = "portfolio_name",
    all.y = TRUE,
    sort = FALSE
  )
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

  out <- list(
    row_results = row_df,
    budget_summary = budget_df,
    instance_scores_wide = eval_res$instance_scores_wide,
    details = eval_res$details,
    portfolios = eval_res$portfolios,
    source_metadata = metadata,
    test_splits = eval_res$test_splits,
    portfolio_mode = portfolio_mode,
    portfolio_col = portfolio_col,
    portfolio_mode_by_method = portfolio_mode_by_method
  )

  rank_summary <- compute_qap_average_rank_by_method(
    instance_scores_wide = eval_res$instance_scores_wide,
    source_metadata = metadata,
    detail_names = names(eval_res$details)
  )
  out$rank_summary <- rank_summary

  if (isTRUE(save_outputs) && !is.null(out_prefix) && nzchar(trimws(out_prefix))) {
    utils::write.csv(row_df, paste0(out_prefix, "_rows.csv"), row.names = FALSE)
    utils::write.csv(budget_df, paste0(out_prefix, "_by_budget.csv"), row.names = FALSE)
    utils::write.csv(eval_res$instance_scores_wide, paste0(out_prefix, "_instance_scores.csv"), row.names = FALSE)
    utils::write.csv(rank_summary$by_method, paste0(out_prefix, "_average_rank_by_method.csv"), row.names = FALSE)
    utils::write.csv(rank_summary$by_instance_method, paste0(out_prefix, "_average_rank_by_instance_method.csv"), row.names = FALSE)
    utils::write.csv(rank_summary$by_portfolio, paste0(out_prefix, "_average_rank_by_portfolio.csv"), row.names = FALSE)
    saveRDS(out, paste0(out_prefix, ".Rds"))
  }

  out
}
