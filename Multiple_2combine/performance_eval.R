# Monte Carlo cost evaluation for portfolio-irace (combine2).
# - Same seed structure as scoring.R: base_seed_i = seed0 + 10000 * i
# - One instance random term per instance (drawn first from seed stream)
# - Config/member noise drawn from the same stream after instance random
# duplicate_independent:
#   TRUE  -> each duplicate member gets its own epsilon draw (for duplicate-benefit analysis)
#   FALSE -> duplicates of same value share one epsilon (strictly aligned with current scoring.R behavior)

parse_port_name = function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(numeric(0))
  s = gsub("[() ]", "", x)
  if (!nzchar(s)) return(numeric(0))
  as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
}

format_port_name_keep_dup = function(vals) {
  vals = as.numeric(vals)
  if (length(vals) == 0L || any(!is.finite(vals))) return(NA_character_)
  vals = sort(vals)
  paste0("(", paste(vals, collapse = ","), ")")
}

make_mc_instance_ids = function(n_mc, instance_start = 1L, even_weight = 2L, odd_weight = 1L) {
  n_mc = as.integer(n_mc)
  if (n_mc <= 0L) return(integer(0))
  even_weight = max(0L, as.integer(even_weight))
  odd_weight = max(0L, as.integer(odd_weight))
  total_weight = even_weight + odd_weight
  if (total_weight <= 0L) stop("At least one of even_weight or odd_weight must be positive.")

  n_even = as.integer(round(n_mc * even_weight / total_weight))
  n_even = max(0L, min(n_mc, n_even))
  n_odd = as.integer(n_mc - n_even)

  start_even = if ((as.integer(instance_start) %% 2L) == 0L) as.integer(instance_start) else as.integer(instance_start + 1L)
  start_odd = if ((as.integer(instance_start) %% 2L) == 1L) as.integer(instance_start) else as.integer(instance_start + 1L)

  even_ids = if (n_even > 0L) seq.int(from = start_even, by = 2L, length.out = n_even) else integer(0)
  odd_ids = if (n_odd > 0L) seq.int(from = start_odd, by = 2L, length.out = n_odd) else integer(0)
  c(even_ids, odd_ids)
}

score_config_one = function(c, i, w = 1, p = 1, instance_random = 0, eps = 0) {
  c_int = as.integer(round(c))
  i_int = as.integer(i)
  penalty = if ((c_int %% 2L) == (i_int %% 2L)) 0 else p
  instance_random + w * abs(c - 50) + penalty + eps
}

score_portfolio_on_instance_mc = function(port_vals,
                                          i,
                                          base_seed,
                                          w = 1,
                                          sigma1 = 0,
                                          sigma2 = 0.1,
                                          p = 10,
                                          duplicate_independent = TRUE,
                                          sort_port = TRUE) {
  port_vals = as.numeric(port_vals)
  if (length(port_vals) == 0L) return(NA_real_)
  if (isTRUE(sort_port)) port_vals = sort(port_vals)

  set.seed(as.integer(base_seed))
  instance_random = rnorm(1, mean = 0, sd = sigma1)

  if (isTRUE(duplicate_independent)) {
    eps = rnorm(length(port_vals), mean = 0, sd = sigma2)
    sc = vapply(seq_along(port_vals), function(k) {
      score_config_one(
        c = port_vals[k], i = i, w = w, p = p,
        instance_random = instance_random, eps = eps[k]
      )
    }, numeric(1))
  } else {
    # Strictly mimic scoring.R: one epsilon per unique config value, then map back.
    cfg_u = sort(unique(port_vals))
    eps_u = rnorm(length(cfg_u), mean = 0, sd = sigma2)
    sc_u = vapply(seq_along(cfg_u), function(k) {
      score_config_one(
        c = cfg_u[k], i = i, w = w, p = p,
        instance_random = instance_random, eps = eps_u[k]
      )
    }, numeric(1))
    sc = sc_u[match(port_vals, cfg_u)]
  }

  min(sc)
}

mc_expected_portfolio_cost = function(port,
                                      n_mc = 5000L,
                                      seed0 = 1L,
                                      instance_start = 1L,
                                      even_weight = 2L,
                                      odd_weight = 1L,
                                      w = 1,
                                      sigma1 = 0,
                                      sigma2 = 0.1,
                                      p = 10,
                                      duplicate_independent = TRUE,
                                      sort_port = TRUE,
                                      return_samples = FALSE) {
  port_vals = if (is.character(port)) parse_port_name(port) else as.numeric(port)
  if (length(port_vals) == 0L || any(!is.finite(port_vals))) {
    stop("Invalid portfolio input. Use numeric vector or '(a,b,...)' string.")
  }

  n_mc = as.integer(n_mc)
  if (n_mc <= 0L) stop("n_mc must be >= 1")

  inst = make_mc_instance_ids(
    n_mc = n_mc,
    instance_start = instance_start,
    even_weight = even_weight,
    odd_weight = odd_weight
  )
  base_seed_vec = as.integer(seed0) + 10000L * inst

  samples = vapply(seq_len(n_mc), function(t) {
    score_portfolio_on_instance_mc(
      port_vals = port_vals,
      i = inst[t],
      base_seed = base_seed_vec[t],
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      duplicate_independent = duplicate_independent,
      sort_port = sort_port
    )
  }, numeric(1))

  est_mean = mean(samples)
  est_sd = stats::sd(samples)
  est_se = est_sd / sqrt(n_mc)

  out = list(
    portfolio = format_port_name_keep_dup(port_vals),
    n_mc = n_mc,
    mean_cost = est_mean,
    sd_cost = est_sd,
    se_cost = est_se,
    ci95_low = est_mean - 1.96 * est_se,
    ci95_high = est_mean + 1.96 * est_se,
    duplicate_independent = isTRUE(duplicate_independent)
  )
  if (isTRUE(return_samples)) out$samples = samples
  out
}

mc_eval_portfolios_cost = function(portfolios,
                                   n_mc = 5000L,
                                   seed0 = 1L,
                                   instance_start = 1L,
                                   even_weight = 2L,
                                   odd_weight = 1L,
                                   w = 1,
                                   sigma1 = 0,
                                   sigma2 = 0.1,
                                   p = 10,
                                   duplicate_independent = TRUE,
                                   sort_port = TRUE) {
  ports = portfolios
  res = lapply(seq_along(ports), function(i) {
    est = mc_expected_portfolio_cost(
      port = ports[[i]],
      n_mc = n_mc,
      seed0 = seed0,
      instance_start = instance_start,
      even_weight = even_weight,
      odd_weight = odd_weight,
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      duplicate_independent = duplicate_independent,
      sort_port = sort_port,
      return_samples = FALSE
    )
    data.frame(
      portfolio_in = as.character(ports[[i]]),
      portfolio = est$portfolio,
      n_mc = est$n_mc,
      mean_cost = est$mean_cost,
      sd_cost = est$sd_cost,
      se_cost = est$se_cost,
      ci95_low = est$ci95_low,
      ci95_high = est$ci95_high,
      duplicate_independent = est$duplicate_independent,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

# Convenience for res_grid-like outputs: evaluate only best_name cost (no gap).
# Required columns: best_name, w, sigma1, sigma2, p, seed_shift
mc_eval_res_grid_best_cost = function(df,
                                      n_mc = 5000L,
                                      seed0_base = 1L,
                                      instance_start = 1L,
                                      even_weight = 2L,
                                      odd_weight = 1L,
                                      best_col = "best_name",
                                      duplicate_independent = TRUE) {
  need_cols = c(best_col, "w", "sigma1", "sigma2", "p", "seed_shift")
  miss = setdiff(need_cols, names(df))
  if (length(miss) > 0L) stop("Missing columns: ", paste(miss, collapse = ", "))

  out = lapply(seq_len(nrow(df)), function(i) {
    est = mc_expected_portfolio_cost(
      port = as.character(df[[best_col]][i]),
      n_mc = n_mc,
      seed0 = as.integer(seed0_base) + as.integer(df$seed_shift[i]),
      instance_start = instance_start,
      even_weight = even_weight,
      odd_weight = odd_weight,
      w = as.numeric(df$w[i]),
      sigma1 = as.numeric(df$sigma1[i]),
      sigma2 = as.numeric(df$sigma2[i]),
      p = as.numeric(df$p[i]),
      duplicate_independent = duplicate_independent,
      return_samples = FALSE
    )
    data.frame(
      mean_best_cost_mc = est$mean_cost,
      sd_best_cost_mc = est$sd_cost,
      se_best_cost_mc = est$se_cost,
      ci95_low_best_cost_mc = est$ci95_low,
      ci95_high_best_cost_mc = est$ci95_high,
      stringsAsFactors = FALSE
    )
  })

  cbind(df, do.call(rbind, out))
}

# Wrapper aligned with old run_overall logic:
# - k_best == 1: evaluate best_name
# - k_best > 1 : evaluate best_k_union (union portfolio of top-k names)
mc_eval_res_grid_cost = function(df,
                                 k_best = 1L,
                                 n_mc = 5000L,
                                 seed0_base = 1L,
                                 instance_start = 1L,
                                 even_weight = 2L,
                                 odd_weight = 1L,
                                 duplicate_independent = TRUE,
                                 port_col = NULL) {
  k_best = as.integer(k_best)
  if (is.null(port_col)) {
    if (k_best <= 1L) {
      port_col = "best_name"
    } else {
      if ("best_k_union" %in% names(df)) {
        port_col = "best_k_union"
      } else if ("best_k_names" %in% names(df)) {
        # fallback: build union from "a;b;c" style names
        df$best_k_union = vapply(as.character(df$best_k_names), function(s) {
          if (is.na(s) || !nzchar(s)) return(NA_character_)
          parts = strsplit(s, ";", fixed = TRUE)[[1]]
          vals = sort(unique(unlist(lapply(parts, parse_port_name))))
          format_port_name_keep_dup(vals)
        }, character(1))
        port_col = "best_k_union"
      } else {
        stop("k_best > 1 but neither best_k_union nor best_k_names exists in df.")
      }
    }
  }
  mc_eval_res_grid_best_cost(
    df = df,
    n_mc = n_mc,
    seed0_base = seed0_base,
    instance_start = instance_start,
    even_weight = even_weight,
    odd_weight = odd_weight,
    best_col = port_col,
    duplicate_independent = duplicate_independent
  )
}

# Summarize evaluated res_grid by budget.
# Expects `cost_col` to exist in df_eval (default from mc_eval_res_grid_best_cost).
summarize_cost_by_budget = function(df_eval,
                                    budget_col = "B_total",
                                    cost_col = "mean_best_cost_mc") {
  if (!(budget_col %in% names(df_eval))) stop("Missing budget column: ", budget_col)
  if (!(cost_col %in% names(df_eval))) stop("Missing cost column: ", cost_col)
  b = as.numeric(df_eval[[budget_col]])
  y = as.numeric(df_eval[[cost_col]])
  ok = is.finite(b) & is.finite(y)
  b = b[ok]
  y = y[ok]
  sp = split(y, b)
  out = do.call(rbind, lapply(names(sp), function(k) {
    z = as.numeric(sp[[k]])
    n = length(z)
    s = stats::sd(z)
    se = if (n > 0) s / sqrt(n) else NA_real_
    data.frame(
      B_total = as.numeric(k),
      n = n,
      mean = mean(z),
      median = stats::median(z),
      min = min(z),
      max = max(z),
      sd = s,
      se = se,
      ci95_low = mean(z) - 1.96 * se,
      ci95_high = mean(z) + 1.96 * se,
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$B_total), , drop = FALSE]
}

# Plot budget-performance curve:
# - min/max band (light)
# - mean +/- SE error bars
# - mean and median lines
plot_cost_by_budget = function(summary_df,
                               x_col = "B_total",
                               mean_col = "mean",
                               median_col = "median",
                               min_col = "min",
                               max_col = "max",
                               se_col = "se",
                               ylab = "MC expected cost",
                               xlab = "Budget (B_total)",
                               main = "Portfolio Performance vs Budget",
                               ylim = NULL,
                               ylim_pad_frac = 0.04) {
  x = as.numeric(summary_df[[x_col]])
  m = as.numeric(summary_df[[mean_col]])
  md = as.numeric(summary_df[[median_col]])
  lo = as.numeric(summary_df[[min_col]])
  hi = as.numeric(summary_df[[max_col]])
  se = as.numeric(summary_df[[se_col]])
  
  ord = order(x)
  x = x[ord]; m = m[ord]; md = md[ord]; lo = lo[ord]; hi = hi[ord]; se = se[ord]
  
  if (is.null(ylim)) {
    y_all = c(lo, hi, md, m - se, m + se)
    y_all = y_all[is.finite(y_all)]
    if (length(y_all) == 0L) {
      ylim = c(0, 1)
    } else {
      yr = range(y_all, na.rm = TRUE)
      span = yr[2] - yr[1]
      pad = if (is.finite(span) && span > 0) span * max(0, as.numeric(ylim_pad_frac)) else 0.5
      ylim = c(yr[1] - pad, yr[2] + pad)
    }
  }
  
  xl = range(x, na.rm = TRUE)
  if (!all(is.finite(xl))) {
    xl = c(0, 1)
  } else if (xl[1] == xl[2]) {
    pad_x = max(1, abs(xl[1]) * 0.05)
    xl = c(xl[1] - pad_x, xl[2] + pad_x)
  }
  
  plot(x, m, type = "n", xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xl)
  
  # min-max ribbon
  polygon(
    x = c(x, rev(x)),
    y = c(lo, rev(hi)),
    col = grDevices::adjustcolor("gray60", alpha.f = 0.20),
    border = NA
  )
  
  # mean +/- SE error bars
  arrows(
    x0 = x, y0 = m - se,
    x1 = x, y1 = m + se,
    angle = 90, code = 3, length = 0.04,
    col = "steelblue3", lwd = 1.2
  )
  
  lines(x, m, type = "b", pch = 16, col = "steelblue4", lwd = 2)
  lines(x, md, type = "b", pch = 17, col = "firebrick3", lwd = 1.5, lty = 2)
  
  legend(
    "topright",
    legend = c("Mean", "Median", "Mean +/- SE", "Min-Max range"),
    col = c("steelblue4", "firebrick3", "steelblue3", "gray60"),
    lty = c(1, 2, 1, 1),
    pch = c(16, 17, NA, 15),
    pt.cex = c(1, 1, NA, 1.5),
    bty = "n"
  )
}

mc_plot_res_grid_cost_by_budget = function(res_grid,
                                           k_best = 1L,
                                           n_mc = 5000L,
                                           seed0_base = 1L,
                                           instance_start = 1L,
                                           even_weight = 2L,
                                           odd_weight = 1L,
                                           duplicate_independent = TRUE,
                                           budget_col = "B_total",
                                           port_col = NULL,
                                           ylab = "MCMC expected cost",
                                           xlab = "Budget (B_total)",
                                           main = "MCMC Cost vs Budget",
                                           ylim = NULL,
                                           ylim_pad_frac = 0.04,
                                           return_data = TRUE) {
  df_eval = mc_eval_res_grid_cost(
    df = res_grid,
    k_best = k_best,
    n_mc = n_mc,
    seed0_base = seed0_base,
    instance_start = instance_start,
    even_weight = even_weight,
    odd_weight = odd_weight,
    duplicate_independent = duplicate_independent,
    port_col = port_col
  )
  sum_budget = summarize_cost_by_budget(
    df_eval = df_eval,
    budget_col = budget_col,
    cost_col = "mean_best_cost_mc"
  )
  plot_cost_by_budget(
    summary_df = sum_budget,
    ylab = ylab,
    xlab = xlab,
    main = main,
    ylim = ylim,
    ylim_pad_frac = ylim_pad_frac
  )
  if (isTRUE(return_data)) {
    return(list(df_eval = df_eval, summary_by_budget = sum_budget))
  }
  invisible(NULL)
}
