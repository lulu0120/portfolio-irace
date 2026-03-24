#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
  library(truncnorm)
}))

script_dir <- local({
  find_scripts_dir <- function(start) {
    if (!is.character(start) || length(start) == 0L || is.na(start) || !nzchar(start)) {
      return(NA_character_)
    }
    current <- normalizePath(start[[1L]], winslash = "/", mustWork = TRUE)
    if (!dir.exists(current)) {
      current <- dirname(current)
    }
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

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    editor_path <- if (!is.null(ctx)) ctx$path else ""
    if (is.character(editor_path) && length(editor_path) > 0L && nzchar(editor_path) && file.exists(editor_path)) {
      resolved <- find_scripts_dir(editor_path[[1L]])
      if (is.character(resolved) && !is.na(resolved)) {
        return(resolved)
      }
    }
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
        stop("Cannot resolve QAP_problem/Scripts. In RStudio, source qap_run_overall.R or set working directory inside the repo.")
      }
    }
  }
})
source(file.path(script_dir, "qap_common.R"), local = .GlobalEnv)
script_dir <- qap_scripts_dir()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

performance_eval_path <- file.path(repo_root, "Multiple_2combine", "performance_eval.R")
if (file.exists(performance_eval_path)) {
  source(performance_eval_path, local = .GlobalEnv)
}

load_qap_portfolio_modules <- function(repo_root) {
  src_dir <- file.path(repo_root, "Multiple_multidim")
  files <- c(
    "initial_and_env.R",
    "statistical_test.R",
    "elites_functions.R",
    "generate_discrete_input.R",
    "within_race.R",
    "wrapper_new.R"
  )
  invisible(lapply(file.path(src_dir, files), function(f) source(f, local = .GlobalEnv)))
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

make_qap_param_config <- function(a_values = 0:7, p_values = 0:9, l_values = 0:5) {
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

get_qap_split_paths <- function(split_seed,
                                script_dir = script_dir) {
  split_dir <- file.path(script_dir, "splits", paste0("seed_", as.integer(split_seed)))
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

parse_qap_port_name <- function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(numeric(0))
  s <- gsub("[() ]", "", x)
  if (!nzchar(s)) return(numeric(0))
  as.numeric(strsplit(s, ",", fixed = TRUE)[[1L]])
}

format_qap_port_name <- function(vals) {
  vals <- sort(unique(as.numeric(vals)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_character_)
  paste0("(", paste(as.integer(vals), collapse = ","), ")")
}

format_qap_port_name_keep_dup <- function(vals) {
  vals <- sort(as.numeric(vals))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_character_)
  paste0("(", paste(as.integer(vals), collapse = ","), ")")
}

get_top_k_qap_names <- function(last_race, k) {
  if (is.null(last_race) || is.null(last_race$alive_names) || length(last_race$alive_names) == 0L) {
    return(character(0))
  }
  k <- max(1L, as.integer(k))
  names0 <- as.character(last_race$alive_names)
  ranks0 <- suppressWarnings(as.numeric(last_race$alive_R))
  if (length(ranks0) == length(names0) && any(is.finite(ranks0))) {
    ord <- order(ranks0, seq_along(names0), na.last = TRUE)
    names0 <- names0[ord]
  }
  names0[seq_len(min(k, length(names0)))]
}

qap_port_name_from_key <- function(state, key) {
  if (!is.null(state) &&
      !is.null(state$ports_all_env) &&
      exists(key, envir = state$ports_all_env, inherits = FALSE)) {
    ticks <- as.numeric(get(key, envir = state$ports_all_env, inherits = FALSE))
    return(format_qap_port_name_keep_dup(cfg_from_tick(ticks)))
  }
  vals <- as.numeric(strsplit(as.character(key), ",", fixed = TRUE)[[1L]])
  format_qap_port_name_keep_dup(cfg_from_tick(vals))
}

rank_qap_history_portfolio_names <- function(last_race, exclude_names = character(0)) {
  if (is.null(last_race) || is.null(last_race$state) || is.null(last_race$state$elite_seen_env)) {
    return(character(0))
  }
  state <- last_race$state
  keys_all <- ls(envir = state$elite_seen_env, all.names = TRUE)
  if (length(keys_all) == 0L) return(character(0))

  names_all <- vapply(keys_all, function(k) qap_port_name_from_key(state, k), character(1))
  keep <- !(names_all %in% exclude_names) & !is.na(names_all) & nzchar(names_all)
  keys_all <- keys_all[keep]
  names_all <- names_all[keep]
  if (length(keys_all) == 0L) return(character(0))

  ports_all <- lapply(keys_all, function(k) parse_qap_port_name(qap_port_name_from_key(state, k)))
  names(ports_all) <- keys_all
  elite_e <- vapply(keys_all, function(k) get_port_e(state, k), integer(1))
  cand_pos <- seq_along(keys_all)
  ordered <- integer(0)
  count_levels <- sort(unique(elite_e), decreasing = TRUE)

  for (count_i in count_levels) {
    pos_i <- cand_pos[elite_e[cand_pos] == count_i]
    if (length(pos_i) <= 1L) {
      ordered <- c(ordered, pos_i)
      next
    }
    tie_keys <- keys_all[pos_i]
    tie_ports <- ports_all[pos_i]
    tie_inst_list <- lapply(tie_keys, function(k) {
      if (exists(k, envir = state$elite_seen_env, inherits = FALSE)) {
        occ_keys <- unique(as.character(get(k, envir = state$elite_seen_env, inherits = FALSE)))
        parsed <- parse_ir_key(occ_keys)
        keep_occ <- is.finite(parsed$i) & is.finite(parsed$r) & is.finite(parsed$u)
        if (!any(keep_occ)) {
          data.frame(i = integer(0), r = integer(0), u = integer(0), stringsAsFactors = FALSE)
        } else {
          parsed[keep_occ, , drop = FALSE]
        }
      } else {
        data.frame(i = integer(0), r = integer(0), u = integer(0), stringsAsFactors = FALSE)
      }
    })
    tie_inst_union <- do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x) > 0L, tie_inst_list))
    R_tie <- rep(Inf, length(pos_i))
    F_tie_mean <- rep(Inf, length(pos_i))
    if (!is.null(tie_inst_union) && nrow(tie_inst_union) > 0L) {
      tie_inst_union <- unique(tie_inst_union)
      F_rows <- vector("list", nrow(tie_inst_union))
      keep_row <- logical(nrow(tie_inst_union))
      for (ii in seq_len(nrow(tie_inst_union))) {
        i0 <- as.integer(tie_inst_union$i[ii])
        r0 <- as.integer(tie_inst_union$r[ii])
        u0 <- as.integer(tie_inst_union$u[ii])
        fvec <- rep(NA_real_, length(pos_i))
        for (pp in seq_along(tie_ports)) {
          cfgs <- as.numeric(tie_ports[[pp]])
          keys_ic <- vapply(cfgs, function(c0) make_icr_key(i0, c0, r = r0, u = u0), character(1))
          has_all <- all(vapply(keys_ic, exists, logical(1), envir = state$score_env, inherits = FALSE))
          if (has_all) {
            sc <- vapply(keys_ic, get, numeric(1), envir = state$score_env, inherits = FALSE)
            fvec[pp] <- min(as.numeric(sc))
          }
        }
        if (all(is.finite(fvec))) {
          F_rows[[ii]] <- fvec
          keep_row[ii] <- TRUE
        }
      }
      if (any(keep_row)) {
        Fmat_tie <- do.call(rbind, F_rows[keep_row])
        R_tie <- colSums2_base(rowRanks_base(Fmat_tie, cols = seq_len(ncol(Fmat_tie))))
        F_tie_mean <- colMeans(Fmat_tie)
      }
    }
    tie_min_id <- vapply(tie_ports, function(P) min(as.integer(P)), integer(1))
    ord_i <- order(R_tie, F_tie_mean, tie_min_id, names_all[pos_i])
    ordered <- c(ordered, pos_i[ord_i])
  }

  names_all[ordered]
}

rank_qap_current_race_config_names <- function(last_race, exclude_names = character(0)) {
  if (is.null(last_race) ||
      is.null(last_race$state) ||
      is.null(last_race$configs_race) ||
      length(last_race$configs_race) == 0L ||
      is.null(last_race$inst_used) ||
      length(last_race$inst_used) == 0L) {
    return(character(0))
  }
  state <- last_race$state
  cfgs <- as.numeric(last_race$configs_race)
  cfg_names <- vapply(cfgs, function(x) format_qap_port_name_keep_dup(x), character(1))
  keep_cfg <- !(cfg_names %in% exclude_names) & !duplicated(cfg_names)
  cfgs <- cfgs[keep_cfg]
  cfg_names <- cfg_names[keep_cfg]
  if (length(cfgs) == 0L) return(character(0))

  inst_now <- as.integer(last_race$inst_used)
  rep_now <- if (!is.null(last_race$inst_repeat_used)) as.integer(last_race$inst_repeat_used) else rep.int(0L, length(inst_now))
  reuse_now <- if (!is.null(last_race$inst_reuse_new_used)) as.integer(last_race$inst_reuse_new_used) else rep.int(0L, length(inst_now))
  n_occ <- min(length(inst_now), length(rep_now), length(reuse_now))
  inst_now <- inst_now[seq_len(n_occ)]
  rep_now <- rep_now[seq_len(n_occ)]
  reuse_now <- reuse_now[seq_len(n_occ)]

  F_rows <- vector("list", n_occ)
  keep_row <- logical(n_occ)
  for (ii in seq_len(n_occ)) {
    i0 <- inst_now[ii]
    r0 <- rep_now[ii]
    u0 <- reuse_now[ii]
    fvec <- vapply(cfgs, function(c0) {
      key0 <- make_icr_key(i0, c0, r = r0, u = u0)
      if (exists(key0, envir = state$score_env, inherits = FALSE)) {
        as.numeric(get(key0, envir = state$score_env, inherits = FALSE))
      } else {
        NA_real_
      }
    }, numeric(1))
    if (all(is.finite(fvec))) {
      F_rows[[ii]] <- fvec
      keep_row[ii] <- TRUE
    }
  }
  if (!any(keep_row)) return(character(0))

  Fmat <- do.call(rbind, F_rows[keep_row])
  R_cfg <- colSums2_base(rowRanks_base(Fmat, cols = seq_len(ncol(Fmat))))
  F_mean <- colMeans(Fmat)
  ord <- order(R_cfg, F_mean, cfgs, cfg_names)
  cfg_names[ord]
}

get_top_k_qap_names_from_result <- function(res, k) {
  k <- max(1L, as.integer(k))
  stop_names <- as.character(res$stop_ranked_names)
  stop_names <- stop_names[!is.na(stop_names) & nzchar(stop_names)]
  selected <- stop_names[seq_len(min(k, length(stop_names)))]
  last_race <- if (length(res$races) > 0L) res$races[[length(res$races)]]$result else NULL

  if (length(selected) < k) {
    hist_names <- rank_qap_history_portfolio_names(last_race, exclude_names = selected)
    if (length(hist_names) > 0L) {
      need <- k - length(selected)
      selected <- c(selected, hist_names[seq_len(min(need, length(hist_names)))])
    }
  }

  if (length(selected) < k) {
    race_cfg_names <- rank_qap_current_race_config_names(last_race, exclude_names = selected)
    if (length(race_cfg_names) > 0L) {
      need <- k - length(selected)
      selected <- c(selected, race_cfg_names[seq_len(min(need, length(race_cfg_names)))])
    }
  }

  if (length(selected) > 0L) {
    return(selected[seq_len(min(k, length(selected)))])
  }
  get_top_k_qap_names(last_race, k = k)
}

decode_qap_port_name <- function(port_name, state, param_config) {
  ids <- parse_qap_port_name(port_name)
  ids <- ids[is.finite(ids)]
  if (length(ids) == 0L) return("<none>")
  paste(vapply(ids, function(id) decode_qap_cfg_label(state, id, param_config), character(1)), collapse = "; ")
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

run_qap_portfolio_once <- function(train_split,
                                   budget,
                                   port_size = 2L,
                                   best_k = 1L,
                                   enable_second_level_test = TRUE,
                                   global_free_rider_elimination = TRUE,
                                   enable_duplicate = TRUE,
                                   seed0 = 1L,
                                   race_seed0 = 700L,
                                   time_limit = 5,
                                   a_values = 0:7,
                                   p_values = 0:9,
                                   l_values = 0:5,
                                   max_T = 120L,
                                   T_first = 5L,
                                   T_each = 1L,
                                   N_iter = 3L,
                                   N_min = 3L,
                                   verbose = TRUE) {
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  load_qap_portfolio_modules(repo_root)

  param_config <- make_qap_param_config(a_values = a_values, p_values = p_values, l_values = l_values)
  instance_table <- read_qap_instance_table(train_split)
  scorer <- make_qap_score_functions(instance_table = instance_table, param_config = param_config, time_limit = time_limit)
  assign("cost_config", scorer$cost_config, envir = .GlobalEnv)
  assign("scores_configs_on_instance", scorer$scores_configs_on_instance, envir = .GlobalEnv)

  N_param <<- length(param_config$params)
  port_size_k <<- as.integer(port_size)
  N_min_default <<- as.integer(N_min)
  T_first <<- as.integer(T_first)
  LB <<- 1
  UB <<- max(length(a_values), length(p_values), length(l_values))

  res <- run_iterated_portfolio_irace(
    param_config = param_config,
    k_portfolio = as.integer(port_size),
    alpha = 0.05,
    max_T = as.integer(max_T),
    min_T_for_test = as.integer(T_first),
    T_each = as.integer(T_each),
    w = NULL,
    sigma1 = 0,
    sigma2 = 0,
    p = NULL,
    seed0 = as.integer(seed0),
    race_seed0 = as.integer(race_seed0),
    B_total = as.integer(budget),
    N_iter = as.integer(N_iter),
    N_min = as.integer(N_min),
    enable_second_level_test = isTRUE(enable_second_level_test),
    global_free_rider_elimination = isTRUE(global_free_rider_elimination),
    enable_duplicate = isTRUE(enable_duplicate),
    stop_if_nj_le = as.integer(best_k),
    verbose = isTRUE(verbose),
    training_instances = instance_table$instance_id
  )

  last_race <- if (length(res$races) > 0L) res$races[[length(res$races)]]$result else NULL
  decode_state <- if (!is.null(res$state)) res$state else if (!is.null(last_race)) last_race$state else NULL
  decoded_best <- if (!is.null(last_race)) {
    best_name_chr <- as.character(res$best_name)
    best_name_chr <- gsub("^\\(|\\)$", "", best_name_chr)
    ids <- suppressWarnings(as.integer(strsplit(best_name_chr, ",", fixed = TRUE)[[1L]]))
    ids <- ids[is.finite(ids)]
    if (length(ids) == 0L) {
      "<none>"
    } else {
      paste(vapply(ids, function(id) decode_qap_cfg_label(decode_state, id, param_config), character(1)), collapse = "; ")
    }
  } else {
    "<none>"
  }

  list(
    result = res,
    best_name_decoded = decoded_best,
    instance_table = instance_table,
    param_config = param_config
  )
}

run_qap_portfolio_grid <- function(split_seeds,
                                   budgets,
                                   port_sizes,
                                   algo_seeds,
                                   best_k = 1L,
                                   enable_second_level_test = TRUE,
                                   global_free_rider_elimination = TRUE,
                                   enable_duplicate = TRUE,
                                   race_seed_offset = 700L,
                                   time_limit = 5,
                                   max_T = 120L,
                                   T_first = 5L,
                                   T_each = 1L,
                                   N_iter = 3L,
                                   N_min = 3L,
                                   debug = FALSE,
                                   debug_eval_hit = TRUE,
                                   debug_softrestart = FALSE,
                                   debug_duplicate = FALSE,
                                   verbose = TRUE) {
  options(DBG = isTRUE(debug))
  options(DBG2 = isTRUE(debug))
  options(DBG_eval_hit = isTRUE(debug_eval_hit))
  options(DBG_softrestart = isTRUE(debug_softrestart))
  options(DBG_duplicate = isTRUE(debug_duplicate))

  grid <- expand.grid(
    split_seed = as.integer(split_seeds),
    budget = as.integer(budgets),
    port_size = as.integer(port_sizes),
    algo_seed = as.integer(algo_seeds),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(grid))
  details <- vector("list", nrow(grid))

  for (idx in seq_len(nrow(grid))) {
    split_seed <- grid$split_seed[[idx]]
    budget <- grid$budget[[idx]]
    port_size <- grid$port_size[[idx]]
    algo_seed <- grid$algo_seed[[idx]]
    race_seed <- as.integer(race_seed_offset + algo_seed)

    split_paths <- get_qap_split_paths(split_seed = split_seed, script_dir = script_dir)
    if (isTRUE(verbose)) {
      cat(sprintf(
        "\n[QAP] split_seed=%d | algo_seed=%d | budget=%d | port_size=%d\n",
        split_seed, algo_seed, budget, port_size
      ))
    }

    one <- run_qap_portfolio_once(
      train_split = split_paths$train_split,
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
      split_seed = split_seed,
      algo_seed = algo_seed,
      race_seed0 = race_seed,
      budget = budget,
      port_size = port_size,
      time_limit = as.numeric(time_limit),
      train_split = normalizePath(split_paths$train_split, winslash = "/", mustWork = TRUE),
      test_split = normalizePath(split_paths$test_split, winslash = "/", mustWork = TRUE),
      train_instances = length(readLines(split_paths$train_split, warn = FALSE)),
      test_instances = length(readLines(split_paths$test_split, warn = FALSE)),
      best_k = as.integer(best_k),
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
  }

  list(
    summary = do.call(rbind, rows),
    details = details
  )
}

default_qap_overall_config <- function() {
  list(
    split_seeds = 11L,
    budgets = c(1000L, 2500L),
    port_sizes = c(1L, 2L, 4L),
    algo_seeds = 1L,
    best_k = 1L,
    enable_second_level_test = TRUE,
    global_free_rider_elimination = TRUE,
    enable_duplicate = TRUE,
    race_seed_offset = 700L,
    time_limit = 5,
    max_T = 120L,
    T_first = 5L,
    T_each = 1L,
    N_iter = 3L,
    N_min = as.integer(max(3L, floor(2 + log2(3)))),
    debug = FALSE,
    debug_eval_hit = TRUE,
    debug_softrestart = FALSE,
    debug_duplicate = FALSE,
    verbose = TRUE,
    name = NULL,
    out_prefix = "qap_run_overall",
    save_outputs = TRUE
  )
}

run_qap_overall_rstudio <- function(split_seeds = 11L,
                                    budgets = c(1000L, 2500L),
                                    port_sizes = c(1L, 2L, 4L),
                                    algo_seeds = 1L,
                                    best_k = 1L,
                                    enable_second_level_test = TRUE,
                                    global_free_rider_elimination = TRUE,
                                    enable_duplicate = TRUE,
                                    race_seed_offset = 700L,
                                    time_limit = 5,
                                    max_T = 120L,
                                    T_first = 5L,
                                    T_each = 1L,
                                    N_iter = 3L,
                                    N_min = as.integer(max(3L, floor(2 + log2(3)))),
                                    debug = FALSE,
                                    debug_eval_hit = TRUE,
                                    debug_softrestart = FALSE,
                                    debug_duplicate = FALSE,
                                    verbose = TRUE,
                                    name = NULL,
                                    out_prefix = "qap_run_overall",
                                    save_outputs = TRUE) {
  if (!is.null(name) && nzchar(as.character(name))) {
    out_prefix <- as.character(name)[[1L]]
  }
  cat("=== QAP Portfolio Run Overall ===\n")
  cat("split_seeds:", paste(as.integer(split_seeds), collapse = ","), "\n")
  cat("budgets:", paste(as.integer(budgets), collapse = ","), "\n")
  cat("port_sizes:", paste(as.integer(port_sizes), collapse = ","), "\n")
  cat("algo_seeds:", paste(as.integer(algo_seeds), collapse = ","), "\n")
  cat("best_k:", as.integer(best_k), "\n")
  cat("time_limit:", as.numeric(time_limit), "\n")
  cat(
    "feature flags:",
    sprintf(
      "second_level_test=%s global_free_rider_elimination=%s duplicate=%s",
      isTRUE(enable_second_level_test),
      isTRUE(global_free_rider_elimination),
      isTRUE(enable_duplicate)
    ),
    "\n"
  )
  cat(
    "debug flags:",
    sprintf(
      "DBG=%s DBG_eval_hit=%s DBG_softrestart=%s DBG_duplicate=%s",
      isTRUE(debug), isTRUE(debug_eval_hit), isTRUE(debug_softrestart), isTRUE(debug_duplicate)
    ),
    "\n"
  )

  res <- run_qap_portfolio_grid(
    split_seeds = as.integer(split_seeds),
    budgets = as.integer(budgets),
    port_sizes = as.integer(port_sizes),
    algo_seeds = as.integer(algo_seeds),
    best_k = as.integer(best_k),
    enable_second_level_test = isTRUE(enable_second_level_test),
    global_free_rider_elimination = isTRUE(global_free_rider_elimination),
    enable_duplicate = isTRUE(enable_duplicate),
    race_seed_offset = as.integer(race_seed_offset),
    time_limit = as.numeric(time_limit),
    max_T = as.integer(max_T),
    T_first = as.integer(T_first),
    T_each = as.integer(T_each),
    N_iter = as.integer(N_iter),
    N_min = as.integer(N_min),
    debug = isTRUE(debug),
    debug_eval_hit = isTRUE(debug_eval_hit),
    debug_softrestart = isTRUE(debug_softrestart),
    debug_duplicate = isTRUE(debug_duplicate),
    verbose = isTRUE(verbose)
  )

  summary_df <- res$summary
  print(summary_df)

  csv_path <- file.path(script_dir, paste0(out_prefix, ".csv"))
  rds_path <- file.path(script_dir, paste0(out_prefix, ".Rds"))

  if (isTRUE(save_outputs)) {
    write.csv(summary_df, csv_path, row.names = FALSE)
    saveRDS(res, rds_path)
    cat("\nSaved summary:", csv_path, "\n")
    cat("Saved details:", rds_path, "\n")
  }

  invisible(res)
}
