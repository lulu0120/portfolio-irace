library(truncnorm)
port_size_k=3
N_param =1
N_min_default = floor(2 + log2(N_param))
T_first = 5
LB = 1
UB = 100
k_best = 1
use_scoring_test_2dim = F
 


options(DBG = F)
options(DBG2 = F)
options(DBG_softrestart = F)
options(DBG_duplicate = F)

repo_root <- local({
  find_repo_root <- function(start) {
    if (!is.character(start) || length(start) == 0L || is.na(start) || !nzchar(start)) {
      return(NA_character_)
    }
    current <- normalizePath(start[[1L]], winslash = "/", mustWork = TRUE)
    if (!dir.exists(current)) current <- dirname(current)
    repeat {
      if (dir.exists(file.path(current, "Multiple_multidim")) &&
          dir.exists(file.path(current, "Multiple_2combine"))) {
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
    editor_path <- if (!is.null(ctx)) ctx$path else ""
    if (is.character(editor_path) && length(editor_path) > 0L && nzchar(editor_path) && file.exists(editor_path)) {
      resolved <- find_repo_root(editor_path[[1L]])
      if (is.character(resolved) && !is.na(resolved)) return(resolved)
    }
  }

  file_arg <- commandArgs(trailingOnly = FALSE)
  file_flag <- grep("^--file=", file_arg, value = TRUE)
  if (length(file_flag) > 0L) {
    resolved <- find_repo_root(sub("^--file=", "", file_flag[[1L]]))
    if (is.character(resolved) && !is.na(resolved)) return(resolved)
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
    resolved <- find_repo_root(tail(ofiles, 1L))
    if (is.character(resolved) && !is.na(resolved)) return(resolved)
  }

  resolved <- find_repo_root(getwd())
  if (is.character(resolved) && !is.na(resolved)) return(resolved)

  stop("Cannot resolve repository root for Multiple_multidim/run_overall.R")
})





rtnorm_trunc = function(n, mean, sd, lower=LB, upper=UB, do_round = T) {
  n = as.integer(n)
  if (n <= 0L) return(numeric(0))

  if (length(mean) == 1L) mean = rep(as.numeric(mean), n)
  if (length(sd) == 1L) sd = rep(as.numeric(sd), n)

  if (length(mean) != n) stop("rtnorm_trunc: length(mean) must be 1 or n")
  if (length(sd) != n) stop("rtnorm_trunc: length(sd) must be 1 or n")

  mean = as.numeric(mean)
  sd = as.numeric(sd)

  x = numeric(n)
  bad_sd = !is.finite(sd) | sd <= 0

  if (any(bad_sd)) {
    x[bad_sd] = pmin(pmax(mean[bad_sd], lower), upper)
  }

  ok = !bad_sd
  if (any(ok)) {
    x[ok] = rtruncnorm(sum(ok), a = lower, b = upper, mean = mean[ok], sd = sd[ok])
  }

  if (do_round) x = round(x)
  x
}




#run the ones 

# ---- load all functions (fast) ----
options(keep.source = FALSE)

src_dir <- file.path(repo_root, "Multiple_multidim")
score_dir <- file.path(repo_root, "Multiple_2combine")
result_dir <- file.path(src_dir, "resultdata")
plot_out_dir <- file.path(repo_root, "1Dctsresult")
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(plot_out_dir)) dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  file.path(src_dir, "initial_and_env.R"),
  file.path(score_dir, "scoring.R"),
  file.path(src_dir, "statistical_test.R"),
  file.path(src_dir, "elites_functions.R"),
  file.path(src_dir, "generate_discrete_input.R"),
  file.path(src_dir, "within_race.R"),
  file.path(src_dir, "wrapper_new.R")
)

invisible(lapply(files, function(f) {
  if (!file.exists(f)) stop("Missing: ", f)
  source(f, local = .GlobalEnv, echo = FALSE, verbose = FALSE, keep.source = FALSE)
}))




param0 = make_param_discrete(
  name = "int",
  values = 1:100
)

make_param_cts <- function(name = "c", lower = 1, upper = 100) {
  stopifnot(is.finite(lower), is.finite(upper), lower < upper)
  list(
    name  = name,
    type  = "cts",
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}

make_param_multi_numeric <- function() {
  list(
    type = "multi",
    params = list(
      list(name = "x_int", type = "int", lower = 1, upper = 100),
      list(name = "x_cts", type = "cts", lower = -10, upper = 10)
    )
  )
}

parse_port_name = function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(numeric(0))
  s = gsub("[() ]", "", x)
  if (!nzchar(s)) return(numeric(0))
  as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
}

decode_best_name_multidim = function(best_name, state) {
  ids = parse_port_name(best_name)
  if (length(ids) == 0L || is.null(state) || is.null(state$cfg_vec_by_id_env)) return(NA_character_)
  parts = vapply(ids, function(id) {
    key = as.character(as.integer(id))
    if (!exists(key, envir = state$cfg_vec_by_id_env, inherits = FALSE)) return("[<missing>]")
    v = as.numeric(get(key, envir = state$cfg_vec_by_id_env, inherits = FALSE))
    if (length(v) == 0L) return("[<empty>]")
    vv = vapply(v, function(x) sprintf("%.4f", as.numeric(x)), character(1))
    paste0("[", paste(vv, collapse = ","), "]")
  }, character(1))
  paste0("(", paste(parts, collapse = ","), ")")
}

best_member_win_counts = function(best_name, last_race) {
  if (is.null(last_race) || is.null(last_race$state) || is.null(last_race$inst_used)) {
    return(NA_character_)
  }
  ids = as.numeric(parse_port_name(best_name))
  if (length(ids) == 0L || any(!is.finite(ids))) return(NA_character_)
  inst = as.integer(last_race$inst_used)
  inst = inst[is.finite(inst)]
  if (length(inst) == 0L) return(paste0("(", paste(rep(0L, length(ids)), collapse = ","), ")"))
  
  ticks = cfg_tick(ids)
  dup_idx = ave(ticks, ticks, FUN = seq_along)
  st = last_race$state
  wins = integer(length(ids))
  
  for (i0 in inst) {
    sc = rep(Inf, length(ids))
    for (k in seq_along(ids)) {
      key_k = if (dup_idx[k] <= 1L) {
        make_ic_key(i0, ids[k])
      } else {
        make_icd_key(i0, ids[k], dup_idx[k])
      }
      if (exists(key_k, envir = st$score_env, inherits = FALSE)) {
        sc[k] = as.numeric(get(key_k, envir = st$score_env, inherits = FALSE))
      } else {
        # Fallback: if duplicate-specific key is missing, reuse base key.
        key_base = make_ic_key(i0, ids[k])
        if (exists(key_base, envir = st$score_env, inherits = FALSE)) {
          sc[k] = as.numeric(get(key_base, envir = st$score_env, inherits = FALSE))
        }
      }
    }
    if (any(is.finite(sc))) {
      wins[which.min(sc)] = wins[which.min(sc)] + 1L
    }
  }
  
  paste0("(", paste(wins, collapse = ","), ")")
}

format_port_name = function(vals) {
  vals = sort(unique(as.numeric(vals)))
  if (length(vals) == 0L) return(NA_character_)
  paste0("(", paste(vals, collapse = ","), ")")
}

format_port_name_keep_dup = function(vals) {
  vals = sort(as.numeric(vals))
  if (length(vals) == 0L || any(!is.finite(vals))) return(NA_character_)
  paste0("(", paste(vals, collapse = ","), ")")
}

default_optimal_name_1d = function(k = max(as.integer(port_size_k), as.integer(k_best))) {
  k = max(1L, as.integer(k))
  if (k >= 3L) "(49.4999,50,50)" else "(49.4999,50)"
}

get_top_k_names = function(last_race, k) {
  if (is.null(last_race) || is.null(last_race$alive_names) || length(last_race$alive_names) == 0L) {
    return(character(0))
  }
  k = max(1L, as.integer(k))
  names0 = as.character(last_race$alive_names)
  ranks0 = suppressWarnings(as.numeric(last_race$alive_R))
  
  if (length(ranks0) == length(names0) && any(is.finite(ranks0))) {
    ord = order(ranks0, seq_along(names0), na.last = TRUE)
    names0 = names0[ord]
  }
  names0[seq_len(min(k, length(names0)))]
}

port_name_from_key = function(state, key) {
  if (!is.null(state) &&
      !is.null(state$ports_all_env) &&
      exists(key, envir = state$ports_all_env, inherits = FALSE)) {
    ticks = as.numeric(get(key, envir = state$ports_all_env, inherits = FALSE))
    return(format_port_name_keep_dup(cfg_from_tick(ticks)))
  }
  vals = as.numeric(strsplit(as.character(key), ",", fixed = TRUE)[[1]])
  format_port_name_keep_dup(cfg_from_tick(vals))
}

rank_history_portfolio_names = function(last_race, exclude_names = character(0)) {
  if (is.null(last_race) || is.null(last_race$state) || is.null(last_race$state$elite_seen_env)) {
    return(character(0))
  }
  state = last_race$state
  keys_all = ls(envir = state$elite_seen_env, all.names = TRUE)
  if (length(keys_all) == 0L) return(character(0))

  names_all = vapply(keys_all, function(k) port_name_from_key(state, k), character(1))
  keep = !(names_all %in% exclude_names) & !is.na(names_all) & nzchar(names_all)
  keys_all = keys_all[keep]
  names_all = names_all[keep]
  if (length(keys_all) == 0L) return(character(0))

  ports_all = lapply(keys_all, function(k) parse_port_name(port_name_from_key(state, k)))
  names(ports_all) = keys_all
  elite_e = vapply(keys_all, function(k) get_port_e(state, k), integer(1))
  cand_pos = seq_along(keys_all)
  ordered = integer(0)
  count_levels = sort(unique(elite_e), decreasing = TRUE)

  for (count_i in count_levels) {
    pos_i = cand_pos[elite_e[cand_pos] == count_i]
    if (length(pos_i) <= 1L) {
      ordered = c(ordered, pos_i)
      next
    }
    tie_keys = keys_all[pos_i]
    tie_ports = ports_all[pos_i]
    tie_inst_list = lapply(tie_keys, function(k) {
      if (exists(k, envir = state$elite_seen_env, inherits = FALSE)) {
        occ_keys = unique(as.character(get(k, envir = state$elite_seen_env, inherits = FALSE)))
        parsed = parse_ir_key(occ_keys)
        keep_occ = is.finite(parsed$i) & is.finite(parsed$r) & is.finite(parsed$u)
        if (!any(keep_occ)) {
          data.frame(i = integer(0), r = integer(0), u = integer(0), stringsAsFactors = FALSE)
        } else {
          parsed[keep_occ, , drop = FALSE]
        }
      } else {
        data.frame(i = integer(0), r = integer(0), u = integer(0), stringsAsFactors = FALSE)
      }
    })
    tie_inst_union = do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x) > 0L, tie_inst_list))
    R_tie = rep(Inf, length(pos_i))
    F_tie_mean = rep(Inf, length(pos_i))
    if (!is.null(tie_inst_union) && nrow(tie_inst_union) > 0L) {
      tie_inst_union = unique(tie_inst_union)
      F_rows = vector("list", nrow(tie_inst_union))
      keep_row = logical(nrow(tie_inst_union))
      for (ii in seq_len(nrow(tie_inst_union))) {
        i0 = as.integer(tie_inst_union$i[ii])
        r0 = as.integer(tie_inst_union$r[ii])
        u0 = as.integer(tie_inst_union$u[ii])
        fvec = rep(NA_real_, length(pos_i))
        for (pp in seq_along(tie_ports)) {
          cfgs = as.numeric(tie_ports[[pp]])
          keys_ic = vapply(cfgs, function(c0) make_icr_key(i0, c0, r = r0, u = u0), character(1))
          has_all = all(vapply(keys_ic, exists, logical(1), envir = state$score_env, inherits = FALSE))
          if (has_all) {
            sc = vapply(keys_ic, get, numeric(1), envir = state$score_env, inherits = FALSE)
            fvec[pp] = min(as.numeric(sc))
          }
        }
        if (all(is.finite(fvec))) {
          F_rows[[ii]] = fvec
          keep_row[ii] = TRUE
        }
      }
      if (any(keep_row)) {
        Fmat_tie = do.call(rbind, F_rows[keep_row])
        R_tie = colSums2_base(rowRanks_base(Fmat_tie, cols = seq_len(ncol(Fmat_tie))))
        F_tie_mean = colMeans(Fmat_tie)
      }
    }
    tie_min_tick = vapply(tie_ports, function(P) min(cfg_tick(P)), integer(1))
    ord_i = order(R_tie, F_tie_mean, tie_min_tick, names_all[pos_i])
    ordered = c(ordered, pos_i[ord_i])
  }

  names_all[ordered]
}

rank_current_race_config_names = function(last_race, exclude_names = character(0)) {
  if (is.null(last_race) ||
      is.null(last_race$state) ||
      is.null(last_race$configs_race) ||
      length(last_race$configs_race) == 0L ||
      is.null(last_race$inst_used) ||
      length(last_race$inst_used) == 0L) {
    return(character(0))
  }
  state = last_race$state
  cfgs = as.numeric(last_race$configs_race)
  cfg_names = vapply(cfgs, function(x) format_port_name_keep_dup(x), character(1))
  keep_cfg = !(cfg_names %in% exclude_names) & !duplicated(cfg_names)
  cfgs = cfgs[keep_cfg]
  cfg_names = cfg_names[keep_cfg]
  if (length(cfgs) == 0L) return(character(0))

  inst_now = as.integer(last_race$inst_used)
  rep_now = if (!is.null(last_race$inst_repeat_used)) as.integer(last_race$inst_repeat_used) else rep.int(0L, length(inst_now))
  reuse_now = if (!is.null(last_race$inst_reuse_new_used)) as.integer(last_race$inst_reuse_new_used) else rep.int(0L, length(inst_now))
  n_occ = min(length(inst_now), length(rep_now), length(reuse_now))
  inst_now = inst_now[seq_len(n_occ)]
  rep_now = rep_now[seq_len(n_occ)]
  reuse_now = reuse_now[seq_len(n_occ)]

  F_rows = vector("list", n_occ)
  keep_row = logical(n_occ)
  for (ii in seq_len(n_occ)) {
    i0 = inst_now[ii]
    r0 = rep_now[ii]
    u0 = reuse_now[ii]
    fvec = vapply(cfgs, function(c0) {
      key0 = make_icr_key(i0, c0, r = r0, u = u0)
      if (exists(key0, envir = state$score_env, inherits = FALSE)) {
        as.numeric(get(key0, envir = state$score_env, inherits = FALSE))
      } else {
        NA_real_
      }
    }, numeric(1))
    if (all(is.finite(fvec))) {
      F_rows[[ii]] = fvec
      keep_row[ii] = TRUE
    }
  }
  if (!any(keep_row)) return(character(0))

  Fmat = do.call(rbind, F_rows[keep_row])
  R_cfg = colSums2_base(rowRanks_base(Fmat, cols = seq_len(ncol(Fmat))))
  F_mean = colMeans(Fmat)
  ord = order(R_cfg, F_mean, cfg_tick(cfgs), cfg_names)
  cfg_names[ord]
}

get_top_k_names_from_result = function(res, k) {
  k = max(1L, as.integer(k))
  stop_names = as.character(res$stop_ranked_names)
  stop_names = stop_names[!is.na(stop_names) & nzchar(stop_names)]
  selected = stop_names[seq_len(min(k, length(stop_names)))]
  last_race = if (length(res$races) > 0L) res$races[[length(res$races)]]$result else NULL

  if (length(selected) < k) {
    hist_names = rank_history_portfolio_names(last_race, exclude_names = selected)
    if (length(hist_names) > 0L) {
      need = k - length(selected)
      selected = c(selected, hist_names[seq_len(min(need, length(hist_names)))])
    }
  }

  if (length(selected) < k) {
    race_cfg_names = rank_current_race_config_names(last_race, exclude_names = selected)
    if (length(race_cfg_names) > 0L) {
      need = k - length(selected)
      selected = c(selected, race_cfg_names[seq_len(min(need, length(race_cfg_names)))])
    }
  }

  if (length(selected) > 0L) {
    return(selected[seq_len(min(k, length(selected)))])
  }

  get_top_k_names(last_race, k = k)
}

deterministic_cost_config = function(c, i, w, p) {
  v = w * abs(c - 50)
  c_int = as.integer(round(c))
  i_int = as.integer(i)
  penalty = if ((c_int %% 2L) == (i_int %% 2L)) 0 else p
  v + penalty
}

sum_portfolio_cost_deterministic = function(port_vals, instances, seed0, w, p) {
  port_vals = as.numeric(port_vals)
  if (length(port_vals) == 0L) return(NA_real_)
  vals = vapply(as.integer(instances), function(i) {
    sc = vapply(port_vals, function(c) {
      deterministic_cost_config(c = c, i = i, w = w, p = p)
    }, numeric(1))
    min(sc)
  }, numeric(1))
  sum(vals)
}

sum_cost_from_port_name_deterministic = function(port_name, instances, seed0, w, p) {
  vals = parse_port_name(port_name)
  if (length(vals) == 0L || any(!is.finite(vals))) return(NA_real_)
  sum_portfolio_cost_deterministic(
    port_vals = vals,
    instances = instances,
    seed0 = seed0,
    w = w,
    p = p
  )
}

cost_gap_vs_optimal = function(best_name,
                               optimal_name = default_optimal_name_1d(),
                               instances = c(1, 2),
                               seed0,
                               w, p) {
  best_vals = parse_port_name(best_name)
  opt_vals = parse_port_name(optimal_name)
  best_cost = sum_portfolio_cost_deterministic(best_vals, instances, seed0, w, p)
  opt_cost = sum_portfolio_cost_deterministic(opt_vals, instances, seed0, w, p)
  list(
    best_cost = best_cost,
    optimal_cost = opt_cost,
    cost_gap = best_cost - opt_cost
  )
}


run_one_iterated_simple = function(settings_iter,
                                   w, sigma1, sigma2, p,
                                   B_total,
                                   seed_shift = 0L,
                                   enable_second_level_test = NULL,
                                   global_free_rider_elimination = NULL,
                                   enable_duplicate = NULL,
                                   optimal_name = default_optimal_name_1d(),
                                   eval_instances = c(1, 2),
                                   best_k = 1L,
                                   k_best = NULL) {
  if (!is.null(k_best)) best_k = k_best
  best_k = max(1L, as.integer(best_k))
  
  overrides = list(
    w = w,
    sigma1 = sigma1,
    sigma2 = sigma2,
    p = p,
    B_total = B_total,
    seed0 = settings_iter$seed0 + as.integer(seed_shift),
    race_seed0 = settings_iter$race_seed0 + as.integer(seed_shift)
  )
  if (!is.null(enable_second_level_test)) {
    overrides$enable_second_level_test = isTRUE(enable_second_level_test)
  }
  if (!is.null(global_free_rider_elimination)) {
    overrides$global_free_rider_elimination = isTRUE(global_free_rider_elimination)
  }
  if (!is.null(enable_duplicate)) {
    overrides$enable_duplicate = isTRUE(enable_duplicate)
  }
  
  args = modifyList(settings_iter, overrides)  # <- 关键：覆盖而不是追加
  res = do.call(run_iterated_portfolio_irace, args)
  seed0_run = settings_iter$seed0 + as.integer(seed_shift)
  best_cost_det = sum_cost_from_port_name_deterministic(
    port_name = res$best_name,
    instances = eval_instances,
    seed0 = seed0_run,
    w = w,
    p = p
  )
  optimal_cost_det = sum_cost_from_port_name_deterministic(
    port_name = optimal_name,
    instances = eval_instances,
    seed0 = seed0_run,
    w = w,
    p = p
  )
  gap = list(
    best_cost = best_cost_det,
    optimal_cost = optimal_cost_det,
    cost_gap = best_cost_det - optimal_cost_det
  )
  
  n_races = length(res$races)
  last_race = if (n_races > 0) res$races[[n_races]]$result else NULL
  alive_count_last = if (!is.null(last_race)) last_race$alive_count else NA_integer_
  alive_names_last = if (!is.null(last_race)) last_race$alive_names else character(0)
  top_k_names_raw = get_top_k_names_from_result(res, k = best_k)
  if (!is.na(res$best_name) && nzchar(res$best_name)) {
    rest = top_k_names_raw[top_k_names_raw != res$best_name]
    top_k_names = c(res$best_name, rest)
  } else {
    top_k_names = top_k_names_raw
  }
  top_k_names = top_k_names[seq_len(min(length(top_k_names), best_k))]
  top_k_union = sort(unique(unlist(lapply(top_k_names, parse_port_name))))
  top_k_union_name = format_port_name(top_k_union)
  
  best_k_sum_cost = NA_real_
  best_k_sum_costgap = NA_real_
  if (length(top_k_union) > 0L) {
    best_k_sum_cost = sum_portfolio_cost_deterministic(
      port_vals = top_k_union,
      instances = eval_instances,
      seed0 = seed0_run,
      w = w,
      p = p
    )
    if (as.integer(best_k) > 1L) {
      best_k_sum_costgap = best_k_sum_cost - gap$optimal_cost
    }
  }
  
  data.frame(
    w = w,
    sigma1 = sigma1,
    sigma2 = sigma2,
    p = p,
    B_total = B_total,
    seed_shift = seed_shift,
    best_name = res$best_name,
    n_races = n_races,
    total_eval_used = res$total_eval_used,
    total_instances_used = res$total_instances_used,
    alive_count = alive_count_last,
    best_cost = gap$best_cost,
    best_k = as.integer(best_k),
    best_k_names = paste(top_k_names, collapse = ";"),
    best_k_union = top_k_union_name,
    best_k_sum_cost = best_k_sum_cost,
    best_k_sum_costgap = best_k_sum_costgap,
    optimal_name = optimal_name,
    optimal_cost = gap$optimal_cost,
    cost_gap = gap$cost_gap,
    best_is_12 = (!is.na(res$best_name)) && res$best_name == optimal_name,
    unique_12 = length(alive_names_last) == 1 && alive_names_last[1] == optimal_name,
    best_name_decoded = if (!is.null(last_race) && !is.null(last_race$state)) decode_best_name_multidim(res$best_name, last_race$state) else NA_character_,
    best_member_win_counts = best_member_win_counts(res$best_name, last_race),
    stringsAsFactors = FALSE
  )
}


test_hyper_iterated_simple = function(w_list, sigma1_list, sigma2_list, p_list,
                                      B_list, seed_list,
                                      settings_iter,
                                      enable_second_level_test = NULL,
                                      global_free_rider_elimination = NULL,
                                      enable_duplicate = NULL,
                                      optimal_name = default_optimal_name_1d(),
                                      eval_instances = c(1, 2),
                                      best_k = 1L,
                                      k_best = NULL,
                                      progress_csv = NULL,
                                      recompute_progress_costs = TRUE) {
  if (!is.null(k_best)) best_k = k_best
  best_k = max(1L, as.integer(best_k))
  
  L = length(w_list)
  stopifnot(
    length(sigma1_list) == L,
    length(sigma2_list) == L,
    length(p_list) == L,
    length(B_list) == L,
    length(seed_list) == L
  )
  
  out = vector("list", L)
  for (i in seq_len(L)) {
    out[[i]] = run_one_iterated_simple(
      settings_iter = settings_iter,
      w = w_list[i],
      sigma1 = sigma1_list[i],
      sigma2 = sigma2_list[i],
      p = p_list[i],
      B_total = B_list[i],
      seed_shift = seed_list[i],
      enable_second_level_test = enable_second_level_test,
      global_free_rider_elimination = global_free_rider_elimination,
      enable_duplicate = enable_duplicate,
      optimal_name = optimal_name,
      eval_instances = eval_instances,
      best_k = best_k
    )

    if (!is.null(progress_csv) && nzchar(progress_csv)) {
      progress_df = do.call(rbind, out[seq_len(i)])
      if (isTRUE(recompute_progress_costs)) {
        progress_df = recompute_costs_from_best(
          progress_df,
          optimal_name = optimal_name,
          instances = eval_instances,
          seed0_base = settings_iter$seed0
        )
      }
      utils::write.csv(progress_df, progress_csv, row.names = FALSE)
    }
  }
  
  do.call(rbind, out)
}

run_grid_iterated_simple = function(param_grid, settings_iter,
                                    enable_second_level_test = NULL,
                                    global_free_rider_elimination = NULL,
                                    enable_duplicate = NULL,
                                    optimal_name = default_optimal_name_1d(),
                                    eval_instances = c(1, 2),
                                    best_k = 1L,
                                    k_best = NULL,
                                    progress_csv = NULL,
                                    recompute_progress_costs = TRUE) {
  if (!is.null(k_best)) best_k = k_best
  best_k = max(1L, as.integer(best_k))
  test_hyper_iterated_simple(
    w_list      = param_grid$w,
    sigma1_list = param_grid$sigma1,
    sigma2_list = param_grid$sigma2,
    p_list      = param_grid$p,
    B_list      = param_grid$B_total,
    seed_list   = param_grid$seed_shift,
    settings_iter = settings_iter,
    enable_second_level_test = enable_second_level_test,
    global_free_rider_elimination = global_free_rider_elimination,
    enable_duplicate = enable_duplicate,
    optimal_name = optimal_name,
    eval_instances = eval_instances,
    best_k = best_k,
    k_best = best_k,
    progress_csv = progress_csv,
    recompute_progress_costs = recompute_progress_costs
  )
}

recompute_costs_from_best = function(df,
                                     optimal_name = default_optimal_name_1d(),
                                     instances = c(1, 2),
                                     seed0_base = 1L) {
  stopifnot(all(c("best_name", "w", "p", "seed_shift") %in% names(df)))
  
  best_cost = vapply(seq_len(nrow(df)), function(i) {
    sum_cost_from_port_name_deterministic(
      port_name = df$best_name[i],
      instances = instances,
      seed0 = as.integer(seed0_base) + as.integer(df$seed_shift[i]),
      w = as.numeric(df$w[i]),
      p = as.numeric(df$p[i])
    )
  }, numeric(1))
  
  optimal_cost = vapply(seq_len(nrow(df)), function(i) {
    sum_cost_from_port_name_deterministic(
      port_name = optimal_name,
      instances = instances,
      seed0 = as.integer(seed0_base) + as.integer(df$seed_shift[i]),
      w = as.numeric(df$w[i]),
      p = as.numeric(df$p[i])
    )
  }, numeric(1))
  
  df$best_cost_recomputed = best_cost
  df$optimal_cost_recomputed = optimal_cost
  df$cost_gap_recomputed = best_cost - optimal_cost
  df
}

set.seed(1) #feed 2:1 even and odd
iterated_settings = list(
  param_config = if (isTRUE(use_scoring_test_2dim)) make_param_multi_numeric() else make_param_cts(name = "c", lower = LB, upper = UB),
  #param_config = make_param_discrete(),
  #generator = function(param, n, seed) generate_configs_discrete(param, n = n, seed = seed),
  k_portfolio = port_size_k,
  alpha = 0.05,
  max_T = 200,
  min_T_for_test = T_first,
  T_each = 1,
  seed0 = 1,
  race_seed0 = 100,
  #training_instances = 1:60,
  training_instances = sample(c(seq(2, 40, by = 2), seq(1, 39, by = 2), seq(42, 80, by = 2))),
  N_iter = N_min_default,
  enable_second_level_test = F,
  global_free_rider_elimination = F,
  enable_duplicate = F,
  weight_type = "uniform",
  domain_update = "keep_all",
  eps_counter = 1,
  verbose = F
)

param_grid = expand.grid(
  w = c(5),
  sigma1 = c( 100),
  sigma2 = c(20),
  p = c(50),
  #B_total = seq(200, 300, by = 100),
  #B_total = seq(400, 1200, by = 100),
  #B_total = c(100,200,300,400,500),
  B_total = c(100,200,300,400,500,600,700,800,900,1000,1100,1200),
  seed_shift = 1:20,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

 
res_grid = run_grid_iterated_simple(
  param_grid,
  iterated_settings,
  k_best = k_best,
  progress_csv = file.path(src_dir, "singlerun_baseline.csv")
)

res_grid$best_k_names
table(res_grid$best_name)

source(file.path(score_dir, "performance_eval.R"),
       local = .GlobalEnv, echo = FALSE, verbose = FALSE)

res_grid$sigma1 <- 0

mc_out = mc_plot_res_grid_cost_by_budget(
  res_grid = res_grid,
  k_best = k_best,
  n_mc = 2000L,
  seed0_base = iterated_settings$seed0,
  instance_start = 1L
)

cat("\n[MCMC summary by budget]\n")
print(mc_out$summary_by_budget)

