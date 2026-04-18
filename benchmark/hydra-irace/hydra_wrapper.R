if (!exists("LB", inherits = TRUE)) LB = 1
if (!exists("UB", inherits = TRUE)) UB = 100
if (!exists("T_first", inherits = TRUE)) T_first = 5
if (!exists("N_param", inherits = TRUE)) N_param = 1
if (!exists("N_min_default", inherits = TRUE)) N_min_default = 2

if (!exists("computeNbConfigurations", mode = "function", inherits = TRUE)) {
  computeNbConfigurations = function(B_j, mu, T_each, j) {
    denom = mu + T_each * min(5L, as.integer(j))
    floor(as.integer(B_j) / denom)
  }
}

if (!exists("get_instances_for_first_test", mode = "function", inherits = TRUE)) {
  get_instances_for_first_test = function(state, race_id, Tnew, race_seed, min_T_for_test) {
    stopifnot(min_T_for_test >= 1L)
    race_id = as.integer(race_id)
    Tnew = as.integer(Tnew)
    race_seed = as.integer(race_seed)
    if (race_id <= 1L) {
      return(seq_len(min_T_for_test))
    }
    q = make_instance_queue(state, race_id = race_id, Tnew = Tnew, seed = race_seed)
    stopifnot(length(q) >= min_T_for_test)
    q[seq_len(min_T_for_test)]
  }
}

compute_hydra_budget_split = function(B_total, n_slots) {
  B_total = as.integer(B_total)
  n_slots = as.integer(n_slots)
  if (n_slots <= 0L) stop("n_slots must be positive.")
  base = B_total %/% n_slots
  rem = B_total %% n_slots
  out = rep(base, n_slots)
  if (rem > 0L) out[seq_len(rem)] = out[seq_len(rem)] + 1L
  out
}

hydra_default_sd0 = function(param_config, lower = LB, upper = UB) {
  param_type = tolower(as.character(param_config$type))
  if (param_type %in% c("multi", "multidim", "multivariate")) {
    ptype = tolower(vapply(param_config$params, function(p) as.character(p$type), character(1)))
    is_num = !(ptype %in% c("categorical", "class", "cat"))
    sd0 = rep(0, length(param_config$params))
    if (any(is_num)) {
      lowers = vapply(param_config$params[is_num], function(p) as.numeric(p$lower), numeric(1))
      uppers = vapply(param_config$params[is_num], function(p) as.numeric(p$upper), numeric(1))
      sd0[is_num] = (uppers - lowers) / 2
    }
    return(sd0)
  }
  (upper - lower) / 2
}

hydra_n_param = function(param_config) {
  param_type = tolower(as.character(param_config$type))
  if (param_type %in% c("multi", "multidim", "multivariate")) {
    return(length(param_config$params))
  }
  1L
}

hydra_init_state = function(state = NULL) {
  if (is.null(state)) state = init_race_state()
  if (is.null(state$slot_score_env)) state$slot_score_env = new.env(parent = emptyenv())
  if (is.null(state$portfolio_slots)) state$portfolio_slots = list()
  if (is.null(state$training_instances)) state$training_instances = integer(0)
  state
}

hydra_slot_key = function(slot_id, i, cfg, dup_order) {
  paste0(
    "slot=", as.integer(slot_id),
    "|i=", as.integer(i),
    "|c=", cfg_tick(cfg),
    "|d=", as.integer(dup_order)
  )
}

hydra_duplicate_order = function(portfolio_slots, cfg) {
  if (length(portfolio_slots) == 0L) return(1L)
  tick = cfg_tick(cfg)
  prev = vapply(portfolio_slots, function(x) cfg_tick(x$cfg), integer(1))
  sum(prev == tick) + 1L
}

hydra_call_scores_configs_on_instance = function(configs, i, base_seed,
                                                 w, sigma1, sigma2, p,
                                                 cat_penalty_notA = 3,
                                                 state = NULL) {
  args = list(
    configs = configs,
    i = i,
    base_seed = base_seed,
    w = w,
    sigma1 = sigma1,
    sigma2 = sigma2,
    p = p
  )
  sc_formals = names(formals(scores_configs_on_instance))
  if ("state" %in% sc_formals) {
    args$state = state
  }
  if ("cat_penalty_notA" %in% sc_formals) {
    args$cat_penalty_notA = cat_penalty_notA
  }
  do.call(scores_configs_on_instance, args)
}

hydra_call_cost_config = function(c, i, w, p, sigma2,
                                  cat_penalty_notA = 3,
                                  instance_random = 0,
                                  eps = NULL,
                                  state = NULL) {
  args = list(
    c = c,
    i = i,
    w = w,
    p = p,
    sigma2 = sigma2,
    instance_random = instance_random,
    eps = eps
  )
  cost_formals = names(formals(cost_config))
  if ("state" %in% cost_formals) {
    args$state = state
  }
  if ("cat_penalty_notA" %in% cost_formals) {
    args$cat_penalty_notA = cat_penalty_notA
  }
  do.call(cost_config, args)
}

hydra_generate_configs = function(param_config, Nj, race_seed, sd0,
                                  lower, upper, round_now, n_param_cur,
                                  slot_id, state) {
  args = list(
    param_config = param_config,
    Nj = Nj,
    required_cfgs = integer(0),
    elite_ranks = numeric(0),
    j = race_seed,
    sd0 = sd0,
    lower = lower,
    upper = upper,
    round = round_now,
    n_param = n_param_cur,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L
  )
  gen_formals = names(formals(generate_configs_numeric_irace_with_sd))
  if ("state" %in% gen_formals) {
    args$state = state
  }
  if ("iter_idx" %in% gen_formals) {
    args$iter_idx = slot_id
  }
  do.call(generate_configs_numeric_irace_with_sd, args)
}

hydra_eval_base_score = function(state, cfg, i, base_seed,
                                 w, sigma1, sigma2, p,
                                 cat_penalty_notA = 3) {
  key = make_ic_key(i, cfg)
  if (exists(key, envir = state$score_env, inherits = FALSE)) {
    return(list(score = as.numeric(get(key, envir = state$score_env, inherits = FALSE)), new_eval = 0L))
  }
  val = hydra_call_scores_configs_on_instance(
    configs = c(cfg),
    i = i,
    base_seed = base_seed,
    w = w,
    sigma1 = sigma1,
    sigma2 = sigma2,
    p = p,
    cat_penalty_notA = cat_penalty_notA,
    state = state
  )[1L]
  cache_set_scores_state(state, keys = key, scores = val)
  list(score = as.numeric(val), new_eval = 1L)
}

hydra_eval_duplicate_score = function(state, cfg, i, dup_order, base_seed,
                                      w, sigma1, sigma2, p,
                                      cat_penalty_notA = 3) {
  if (dup_order <= 1L) {
    return(hydra_eval_base_score(
      state = state,
      cfg = cfg,
      i = i,
      base_seed = base_seed,
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      cat_penalty_notA = cat_penalty_notA
    ))
  }

  key_icd = make_icd_key(i, cfg, dup_order)
  if (exists(key_icd, envir = state$score_env, inherits = FALSE)) {
    return(list(score = as.numeric(get(key_icd, envir = state$score_env, inherits = FALSE)), new_eval = 0L))
  }

  key_inst = paste0("i=", as.integer(i), "|inst_rand")
  if (exists(key_inst, envir = state$score_env, inherits = FALSE)) {
    inst_rand = as.numeric(get(key_inst, envir = state$score_env, inherits = FALSE))
  } else {
    set.seed(as.integer(base_seed))
    inst_rand = rnorm(1, mean = 0, sd = sigma1)
    assign(key_inst, as.numeric(inst_rand), envir = state$score_env)
  }

  cfg_seed = abs(as.integer(cfg_tick(cfg))) %% 1000000L
  set.seed(as.integer(base_seed + 100L * dup_order + cfg_seed))
  eps = rnorm(1, mean = 0, sd = sigma2)
  val = hydra_call_cost_config(
    c = cfg,
    i = i,
    w = w,
    p = p,
    sigma2 = sigma2,
    cat_penalty_notA = cat_penalty_notA,
    instance_random = inst_rand,
    eps = eps,
    state = state
  )
  assign(key_icd, as.numeric(val), envir = state$score_env)
  list(score = as.numeric(val), new_eval = 1L)
}

hydra_eval_slot_score = function(state, cfg, i, dup_order, slot_id, base_seed,
                                 w, sigma1, sigma2, p,
                                 cat_penalty_notA = 3) {
  slot_key = hydra_slot_key(slot_id, i, cfg, dup_order)
  if (exists(slot_key, envir = state$slot_score_env, inherits = FALSE)) {
    return(list(score = as.numeric(get(slot_key, envir = state$slot_score_env, inherits = FALSE)), new_eval = 0L))
  }

  out = hydra_eval_duplicate_score(
    state = state,
    cfg = cfg,
    i = i,
    dup_order = dup_order,
    base_seed = base_seed,
    w = w,
    sigma1 = sigma1,
    sigma2 = sigma2,
    p = p,
    cat_penalty_notA = cat_penalty_notA
  )
  assign(slot_key, as.numeric(out$score), envir = state$slot_score_env)
  out
}

hydra_eval_fixed_portfolio_score = function(state, i, seed0,
                                            w, sigma1, sigma2, p,
                                            cat_penalty_notA = 3) {
  if (length(state$portfolio_slots) == 0L) {
    return(list(score = Inf, new_eval = 0L))
  }

  base_seed = as.integer(seed0 + 10000L * as.integer(i))
  scores = numeric(length(state$portfolio_slots))
  eval_added = 0L
  for (idx in seq_along(state$portfolio_slots)) {
    slot = state$portfolio_slots[[idx]]
    out = hydra_eval_slot_score(
      state = state,
      cfg = slot$cfg,
      i = i,
      dup_order = slot$dup_order,
      slot_id = slot$slot_id,
      base_seed = base_seed,
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      cat_penalty_notA = cat_penalty_notA
    )
    scores[idx] = out$score
    eval_added = eval_added + out$new_eval
  }
  list(score = min(scores), new_eval = eval_added)
}

hydra_portfolio_score_from_cache = function(state, i) {
  if (length(state$portfolio_slots) == 0L) return(Inf)
  vals = vapply(state$portfolio_slots, function(slot) {
    slot_key = hydra_slot_key(slot$slot_id, i, slot$cfg, slot$dup_order)
    if (!exists(slot_key, envir = state$slot_score_env, inherits = FALSE)) return(NA_real_)
    as.numeric(get(slot_key, envir = state$slot_score_env, inherits = FALSE))
  }, numeric(1))
  vals = vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_real_)
  min(vals)
}

hydra_count_missing_fixed_slot_scores = function(state, instances) {
  instances = as.integer(instances)
  if (length(state$portfolio_slots) == 0L) return(0L)
  sum(vapply(instances, function(i) {
    sum(vapply(state$portfolio_slots, function(slot) {
      slot_key = hydra_slot_key(slot$slot_id, i, slot$cfg, slot$dup_order)
      !exists(slot_key, envir = state$slot_score_env, inherits = FALSE)
    }, logical(1)))
  }, integer(1)))
}

run_hydra_slot_irace = function(param_config,
                                state,
                                slot_id,
                                slot_budget,
                                training_instances,
                                alpha = 0.05,
                                min_T_for_test = T_first,
                                T_each = 1,
                                N_min = N_min_default,
                                w = 1,
                                sigma1 = 1,
                                sigma2 = 0.1,
                                p = 10,
                                seed0 = 1,
                                race_seed0 = 999,
                                lower = LB,
                                upper = UB,
                                cat_penalty_notA = 3,
                                verbose = TRUE) {
  state = hydra_init_state(state)
  slot_budget = as.integer(slot_budget)
  if (slot_budget <= 0L) stop("slot_budget must be positive.")

  sd0 = hydra_default_sd0(param_config, lower = lower, upper = upper)
  n_param_cur = hydra_n_param(param_config)
  n_configs = computeNbConfigurations(
    B_j = slot_budget,
    mu = min_T_for_test,
    T_each = T_each,
    j = slot_id
  )
  n_configs = as.integer(n_configs)

  param_type = tolower(as.character(param_config$type))
  if (param_type %in% c("discrete", "class", "categorical")) {
    n_configs = min(n_configs, length(param_config$values))
  }
  if (n_configs < 1L) {
    stop(sprintf("Hydra slot %d cannot start: Nj=%d < 1", slot_id, n_configs))
  }

  race_seed = as.integer(race_seed0 + slot_id - 1L)
  inst_first_test = as.integer(training_instances[seq_len(min(length(training_instances), min_T_for_test))])
  if (length(inst_first_test) < min_T_for_test) {
    stop("training_instances length is smaller than min_T_for_test.")
  }
  need_fixed = hydra_count_missing_fixed_slot_scores(state, inst_first_test)
  need_new_upper = as.integer(min_T_for_test) * as.integer(n_configs)
  need_eval_to_test = as.integer(need_fixed + need_new_upper)
  if (need_eval_to_test > slot_budget) {
    stop(sprintf(
      "Hydra slot %d cannot reach first test under budget: need=%d, budget=%d",
      slot_id, need_eval_to_test, slot_budget
    ))
  }

  cfgs = hydra_generate_configs(
    param_config = param_config,
    Nj = n_configs,
    race_seed = race_seed,
    sd0 = sd0,
    lower = lower,
    upper = upper,
    round_now = !(param_type %in% c("cts", "numeric", "real", "continuous")),
    n_param_cur = n_param_cur,
    slot_id = slot_id,
    state = state
  )$configs

  active_idx = seq_along(cfgs)
  history = matrix(NA_real_, nrow = length(training_instances), ncol = length(cfgs))
  eval_used = 0L
  used_T = 0L
  dup_orders = vapply(cfgs, function(cfg) hydra_duplicate_order(state$portfolio_slots, cfg), integer(1))
  used_instances = integer(0)

  for (t in seq_along(training_instances)) {
    i_t = as.integer(training_instances[t])
    used_T = t
    used_instances = c(used_instances, i_t)
    fixed_now = hydra_eval_fixed_portfolio_score(
      state = state,
      i = i_t,
      seed0 = seed0,
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      cat_penalty_notA = cat_penalty_notA
    )
    best_prev = fixed_now$score
    active_cfgs = cfgs[active_idx]
    active_dup = dup_orders[active_idx]

    need_eval = fixed_now$new_eval
    for (m in seq_along(active_cfgs)) {
      slot_key = hydra_slot_key(slot_id, i_t, active_cfgs[m], active_dup[m])
      if (!exists(slot_key, envir = state$slot_score_env, inherits = FALSE)) {
        need_eval = need_eval + 1L
      }
    }

    if (eval_used + need_eval > slot_budget) {
      used_T = t - 1L
      if (verbose) {
        cat(sprintf(
          "Hydra slot %d stops before instance %d: need=%d, used=%d, budget=%d\n",
          slot_id, i_t, need_eval, eval_used, slot_budget
        ))
      }
      break
    }

    vals_t = numeric(length(active_cfgs))
    for (m in seq_along(active_cfgs)) {
      base_seed = as.integer(seed0 + 10000L * i_t)
      out = hydra_eval_slot_score(
        state = state,
        cfg = active_cfgs[m],
        i = i_t,
        dup_order = active_dup[m],
        slot_id = slot_id,
        base_seed = base_seed,
        w = w,
        sigma1 = sigma1,
        sigma2 = sigma2,
        p = p,
        cat_penalty_notA = cat_penalty_notA
      )
      vals_t[m] = min(best_prev, out$score)
      eval_used = eval_used + out$new_eval
    }
    history[t, active_idx] = vals_t
    update_elite_seen(
      state,
      ports = as.list(active_cfgs),
      instances = i_t,
      sort_members = TRUE
    )

    if (verbose) {
      cat(sprintf(
        "Hydra slot=%d t=%d/%d i=%d alive=%d eval_used=%d\n",
        slot_id, t, length(training_instances), i_t, length(active_idx), eval_used
      ))
    }

    if (t < min_T_for_test) next
    if (length(active_idx) <= N_min) break
    if (((t - min_T_for_test) %% T_each) != 0L) next

    out_test = aux_friedman_irace(
      results = history[seq_len(t), active_idx, drop = FALSE],
      conf.level = 1 - alpha
    )
    keep_idx = which(out_test$alive)
    if (length(keep_idx) < length(active_idx)) {
      active_idx = active_idx[keep_idx]
    }
  }

  if (used_T <= 0L) stop("Hydra slot race ended before any instance was evaluated.")

  final_mat = history[seq_len(used_T), active_idx, drop = FALSE]
  if (ncol(final_mat) == 1L) {
    R_final = 1
  } else {
    R_final = colSums2_base(
      rowRanks_base(final_mat, cols = seq_len(ncol(final_mat)))
    )
  }
  if (length(active_idx) > 0L && used_T > 0L) {
    update_elite_seen(
      state,
      ports = as.list(cfgs[active_idx]),
      instances = used_instances,
      sort_members = TRUE
    )
  }

  elite_keys_final = vapply(as.list(cfgs[active_idx]), make_port_key, character(1), sort_members = TRUE)
  elite_e = vapply(elite_keys_final, function(k) get_port_e(state, k), integer(1))
  max_e = max(elite_e)
  cand_pos = which(elite_e == max_e)

  if (length(cand_pos) > 1L) {
    tie_keys = elite_keys_final[cand_pos]
    tie_cfgs = cfgs[active_idx[cand_pos]]

    tie_inst_list = lapply(tie_keys, function(k) {
      if (exists(k, envir = state$elite_seen_env, inherits = FALSE)) {
        as.integer(get(k, envir = state$elite_seen_env, inherits = FALSE))
      } else integer(0)
    })
    tie_inst_union = sort(unique(unlist(tie_inst_list)))

    R_tie = rep(Inf, length(cand_pos))
    if (length(tie_inst_union) > 0L) {
      F_rows = vector("list", length(tie_inst_union))
      keep_row = logical(length(tie_inst_union))

      for (ii in seq_along(tie_inst_union)) {
        i0 = as.integer(tie_inst_union[ii])
        fvec = rep(NA_real_, length(cand_pos))

        for (pp in seq_along(tie_cfgs)) {
          key_ic = make_ic_key(i0, tie_cfgs[pp])
          if (exists(key_ic, envir = state$score_env, inherits = FALSE)) {
            fvec[pp] = as.numeric(get(key_ic, envir = state$score_env, inherits = FALSE))
          }
        }

        if (all(is.finite(fvec))) {
          F_rows[[ii]] = fvec
          keep_row[ii] = TRUE
        }
      }

      if (any(keep_row)) {
        Fmat_tie = do.call(rbind, F_rows[keep_row])
        R_tie = colSums2_base(
          rowRanks_base(Fmat_tie, cols = seq_len(ncol(Fmat_tie)))
        )
      }
    }

    F_tie_mean = rep(Inf, length(cand_pos))
    if (exists("Fmat_tie") && is.matrix(Fmat_tie) && nrow(Fmat_tie) > 0L) {
      F_tie_mean = colMeans(Fmat_tie)
    }

    if (all(is.infinite(R_tie)) && length(R_final) == length(active_idx) && any(is.finite(R_final[cand_pos]))) {
      r_sub = R_final[cand_pos]
      r_sub[!is.finite(r_sub)] = Inf
      cand_pos = cand_pos[order(r_sub, F_tie_mean, cand_pos)]
    } else {
      cand_pos = cand_pos[order(R_tie, F_tie_mean, cand_pos)]
    }
  }

  best_local = cand_pos[1L]
  best_global_idx = active_idx[best_local]
  best_cfg = cfgs[best_global_idx]
  best_dup_order = dup_orders[best_global_idx]
  winner_scores = history[seq_len(used_T), best_global_idx]

  list(
    best_cfg = best_cfg,
    best_dup_order = best_dup_order,
    n_configs = n_configs,
    configs = cfgs,
    active_final = cfgs[active_idx],
    eval_used = eval_used,
    used_T = used_T,
    used_instances = used_instances,
    history = history[seq_len(max(1L, used_T)), , drop = FALSE],
    winner_scores = winner_scores,
    slot_budget = slot_budget,
    slot_objective = mean(winner_scores)
  )
}

run_hydra_wrapper_irace = function(param_config,
                                   target_portfolio_size,
                                   B_total,
                                   training_instances = NULL,
                                   alpha = 0.05,
                                   min_T_for_test = T_first,
                                   T_each = 1,
                                   N_min = N_min_default,
                                   w = 1,
                                   sigma1 = 1,
                                   sigma2 = 0.1,
                                   p = 10,
                                   seed0 = 1,
                                   race_seed0 = 999,
                                   lower = LB,
                                   upper = UB,
                                   cat_penalty_notA = 3,
                                   state = NULL,
                                   verbose = TRUE) {
  target_portfolio_size = as.integer(target_portfolio_size)
  if (target_portfolio_size <= 0L) stop("target_portfolio_size must be positive.")

  state = hydra_init_state(state)
  state$N_iter_total = target_portfolio_size
  if (is.null(training_instances)) {
    training_instances = seq_len(as.integer(min_T_for_test))
  }
  training_instances = as.integer(training_instances)
  training_instances = training_instances[is.finite(training_instances)]
  if (length(training_instances) == 0L) stop("training_instances must not be empty.")
  state$training_instances = training_instances

  budgets = compute_hydra_budget_split(B_total, target_portfolio_size)
  slots = vector("list", target_portfolio_size)

  for (k in seq_len(target_portfolio_size)) {
    if (verbose) {
      cat(sprintf(
        "\n=== Hydra iteration %d/%d | slot_budget=%d | instances=%d ===\n",
        k, target_portfolio_size, budgets[k], length(training_instances)
      ))
    }

    slot_res = run_hydra_slot_irace(
      param_config = param_config,
      state = state,
      slot_id = k,
      slot_budget = budgets[k],
      training_instances = training_instances,
      alpha = alpha,
      min_T_for_test = min_T_for_test,
      T_each = T_each,
      N_min = N_min,
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      seed0 = seed0,
      race_seed0 = race_seed0,
      lower = lower,
      upper = upper,
      cat_penalty_notA = cat_penalty_notA,
      verbose = verbose
    )

    state$portfolio_slots[[k]] = list(
      slot_id = as.integer(k),
      cfg = slot_res$best_cfg,
      dup_order = as.integer(slot_res$best_dup_order)
    )

    slot_res$portfolio_after = vapply(state$portfolio_slots, function(x) as.character(x$cfg), character(1))
    slot_res$best_so_far_after = vapply(slot_res$used_instances, function(i) {
      hydra_portfolio_score_from_cache(state, i)
    }, numeric(1))
    slots[[k]] = slot_res
  }

  portfolio_cfgs = vapply(state$portfolio_slots, function(x) x$cfg, numeric(1))
  portfolio_dup_orders = vapply(state$portfolio_slots, function(x) x$dup_order, integer(1))
  evaluated_instances = sort(unique(unlist(lapply(slots, function(x) x$used_instances))))
  best_so_far = if (length(evaluated_instances) > 0L) {
    vapply(evaluated_instances, function(i) hydra_portfolio_score_from_cache(state, i), numeric(1))
  } else {
    numeric(0)
  }

  list(
    portfolio = portfolio_cfgs,
    duplicate_orders = portfolio_dup_orders,
    best_so_far = best_so_far,
    evaluated_instances = evaluated_instances,
    training_instances = training_instances,
    slot_budgets = budgets,
    slot_results = slots,
    state = state,
    portfolio_objective = if (length(best_so_far) > 0L) mean(best_so_far) else NA_real_
  )
}

hydra_parse_single_cfg = function(best_name) {
  if (length(best_name) != 1L || is.na(best_name) || !nzchar(best_name)) return(NA_real_)
  s = gsub("[() ]", "", as.character(best_name))
  vals = suppressWarnings(as.numeric(strsplit(s, ",", fixed = TRUE)[[1L]]))
  vals = vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_real_)
  vals[1L]
}

hydra_next_dup_order = function(slots, cfg) {
  if (length(slots) == 0L) return(1L)
  ticks = vapply(slots, function(x) cfg_tick(x$cfg), integer(1))
  sum(ticks == cfg_tick(cfg)) + 1L
}

# Override the earlier prototype with a pure outer wrapper that delegates each
# Hydra iteration to the original solver (`run_iterated_portfolio_irace`).
run_hydra_wrapper_irace = function(param_config,
                                   target_portfolio_size,
                                   B_total,
                                   training_instances = NULL,
                                   alpha = 0.05,
                                   min_T_for_test = T_first,
                                   T_each = 1,
                                   N_min = N_min_default,
                                   w = 1,
                                   sigma1 = 1,
                                   sigma2 = 0.1,
                                   p = 10,
                                   seed0 = 1,
                                   race_seed0 = 999,
                                   lower = LB,
                                   upper = UB,
                                   state = NULL,
                                   verbose = TRUE,
                                   ...) {
  target_portfolio_size = as.integer(target_portfolio_size)
  if (target_portfolio_size <= 0L) stop("target_portfolio_size must be positive.")

  budgets = compute_hydra_budget_split(B_total, target_portfolio_size)
  selected_slots = list()
  outer_results = vector("list", target_portfolio_size)

  for (k in seq_len(target_portfolio_size)) {
    args = list(
      param_config = param_config,
      k_portfolio = 1L,
      alpha = alpha,
      state = state,
      min_T_for_test = min_T_for_test,
      T_each = T_each,
      w = w,
      sigma1 = sigma1,
      sigma2 = sigma2,
      p = p,
      seed0 = seed0,
      race_seed0 = race_seed0,
      B_total = budgets[k],
      N_min = N_min,
      training_instances = training_instances,
      hydra_fixed_slots = selected_slots,
      verbose = verbose,
      lower = lower,
      upper = upper
    )
    dots = list(...)
    if (length(dots) > 0L) args = modifyList(args, dots)

    res = do.call(run_iterated_portfolio_irace, args)
    state = res$cache_env
    if (is.null(state) && length(res$races) > 0L) {
      state = res$races[[length(res$races)]]$result$state
    }
    if (is.null(state)) {
      stop("Hydra wrapper could not recover solver state from run_iterated_portfolio_irace().")
    }

    cfg = hydra_parse_single_cfg(res$best_name)
    if (!is.finite(cfg)) {
      stop("Hydra wrapper expected a single-config best_name but got: ", res$best_name)
    }

    selected_slots[[k]] = list(
      slot_id = as.integer(k),
      cfg = as.numeric(cfg),
      dup_order = hydra_next_dup_order(selected_slots, cfg)
    )

    outer_results[[k]] = list(
      hydra_iteration = k,
      slot_budget = budgets[k],
      best_name = res$best_name,
      selected_cfg = cfg,
      solver_result = res
    )
  }

  list(
    portfolio = vapply(selected_slots, function(x) x$cfg, numeric(1)),
    duplicate_orders = vapply(selected_slots, function(x) x$dup_order, integer(1)),
    slot_budgets = budgets,
    slot_results = outer_results,
    selected_slots = selected_slots,
    state = state
  )
}

hydra_parse_port_name = function(x) {
  if (length(x) != 1L || is.na(x) || !nzchar(x)) return(numeric(0))
  s = gsub("[() ]", "", x)
  if (!nzchar(s)) return(numeric(0))
  as.numeric(strsplit(s, ",", fixed = TRUE)[[1]])
}

hydra_format_port_name_keep_dup = function(vals) {
  vals = as.numeric(vals)
  if (length(vals) == 0L || any(!is.finite(vals))) return(NA_character_)
  vals = sort(vals)
  paste0("(", paste(vals, collapse = ","), ")")
}

hydra_default_optimal_name_1d = function(k = 2L) {
  k = max(1L, as.integer(k))
  if (k >= 3L) "(49.4999,50,50)" else "(49.4999,50)"
}

hydra_deterministic_cost_config = function(c, i, w, p) {
  v = w * abs(c - 50)
  c_int = as.integer(round(c))
  i_int = as.integer(i)
  penalty = if ((c_int %% 2L) == (i_int %% 2L)) 0 else p
  v + penalty
}

hydra_sum_portfolio_cost_deterministic = function(port_vals, instances, seed0, w, p) {
  port_vals = as.numeric(port_vals)
  if (length(port_vals) == 0L) return(NA_real_)
  vals = vapply(as.integer(instances), function(i) {
    sc = vapply(port_vals, function(c) {
      hydra_deterministic_cost_config(c = c, i = i, w = w, p = p)
    }, numeric(1))
    min(sc)
  }, numeric(1))
  sum(vals)
}

hydra_sum_cost_from_port_name_deterministic = function(port_name, instances, seed0, w, p) {
  vals = hydra_parse_port_name(port_name)
  if (length(vals) == 0L || any(!is.finite(vals))) return(NA_real_)
  hydra_sum_portfolio_cost_deterministic(
    port_vals = vals,
    instances = instances,
    seed0 = seed0,
    w = w,
    p = p
  )
}

run_one_hydra_simple = function(settings_hydra,
                                w, sigma1, sigma2, p,
                                B_total,
                                seed_shift = 0L,
                                optimal_name = hydra_default_optimal_name_1d(),
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
    seed0 = settings_hydra$seed0 + as.integer(seed_shift),
    race_seed0 = settings_hydra$race_seed0 + as.integer(seed_shift)
  )
  args = modifyList(settings_hydra, overrides)
  res = do.call(run_hydra_wrapper_irace, args)

  seed0_run = settings_hydra$seed0 + as.integer(seed_shift)
  best_name = hydra_format_port_name_keep_dup(res$portfolio)
  n_races = length(res$slot_results)
  total_eval_used = sum(vapply(res$slot_results, function(x) x$solver_result$total_eval_used, numeric(1)))
  total_instances_used = sum(vapply(res$slot_results, function(x) x$solver_result$total_instances_used, numeric(1)))
  last_solver = if (n_races > 0L) res$slot_results[[n_races]]$solver_result else NULL
  last_race_result = if (!is.null(last_solver) && length(last_solver$races) > 0L) {
    last_solver$races[[length(last_solver$races)]]$result
  } else NULL
  alive_count_last = if (!is.null(last_race_result)) last_race_result$alive_count else NA_integer_

  best_cost_det = hydra_sum_cost_from_port_name_deterministic(
    port_name = best_name,
    instances = eval_instances,
    seed0 = seed0_run,
    w = w,
    p = p
  )
  optimal_cost_det = hydra_sum_cost_from_port_name_deterministic(
    port_name = optimal_name,
    instances = eval_instances,
    seed0 = seed0_run,
    w = w,
    p = p
  )
  cost_gap = best_cost_det - optimal_cost_det

  data.frame(
    w = w,
    sigma1 = sigma1,
    sigma2 = sigma2,
    p = p,
    B_total = B_total,
    seed_shift = seed_shift,
    best_name = best_name,
    n_races = n_races,
    total_eval_used = total_eval_used,
    total_instances_used = total_instances_used,
    alive_count = alive_count_last,
    best_cost = best_cost_det,
    best_k = as.integer(best_k),
    best_k_names = best_name,
    best_k_union = best_name,
    best_k_sum_cost = best_cost_det,
    best_k_sum_costgap = if (best_k > 1L) cost_gap else NA_real_,
    optimal_name = optimal_name,
    optimal_cost = optimal_cost_det,
    cost_gap = cost_gap,
    best_is_12 = (!is.na(best_name)) && best_name == optimal_name,
    unique_12 = (!is.na(best_name)) && best_name == optimal_name,
    best_name_decoded = NA_character_,
    best_member_win_counts = NA_character_,
    stringsAsFactors = FALSE
  )
}

test_hyper_hydra_simple = function(w_list, sigma1_list, sigma2_list, p_list,
                                   B_list, seed_list,
                                   settings_hydra,
                                   optimal_name = hydra_default_optimal_name_1d(),
                                   eval_instances = c(1, 2),
                                   best_k = 1L,
                                   k_best = NULL) {
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
    out[[i]] = run_one_hydra_simple(
      settings_hydra = settings_hydra,
      w = w_list[i],
      sigma1 = sigma1_list[i],
      sigma2 = sigma2_list[i],
      p = p_list[i],
      B_total = B_list[i],
      seed_shift = seed_list[i],
      optimal_name = optimal_name,
      eval_instances = eval_instances,
      best_k = best_k
    )
  }
  do.call(rbind, out)
}

run_grid_hydra_simple = function(param_grid, settings_hydra,
                                 optimal_name = hydra_default_optimal_name_1d(),
                                 eval_instances = c(1, 2),
                                 best_k = 1L,
                                 k_best = NULL) {
  if (!is.null(k_best)) best_k = k_best
  best_k = max(1L, as.integer(best_k))
  test_hyper_hydra_simple(
    w_list = param_grid$w,
    sigma1_list = param_grid$sigma1,
    sigma2_list = param_grid$sigma2,
    p_list = param_grid$p,
    B_list = param_grid$B_total,
    seed_list = param_grid$seed_shift,
    settings_hydra = settings_hydra,
    optimal_name = optimal_name,
    eval_instances = eval_instances,
    best_k = best_k,
    k_best = best_k
  )
}
