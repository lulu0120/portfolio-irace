
#what we will see for first test
get_instances_for_first_test = function(state, race_id, Tnew, race_seed, min_T_for_test,
                                        training_instances = NULL) {
  stopifnot(min_T_for_test >= 1L)
  race_id = as.integer(race_id)
  Tnew = as.integer(Tnew)
  race_seed = as.integer(race_seed)
  
  if (race_id <= 1L) {
    if (is.null(training_instances) || length(training_instances) == 0L) {
      return(list(instances = seq_len(min_T_for_test), repeats = rep.int(0L, min_T_for_test), reuse_new_ids = rep.int(0L, min_T_for_test)))
    }
    base_pool = as.integer(training_instances)
    take = min(min_T_for_test, length(base_pool))
    inst = base_pool[seq_len(take)]
    repv = rep.int(0L, take)
    reusev = rep.int(0L, take)
    if (take < min_T_for_test) {
      rest = min_T_for_test - take
      counter = setNames(vapply(unique(base_pool), function(i) get_instance_repeat_count(state, i), integer(1)), as.character(unique(base_pool)))
      set.seed(race_seed)
      extra = sample(base_pool, size = rest, replace = TRUE)
      extra_rep = integer(rest)
      extra_reuse = integer(rest)
      for (idx in seq_along(extra)) {
        key = as.character(extra[idx])
        counter[[key]] = as.integer(counter[[key]] + 1L)
        extra_rep[idx] = as.integer(counter[[key]])
        extra_reuse[idx] = next_reuse_new_id(state)
      }
      inst = c(inst, extra)
      repv = c(repv, extra_rep)
      reusev = c(reusev, extra_reuse)
    }
    return(list(instances = inst, repeats = repv, reuse_new_ids = reusev))
  }
  
  q = make_instance_queue(
    state,
    race_id = race_id,
    Tnew = Tnew,
    seed = race_seed,
    training_instances = training_instances,
    min_required = min_T_for_test
  )
   
  stopifnot(length(q$instances) >= min_T_for_test)
  list(
    instances = q$instances[seq_len(min_T_for_test)],
    repeats = q$repeats[seq_len(min_T_for_test)],
    reuse_new_ids = q$reuse_new_ids[seq_len(min_T_for_test)]
  )
}

#what was seen
count_missing_ic_pairs = function(state, instances, configs, repeats = NULL, reuse_new_ids = NULL) {
  instances = as.integer(instances)
  if (is.null(repeats)) repeats = rep.int(0L, length(instances))
  repeats = as.integer(repeats)
  if (is.null(reuse_new_ids)) reuse_new_ids = rep.int(0L, length(instances))
  reuse_new_ids = as.integer(reuse_new_ids)
  configs_tick = cfg_tick(configs)
  sum(vapply(seq_along(instances), function(idx) {
    i = instances[idx]
    r = repeats[idx]
    u = reuse_new_ids[idx]
    keys = vapply(configs_tick, function(ct) paste0("i=", as.integer(i), "|r=", as.integer(r), "|u=", as.integer(u), "|c=", as.integer(ct)), character(1))
    sum(!vapply(keys, exists, logical(1), envir = state$score_env, inherits = FALSE))
  }, integer(1)))
}

computeNbConfigurations = function(B_j, mu, T_each, j) {
  denom = mu + T_each * min(5L, j)
  Nj = floor(B_j / denom) 
}

rank_stop_elites = function(last_race) {
  if (is.null(last_race) ||
      is.null(last_race$alive_portfolios) ||
      length(last_race$alive_portfolios) == 0L ||
      is.null(last_race$state)) {
    return(list(
      ordered_names = character(0),
      seen_counts = integer(0)
    ))
  }
  state = last_race$state
  ports = last_race$alive_portfolios
  names0 = as.character(last_race$alive_names)
  keys0 = vapply(ports, function(P) make_port_key(P, sort_members = TRUE), character(1))
  elite_e = vapply(keys0, function(k) get_port_e(state, k), integer(1))
  cand_pos = seq_along(keys0)
  ordered = integer(0)
  count_levels = sort(unique(elite_e), decreasing = TRUE)
  
  for (count_i in count_levels) {
    pos_i = cand_pos[elite_e[cand_pos] == count_i]
    if (length(pos_i) <= 1L) {
      ordered = c(ordered, pos_i)
      next
    }
    tie_keys = keys0[pos_i]
    tie_ports = ports[pos_i]
    tie_inst_list = lapply(tie_keys, function(k) {
      if (exists(k, envir = state$elite_seen_env, inherits = FALSE)) {
        occ_keys = unique(as.character(get(k, envir = state$elite_seen_env, inherits = FALSE)))
        parsed = parse_ir_key(occ_keys)
        keep = is.finite(parsed$i) & is.finite(parsed$r) & is.finite(parsed$u)
        if (!any(keep)) {
          data.frame(i = integer(0), r = integer(0), u = integer(0), stringsAsFactors = FALSE)
        } else {
          parsed[keep, , drop = FALSE]
        }
      } else integer(0)
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
        R_tie = colSums2_base(
          rowRanks_base(Fmat_tie, cols = seq_len(ncol(Fmat_tie)))
        )
        F_tie_mean = colMeans(Fmat_tie)
      }
    }
    tie_min_tick = vapply(tie_ports, function(P) min(cfg_tick(P)), integer(1))
    ranks_now = suppressWarnings(as.numeric(last_race$alive_R[pos_i]))
    if (all(is.infinite(R_tie)) && any(is.finite(ranks_now))) {
      ranks_now[!is.finite(ranks_now)] = Inf
      ord_i = order(ranks_now, F_tie_mean, tie_min_tick, names0[pos_i])
    } else {
      ord_i = order(R_tie, F_tie_mean, tie_min_tick, names0[pos_i])
    }
    ordered = c(ordered, pos_i[ord_i])
  }
  
  list(
    ordered_names = names0[ordered],
    seen_counts = as.integer(elite_e[ordered])
  )
}

 
 


run_iterated_portfolio_irace = function(
    param_config,
    k_portfolio = port_size_k,
    alpha = 0.05,
    max_T = 200,
    state = NULL,
    min_T_for_test = T_first,
    T_each = 1,
    w = 1, sigma1 = 1, sigma2 = 0.1, p = 10,
    seed0 = 1,
    race_seed0 = 999,
    B_total,
    
    N_iter = N_min_default,
    weight_type = c("uniform", "portfolio_count", "best_count", "keep"), #目前不用这两个
    domain_update = c("alive_union", "keep_all"),
    eps_counter = 1,
    N_min = N_min_default,
    enable_second_level_test = TRUE,
    global_free_rider_elimination = TRUE,
    enable_duplicate = TRUE,
    hydra_fixed_slots = NULL,
    reuse_elite_portfolios = TRUE,
    stop_if_nj_le = NULL,
    verbose = TRUE,
    training_instances = NULL,
    lower = LB,
    upper = UB
) {
  cat(sprintf(
    "start race | seed=%d | budget=%d | portsize=%d\n",
    as.integer(seed0), as.integer(B_total), as.integer(k_portfolio)
  ))
  if (is.null(state)) state = init_race_state()
  state$N_iter_total = as.integer(N_iter)
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
  } else {
    sd0 = (upper - lower) / 2
  }
  #weight_type   = match.arg(weight_type)
  #domain_update = match.arg(domain_update)
  
  races = list()
  
  param_cur   = param_config
  
  if (param_type %in% c("discrete","class","categorical")) {
    if (verbose) cat("enter categorical mode (multidim generator) \n")
    full_values = NULL
    seen_cfgs   = numeric(0)
    cnt_state = NULL
  } else if (param_type %in% c("multi", "multidim", "multivariate")) {
    if (verbose) cat("enter multidim mode \n")
    full_values = NULL
    seen_cfgs   = numeric(0)
    cnt_state = NULL
  } else {
    if (verbose) cat("enter cts mode \n")
    full_values = NULL
    seen_cfgs   = numeric(0)
    cnt_state = NULL   # cts 模式：先不做 counter/domain update
  }
  
  
  total_eval_used      = 0L
  total_instances_used = 0L
  
 
  
  fixed_portfolios = NULL
  res_prev = NULL
  j = 0L
  
  repeat {
    
    remaining_budget = B_total - total_eval_used
    if (remaining_budget <= 0) {
      if (verbose) cat("Stop: total budget exhausted\n")
      break
    }
    
    j = j + 1L
    
    n_races_left = max(1L, as.integer(N_iter) - j + 1L) #超过estimated后，每个race都可以有所有的预算
    
    B_j = floor(remaining_budget / n_races_left)
    if (B_j <= 0) {
      if (verbose) cat("Stop: Bj <= 0\n")
      break
    }
    
    n_configs = computeNbConfigurations(
      B_j   = B_j,
      mu    = min_T_for_test,
      T_each = T_each,
      j     = j
    )
    
    if (param_type %in% c("discrete","class","categorical")) {
      n_configs = min(as.integer(n_configs), length(param_cur$values))
    } else {
      n_configs = as.integer(n_configs)
    }
    
    if (n_configs <= k_portfolio) {
      if (verbose) cat("Stop: n_configs <= k_portfolio\n")
      break
    }

    if (!is.null(stop_if_nj_le) &&
        is.finite(stop_if_nj_le) &&
        as.integer(n_configs) <= as.integer(stop_if_nj_le)) {
      if (verbose) cat(sprintf("Stop: Nj=%d <= best_k=%d\n", as.integer(n_configs), as.integer(stop_if_nj_le)))
      break
    }
    
    race_seed = as.integer(race_seed0 + j - 1L)
    
    if (verbose) {
      cat(sprintf("\n=== Race %d | Bj=%d | Nj=%d ===\n",
                  j, B_j, n_configs))
    }
    
    # hard gate: Nj must cover required cfgs (if reuse_elite_portfolios)
    
    required_cfgs_now = if (reuse_elite_portfolios &&
                            !is.null(fixed_portfolios) &&
                            length(fixed_portfolios) > 0L) {
      # fixed_portfolios 是 list of numeric configs
      cfg_from_tick(sort(unique(unlist(lapply(fixed_portfolios, cfg_tick)))))  # 返回 numeric，但 identity 不丢
    } else {
      numeric(0)
    }
    
    
    if (length(required_cfgs_now) > 0L && n_configs < length(required_cfgs_now)) {
      if (verbose) {
        cat(sprintf(
          "Stop: Nj=%d < |required_cfgs|=%d, cannot include all required configs\n",
          n_configs, length(required_cfgs_now)
        ))
      }
      break
    }
 
    
    param_before  = param_cur
    cnt_before    = cnt_state
    fixed_before  = fixed_portfolios
 
    # BUDGET GATE: must be able to reach FIRST test in this race
    
    Nj    = as.integer(n_configs)
    Tneed = as.integer(min_T_for_test)
    
    inst_first_test = get_instances_for_first_test(
      state = state,
      race_id = j,
      Tnew = 1L,
      race_seed = race_seed,
      min_T_for_test = Tneed,
      training_instances = training_instances
    )
    
 
    # (A) count what was not seen by elites
    need_req = 0L
    if (length(required_cfgs_now) > 0L) {
      need_req = count_missing_ic_pairs(
        state = state,
        instances = inst_first_test$instances,
        configs   = required_cfgs_now
        ,
        repeats = inst_first_test$repeats,
        reuse_new_ids = inst_first_test$reuse_new_ids
      )
    }
    
    # (B)  remaining unknown 
    n_new_cfgs = max(0L, Nj - length(required_cfgs_now))
    need_new_upper = Tneed * n_new_cfgs
    
    need_hydra_fixed = 0L
    if (!is.null(hydra_fixed_slots) && length(hydra_fixed_slots) > 0L) {
      need_hydra_fixed = md_hydra_count_missing_fixed_slot_scores(
        state = state,
        instances = inst_first_test$instances,
        hydra_fixed_slots = hydra_fixed_slots,
        repeats = inst_first_test$repeats,
        reuse_new_ids = inst_first_test$reuse_new_ids
      )
    }
    
    need_eval_to_test = as.integer(need_req + need_new_upper + need_hydra_fixed)
    
    #gate for budget
    if (need_eval_to_test > remaining_budget) { 
      if (verbose) cat(sprintf(
        "Stop: remaining budget insufficient to reach first test (need=%d, remaining=%d, race=%d)\n",
        need_eval_to_test, remaining_budget, j
      ))
      break
    }
    # at least try to finish this race
    if (need_eval_to_test > B_j) {
      if (verbose) cat(sprintf(
        "Adjust: Bj -> remaining to finish this race (need=%d, Bj=%d, remaining=%d, race=%d)\n",
        need_eval_to_test, B_j, remaining_budget, j
      ))
      B_j = remaining_budget
    }

    if (need_eval_to_test > B_j) stop("Internal error: should not happen.")
 
    res = run_portfolio_race(
      param_config = param_cur,
      n_configs = Nj,
      elite_port_rankpos = if (!is.null(res_prev)) res_prev$alive_rankpos else NULL,
      sd0 = sd0,
      k_portfolio = k_portfolio,
      alpha = alpha,
      max_T = max_T,
      state   = state,
      race_id = j,
      Tnew    = 1L,
      min_T_for_test = min_T_for_test,
      T_each = T_each,
      max_eval = B_j,
      stop_when_budget_cannot_cover_next_instance = TRUE,
      w = w, sigma1 = sigma1, sigma2 = sigma2, p = p,
      seed0 = seed0,
      race_seed = race_seed,
      fixed_portfolios = if (reuse_elite_portfolios) fixed_portfolios else NULL,
      hydra_fixed_slots = hydra_fixed_slots,
      N_min = N_min,
      enable_second_level_test = enable_second_level_test,
      global_free_rider_elimination = global_free_rider_elimination,
      enable_duplicate = enable_duplicate,
      verbose = verbose,
      training_instances = training_instances
    )
    
    state = res$state
    sd0 = res$sd0_next
    res_prev = res

    if (isTRUE(enable_duplicate)) {
      gate_no_freerider = isTRUE(res$best_group_no_freerider)
      dup_keep_tick = sort(unique(as.integer(res$best_group_dup_keep_tick)))
      marked_tick = sort(unique(as.integer(res$allow_mark_tick)))
      must_tick_out = if (length(res$must_dup_tick) > 0L) as.integer(res$must_dup_tick[1L]) else integer(0)
      allow_cnt = state$allow_dup_count
      if (length(allow_cnt) > 0L) {
        allow_cnt = as.integer(allow_cnt)
        names(allow_cnt) = names(state$allow_dup_count)
      } else {
        allow_cnt = integer(0)
      }
      if (length(marked_tick) > 0L) {
        for (tt in marked_tick) {
          key = as.character(tt)
          prev = if (key %in% names(allow_cnt)) as.integer(allow_cnt[[key]]) else 1L
          allow_cnt[[key]] = min(as.integer(k_portfolio), prev + 1L)
        }
      }
      state$allow_dup_count = allow_cnt
      state$allow_dup_tick_active = marked_tick
      state$must_dup_tick = must_tick_out
    } else {
      marked_tick = integer(0)
      gate_no_freerider = FALSE
      dup_keep_tick = integer(0)
      state$allow_dup_count = integer(0)
      state$allow_dup_tick_active = integer(0)
      state$must_dup_tick = integer(0)
    }
    if (isTRUE(getOption("DBG2", FALSE))) {
      cat("\n[allow_duplicate_update]\n")
      cat("  race_id =", j, "\n")
      cat("  marked_ticks =", if (length(marked_tick) == 0L) "<none>" else paste(marked_tick, collapse = ","), "\n")
      cat("  allow_dup_tick =", if (length(state$allow_dup_tick_active) == 0L) "<none>" else paste(state$allow_dup_tick_active, collapse = ","), "\n")
      cat("  must_dup_tick =", if (length(state$must_dup_tick) == 0L) "<none>" else paste(state$must_dup_tick, collapse = ","), "\n")
      cat("  gate_no_freerider =", gate_no_freerider, "\n")
      cat("  gate_dup_keep_tick =", if (length(dup_keep_tick) == 0L) "<none>" else paste(dup_keep_tick, collapse = ","), "\n")
    }
    
    total_eval_used      = total_eval_used + res$eval_used
    total_instances_used = total_instances_used + res$used_instances
    
    if (!(param_type %in% c("discrete","class","categorical"))) {
      # cts: 不做 domain shrinking / weights update（保持 range）
      # param_cur 保持不变
    }
    
    
    if (reuse_elite_portfolios) {
      fixed_portfolios = res$alive_portfolios
    }
    
    races[[j]] = list(
      race_id = j,
      B_j = B_j,
      n_configs = Nj,
      
      param_before = param_before,
      weights_before = param_before$weights,
      cnt_state_before = cnt_before,
      fixed_portfolios_before = fixed_before,
      
      result = res,
      
      param_after = param_cur,
      weights_after = param_cur$weights,
      cnt_state_after = cnt_state,
      fixed_portfolios_after = fixed_portfolios,
      
      total_eval_used = total_eval_used,
      total_instances_used = total_instances_used
    )
  }
  
  n_done = length(races)
  last_race = if (n_done == 0L) NULL else races[[n_done]]$result
  ranked_stop = rank_stop_elites(last_race)
  stop_ranked_names = as.character(ranked_stop$ordered_names)
  best_last = if (length(stop_ranked_names) > 0L) {
    stop_ranked_names[[1L]]
  } else if (n_done > 0L) {
    races[[n_done]]$result$best_name
  } else {
    NA_character_
  }
  
  list(
    best_name = best_last,
    stop_ranked_names = stop_ranked_names,
    stop_ranked_seen_counts = as.integer(ranked_stop$seen_counts),
    races = races,
    state = state,
    final_param = param_cur,
    cnt_state = cnt_state,
    total_eval_used = total_eval_used,
    total_instances_used = total_instances_used
  )
}
