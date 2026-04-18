hydra_slot_key_local = function(slot_id, i, cfg, dup_order, repeat_id = 0L, reuse_new_id = 0L) {
  paste0(
    "hslot=", as.integer(slot_id),
    "|i=", as.integer(i),
    "|r=", as.integer(repeat_id),
    "|u=", as.integer(reuse_new_id),
    "|c=", cfg_tick(cfg),
    "|d=", as.integer(dup_order)
  )
}

hydra_call_scores_configs_local = function(configs, i, base_seed,
                                           w, sigma1, sigma2, p,
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
  if ("state" %in% sc_formals) args$state = state
  do.call(scores_configs_on_instance, args)
}

hydra_call_cost_local = function(c, i, w, p, sigma2,
                                 instance_random = 0,
                                 eps = NULL,
                                 duplicate_tick = NULL,
                                 seed_shift = NULL,
                                 state = NULL) {
  if (is.null(eps) && !is.null(duplicate_tick) && !is.null(seed_shift)) {
    # Match the duplicate handling used in the non-Hydra portfolio race:
    # each duplicate slot gets a deterministic seed offset, but the inner
    # scoring function still receives an ordinary eps draw.
    dup_order = as.integer(duplicate_tick) + 1L
    tick_c = as.integer(cfg_tick(c))
    set.seed(as.integer(seed_shift) + 100L * dup_order + abs(tick_c) %% 1000000L)
    eps = rnorm(1, mean = 0, sd = sigma2)
  }
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
  if ("duplicate_tick" %in% cost_formals && !is.null(duplicate_tick)) args$duplicate_tick = duplicate_tick
  if ("seed_shift" %in% cost_formals && !is.null(seed_shift)) args$seed_shift = seed_shift
  if ("state" %in% cost_formals) args$state = state
  do.call(cost_config, args)
}

md_hydra_count_missing_fixed_slot_scores = function(state, instances, hydra_fixed_slots, repeats = NULL, reuse_new_ids = NULL) {
  instances = as.integer(instances)
  if (is.null(repeats)) repeats = rep.int(0L, length(instances))
  repeats = as.integer(repeats)
  if (is.null(reuse_new_ids)) reuse_new_ids = rep.int(0L, length(instances))
  reuse_new_ids = as.integer(reuse_new_ids)
  if (length(hydra_fixed_slots) == 0L) return(0L)
  sum(vapply(seq_along(instances), function(idx) {
    i = instances[idx]
    r = repeats[idx]
    u = reuse_new_ids[idx]
    sum(vapply(hydra_fixed_slots, function(slot) {
      key = hydra_slot_key_local(slot$slot_id, i, slot$cfg, slot$dup_order, repeat_id = r, reuse_new_id = u)
      !exists(key, envir = state$score_env, inherits = FALSE)
    }, logical(1)))
  }, integer(1)))
}

md_hydra_eval_fixed_portfolio_score = function(state, hydra_fixed_slots, i_t, base_seed_t,
                                            w, sigma1, sigma2, p, repeat_id = 0L, reuse_new_id = 0L) {
  if (length(hydra_fixed_slots) == 0L) return(list(score = Inf, new_eval = 0L))

  key_inst = paste0("i=", as.integer(i_t), "|r=", as.integer(repeat_id), "|u=", as.integer(reuse_new_id), "|inst_rand")
  if (exists(key_inst, envir = state$score_env, inherits = FALSE)) {
    inst_rand = as.numeric(get(key_inst, envir = state$score_env, inherits = FALSE))
  } else {
    set.seed(as.integer(base_seed_t))
    inst_rand = rnorm(1, mean = 0, sd = sigma1)
    assign(key_inst, as.numeric(inst_rand), envir = state$score_env)
  }

  scores = numeric(length(hydra_fixed_slots))
  eval_added = 0L
  for (idx in seq_along(hydra_fixed_slots)) {
    slot = hydra_fixed_slots[[idx]]
    key = hydra_slot_key_local(slot$slot_id, i_t, slot$cfg, slot$dup_order, repeat_id = repeat_id, reuse_new_id = reuse_new_id)
    if (exists(key, envir = state$score_env, inherits = FALSE)) {
      scores[idx] = as.numeric(get(key, envir = state$score_env, inherits = FALSE))
      next
    }

    if (as.integer(slot$dup_order) <= 1L) {
      out = hydra_call_scores_configs_local(
        configs = c(slot$cfg),
        i = i_t,
        base_seed = base_seed_t,
        w = w,
        sigma1 = sigma1,
        sigma2 = sigma2,
        p = p,
        state = state
      )[1L]
      assign(key, as.numeric(out), envir = state$score_env)
      scores[idx] = as.numeric(out)
      eval_added = eval_added + 1L
    } else {
      key_icd = make_icd_key(i_t, slot$cfg, slot$dup_order, r = repeat_id, u = reuse_new_id)
      if (exists(key_icd, envir = state$score_env, inherits = FALSE)) {
        val = as.numeric(get(key_icd, envir = state$score_env, inherits = FALSE))
      } else {
        val = hydra_call_cost_local(
          c = slot$cfg,
          i = i_t,
          w = w,
          p = p,
          sigma2 = sigma2,
          instance_random = inst_rand,
          duplicate_tick = as.integer(slot$dup_order - 1L),
          seed_shift = as.integer(base_seed_t),
          state = state
        )
        assign(key_icd, as.numeric(val), envir = state$score_env)
      }
      assign(key, as.numeric(val), envir = state$score_env)
      scores[idx] = as.numeric(val)
      eval_added = eval_added + 1L
    }
  }
  list(score = min(scores), new_eval = eval_added)
}

run_portfolio_race = function(param_config,
                              n_configs,
                              k_portfolio = port_size_k,
                              elite_port_rankpos = NULL,
                              sd0 = NULL,
                              alpha = 0.05,
                              max_T = 200,
                              min_T_for_test = T_first,
                              
                              # ---- cache/state ----
 
                              state = NULL,              # NEW: race state (score_env + port_seen_env)
                              race_id = 1L,              # NEW: 第几轮 race（>=2 才启动 instance schedule）
                              Tnew = 1L,                 # NEW: 每个新 race prepend 的新 instance 数
                              
                              # ---- test schedule ----
                              T_each = 1,
                              T_minnodrop = min_T_for_test + 1,
                              Tmaxnodrop = 2,
                              
                              # ---- budget ----
                              max_eval = Inf,
                              stop_when_budget_cannot_cover_next_instance = TRUE,
                              
                              # ---- toy cost ----
                              w = 1, sigma1 = 1, sigma2 = 0.1, p = 10,
                              seed0 = 1,
                              
                              # ---- elites ----
                              fixed_portfolios = NULL,   # 上轮 elite portfolios（用于：保留 & e）
                              hydra_fixed_slots = NULL,  # fixed historical slots from outer Hydra wrapper
                              race_seed = 999,
                              N_min = N_min_default,
                              enable_second_level_test = TRUE,
                              global_free_rider_elimination = TRUE,
                              enable_duplicate = TRUE,
                              verbose = TRUE,
                              training_instances = NULL) {
  build_portfolios_multiset = function(configs_vec, k_portfolio) {
    vals = as.numeric(configs_vec)
    keys = as.character(cfg_tick(vals))
    ord = order(vals)
    vals = vals[ord]
    keys = keys[ord]
    keys_u = unique(keys)
    reps = vapply(keys_u, function(k) vals[match(k, keys)], numeric(1))
    caps = as.integer(tabulate(match(keys, keys_u), nbins = length(keys_u)))
    n_u = length(reps)
    if (n_u == 0L || k_portfolio <= 0L) return(list())
    out = list()
    cur = integer(n_u)
    rec = function(pos, remain) {
      if (pos > n_u) {
        if (remain == 0L) {
          out[[length(out) + 1L]] <<- as.numeric(rep(reps, times = cur))
        }
        return(invisible(NULL))
      }
      take_max = min(caps[pos], remain)
      for (take in 0:take_max) {
        cur[pos] <<- as.integer(take)
        rec(pos + 1L, remain - take)
      }
      cur[pos] <<- 0L
      invisible(NULL)
    }
    rec(1L, as.integer(k_portfolio))
    out
  }

  portfolio_scores_with_duplicate = function(portfolios, cfg_scores, i_t, base_seed_t,
                                             repeat_id, reuse_new_id, state, w, sigma1, sigma2, p) {
    # Keep per-instance shared random term identical to scores_configs_on_instance().
    key_inst = paste0("i=", as.integer(i_t), "|r=", as.integer(repeat_id), "|u=", as.integer(reuse_new_id), "|inst_rand")
    if (exists(key_inst, envir = state$score_env, inherits = FALSE)) {
      inst_rand = as.numeric(get(key_inst, envir = state$score_env, inherits = FALSE))
    } else {
      set.seed(as.integer(base_seed_t))
      inst_rand = rnorm(1, mean = 0, sd = sigma1)
      assign(key_inst, as.numeric(inst_rand), envir = state$score_env)
    }
    
    vapply(portfolios, function(P) {
      cfg_vals = as.numeric(P)
      cfg_ticks = as.integer(cfg_tick(cfg_vals))
      dup_idx = ave(cfg_ticks, cfg_ticks, FUN = seq_along)
      sc = numeric(length(P))
      for (k in seq_along(P)) {
        cfg_k = cfg_vals[k]
        tick_k = cfg_ticks[k]
        d = as.integer(dup_idx[k])
        if (d <= 1L) {
          sc[k] = cfg_scores[as.character(tick_k)]
        } else {
          key_icd = make_icd_key(i_t, cfg_k, d, r = repeat_id, u = reuse_new_id)
          if (exists(key_icd, envir = state$score_env, inherits = FALSE)) {
            sc[k] = as.numeric(get(key_icd, envir = state$score_env, inherits = FALSE))
          } else {
            val = hydra_call_cost_local(
              c = cfg_k,
              i = i_t,
              w = w,
              p = p,
              sigma2 = sigma2,
              instance_random = inst_rand,
              duplicate_tick = as.integer(d - 1L),
              seed_shift = as.integer(base_seed_t),
              state = state
            )
            assign(key_icd, as.numeric(val), envir = state$score_env)
            sc[k] = as.numeric(val)
          }
        }
      }
      min(sc)
    }, numeric(1))
  }
  
  if (is.null(state)) state = init_race_state()
  if (is.null(state$score_env))     state$score_env     = new.env(parent = emptyenv())
  if (is.null(state$port_seen_env)) state$port_seen_env = new.env(parent = emptyenv())
  if (is.null(state$elite_seen_env)) state$elite_seen_env = new.env(parent = emptyenv())
  if (is.null(state$ports_all_env)) state$ports_all_env = new.env(parent = emptyenv())
  if (is.null(state$sd_env))        state$sd_env        = new.env(parent = emptyenv())   
  if (is.null(state$sd_env)) state$sd_env = new.env(parent = emptyenv())
  if (is.null(state$allow_dup_count)) state$allow_dup_count = integer(0)
  if (is.null(state$allow_dup_tick_active)) state$allow_dup_tick_active = integer(0)
  if (is.null(state$must_dup_tick)) state$must_dup_tick = integer(0)

  
  # ---------- generate configs for this race (numeric-irace style) ----------
  #required_cfgs = integer(0)
  #if (!is.null(fixed_portfolios) && length(fixed_portfolios) > 0) {
  #  required_cfgs = sort(unique(as.integer(unlist(fixed_portfolios))))
  #}
  
  elite_cfg_info = config_ranks_from_elite_portfolios(
    elite_portfolios = fixed_portfolios,
    elite_port_rankpos = elite_port_rankpos
  )
  if (isTRUE(verbose) &&
      as.integer(race_id) > 1L &&
      !is.null(fixed_portfolios) &&
      length(fixed_portfolios) > 0L) {
    elite_labels = vapply(
      fixed_portfolios,
      function(P) paste0("(", paste(as.numeric(P), collapse = ","), ")"),
      character(1)
    )
    cat("[race_elite_set] ")
    cat(sprintf("race_id=%d | n_elite=%d\n", as.integer(race_id), length(fixed_portfolios)))
    cat("  elites=", paste(elite_labels, collapse = " | "), "\n", sep = "")
    if (!is.null(elite_port_rankpos) && length(elite_port_rankpos) == length(fixed_portfolios)) {
      cat("  elite_rankpos=", paste(as.numeric(elite_port_rankpos), collapse = ","), "\n", sep = "")
    }
  }
  DBG <- isTRUE(getOption("DBG", FALSE))
  DBG2 <- isTRUE(getOption("DBG2", FALSE))
  if (DBG) {
    cat("\n[run_portfolio_race]\n")
    cat("  race_id =", race_id, "\n")
    cat("  n_configs =", n_configs, "\n")
    cat("  fixed_portfolios_len =", if (is.null(fixed_portfolios)) 0 else length(fixed_portfolios), "\n")
    cat("  elite_port_rankpos_len =", if (is.null(elite_port_rankpos)) 0 else length(elite_port_rankpos), "\n")
    if (!is.null(elite_port_rankpos)) cat("  elite_port_rankpos_hasNA =", any(is.na(elite_port_rankpos)), "\n")
    #cat("  required_cfgs_from_fixed =", paste(required_cfgs, collapse=","), "\n")
  }
  type_now = tolower(as.character(param_config$type))
  round_now = ! (type_now %in% c("cts", "numeric", "real", "continuous"))
  allow_dup_tick_active = if (!isTRUE(enable_duplicate) || as.integer(race_id) <= 1L) integer(0) else as.integer(state$allow_dup_tick_active)
  allow_dup_count = state$allow_dup_count
  must_dup_tick = if (!isTRUE(enable_duplicate) || as.integer(race_id) <= 1L) integer(0) else as.integer(state$must_dup_tick)
  DBG_duplicate <- isTRUE(getOption("DBG_duplicate", FALSE))
  if (DBG_duplicate) {
    cnt_show = integer(0)
    if (length(allow_dup_tick_active) > 0L && length(allow_dup_count) > 0L) {
      keys = as.character(allow_dup_tick_active)
      cnt_show = as.integer(allow_dup_count[keys])
      cnt_show[is.na(cnt_show)] = 1L
      names(cnt_show) = keys
    }
    cat("\n[duplicate_race_start]\n")
    cat("  race_id =", race_id, "\n")
    cat("  allow_dup_tick_active =",
        if (length(allow_dup_tick_active) == 0L) "<none>" else paste(allow_dup_tick_active, collapse = ","), "\n")
    cat("  allow_dup_cap =",
        if (length(cnt_show) == 0L) "<none>" else paste(sprintf("%s:%d", names(cnt_show), cnt_show), collapse = ", "), "\n")
    cat("  must_dup_tick =",
        if (length(must_dup_tick) == 0L) "<none>" else paste(must_dup_tick, collapse = ","), "\n")
  }
  
  gen = generate_configs_numeric_irace_with_sd(
    param_config = param_config,          # 让 wrapper 自己分发 discrete/cts 
    #比如param_config = list(type="discrete", values=1:100)或者param_config = list(type="cts", lower=1, upper=100)


    Nj = n_configs,
    required_cfgs = elite_cfg_info$cfgs,
    elite_ranks = elite_cfg_info$ranks,
    j = race_seed,
    sd0 = sd0,
    lower = LB,
    upper = UB,
    round = round_now,
    n_param = N_param,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = allow_dup_tick_active,
    allow_dup_count = allow_dup_count,
    must_dup_tick = must_dup_tick,
    k_portfolio = k_portfolio,
    state = state,
    iter_idx = race_id
  )
  
  
  configs_race = gen$configs
  portfolios = build_portfolios_multiset(configs_race, k_portfolio)
  
  DBG <- isTRUE(getOption("DBG", FALSE))
  if (DBG) {
    cat("\n[run_portfolio_race]\n")
    cat("  configs_race =", paste(configs_race, collapse = ","), "\n")
  }
  sd0_next = gen$sd0_next
  
  if (length(portfolios) == 0L) stop("Not enough configs to form portfolios")
  stopifnot(length(configs_race) == n_configs)
  
  # ---------- build portfolios ----------
  if (!is.null(fixed_portfolios) && length(fixed_portfolios) > 0) {
    fixed_norm = lapply(fixed_portfolios, function(P) sort(cfg_trunc4(P)))
    fixed_norm = Filter(function(P) all(P %in% configs_race) && length(P) == k_portfolio, fixed_norm)

    allP = c(fixed_norm, portfolios)
    keys_allP = vapply(allP, function(P) make_port_key(P, sort_members = TRUE), character(1))
    portfolios = allP[!duplicated(keys_allP)]
  }
  

  
  port_name = function(P) paste0("(", paste(P, collapse = ","), ")")
  P_names_all = vapply(portfolios, port_name, character(1))
  
  # portfolio keys (用于 port_seen_env & elite gating)
  P_keys_all = vapply(portfolios, function(P) make_port_key(P, sort_members = TRUE), character(1))
  
  # elite keys：本 race 开始时的 fixed_portfolios
  elite_keys = if (!is.null(fixed_portfolios) && length(fixed_portfolios) > 0) {
    vapply(fixed_portfolios, function(P) make_port_key(P, sort_members = TRUE), character(1))
  } else character(0)
  
  # fallback ranking from previous race elites (used when no effective test is run)
  prev_elite_rank_map = numeric(0)
  if (!is.null(fixed_portfolios) && length(fixed_portfolios) > 0 &&
      !is.null(elite_port_rankpos) && length(elite_port_rankpos) == length(fixed_portfolios)) {
    prev_keys = vapply(fixed_portfolios, function(P) make_port_key(P, sort_members = TRUE), character(1))
    prev_ranks = as.numeric(elite_port_rankpos)
    ok = is.finite(prev_ranks)
    if (any(ok)) {
      prev_split = split(prev_ranks[ok], prev_keys[ok])
      prev_elite_rank_map = vapply(prev_split, min, numeric(1))
    }
  }
  
  # elite 保护阈值：max_e + Tnew
  elite_threshold = compute_elite_thresholds(state, elite_keys, Tnew = Tnew)
  
  # ---------- instance schedule (race>=2) ----------
  instance_queue = make_instance_queue(
    state,
    race_id = as.integer(race_id),
    Tnew = as.integer(Tnew),
    seed = as.integer(race_seed),
    training_instances = training_instances,
    min_required = min_T_for_test
  )
  q_pos = 1L
  
  # 当前 race 开始后的 acquire 阶段：只有初始 queue 跑完后才会进入。
  acquire_occurrence = function() {
    if (is.null(training_instances) || length(training_instances) == 0L) {
      seen_old = get_seen_instances_from_port_env(state$port_seen_env)
      max_old = if (length(seen_old) == 0L) 0L else max(seen_old)
      return(list(i = as.integer(max_old + 1L), r = 0L, u = 0L, kind = "new_base"))
    }
    occ = make_instance_queue(
      state,
      race_id = max(2L, as.integer(race_id)),
      Tnew = 1L,
      seed = as.integer(race_seed + 2000L + used_T + 1L),
      training_instances = training_instances,
      min_required = 1L
    )
    list(
      i = as.integer(occ$instances[1L]),
      r = as.integer(occ$repeats[1L]),
      u = as.integer(occ$reuse_new_ids[1L]),
      kind = as.character(occ$kinds[1L])
    )
  }
  decode_occurrence_keys = function(x) {
    parse_ir_key(x)
  }
  
  # ---------- racing loop ----------
  active_idx = seq_along(portfolios)
  n_ports_all = length(portfolios)
  history_F = matrix(NA_real_, nrow = max_T, ncol = n_ports_all)
  
  used_T = 0L
  eval_used = 0L
  best_cfg_each_t = numeric(max_T)
  best_port_name_each_t = character(max_T)
  
  inst_used = integer(0)  # 本 race 里实际跑过的 base instance id
  inst_repeat_used = integer(0)
  inst_reuse_new_used = integer(0)
  
  no_drop_streak = 0L
  
  for (t in seq_len(max_T)) {
    if (length(active_idx) == 0L) {
      used_T = t - 1L
      if (verbose) cat(sprintf("Stop: no active portfolios remain before t=%d\n", t))
      break
    }
    
    # ---- choose instance id i_t ----
    if (race_id <= 1L) {
      if (!is.null(training_instances) && length(training_instances) > 0L) {
        if (t <= length(training_instances)) {
          occ_t = list(i = as.integer(training_instances[t]), r = 0L, u = 0L, kind = "new_base")
        } else {
          if (q_pos <= length(instance_queue$instances)) {
            occ_t = list(
              i = as.integer(instance_queue$instances[q_pos]),
              r = as.integer(instance_queue$repeats[q_pos]),
              u = as.integer(instance_queue$reuse_new_ids[q_pos]),
              kind = as.character(instance_queue$kinds[q_pos])
            )
            q_pos = q_pos + 1L
          } else {
            occ_t = acquire_occurrence()
          }
        }
      } else {
        occ_t = list(i = as.integer(t), r = 0L, u = 0L, kind = "new_base")
      }
    } else {
      if (q_pos <= length(instance_queue$instances)) {
        occ_t = list(
          i = as.integer(instance_queue$instances[q_pos]),
          r = as.integer(instance_queue$repeats[q_pos]),
          u = as.integer(instance_queue$reuse_new_ids[q_pos]),
          kind = as.character(instance_queue$kinds[q_pos])
        )
        q_pos = q_pos + 1L
      } else {
        occ_t = acquire_occurrence()
      }
    }
    
    i_t = as.integer(occ_t$i)
    repeat_t = as.integer(occ_t$r)
    reuse_new_t = as.integer(occ_t$u)
    inst_used = c(inst_used, i_t)
    inst_repeat_used = c(inst_repeat_used, repeat_t)
    inst_reuse_new_used = c(inst_reuse_new_used, reuse_new_t)

    base_seed_t = seed0 + 10000L * as.integer(i_t) + 100L * as.integer(repeat_t) + 1000000L * as.integer(reuse_new_t)
    used_T = t
    
    active_ports = portfolios[active_idx]
    active_keys  = P_keys_all[active_idx]
    
    active_configs = sort(unique(unlist(active_ports)))
    
    # ---- lookup / evaluate configs (cached by instance+config) ----
    got = cache_get_scores_state(state, i = i_t, configs = active_configs, repeat_id = repeat_t, reuse_new_id = reuse_new_t)
    cfg_scores = got$scores
    hit = got$hit
    keys = got$keys
    
    todo_configs = active_configs[!hit]
    cache_hit_count = sum(hit)
    cache_miss_count = length(todo_configs)
    need_eval = cache_miss_count
    if (!is.null(hydra_fixed_slots) && length(hydra_fixed_slots) > 0L) {
      need_eval = need_eval + md_hydra_count_missing_fixed_slot_scores(
        state = state,
        instances = i_t,
        hydra_fixed_slots = hydra_fixed_slots,
        repeats = repeat_t,
        reuse_new_ids = reuse_new_t
      )
    }
    
    if (stop_when_budget_cannot_cover_next_instance && (eval_used + need_eval > max_eval)) {
      used_T = t - 1L
      if (verbose) cat(sprintf(
        "Stop: budget not enough for next instance (eval_used=%d, need_new=%d, max_eval=%d)\n",
        eval_used, need_eval, max_eval
      ))
      break
    }
    
    if (need_eval > 0) {
      new_scores = hydra_call_scores_configs_local(
        configs = todo_configs,
        i = i_t,
        base_seed = base_seed_t,
        w = w,
        sigma1 = sigma1,
        sigma2 = sigma2,
        p = p,
        state = state
      )
      cfg_scores[!hit] = as.numeric(new_scores)
      cache_set_scores_state(state, keys = keys[!hit], scores = cfg_scores[!hit])
      eval_used = eval_used + need_eval
    }

    if (isTRUE(verbose) && isTRUE(getOption("DBG_eval_hit", FALSE))) {
      cat(sprintf(
        "  eval_hit: i=%d r=%d u=%d hit=%d miss=%d new_eval=%d total_eval=%d\n",
        i_t, repeat_t, reuse_new_t, cache_hit_count, cache_miss_count, need_eval, eval_used
      ))
    }
    
    names(cfg_scores) = as.character(cfg_tick(active_configs))

    hydra_fixed_score_t = NULL
    if (!is.null(hydra_fixed_slots) && length(hydra_fixed_slots) > 0L) {
      hydra_fixed_score_t = md_hydra_eval_fixed_portfolio_score(
        state = state,
        hydra_fixed_slots = hydra_fixed_slots,
        i_t = i_t,
        base_seed_t = base_seed_t,
        w = w,
        sigma1 = sigma1,
        sigma2 = sigma2,
        p = p,
        repeat_id = repeat_t,
        reuse_new_id = reuse_new_t
      )
      eval_used = eval_used + hydra_fixed_score_t$new_eval
    }
    
    
    # ---- compute portfolio scores ----
    has_intra_port_dup = any(vapply(active_ports, function(P) any(duplicated(cfg_tick(P))), logical(1)))
    if (isTRUE(enable_duplicate) && has_intra_port_dup) {
      F_t = portfolio_scores_with_duplicate(
        portfolios = active_ports,
        cfg_scores = cfg_scores,
        i_t = i_t,
        base_seed_t = base_seed_t,
        repeat_id = repeat_t,
        reuse_new_id = reuse_new_t,
        state = state,
        w = w,
        sigma1 = sigma1,
        sigma2 = sigma2,
        p = p
      )
    } else {
      F_t = scores_portfolios_on_instance(active_ports, cfg_scores)
    }
    if (!is.null(hydra_fixed_score_t)) {
      F_t = pmin(F_t, hydra_fixed_score_t$score)
    }
    history_F[t, active_idx] = F_t
    
    # ---- update all-seen instances for scheduling/cache ----
    update_port_seen(state, ports = active_ports, i = i_t, repeat_id = repeat_t, reuse_new_id = reuse_new_t, sort_members = TRUE)
    # ---- update e-count only for current race's carried elites ----
    if (length(elite_keys) > 0L) {
      elite_now = active_ports[active_keys %in% elite_keys]
      if (length(elite_now) > 0L) {
        update_elite_seen(state, ports = elite_now, instances = i_t, repeats = repeat_t, reuse_new_ids = reuse_new_t, sort_members = TRUE)
      }
    }
    
    # ---- track per-iteration best ----
    j_best_t = which.min(F_t)
    Pbest = active_ports[[j_best_t]]
    Pbest_tick = cfg_tick(Pbest)
    best_cfg_each_t[t] = Pbest[ which.min(cfg_scores[as.character(Pbest_tick)]) ]
    best_port_name_each_t[t] = P_names_all[active_idx][j_best_t]
    
    if (verbose) {
      if (reuse_new_t > 0L) {
        cat(sprintf("t=%d | i=%d | repeat=%d | reuse_new=%d | alive=%d | eval_used=%d\n",
                    t, i_t, repeat_t, reuse_new_t, length(active_idx), eval_used))
      } else if (repeat_t > 0L) {
        cat(sprintf("t=%d | i=%d | repeat=%d | alive=%d | eval_used=%d\n",
                    t, i_t, repeat_t, length(active_idx), eval_used))
      } else {
        cat(sprintf("t=%d | i=%d | alive=%d | eval_used=%d\n",
                    t, i_t, length(active_idx), eval_used))
      }
    }
    
    # ---- stopping / testing schedule ----
    # Do not stop by N_min before reaching the first testing point.
    if (t < min_T_for_test) next
    if (length(active_idx) <= N_min) break
    if (((t - min_T_for_test) %% T_each) != 0) next
    
    # ---- statistical test on current alive ----
    Fmat = history_F[1:t, active_idx, drop = FALSE]
    
    out_test = aux_friedman_irace(
      results = Fmat,
      conf.level = 1 - alpha
    )
    out_alive = apply_elite_protection_per_elite(
      out_test_alive = out_test$alive,
      port_keys_alive_order = active_keys,
      elite_thresholds = elite_threshold,
      state = state
    )
    keep_cols = which(out_alive)
    if (length(keep_cols) == 0L) {
      ranks_now = as.numeric(out_test$ranks)
      if (length(ranks_now) == length(active_idx) && any(is.finite(ranks_now))) {
        keep_cols = which.min(replace(ranks_now, !is.finite(ranks_now), Inf))
      } else if (length(F_t) > 0L && any(is.finite(F_t))) {
        keep_cols = which.min(replace(as.numeric(F_t), !is.finite(F_t), Inf))
      } else {
        keep_cols = 1L
      }
      out_alive[] = FALSE
      out_alive[keep_cols] = TRUE
      if (DBG2) {
        cat("\n[test_keep_one_fallback]\n")
        cat("  race_id =", race_id, "| t =", t, "\n")
        cat("  reason = no survivor after primary test\n")
        cat("  kept_key =", active_keys[keep_cols], "\n")
      }
    }
    
    # Free-rider filter for k>1:
    # among survivors with identical portfolio performance on all past instances,
    # compare only the non-common members (free riders) and drop portfolios whose
    # free riders are statistically worse.
    if (isTRUE(enable_second_level_test) && isTRUE(global_free_rider_elimination) &&
        k_portfolio > 1L && length(keep_cols) > 1L) {
      F_keep = Fmat[, keep_cols, drop = FALSE]
      if (!is.matrix(F_keep)) F_keep = matrix(F_keep, ncol = length(keep_cols))
      perf_key = apply(F_keep, 2, function(v) paste(format(signif(v, 12), scientific = FALSE, trim = TRUE), collapse = "|"))
      perf_groups = split(seq_along(keep_cols), perf_key)
      if (length(perf_groups) > 0L) {
        inst_now = inst_used[seq_len(min(length(inst_used), t))]
        rep_now = inst_repeat_used[seq_len(min(length(inst_repeat_used), t))]
        reuse_now = inst_reuse_new_used[seq_len(min(length(inst_reuse_new_used), t))]
        tol_core = 1e-12
        grp_info = vector("list", length(perf_groups))
        core_all = integer(0)
        free_union_global = integer(0)
        
        gi = 0L
        for (g in perf_groups) {
          gi = gi + 1L
          ports_g = active_ports[g]
          keys_g = active_keys[g]
          ticks_list = lapply(ports_g, cfg_tick)
          
          core_best = integer(0)
          for (kk in seq_along(inst_now)) {
            i0 = as.integer(inst_now[kk])
            r0 = as.integer(rep_now[kk])
            u0 = as.integer(reuse_now[kk])
            for (pp in seq_along(ticks_list)) {
              tt = as.integer(ticks_list[[pp]])
              keys_ic = vapply(tt, function(ct) make_icr_key(i0, cfg_from_tick(ct), r = r0, u = u0), character(1))
              if (!all(vapply(keys_ic, exists, logical(1), envir = state$score_env, inherits = FALSE))) next
              sc = vapply(keys_ic, get, numeric(1), envir = state$score_env, inherits = FALSE)
              smin = min(sc)
              winners = tt[abs(sc - smin) <= tol_core]
              core_best = union(core_best, winners)
            }
          }
          
          free_list = lapply(ticks_list, function(tt) unique(setdiff(tt, core_best)))
          free_union_group = sort(unique(unlist(free_list)))
          
          grp_info[[gi]] = list(
            idx = g,
            keys = keys_g,
            ticks_list = ticks_list,
            core_best = core_best,
            free_list = free_list,
            free_union_group = free_union_group
          )
          core_all = union(core_all, core_best)
          free_union_global = union(free_union_global, free_union_group)
        }
        
        free_union_global = sort(as.integer(free_union_global))
        free_noncore_global = sort(setdiff(free_union_global, core_all))
        
        drop_in_keep = integer(0)
        free_drop_global = integer(0)
        groups_affected = integer(0)
        alive_before = length(keep_cols)
        
        if (length(free_noncore_global) > 1L) {
          S = matrix(NA_real_, nrow = length(inst_now), ncol = length(free_noncore_global))
          for (ri in seq_along(inst_now)) {
            i0 = as.integer(inst_now[ri])
            r0 = as.integer(rep_now[ri])
            u0 = as.integer(reuse_now[ri])
            for (ci in seq_along(free_noncore_global)) {
              c0 = cfg_from_tick(free_noncore_global[ci])
              key_ic = make_icr_key(i0, c0, r = r0, u = u0)
              if (exists(key_ic, envir = state$score_env, inherits = FALSE)) {
                S[ri, ci] = get(key_ic, envir = state$score_env, inherits = FALSE)
              }
            }
          }
          
          row_ok = apply(S, 1, function(r) all(is.finite(r)))
          if (sum(row_ok) >= 2L) {
            S = S[row_ok, , drop = FALSE]
            out_free = aux_friedman_irace(
              results = S,
              conf.level = 1 - alpha
            )
            free_drop_global = free_noncore_global[!out_free$alive]
            
            if (length(free_drop_global) > 0L) {
              for (gi in seq_along(grp_info)) {
                ginfo = grp_info[[gi]]
                g = ginfo$idx
                # eliminate any portfolio containing rejected non-core freerider.
                g_drop = which(vapply(ginfo$ticks_list, function(tt) {
                  any(tt %in% free_drop_global)
                }, logical(1)))
                if (length(g_drop) > 0L) {
                  drop_in_keep = c(drop_in_keep, g[g_drop])
                  groups_affected = c(groups_affected, gi)
                }
              }
              
              # Protection: keep at least one representative per performance group.
              if (length(drop_in_keep) > 0L) {
                drop_in_keep = unique(drop_in_keep)
                keep_after = setdiff(seq_along(keep_cols), drop_in_keep)
                
                ranks_keep = out_test$ranks[keep_cols]
                for (gi in seq_along(grp_info)) {
                  g = grp_info[[gi]]$idx
                  if (!any(g %in% keep_after)) {
                    g_ranks = ranks_keep[g]
                    if (any(is.finite(g_ranks))) {
                      pick = g[which.min(g_ranks)]
                    } else {
                      pick = g[1L]
                    }
                    drop_in_keep = setdiff(drop_in_keep, pick)
                    keep_after = union(keep_after, pick)
                  }
                }
              }
            }
          }
        }
        
        if (length(drop_in_keep) > 0L) {
          out_alive[keep_cols[unique(drop_in_keep)]] = FALSE
          keep_cols = which(out_alive)
          if (length(keep_cols) == 0L) {
            ranks_now = as.numeric(out_test$ranks)
            if (length(ranks_now) == length(active_idx) && any(is.finite(ranks_now))) {
              keep_cols = which.min(replace(ranks_now, !is.finite(ranks_now), Inf))
            } else if (length(F_t) > 0L && any(is.finite(F_t))) {
              keep_cols = which.min(replace(as.numeric(F_t), !is.finite(F_t), Inf))
            } else {
              keep_cols = 1L
            }
            out_alive[] = FALSE
            out_alive[keep_cols] = TRUE
            if (DBG2) {
              cat("\n[test_keep_one_fallback]\n")
              cat("  race_id =", race_id, "| t =", t, "\n")
              cat("  reason = no survivor after second-level test\n")
              cat("  kept_key =", active_keys[keep_cols], "\n")
            }
          }
        }
        
        if (DBG2) {
          cat("\n[global_free_rider_filter]\n")
          cat("  race_id =", race_id, "| t =", t, "\n")
          cat("  groups_n =", length(grp_info), "\n")
          cat("  core_all_n =", length(core_all), "\n")
          cat("  free_global_n =", length(free_noncore_global), "\n")
          cat("  free_drop_global =", if (length(free_drop_global) == 0L) "<none>" else paste(free_drop_global, collapse = ","), "\n")
          cat("  groups_affected_n =", length(unique(groups_affected)), "\n")
          cat("  alive_before =", alive_before, "| alive_after =", length(keep_cols), "\n")
        }
      }
    }
    
    dropped_any_effective = length(keep_cols) < length(active_idx)
    
    if (!dropped_any_effective) {
      if (t >= T_minnodrop) {
        no_drop_streak = no_drop_streak + 1L
        if (no_drop_streak >= Tmaxnodrop) {
          if (verbose) cat(sprintf("Stop: %d consecutive tests without discard (t=%d)\n",
                                   Tmaxnodrop, t))
          break
        }
      }
      next
    }
    
    no_drop_streak = 0L
    
    # apply elimination
    active_idx = active_idx[keep_cols]
  }
  
  if (used_T >= min_T_for_test && used_T > 0L) {
    Fmat_final = history_F[1:used_T, active_idx, drop = FALSE]
    R_final = colSums2_base(
      rowRanks_base(Fmat_final, cols = seq_len(ncol(Fmat_final)))
    )
    
    ord_all = order(R_final)
    rankpos_all = rank(R_final, ties.method = "average")
    
    N_elite = min(length(active_idx), N_min)
    keep_pos = ord_all[seq_len(N_elite)]
    
    active_idx_elite = active_idx[keep_pos]
    R_elite = R_final[keep_pos]
    rankpos_elite = rankpos_all[keep_pos]
    
    active_idx = active_idx_elite
    R_final = R_elite
    
    best_idx = which.min(R_final)
    best_portfolio_name = P_names_all[active_idx][best_idx]
    
  } else {
    N_elite = min(length(active_idx), N_min)
    
    # If we break before enough tests, choose by previous elite ranking when available.
    if (N_elite > 0L && length(prev_elite_rank_map) > 0L) {
      keys_now = P_keys_all[active_idx]
      prev_rank_now = prev_elite_rank_map[keys_now]
      prev_rank_now[is.na(prev_rank_now)] = Inf
      ord_fallback = order(prev_rank_now, seq_along(active_idx))
      active_idx = active_idx[ord_fallback]
    }
    
    active_idx = active_idx[seq_len(N_elite)]
    best_idx = 1L
    best_portfolio_name = P_names_all[active_idx][best_idx]
    
    if (length(active_idx) > 0L && length(prev_elite_rank_map) > 0L) {
      ranks_now = prev_elite_rank_map[P_keys_all[active_idx]]
      if (all(is.na(ranks_now))) {
        R_final = rep(NA_real_, length(active_idx))
        rankpos_elite = rep(NA_real_, length(active_idx))
      } else {
        miss = is.na(ranks_now)
        if (any(miss)) {
          max_rank = max(ranks_now[!miss], na.rm = TRUE)
          ranks_now[miss] = max_rank + seq_len(sum(miss))
        }
        R_final = as.numeric(ranks_now)
        rankpos_elite = rank(R_final, ties.method = "average")
      }
    } else {
      R_final = rep(NA_real_, length(active_idx))
      rankpos_elite = rep(NA_real_, length(active_idx))
    }
    
    if (DBG2 && used_T < min_T_for_test) {
      cat("\n[early_terminate_fallback]\n")
      cat("  race_id =", race_id,
          "| used_T =", used_T,
          "| min_T_for_test =", min_T_for_test, "\n")
      
      if (length(prev_elite_rank_map) > 0L) {
        ord_prev = order(as.numeric(prev_elite_rank_map), names(prev_elite_rank_map))
        prev_keys = names(prev_elite_rank_map)[ord_prev]
        prev_vals = as.numeric(prev_elite_rank_map)[ord_prev]
        cat("  prev_elite_rank_map =",
            paste(sprintf("%s:%g", prev_keys, prev_vals), collapse = ", "),
            "\n")
      } else {
        cat("  prev_elite_rank_map = <empty>\n")
      }
      
      cat("  fallback_best =", best_portfolio_name, "\n")
    }
  }
  
  # Count current race instances only for this race's final elites.
  if (length(active_idx) > 0L && used_T > 0L) {
    final_elites = portfolios[active_idx]
    update_elite_seen(
      state,
      ports = final_elites,
      instances = inst_used[seq_len(min(length(inst_used), used_T))],
      repeats = inst_repeat_used[seq_len(min(length(inst_repeat_used), used_T))],
      reuse_new_ids = inst_reuse_new_used[seq_len(min(length(inst_reuse_new_used), used_T))],
      sort_members = TRUE
    )
  }
  
  # Final best must be one of elites, and among elites with maximal seen instances (e).
  if (length(active_idx) > 0L) {
    elite_keys_final = P_keys_all[active_idx]
    elite_e = vapply(elite_keys_final, function(k) get_port_e(state, k), integer(1))
    max_e = max(elite_e)
    cand_pos = which(elite_e == max_e)
    
    if (length(cand_pos) > 1L) {
      # Tie on e: compare candidates using cumulative elite-seen instances,
      # not only this race's R_final.
      tie_keys = elite_keys_final[cand_pos]
      tie_ports = portfolios[active_idx[cand_pos]]
      
      tie_inst_list = lapply(tie_keys, function(k) {
        if (exists(k, envir = state$elite_seen_env, inherits = FALSE)) {
          decode_occurrence_keys(get(k, envir = state$elite_seen_env, inherits = FALSE))
        } else integer(0)
      })
      tie_inst_union = do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x) > 0L, tie_inst_list))
      
      R_tie = rep(Inf, length(cand_pos))
      if (!is.null(tie_inst_union) && nrow(tie_inst_union) > 0L) {
        tie_inst_union = unique(tie_inst_union)
        F_rows = vector("list", nrow(tie_inst_union))
        keep_row = logical(nrow(tie_inst_union))
        
        for (ii in seq_len(nrow(tie_inst_union))) {
          i0 = as.integer(tie_inst_union$i[ii])
          r0 = as.integer(tie_inst_union$r[ii])
          u0 = as.integer(tie_inst_union$u[ii])
          fvec = rep(NA_real_, length(cand_pos))
          
          for (pp in seq_along(tie_ports)) {
            Pp = tie_ports[[pp]]
            cfgs = as.numeric(Pp)
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
        }
      }
      
      F_tie_sum = rep(Inf, length(cand_pos))
      F_tie_mean = rep(Inf, length(cand_pos))
      if (exists("Fmat_tie") && is.matrix(Fmat_tie) && nrow(Fmat_tie) > 0L) {
        F_tie_sum = colSums(Fmat_tie)
        F_tie_mean = colMeans(Fmat_tie)
      }
      tie_min_tick = vapply(tie_ports, function(P) min(cfg_tick(P)), integer(1))
      
      if (all(is.infinite(R_tie)) && length(R_final) == length(active_idx) && any(is.finite(R_final[cand_pos]))) {
        r_sub = R_final[cand_pos]
        r_sub[!is.finite(r_sub)] = Inf
        cand_pos = cand_pos[order(r_sub, F_tie_mean, tie_min_tick, cand_pos)]
      } else {
        cand_pos = cand_pos[order(R_tie, F_tie_mean, tie_min_tick, cand_pos)]
      }
      
      if (DBG2) {
        cat("  tie_on_e_keys =", paste(tie_keys, collapse = ","), "\n")
        cat("  tie_on_e_union_instances_n =", if (is.null(tie_inst_union)) 0L else nrow(tie_inst_union), "\n")
        cat("  tie_on_e_R =", paste(R_tie, collapse = ","), "\n")
        cat("  tie_on_e_mean =", paste(F_tie_mean, collapse = ","), "\n")
        cat("  tie_on_e_Fsum =", paste(F_tie_sum, collapse = ","), "\n")
      }
    }
    
    best_idx = cand_pos[1L]
    best_portfolio_name = P_names_all[active_idx][best_idx]
    
    if (DBG2) {
      cat("\n[best_from_max_e_elite]\n")
      cat("  race_id =", race_id, "\n")
      cat("  elite_keys =", paste(elite_keys_final, collapse = ","), "\n")
      cat("  elite_e =", paste(elite_e, collapse = ","), "\n")
      cat("  max_e =", max_e, "\n")
      cat("  best_name =", best_portfolio_name, "\n")
    }
  } else {
    best_portfolio_name = NA_character_
  }

  is_md = exists("is_multidim_param", mode = "function") && isTRUE(is_multidim_param(param_config))
  allow_mark_tick = integer(0)
  must_dup_tick_out = integer(0)
  best_group_no_freerider = FALSE
  best_group_dup_keep_tick = integer(0)
  core_best_out = integer(0)
  if (length(active_idx) > 0L && used_T > 0L) {
    Fmat_alive = history_F[1:used_T, active_idx, drop = FALSE]
    alive_ports = portfolios[active_idx]
    alive_names = P_names_all[active_idx]
    best_pos = which(alive_names == best_portfolio_name)
    if (length(best_pos) > 0L && is.matrix(Fmat_alive) && ncol(Fmat_alive) == length(active_idx)) {
      best_pos = best_pos[1L]
      perf_key = apply(Fmat_alive, 2, function(v) paste(format(signif(v, 12), scientific = FALSE, trim = TRUE), collapse = "|"))
      grp = which(perf_key == perf_key[best_pos])
      inst_now = inst_used[seq_len(min(length(inst_used), used_T))]
      rep_now = inst_repeat_used[seq_len(min(length(inst_repeat_used), used_T))]
      reuse_now = inst_reuse_new_used[seq_len(min(length(inst_reuse_new_used), used_T))]
      tol_core = 1e-12
      core_best = integer(0)
      for (kk in seq_along(inst_now)) {
        i0 = as.integer(inst_now[kk])
        r0 = as.integer(rep_now[kk])
        u0 = as.integer(reuse_now[kk])
        for (pp in grp) {
          Pt = if (is_md) as.integer(alive_ports[[pp]]) else cfg_tick(alive_ports[[pp]])
          keys_ic = vapply(Pt, function(ct) {
            c_key = if (is_md) as.numeric(ct) else cfg_from_tick(ct)
            make_icr_key(i0, c_key, r = r0, u = u0)
          }, character(1))
          if (!all(vapply(keys_ic, exists, logical(1), envir = state$score_env, inherits = FALSE))) next
          sc = vapply(keys_ic, get, numeric(1), envir = state$score_env, inherits = FALSE)
          smin = min(sc)
          core_best = union(core_best, Pt[abs(sc - smin) <= tol_core])
        }
      }
      core_best_out = sort(as.integer(core_best))
      grp_ticks = lapply(grp, function(pp) {
        if (is_md) as.integer(alive_ports[[pp]]) else cfg_tick(alive_ports[[pp]])
      })
      grp_union_tick = sort(unique(as.integer(unlist(grp_ticks))))
      free_union_best = sort(setdiff(grp_union_tick, core_best))
      best_group_no_freerider = (length(free_union_best) == 0L)

      best_port_ticks = if (is_md) as.integer(alive_ports[[best_pos]]) else cfg_tick(alive_ports[[best_pos]])
      dup_tab_best = table(as.integer(best_port_ticks))
      best_group_dup_keep_tick = sort(as.integer(names(dup_tab_best)[dup_tab_best >= 2L]))
      # Duplicate candidates should always come from the configs that
      # actually represent the best-performance group, even when the
      # group contains no free rider and even when only one portfolio
      # remains in the group.
      allow_mark_tick = sort(as.integer(core_best))
    }
    wins = if (is_md) as.integer(best_cfg_each_t[seq_len(used_T)]) else cfg_tick(best_cfg_each_t[seq_len(used_T)])
    wins = wins[is.finite(wins)]
    if (length(wins) > 0L && length(allow_mark_tick) > 0L) {
      tab = table(as.character(wins))
      max_n = max(tab)
      cand = as.integer(names(tab)[tab == max_n])
      cand = intersect(cand, allow_mark_tick)
      if (length(cand) == 0L) cand = allow_mark_tick
      if (length(cand) > 1L) {
        inst_now = inst_used[seq_len(min(length(inst_used), used_T))]
        rep_now = inst_repeat_used[seq_len(min(length(inst_repeat_used), used_T))]
        reuse_now = inst_reuse_new_used[seq_len(min(length(inst_reuse_new_used), used_T))]
        mscore = vapply(cand, function(ct) {
          keys = vapply(seq_along(inst_now), function(kk) {
            i0 = as.integer(inst_now[kk])
            r0 = as.integer(rep_now[kk])
            u0 = as.integer(reuse_now[kk])
            c_key = if (is_md) as.numeric(ct) else cfg_from_tick(ct)
            make_icr_key(i0, c_key, r = r0, u = u0)
          }, character(1))
          ok = vapply(keys, exists, logical(1), envir = state$score_env, inherits = FALSE)
          if (!any(ok)) return(Inf)
          mean(vapply(keys[ok], get, numeric(1), envir = state$score_env, inherits = FALSE))
        }, numeric(1))
        cand = cand[order(mscore, cand)]
      }
      must_dup_tick_out = cand[1L]
    }
  }
  if (!isTRUE(enable_duplicate)) {
    allow_mark_tick = integer(0)
    must_dup_tick_out = integer(0)
  }
  
  list(
    configs_race = configs_race,
    alive_portfolios = portfolios[active_idx],
    alive_names = P_names_all[active_idx],
    best_name = best_portfolio_name,
    alive_count = length(active_idx),
    
    alive_R = R_final,
    alive_rankpos = rank(R_final, ties.method="average"),
    alive_order = order(R_final),
    
    sd0_next = sd0_next,
    
    used_instances = used_T,
    inst_used = inst_used[seq_len(min(length(inst_used), used_T))],  # 真实 instance ids
    inst_repeat_used = inst_repeat_used[seq_len(min(length(inst_repeat_used), used_T))],
    inst_reuse_new_used = inst_reuse_new_used[seq_len(min(length(inst_reuse_new_used), used_T))],
    eval_used = eval_used,
    
    best_port_name_each_t = best_port_name_each_t[1:max(0L, used_T)],
    best_cfg_each_t = best_cfg_each_t[1:max(0L, used_T)],
    
 
    state = state,
    
    elite_keys = elite_keys,
    elite_threshold = elite_threshold,
    core_best = core_best_out,
    allow_mark_tick = allow_mark_tick,
    must_dup_tick = must_dup_tick_out,
    best_group_no_freerider = best_group_no_freerider,
    best_group_dup_keep_tick = best_group_dup_keep_tick
  )
}
