
draw_one_child_tick = function(parent_tick, sd0, lower, upper, round) {
  parent = cfg_from_tick(parent_tick)
  sd_p   = get_cfg_sd(state, parent_tick, sd0)
  
  x = rtnorm_trunc(
    n = 1, mean = parent, sd = sd_p,
    lower = lower, upper = upper,
    do_round = round
  )
  cfg_tick(cfg_trunc4(x))
}


config_ranks_from_elite_portfolios = function(elite_portfolios, elite_port_rankpos) {

  if (is.null(elite_portfolios) || length(elite_portfolios) == 0L ||
      is.null(elite_port_rankpos) || length(elite_port_rankpos) == 0L) {
    return(list(cfgs = integer(0), ranks = integer(0)))
  }
  stopifnot(length(elite_portfolios) == length(elite_port_rankpos))
  
  elite_portfolios_tick = lapply(elite_portfolios, cfg_tick)
  cfgs = sort(unique(unlist(elite_portfolios_tick)))

  if (length(cfgs) == 0L) return(list(cfgs = integer(0), ranks = integer(0)))
  best_portpos = vapply(cfgs, function(c) {
    idx = vapply(elite_portfolios_tick, function(Pt) c %in% Pt, logical(1))
    min(elite_port_rankpos[idx])
  }, numeric(1))
  
  ranks = rank(best_portpos, ties.method = "first")
  
  if (DBG) {
    cat("\n[config_ranks_from_elite_portfolios]\n")
    cat("  cfgs =", paste(cfgs, collapse = ","), "\n")
    cat("  ranks =", paste(ranks, collapse = ","), "\n")
  }
  list(cfgs = cfg_from_tick(cfgs), ranks = as.integer(ranks))
}

parent_probs_from_ranks = function(ranks) {  
  r = as.integer(ranks)
  N = length(r)
  stopifnot(all(r >= 1L), all(r <= N))
  w = N - r + 1L
  w / sum(w)
}

build_dup_cap = function(all_tick,
                         allow_dup_tick = integer(0),
                         allow_dup_count = integer(0),
                         k_portfolio = 1L,
                         must_dup_tick = integer(0)) {
  all_tick = as.integer(all_tick)
  k_portfolio = as.integer(k_portfolio)
  cap = setNames(rep(1L, length(all_tick)), as.character(all_tick))
  if (length(allow_dup_tick) > 0L) {
    allow_dup_tick = as.integer(unique(allow_dup_tick))
    cnt_map = integer(0)
    if (length(allow_dup_count) > 0L) {
      cnt_map = as.integer(allow_dup_count)
      names(cnt_map) = as.character(names(allow_dup_count))
    }
    for (tt in allow_dup_tick) {
      key = as.character(tt)
      base = if (key %in% names(cnt_map)) as.integer(cnt_map[[key]]) else 1L
      cap[key] = min(k_portfolio, max(1L, base))
    }
  }
  if (length(must_dup_tick) > 0L) {
    key = as.character(as.integer(must_dup_tick[1L]))
    if (key %in% names(cap)) cap[key] = max(cap[key], 2L)
  }
  cap
}

append_with_cap = function(chosen_tick, cand_tick, cap_map, target_n) {
  if (length(chosen_tick) >= target_n || length(cand_tick) == 0L) return(chosen_tick)
  cnt = table(as.character(chosen_tick))
  for (tt in as.integer(cand_tick)) {
    if (length(chosen_tick) >= target_n) break
    key = as.character(tt)
    cap_t = if (key %in% names(cap_map)) as.integer(cap_map[[key]]) else 1L
    cur_t = if (key %in% names(cnt)) as.integer(cnt[[key]]) else 0L
    if (cur_t < cap_t) {
      chosen_tick = c(chosen_tick, tt)
      cnt[[key]] = cur_t + 1L
    }
  }
  chosen_tick
}

sd_prob_shrink = function(sd0, N_new, n_param = N_param) {
  if (N_new <= 0) return(sd0)
  sd0 * (1 / N_new)^(1 / n_param)
}

bringback_sd_irace = function(sd_now, N_new, lower = LB, upper = UB, n_param = N_param) {
  N_new = as.integer(N_new)
  if (N_new <= 0L) return(as.numeric(sd_now))
  
  term1 = as.numeric(sd_now) * (N_new)^(2 / n_param)
  term2 = ((upper - lower) / 2) * (1 / N_new)^(1 / n_param)
  min(term1, term2)
}

duplicate_group_idx_tick = function(children_tick) {
  dup_flag = duplicated(children_tick) | duplicated(children_tick, fromLast = TRUE)
  which(dup_flag)
}

sobol_u_1d = function(n, seed = 1L) {
  n = as.integer(n)
  if (n <= 0L) return(numeric(0))
  
  seed = suppressWarnings(as.integer(seed))
  if (length(seed) != 1L || is.na(seed)) seed = 1L
  
  if (requireNamespace("randtoolbox", quietly = TRUE)) {
    set.seed(seed)
    u = as.numeric(randtoolbox::sobol(
      n = n, dim = 1, scrambling = 1, init = TRUE, normal = FALSE
    ))
    return(pmin(pmax(u, 0), 1 - .Machine$double.eps))
  }
  
  # Fallback: base-2 van der Corput with deterministic offset by seed.
  vdc_base2 = function(i) {
    x = 0
    denom = 1
    while (i > 0L) {
      denom = denom * 2
      x = x + (i %% 2L) / denom
      i = i %/% 2L
    }
    x
  }
  start = abs(seed) %% 1048576L
  u = vapply(seq_len(n), function(k) vdc_base2(start + k), numeric(1))
  pmin(pmax(u, 0), 1 - .Machine$double.eps)
}

sobol_select_from_pool_tick = function(pool_tick, n, seed = 1L, include_bounds = TRUE) {
  pool_tick = sort(unique(as.integer(pool_tick)))
  n = as.integer(n)
  if (n <= 0L) return(integer(0))
  if (length(pool_tick) < n) stop("Pool too small for requested Sobol selection size.")
  
  picked = integer(0)
  if (isTRUE(include_bounds) && n >= 2L) {
    picked = c(pool_tick[1L], pool_tick[length(pool_tick)])
  }
  
  need = n - length(picked)
  if (need > 0L) {
    u = sobol_u_1d(max(need * 4L, need + 8L), seed = seed)
    idx = pmin(length(pool_tick), pmax(1L, floor(u * length(pool_tick)) + 1L))
    cand = pool_tick[idx]
    cand = setdiff(unique(cand), picked)
    if (length(cand) > 0L) picked = c(picked, cand[seq_len(min(need, length(cand)))])
  }
  
  if (length(picked) < n) {
    remaining = setdiff(pool_tick, picked)
    need = n - length(picked)
    if (need > 0L) {
      # Deterministic spread fill for any leftover.
      idx_fill = unique(floor(seq(1, length(remaining), length.out = need)))
      idx_fill = pmin(pmax(idx_fill, 1L), length(remaining))
      picked = c(picked, remaining[idx_fill])
    }
  }
  
  if (length(picked) < n) {
    remaining = setdiff(pool_tick, picked)
    picked = c(picked, remaining[seq_len(n - length(picked))])
  }
  
  sort(unique(as.integer(picked)))[seq_len(n)]
}

is_multidim_param = function(param_config) {
  type = tolower(as.character(param_config$type))
  identical(type, "multi") && !is.null(param_config$params) && length(param_config$params) > 0L
}

md_param_info = function(param_config) {
  stopifnot(is_multidim_param(param_config))
  ps = param_config$params
  D = length(ps)
  ptype = tolower(vapply(ps, function(p) as.character(p$type), character(1)))
  is_cat = ptype %in% c("categorical", "class", "cat")
  is_num = !is_cat
  
  lower = rep(NA_real_, D)
  upper = rep(NA_real_, D)
  if (any(is_num)) {
    lower[is_num] = vapply(ps[is_num], function(p) as.numeric(p$lower), numeric(1))
    upper[is_num] = vapply(ps[is_num], function(p) as.numeric(p$upper), numeric(1))
  }
  round_flag = is_num & (ptype %in% c("int", "integer", "discrete"))
  
  cat_levels = vector("list", D)
  cat_labels = vector("list", D)
  if (any(is_cat)) {
    for (d in which(is_cat)) {
      vals = ps[[d]]$values
      if (is.null(vals)) vals = ps[[d]]$levels
      if (is.null(vals) || length(vals) == 0L) {
        stop("Categorical param needs non-empty `values` or `levels` at dim ", d)
      }
      vals = as.character(vals)
      cat_labels[[d]] = vals
      cat_levels[[d]] = seq_along(vals) # internal code: 1..K
    }
  }
  
  list(
    D = D,
    lower = lower,
    upper = upper,
    ptype = ptype,
    round_flag = round_flag,
    is_cat = is_cat,
    is_num = is_num,
    num_idx = which(is_num),
    cat_idx = which(is_cat),
    cat_levels = cat_levels,
    cat_labels = cat_labels
  )
}

balanced_cat_codes = function(n, k, seed = 1L) {
  n = as.integer(n); k = as.integer(k)
  if (n <= 0L || k <= 0L) return(integer(0))
  set.seed(as.integer(seed))
  out = integer(0)
  while (length(out) < n) {
    out = c(out, sample.int(k, size = k, replace = FALSE))
  }
  out[seq_len(n)]
}

md_sample_initial_sobol = function(Nj, info, seed_j) {
  Nj = as.integer(Nj)
  D = info$D
  if (Nj <= 0L) return(matrix(numeric(0), nrow = 0L, ncol = D))
  
  seed_j = suppressWarnings(as.integer(seed_j))
  if (length(seed_j) != 1L || is.na(seed_j)) seed_j = 1L
  
  out = matrix(NA_real_, nrow = Nj, ncol = D)
  if (length(info$num_idx) > 0L) {
    dim_num = length(info$num_idx)
    if (requireNamespace("randtoolbox", quietly = TRUE)) {
      set.seed(seed_j)
      U = randtoolbox::sobol(
        n = Nj, dim = dim_num, scrambling = 1, init = TRUE, normal = FALSE
      )
      U = as.matrix(U)
    } else {
      set.seed(seed_j)
      U = matrix(stats::runif(Nj * dim_num), nrow = Nj, ncol = dim_num)
    }
    U = pmin(pmax(U, 0), 1 - .Machine$double.eps)
    for (k in seq_along(info$num_idx)) {
      d = info$num_idx[k]
      x = info$lower[d] + U[, k] * (info$upper[d] - info$lower[d])
      if (isTRUE(info$round_flag[d])) x = round(x)
      out[, d] = cfg_trunc4(x, digits = CFG_DIGITS)
    }
  }
  if (length(info$cat_idx) > 0L) {
    for (d in info$cat_idx) {
      K = length(info$cat_levels[[d]])
      out[, d] = as.numeric(balanced_cat_codes(Nj, K, seed = seed_j + 1000L * d))
    }
  }
  out
}

md_cfg_key_local = function(x_vec) {
  paste(cfg_tick(as.numeric(x_vec)), collapse = "|")
}

md_state_get_id = function(state, x_vec) {
  if (is.null(state) || is.null(state$cfg_id_by_key_env)) return(NA_integer_)
  key = md_cfg_key_local(x_vec)
  if (!exists(key, envir = state$cfg_id_by_key_env, inherits = FALSE)) return(NA_integer_)
  as.integer(get(key, envir = state$cfg_id_by_key_env, inherits = FALSE))
}

md_state_register = function(state, x_vec) {
  md_state_get_or_create_id(state, x_vec)
}

md_state_get_vec_many = function(state, ids) {
  ids = as.integer(ids)
  mats = lapply(ids, function(id) md_state_get_vec(state, id))
  if (length(mats) == 0L) return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  do.call(rbind, mats)
}

generate_configs_numeric_irace_multidim = function(
    param_config,
    Nj,
    required_cfgs = integer(0),
    elite_ranks,
    j,
    sd0 = NULL,
    n_param = N_param,
    max_tries = 300L,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L,
    state = NULL,
    iter_idx = 1L
) {
  if (is.null(state)) stop("multidim generation requires state.")
  seed_j = suppressWarnings(as.integer(j))
  if (length(seed_j) == 1L && !is.na(seed_j)) set.seed(seed_j)
  DBG_softrestart = isTRUE(getOption("DBG_softrestart", FALSE))
  if (is.null(state$cat_prob_env)) state$cat_prob_env = new.env(parent = emptyenv())
  iter_idx = as.integer(iter_idx)
  if (!is.finite(iter_idx) || iter_idx < 1L) iter_idx = 1L
  N_iter_total = suppressWarnings(as.integer(state$N_iter_total))
  if (length(N_iter_total) != 1L || !is.finite(N_iter_total) || N_iter_total <= 0L) {
    N_iter_total = max(1L, iter_idx)
  }
  beta_iter = (iter_idx - 1) / N_iter_total
  beta_iter = max(0, min(1, beta_iter))
  info = md_param_info(param_config)
  D = info$D
  if (is.null(sd0)) {
    sd0 = rep(0, D)
    if (length(info$num_idx) > 0L) {
      sd0[info$num_idx] = (info$upper[info$num_idx] - info$lower[info$num_idx]) / 2
    }
  }
  sd0 = as.numeric(rep_len(sd0, D))
  
  required_ids = as.integer(required_cfgs)
  required_ids = required_ids[is.finite(required_ids)]
  required_ids = unique(required_ids)
  elite_ids = required_ids
  
  if (length(elite_ids) == 0L) {
    init_mat = md_sample_initial_sobol(Nj = Nj, info = info, seed_j = seed_j)
    ids = vapply(seq_len(nrow(init_mat)), function(i) md_state_register(state, init_mat[i, ]), integer(1))
    return(as.numeric(ids))
  }
  
  if (length(elite_ids) != length(elite_ranks)) stop("elite_cfgs and elite_ranks length mismatch.")
  if (Nj < length(required_ids)) stop("Nj < length(required_cfgs).")
  elite_mat = md_state_get_vec_many(state, elite_ids)
  probs = parent_probs_from_ranks(elite_ranks)
  
  N_elite = length(required_ids)
  N_need_must = if (length(must_dup_tick) > 0L && Nj > N_elite) 1L else 0L
  N_new = Nj - N_elite - N_need_must
  if (N_new <= 0L) {
    base_ids = required_ids
    if (N_need_must > 0L) base_ids = c(base_ids, as.integer(must_dup_tick[1L]))
    return(as.numeric(base_ids))
  }
  
  sd_j = sd_prob_shrink(sd0, N_new, n_param = max(1L, D))
  sd_state = setNames(lapply(seq_along(elite_ids), function(i) as.numeric(sd_j)), as.character(elite_ids))
  soft_restart_used = FALSE
  
  all_ids_known = as.integer(ls(state$cfg_vec_by_id_env))
  cap_map = build_dup_cap(
    all_tick = all_ids_known,
    allow_dup_tick = allow_dup_tick,
    allow_dup_count = allow_dup_count,
    k_portfolio = k_portfolio,
    must_dup_tick = must_dup_tick
  )
  
  cat_prob_key = function(parent_id, d) {
    paste0("pid=", as.integer(parent_id), "|d=", as.integer(d))
  }
  get_cat_prob = function(parent_id, d, parent_code) {
    K = length(info$cat_levels[[d]])
    key = cat_prob_key(parent_id, d)
    if (exists(key, envir = state$cat_prob_env, inherits = FALSE)) {
      p = as.numeric(get(key, envir = state$cat_prob_env, inherits = FALSE))
    } else {
      p = rep(1 / K, K)
    }
    # irace categorical update before sampling (j>1)
    if (iter_idx > 1L) {
      p = p * (1 - beta_iter)
      idx = as.integer(round(parent_code))
      idx = pmax(1L, pmin(K, idx))
      p[idx] = p[idx] + beta_iter
    }
    if (!all(is.finite(p)) || sum(p) <= 0) p = rep(1 / K, K)
    p = p / sum(p)
    assign(key, p, envir = state$cat_prob_env)
    p
  }
  
  sample_child_from_parent = function(parent_id) {
    pvec = md_state_get_vec(state, parent_id)
    if (is.null(pvec)) stop("Missing parent vector in state for id=", parent_id)
    sd_vec = as.numeric(sd_state[[as.character(parent_id)]])
    if (length(sd_vec) != D) sd_vec = rep_len(sd0, D)
    x = numeric(D)
    for (d in seq_len(D)) {
      if (isTRUE(info$is_cat[d])) {
        p = get_cat_prob(parent_id, d, parent_code = pvec[d])
        x[d] = as.numeric(sample.int(length(p), size = 1L, replace = TRUE, prob = p))
      } else {
        x[d] = rtnorm_trunc(
          n = 1,
          mean = pvec[d],
          sd = sd_vec[d],
          lower = info$lower[d],
          upper = info$upper[d],
          do_round = isTRUE(info$round_flag[d])
        )
        if (info$round_flag[d]) x[d] = round(x[d])
        x[d] = cfg_trunc4(x[d], digits = CFG_DIGITS)
      }
    }
    x
  }
  
  draw_batch_ids = function(m, allow_soft_restart = TRUE) {
    sample_once = function() {
      if (length(elite_ids) == 1L) {
        parents = rep(as.integer(elite_ids[1L]), m)
      } else {
        parents = as.integer(sample(elite_ids, size = m, replace = TRUE, prob = probs))
      }
      child_mat = t(vapply(parents, function(pid) sample_child_from_parent(pid), numeric(D)))
      child_keys = vapply(seq_len(nrow(child_mat)), function(i) md_cfg_key_local(child_mat[i, ]), character(1))
      child_ids = vapply(seq_len(nrow(child_mat)), function(i) md_state_register(state, child_mat[i, ]), integer(1))
      list(child_ids = as.integer(child_ids), parents = parents, child_keys = child_keys)
    }
    
    draw = sample_once()
    child_ids = draw$child_ids
    parents = draw$parents
    child_keys = draw$child_keys
    
    if (isTRUE(soft_restart) && allow_soft_restart) {
      restart_round = 0L
      while (restart_round < as.integer(soft_restart_max_times)) {
        dup_idx = which(duplicated(child_keys) | duplicated(child_keys, fromLast = TRUE))
        if (length(dup_idx) == 0L) break
        parents_to_restart = unique(as.character(parents[dup_idx]))
        if (length(parents_to_restart) == 0L) break
        if (DBG_softrestart) {
          cat("\n[soft_restart_triggered]\n")
          cat("  mode = multidim | race_seed =", j, "\n")
          cat("  parents_id =", paste(parents_to_restart, collapse = ","), "\n")
        }
        for (pid in parents_to_restart) {
          sd_prev = as.numeric(sd_state[[pid]])
          sd_new = numeric(D)
          for (d in seq_len(D)) {
            if (isTRUE(info$is_cat[d])) {
              key = cat_prob_key(as.integer(pid), d)
              K = length(info$cat_levels[[d]])
              if (exists(key, envir = state$cat_prob_env, inherits = FALSE)) {
                p = as.numeric(get(key, envir = state$cat_prob_env, inherits = FALSE))
              } else {
                p = rep(1 / K, K)
              }
              if (!all(is.finite(p)) || sum(p) <= 0) p = rep(1 / K, K)
              p = p / sum(p)
              p = 0.9 * p + 0.1 * max(p)
              p = p / sum(p)
              assign(key, p, envir = state$cat_prob_env)
              sd_new[d] = sd_prev[d]
            } else {
              sd_new[d] = bringback_sd_irace(
                sd_now = sd_prev[d],
                N_new = N_new,
                lower = info$lower[d],
                upper = info$upper[d],
                n_param = max(1L, D)
              )
            }
          }
          sd_state[[pid]] = sd_new
        }
        draw = sample_once()
        child_ids = draw$child_ids
        parents = draw$parents
        child_keys = draw$child_keys
        restart_round = restart_round + 1L
      }
      if (restart_round > 0L) soft_restart_used <<- TRUE
    }
    child_ids
  }
  
  chosen_ids = as.integer(required_ids)
  if (N_need_must > 0L) {
    chosen_ids = append_with_cap(chosen_ids, as.integer(must_dup_tick[1L]), cap_map = cap_map, target_n = Nj)
  }
  
  tries = 0L
  while (length(chosen_ids) < Nj && tries < max_tries) {
    need = Nj - length(chosen_ids)
    cand_ids = draw_batch_ids(need, allow_soft_restart = !soft_restart_used)
    all_known_now = sort(unique(as.integer(ls(state$cfg_vec_by_id_env))))
    cap_map = build_dup_cap(
      all_tick = all_known_now,
      allow_dup_tick = allow_dup_tick,
      allow_dup_count = allow_dup_count,
      k_portfolio = k_portfolio,
      must_dup_tick = must_dup_tick
    )
    chosen_ids = append_with_cap(chosen_ids, as.integer(cand_ids), cap_map = cap_map, target_n = Nj)
    tries = tries + 1L
  }
  
  if (length(chosen_ids) < Nj) {
    all_ids = sort(unique(as.integer(ls(state$cfg_vec_by_id_env))))
    if (length(all_ids) > 0L) {
      chosen_ids = append_with_cap(
        chosen_ids,
        cand_tick = sample(rep(all_ids, each = 2L), size = 2L * length(all_ids), replace = FALSE),
        cap_map = cap_map,
        target_n = Nj
      )
    }
  }
  
  if (length(chosen_ids) < Nj) {
    extra_try = 0L
    while (length(chosen_ids) < Nj && extra_try < max_tries) {
      x_new = md_sample_initial_sobol(
        Nj = max(4L, Nj - length(chosen_ids)),
        info = info,
        seed_j = as.integer(seed_j + 1000L + extra_try)
      )
      id_new = vapply(seq_len(nrow(x_new)), function(i) md_state_register(state, x_new[i, ]), integer(1))
      all_ids = sort(unique(as.integer(ls(state$cfg_vec_by_id_env))))
      cap_map = build_dup_cap(
        all_tick = all_ids,
        allow_dup_tick = allow_dup_tick,
        allow_dup_count = allow_dup_count,
        k_portfolio = k_portfolio,
        must_dup_tick = must_dup_tick
      )
      chosen_ids = append_with_cap(chosen_ids, as.integer(id_new), cap_map = cap_map, target_n = Nj)
      extra_try = extra_try + 1L
    }
  }
  if (length(chosen_ids) < Nj) {
    # Final fail-safe: fill the remainder with fresh unique samples, ignoring
    # duplicate-cap constraints to avoid hard failure from cap interactions.
    remain = Nj - length(chosen_ids)
    hard_try = 0L
    while (remain > 0L && hard_try < max_tries) {
      x_new = md_sample_initial_sobol(
        Nj = max(8L, remain * 2L),
        info = info,
        seed_j = as.integer(seed_j + 5000L + hard_try)
      )
      id_new = vapply(seq_len(nrow(x_new)), function(i) md_state_register(state, x_new[i, ]), integer(1))
      id_new = setdiff(unique(as.integer(id_new)), as.integer(chosen_ids))
      if (length(id_new) > 0L) {
        take_n = min(remain, length(id_new))
        chosen_ids = c(chosen_ids, id_new[seq_len(take_n)])
        remain = Nj - length(chosen_ids)
      }
      hard_try = hard_try + 1L
    }
  }
  
  if (length(chosen_ids) < Nj) stop("Cannot fill Nj in multidim generator.")
  as.numeric(chosen_ids)
}

wrap_scalar_categorical_as_multidim = function(param_config) {
  vals = param_config$values
  if (is.null(vals)) vals = param_config$levels
  if (is.null(vals) || length(vals) == 0L) {
    stop("Scalar categorical param needs non-empty `values` or `levels`.")
  }
  list(
    type = "multi",
    params = list(
      list(
        name = if (!is.null(param_config$name)) as.character(param_config$name) else "x1",
        type = "categorical",
        values = as.character(vals)
      )
    )
  )
}

scalar_categorical_register_ids = function(state, values) {
  vals_chr = as.character(values)
  ids = vapply(seq_along(vals_chr), function(idx) {
    md_state_register(state, as.numeric(idx))
  }, integer(1))
  names(ids) = vals_chr
  ids
}

scalar_categorical_values_to_ids = function(values, value_to_id) {
  if (length(values) == 0L) return(integer(0))
  vals_chr = as.character(values)
  ids = unname(as.integer(value_to_id[vals_chr]))
  ids[is.finite(ids)]
}

scalar_categorical_ids_to_values = function(ids, state, values) {
  ids = as.integer(ids)
  vals_chr = as.character(values)
  out = vapply(ids, function(id) {
    vec = md_state_get_vec(state, id)
    if (is.null(vec) || length(vec) == 0L || !is.finite(vec[[1L]])) return(NA_character_)
    idx = as.integer(round(vec[[1L]]))
    idx = max(1L, min(length(vals_chr), idx))
    vals_chr[[idx]]
  }, character(1))
  if (is.numeric(values)) {
    as.numeric(out)
  } else {
    out
  }
}

generate_configs_numeric_irace_scalar_categorical_multidim = function(
    param_config,
    Nj,
    required_cfgs = integer(0),
    elite_ranks,
    j,
    sd0 = NULL,
    n_param = 1L,
    max_tries = 200L,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L,
    state = NULL,
    iter_idx = 1L
) {
  if (is.null(state)) stop("scalar categorical multidim adapter requires state.")
  raw_values = if (!is.null(param_config$values)) param_config$values else param_config$levels
  wrapped = wrap_scalar_categorical_as_multidim(param_config)
  value_to_id = scalar_categorical_register_ids(state, raw_values)

  required_ids = scalar_categorical_values_to_ids(required_cfgs, value_to_id)
  allow_dup_ids = scalar_categorical_values_to_ids(allow_dup_tick, value_to_id)
  must_dup_ids = scalar_categorical_values_to_ids(must_dup_tick, value_to_id)

  allow_dup_count_ids = integer(0)
  if (length(allow_dup_count) > 0L) {
    raw_names = names(allow_dup_count)
    if (!is.null(raw_names) && length(raw_names) > 0L) {
      mapped_names = as.character(unname(value_to_id[as.character(raw_names)]))
      keep = !is.na(mapped_names) & nzchar(mapped_names)
      allow_dup_count_ids = as.integer(allow_dup_count[keep])
      names(allow_dup_count_ids) = mapped_names[keep]
    }
  }

  chosen_ids = generate_configs_numeric_irace_multidim(
    param_config = wrapped,
    Nj = Nj,
    required_cfgs = required_ids,
    elite_ranks = elite_ranks,
    j = j,
    sd0 = if (is.null(sd0)) 0 else sd0,
    n_param = max(1L, as.integer(n_param)),
    max_tries = max_tries,
    soft_restart = soft_restart,
    soft_restart_max_times = soft_restart_max_times,
    allow_dup_tick = allow_dup_ids,
    allow_dup_count = allow_dup_count_ids,
    must_dup_tick = must_dup_ids,
    k_portfolio = k_portfolio,
    state = state,
    iter_idx = iter_idx
  )

  scalar_categorical_ids_to_values(chosen_ids, state = state, values = raw_values)
}

# ----------------------------- helpers -----------------------------
# Expect param_config like:
#  - discrete: list(type="discrete", values=1:100)
#  - cts:      list(type="cts", lower=1, upper=100)
#
# cfg_tick/cfg_from_tick/cfg_trunc4/CFG_DIGITS already defined.

# ------------------- implementation: discrete pool -----------------
generate_configs_numeric_irace_discrete = function(
    configs,                    # numeric pool (e.g., 1:20)
    Nj,
    required_cfgs = integer(0),  # alive portfolios' member configs (numeric or tick ok)
    elite_ranks,
    j,
    sd0 = NULL,
    lower = LB,
    upper = UB,
    round = FALSE,
    n_param = N_param,
    max_tries = 200L,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L
) {
  seed_j = suppressWarnings(as.integer(j))
  if (length(seed_j) == 1L && !is.na(seed_j)) set.seed(seed_j)
  DBG_softrestart = isTRUE(getOption("DBG_softrestart", FALSE))
  
  if (is.null(sd0)) sd0 = (upper - lower) / 2
  
  # discrete pool -> tick identity
  configs_tick  = sort(unique(cfg_tick(configs)))
  required_tick = sort(unique(cfg_tick(required_cfgs)))
  elite_tick    = required_tick
  
  if (length(elite_tick) == 0L) {
    if (length(configs_tick) < Nj) stop("Discrete pool too small for Nj unique configs.")
    initial_tick = sobol_select_from_pool_tick(
      pool_tick = configs_tick,
      n = Nj,
      seed = seed_j,
      include_bounds = TRUE
    )
    return(sort(cfg_from_tick(initial_tick)))
  }
  
  if (length(elite_tick) != length(elite_ranks)) stop("elite_cfgs and elite_ranks length mismatch.")
  if (Nj < length(required_tick)) stop("Nj < length(required_cfgs).")
  
  probs = parent_probs_from_ranks(elite_ranks)
  
  N_elite = length(required_tick)
  N_need_must = if (length(must_dup_tick) > 0L && Nj > N_elite) 1L else 0L
  N_new   = Nj - N_elite - N_need_must
  if (N_new <= 0L) {
    base_tick = as.integer(required_tick)
    if (N_need_must > 0L) base_tick = c(base_tick, as.integer(must_dup_tick[1L]))
    return(cfg_from_tick(base_tick))
  }
  
  sd_j = sd_prob_shrink(sd0, N_new, n_param)
  sd_state = setNames(rep(sd_j, length(elite_tick)), as.character(elite_tick))
  soft_restart_used = FALSE
  cap_map = build_dup_cap(
    all_tick = configs_tick,
    allow_dup_tick = allow_dup_tick,
    allow_dup_count = allow_dup_count,
    k_portfolio = k_portfolio,
    must_dup_tick = must_dup_tick
  )
  
  draw_batch_tick = function(m, allow_soft_restart = TRUE) {
    sample_once = function() {
      if (length(elite_tick) == 1L) {
        parents_tick = rep(elite_tick, m)
      } else {
        parents_tick = sample(elite_tick, size = m, replace = TRUE, prob = probs)
      }
      
      parents = cfg_from_tick(parents_tick)
      sd_vec = as.numeric(sd_state[as.character(parents_tick)])
      x = rtnorm_trunc(
        n = m,
        mean = parents,
        sd = sd_vec,
        lower = lower,
        upper = upper,
        do_round = round
      )
      x = cfg_trunc4(x, digits = CFG_DIGITS)
      list(children_tick = cfg_tick(x), parents_tick = as.integer(parents_tick))
    }
    
    draw = sample_once()
    children_tick = draw$children_tick
    parents_tick  = draw$parents_tick
    
    if (isTRUE(soft_restart) && allow_soft_restart) {
      restart_round = 0L
      while (restart_round < as.integer(soft_restart_max_times)) {
        dup_idx = duplicate_group_idx_tick(children_tick)
        if (length(dup_idx) == 0L) break
        
        parents_to_restart = unique(as.character(parents_tick[dup_idx]))
        if (length(parents_to_restart) == 0L) break
        
        if (DBG_softrestart) {
          p_tick = as.integer(parents_to_restart)
          p_cfg = cfg_from_tick(p_tick)
          cat("\n[soft_restart_triggered]\n")
          cat("  mode = discrete | race_seed =", j, "\n")
          cat("  parents_tick =", paste(p_tick, collapse = ","), "\n")
          cat("  parents_cfg =", paste(p_cfg, collapse = ","), "\n")
        }
        
        for (cfg in parents_to_restart) {
          sd_state[cfg] = bringback_sd_irace(
            sd_now = sd_state[cfg],
            N_new = N_new,
            lower = lower,
            upper = upper,
            n_param = n_param
          )
        }
        
        draw = sample_once()
        children_tick = draw$children_tick
        parents_tick  = draw$parents_tick
        restart_round = restart_round + 1L
      }
      
      if (restart_round > 0L) soft_restart_used <<- TRUE
    }
    
    children_tick
  }
  
  chosen_tick = as.integer(required_tick)
  if (N_need_must > 0L) {
    chosen_tick = append_with_cap(
      chosen_tick = chosen_tick,
      cand_tick = as.integer(must_dup_tick[1L]),
      cap_map = cap_map,
      target_n = Nj
    )
  }
  tries = 0L
  
  while (length(chosen_tick) < Nj && tries < max_tries) {
    need = Nj - length(chosen_tick)
    cand_tick = draw_batch_tick(need, allow_soft_restart = !soft_restart_used)
    chosen_tick = append_with_cap(
      chosen_tick = chosen_tick,
      cand_tick = as.integer(cand_tick),
      cap_map = cap_map,
      target_n = Nj
    )
    tries = tries + 1L
  }
  
  if (length(chosen_tick) < Nj) {
    need = Nj - length(chosen_tick)
    pool = as.integer(configs_tick)
    if (length(pool) > 0L) {
      reps = rep(pool, each = 2L)
      chosen_tick = append_with_cap(
        chosen_tick = chosen_tick,
        cand_tick = sample(reps, size = length(reps), replace = FALSE),
        cap_map = cap_map,
        target_n = Nj
      )
      need = Nj - length(chosen_tick)
    }
    if (need > 0L) stop("Cannot fill Nj from discrete pool under duplicate constraints.")
  }
  
  cfg_from_tick(chosen_tick)
}

# ------------------- implementation: continuous range ----------------
generate_configs_numeric_irace_cts = function(
    lower,
    upper,
    Nj,
    required_cfgs = integer(0),
    elite_ranks,
    j,
    sd0 = NULL,
    round = FALSE,
    n_param = N_param,
    max_tries = 200L,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L
) {
  seed_j = suppressWarnings(as.integer(j))
  if (length(seed_j) == 1L && !is.na(seed_j)) set.seed(seed_j)
  DBG_softrestart = isTRUE(getOption("DBG_softrestart", FALSE))
  
  if (is.null(sd0)) sd0 = (upper - lower) / 2
  
  # build tick range (4-digit grid)
  lo_tick = cfg_tick(lower)
  hi_tick = cfg_tick(upper)
  if (hi_tick < lo_tick) stop("upper < lower")
  configs_tick = seq.int(lo_tick, hi_tick)
  
  required_tick = sort(unique(cfg_tick(required_cfgs)))
  elite_tick    = required_tick
  
  if (length(elite_tick) == 0L) {
    if (length(configs_tick) < Nj) stop("CTS range too small for Nj unique configs.")
    initial_tick = sobol_select_from_pool_tick(
      pool_tick = configs_tick,
      n = Nj,
      seed = seed_j,
      include_bounds = TRUE
    )
    return(sort(cfg_from_tick(initial_tick)))
  }
  
  if (length(elite_tick) != length(elite_ranks)) stop("elite_cfgs and elite_ranks length mismatch.")
  if (Nj < length(required_tick)) stop("Nj < length(required_cfgs).")
  
  probs = parent_probs_from_ranks(elite_ranks)
  
  N_elite = length(required_tick)
  N_need_must = if (length(must_dup_tick) > 0L && Nj > N_elite) 1L else 0L
  N_new   = Nj - N_elite - N_need_must
  if (N_new <= 0L) {
    base_tick = as.integer(required_tick)
    if (N_need_must > 0L) base_tick = c(base_tick, as.integer(must_dup_tick[1L]))
    return(cfg_from_tick(base_tick))
  }
  
  sd_j = sd_prob_shrink(sd0, N_new, n_param)
  sd_state = setNames(rep(sd_j, length(elite_tick)), as.character(elite_tick))
  soft_restart_used = FALSE
  cap_map = build_dup_cap(
    all_tick = configs_tick,
    allow_dup_tick = allow_dup_tick,
    allow_dup_count = allow_dup_count,
    k_portfolio = k_portfolio,
    must_dup_tick = must_dup_tick
  )
  
  draw_batch_tick = function(m, allow_soft_restart = TRUE) {
    sample_once = function() {
      if (length(elite_tick) == 1L) {
        parents_tick = rep(elite_tick, m)
      } else {
        parents_tick = sample(elite_tick, size = m, replace = TRUE, prob = probs)
      }
      
      parents = cfg_from_tick(parents_tick)
      sd_vec = as.numeric(sd_state[as.character(parents_tick)])
      x = rtnorm_trunc(
        n = m,
        mean = parents,
        sd = sd_vec,
        lower = lower,
        upper = upper,
        do_round = round
      )
      x = cfg_trunc4(x, digits = CFG_DIGITS)
      list(children_tick = cfg_tick(x), parents_tick = as.integer(parents_tick))
    }
    
    draw = sample_once()
    children_tick = draw$children_tick
    parents_tick  = draw$parents_tick
    
    if (isTRUE(soft_restart) && allow_soft_restart) {
      restart_round = 0L
      while (restart_round < as.integer(soft_restart_max_times)) {
        dup_idx = duplicate_group_idx_tick(children_tick)
        if (length(dup_idx) == 0L) break
        
        parents_to_restart = unique(as.character(parents_tick[dup_idx]))
        if (length(parents_to_restart) == 0L) break
        
        if (DBG_softrestart) {
          p_tick = as.integer(parents_to_restart)
          p_cfg = cfg_from_tick(p_tick)
          cat("\n[soft_restart_triggered]\n")
          cat("  mode = cts | race_seed =", j, "\n")
          cat("  parents_tick =", paste(p_tick, collapse = ","), "\n")
          cat("  parents_cfg =", paste(p_cfg, collapse = ","), "\n")
        }
        
        for (cfg in parents_to_restart) {
          sd_state[cfg] = bringback_sd_irace(
            sd_now = sd_state[cfg],
            N_new = N_new,
            lower = lower,
            upper = upper,
            n_param = n_param
          )
        }
        
        draw = sample_once()
        children_tick = draw$children_tick
        parents_tick  = draw$parents_tick
        restart_round = restart_round + 1L
      }
      
      if (restart_round > 0L) soft_restart_used <<- TRUE
    }
    
    children_tick
  }
  
  chosen_tick = as.integer(required_tick)
  if (N_need_must > 0L) {
    chosen_tick = append_with_cap(
      chosen_tick = chosen_tick,
      cand_tick = as.integer(must_dup_tick[1L]),
      cap_map = cap_map,
      target_n = Nj
    )
  }
  tries = 0L
  
  while (length(chosen_tick) < Nj && tries < max_tries) {
    need = Nj - length(chosen_tick)
    cand_tick = draw_batch_tick(need, allow_soft_restart = !soft_restart_used)
    chosen_tick = append_with_cap(
      chosen_tick = chosen_tick,
      cand_tick = as.integer(cand_tick),
      cap_map = cap_map,
      target_n = Nj
    )
    tries = tries + 1L
  }
  
  if (length(chosen_tick) < Nj) {
    need = Nj - length(chosen_tick)
    reps = rep(as.integer(configs_tick), each = 2L)
    chosen_tick = append_with_cap(
      chosen_tick = chosen_tick,
      cand_tick = sample(reps, size = length(reps), replace = FALSE),
      cap_map = cap_map,
      target_n = Nj
    )
    if (length(chosen_tick) < Nj) {
      stop("Cannot fill Nj with current duplicate constraints in [lower, upper].")
    }
  }
  
  cfg_from_tick(chosen_tick)
}

# ----------------------------- wrapper -----------------------------
generate_configs_numeric_irace = function(
    param_config,                # list(type=..., values=... OR lower/upper=...)
    Nj,
    required_cfgs = integer(0),
    elite_ranks,
    j,
    sd0 = NULL,
    lower = LB,
    upper = UB,
    round = FALSE,
    n_param = N_param,
    max_tries = 200L,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L,
    state = NULL,
    iter_idx = 1L
) {
  type = tolower(as.character(param_config$type))
  
  if (type %in% c("multi", "multidim", "multivariate")) {
    return(generate_configs_numeric_irace_multidim(
      param_config = param_config,
      Nj = Nj,
      required_cfgs = required_cfgs,
      elite_ranks = elite_ranks,
      j = j,
      sd0 = sd0,
      n_param = n_param,
      max_tries = max_tries,
      soft_restart = soft_restart,
      soft_restart_max_times = soft_restart_max_times,
      allow_dup_tick = allow_dup_tick,
      allow_dup_count = allow_dup_count,
      must_dup_tick = must_dup_tick,
      k_portfolio = k_portfolio,
      state = state,
      iter_idx = iter_idx
    ))
  }
  
  if (type %in% c("discrete", "class", "categorical")) {
    return(generate_configs_numeric_irace_scalar_categorical_multidim(
      param_config = param_config,
      Nj = Nj,
      required_cfgs = required_cfgs,
      elite_ranks = elite_ranks,
      j = j,
      sd0 = sd0,
      n_param = n_param,
      max_tries = max_tries,
      soft_restart = soft_restart,
      soft_restart_max_times = soft_restart_max_times,
      allow_dup_tick = allow_dup_tick,
      allow_dup_count = allow_dup_count,
      must_dup_tick = must_dup_tick,
      k_portfolio = k_portfolio,
      state = state,
      iter_idx = iter_idx
    ))
  }
  
  if (type %in% c("cts", "numeric", "real", "continuous")) {
    lo = if (!is.null(param_config$lower)) param_config$lower else lower
    hi = if (!is.null(param_config$upper)) param_config$upper else upper
    return(generate_configs_numeric_irace_cts(
      lower = lo,
      upper = hi,
      Nj = Nj,
      required_cfgs = required_cfgs,
      elite_ranks = elite_ranks,
      j = j,
      sd0 = sd0,
      round = round,
      n_param = n_param,
      max_tries = max_tries,
      soft_restart = soft_restart,
      soft_restart_max_times = soft_restart_max_times,
      allow_dup_tick = allow_dup_tick,
      allow_dup_count = allow_dup_count,
      must_dup_tick = must_dup_tick,
      k_portfolio = k_portfolio
    ))
  }
  
  stop("Unknown param_config$type: ", param_config$type)
}



# Example of usage:
#   param_config = list(type = "discrete", values = 1:100)
# configs_race = generate_configs_numeric_irace(
#   param_config = param_config,
#   Nj = n_configs,
#   required_cfgs = elite_cfg_info$cfgs,
#   elite_ranks = elite_cfg_info$ranks,
#   j = race_seed,
#   sd0 = sd0
# )
# param_config = list(type = "cts", lower = 1, upper = 100)
# configs_race = generate_configs_numeric_irace(
#   param_config = param_config,
#   Nj = n_configs,
#   required_cfgs = elite_cfg_info$cfgs,
#   elite_ranks = elite_cfg_info$ranks,
#   j = race_seed,
#   sd0 = sd0
# )


scores_portfolios_on_instance = function(portfolios, cfg_scores) {
  vapply(portfolios, function(P) {
    Pt = cfg_tick(P)
    min(cfg_scores[as.character(Pt)])
    
  }, numeric(1))
}
update_sd0_irace = function(sd0, N_new, n_param = N_param) {
  if (N_new <= 0L) return(sd0)
  sd_prob_shrink(sd0, N_new = N_new, n_param = n_param)
}

generate_configs_numeric_irace_with_sd = function(
    param_config,
    Nj,
    required_cfgs = integer(0),
    elite_ranks,
    j,
    sd0,
    lower = LB,
    upper = UB,
    round = FALSE,
    n_param = N_param,
    soft_restart = TRUE,
    soft_restart_max_times = 1L,
    allow_dup_tick = integer(0),
    allow_dup_count = integer(0),
    must_dup_tick = integer(0),
    k_portfolio = 1L,
    state = NULL,
    iter_idx = 1L
) {
  DBG = isTRUE(getOption("DBG", FALSE))
  if (DBG) {
    cat("\n[generate_configs_numeric_irace_with_sd]\n")
    cat("  Nj =", Nj, "\n")
    cat("  required_cfgs_len =", length(required_cfgs), "\n")
    cat("  elite_ranks_len =", length(elite_ranks), "\n")
  }
  
  required_tick = sort(unique(cfg_tick(required_cfgs)))
  N_new = as.integer(Nj - length(required_tick))
  
  chosen = generate_configs_numeric_irace(
    param_config = param_config,
    Nj = Nj,
    required_cfgs = required_cfgs,
    elite_ranks = elite_ranks,
    j = j,
    sd0 = sd0,
    lower = lower,
    upper = upper,
    round = round,
    n_param = n_param,
    soft_restart = soft_restart,
    soft_restart_max_times = soft_restart_max_times,
    allow_dup_tick = allow_dup_tick,
    allow_dup_count = allow_dup_count,
    must_dup_tick = must_dup_tick,
    k_portfolio = k_portfolio,
    state = state,
    iter_idx = iter_idx
  )
  
  sd0_next = update_sd0_irace(sd0, N_new = N_new, n_param = n_param)
  list(configs = chosen, sd0_next = sd0_next)
}
