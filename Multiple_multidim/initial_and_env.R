
CFG_DIGITS = 4L
CFG_SCALE  = 10^CFG_DIGITS   # 1e4

#######key 必须是int不然容易出floating的问题

# scale up and floor down to match the irace scale
cfg_tick = function(x, scale = CFG_SCALE) {
  as.integer(floor(as.numeric(x) * scale))
}

cfg_trunc4 = function(x, digits = CFG_DIGITS) {
  scale = 10^digits
  floor(as.numeric(x) * scale) / scale
}


cfg_from_tick = function(tick, scale = CFG_SCALE) {
  as.numeric(tick) / as.numeric(scale)
}


#use tick over floating number
make_ir_key = function(i, r = 0L, u = 0L) {
  paste0("i=", as.integer(i), "|r=", as.integer(r), "|u=", as.integer(u))
}

parse_ir_key = function(key) {
  toks = regmatches(as.character(key), gregexpr("[0-9]+", as.character(key)))
  ii = vapply(toks, function(x) if (length(x) >= 1L) as.integer(x[1L]) else NA_integer_, integer(1))
  rr = vapply(toks, function(x) if (length(x) >= 2L) as.integer(x[2L]) else 0L, integer(1))
  uu = vapply(toks, function(x) if (length(x) >= 3L) as.integer(x[3L]) else 0L, integer(1))
  data.frame(i = ii, r = rr, u = uu, stringsAsFactors = FALSE)
}

make_icr_key = function(i, c, r = 0L, u = 0L) {
  paste0("i=", as.integer(i), "|r=", as.integer(r), "|u=", as.integer(u), "|c=", cfg_tick(c))
}

make_ic_key = function(i, c) {
  make_icr_key(i = i, c = c, r = 0L, u = 0L)
}

# key for duplicate-aware cache: same instance/config but different duplicate index.
make_icd_key = function(i, c, d, r = 0L, u = 0L) {
  paste0("i=", as.integer(i), "|r=", as.integer(r), "|u=", as.integer(u), "|c=", cfg_tick(c), "|d=", as.integer(d))
}



cache_get_scores_state = function(state, i, configs, repeat_id = 0L, reuse_new_id = 0L) {
  keys = vapply(configs, function(c) make_icr_key(i = i, c = c, r = repeat_id, u = reuse_new_id), character(1))
  hit  = vapply(keys, exists, logical(1), envir = state$score_env, inherits = FALSE)
  scores = rep(NA_real_, length(configs))
  if (any(hit)) {
    scores[hit] = vapply(keys[hit], get, numeric(1), envir = state$score_env, inherits = FALSE)
  }
  list(scores = scores, hit = hit, keys = keys)
}

cache_set_scores_state = function(state, keys, scores) {
  for (k in seq_along(keys)) assign(keys[k], as.numeric(scores[k]), envir = state$score_env)
  invisible(state)
}


make_port_key = function(P, sort_members = TRUE, sep = ",") {
  P = cfg_tick(P)                 # tick for identity
  if (sort_members) P = sort(P)
  paste(P, collapse = sep)
}

init_race_state = function() {
  list(
    score_env    = new.env(parent = emptyenv()),   
    port_seen_env = new.env(parent = emptyenv()), # 	port_seen_env[[key]] = integer vector;save the seen for key
    elite_seen_env = new.env(parent = emptyenv()), # elite-only seen instances for e counting
    ports_all_env = new.env(parent = emptyenv()),
    sd_env        = new.env(parent = emptyenv()),
    cat_prob_env  = new.env(parent = emptyenv()),
    cfg_id_by_key_env = new.env(parent = emptyenv()),   # multidim: key(vec ticks) -> id
    cfg_vec_by_id_env = new.env(parent = emptyenv()),   # multidim: id -> numeric vector
    next_cfg_id = 1L,
    inst_repeat_env = new.env(parent = emptyenv()),
    reuse_new_counter_env = new.env(parent = emptyenv()),
    allow_dup_count = integer(0),      # named by cfg tick; stored value is max allowed copies
    allow_dup_tick_active = integer(0),# only this set can duplicate in next race
    must_dup_tick = integer(0)         # one tick to force one extra copy in next race
  )
}

next_reuse_new_id = function(state) {
  if (is.null(state$reuse_new_counter_env)) state$reuse_new_counter_env = new.env(parent = emptyenv())
  key = "counter"
  cur = if (exists(key, envir = state$reuse_new_counter_env, inherits = FALSE)) {
    as.integer(get(key, envir = state$reuse_new_counter_env, inherits = FALSE))
  } else 0L
  cur = as.integer(cur + 1L)
  assign(key, cur, envir = state$reuse_new_counter_env)
  cur
}

md_cfg_key = function(x_vec, scale = CFG_SCALE) {
  paste(cfg_tick(as.numeric(x_vec), scale = scale), collapse = "|")
}

md_state_get_vec = function(state, id) {
  key = as.character(as.integer(id))
  if (!exists(key, envir = state$cfg_vec_by_id_env, inherits = FALSE)) return(NULL)
  as.numeric(get(key, envir = state$cfg_vec_by_id_env, inherits = FALSE))
}

md_state_get_or_create_id = function(state, x_vec) {
  if (is.null(state$cfg_id_by_key_env)) state$cfg_id_by_key_env = new.env(parent = emptyenv())
  if (is.null(state$cfg_vec_by_id_env)) state$cfg_vec_by_id_env = new.env(parent = emptyenv())
  x_vec = as.numeric(x_vec)
  key = md_cfg_key(x_vec)
  if (exists(key, envir = state$cfg_id_by_key_env, inherits = FALSE)) {
    return(as.integer(get(key, envir = state$cfg_id_by_key_env, inherits = FALSE)))
  }
  # Use env-backed cardinality as the ID source so ID allocation is stable
  # even though `state` is passed by value as a list.
  id = as.integer(length(ls(envir = state$cfg_vec_by_id_env, all.names = TRUE)) + 1L)
  assign(key, id, envir = state$cfg_id_by_key_env)
  assign(as.character(id), x_vec, envir = state$cfg_vec_by_id_env)
  id
}

update_port_seen = function(state, ports, i, repeat_id = 0L, reuse_new_id = 0L, sort_members = TRUE) {
  i = as.integer(i) 
  repeat_id = as.integer(repeat_id)
  reuse_new_id = as.integer(reuse_new_id)
  i_key = make_ir_key(i = i, r = repeat_id, u = reuse_new_id)
  
  for (P in ports) {
    key = make_port_key(P, sort_members = sort_members)
    
    if (!exists(key, envir = state$ports_all_env, inherits = FALSE)) {
      assign(key, sort(cfg_tick(P)), envir = state$ports_all_env)
    }
    
    old = if (exists(key, envir = state$port_seen_env, inherits = FALSE)) {
      get(key, envir = state$port_seen_env, inherits = FALSE)
    } else integer(0)
    
    assign(key, unique(c(old, i_key)), envir = state$port_seen_env)
  }
  base_key = as.character(i)
  prev_rep = if (exists(base_key, envir = state$inst_repeat_env, inherits = FALSE)) {
    as.integer(get(base_key, envir = state$inst_repeat_env, inherits = FALSE))
  } else -1L
  if (repeat_id > prev_rep) {
    assign(base_key, repeat_id, envir = state$inst_repeat_env)
  }
  invisible(state)
}

get_port_e = function(state, P_or_key, sort_members = TRUE) {  
  key = if (is.character(P_or_key)) P_or_key else make_port_key(P_or_key, sort_members)
  if (is.null(state$elite_seen_env)) return(0L)
  if (!exists(key, envir = state$elite_seen_env, inherits = FALSE)) return(0L)
  length(get(key, envir = state$elite_seen_env, inherits = FALSE))
}

update_elite_seen = function(state, ports, instances, repeats = NULL, reuse_new_ids = NULL, sort_members = TRUE) {
  if (is.null(state$elite_seen_env)) state$elite_seen_env = new.env(parent = emptyenv())
  inst = as.integer(instances)
  if (is.null(repeats)) repeats = rep.int(0L, length(inst))
  repeats = as.integer(repeats)
  if (is.null(reuse_new_ids)) reuse_new_ids = rep.int(0L, length(inst))
  reuse_new_ids = as.integer(reuse_new_ids)
  keep = is.finite(inst) & is.finite(repeats) & is.finite(reuse_new_ids)
  inst = inst[keep]
  repeats = repeats[keep]
  reuse_new_ids = reuse_new_ids[keep]
  occ = unique(vapply(seq_along(inst), function(k) make_ir_key(inst[k], repeats[k], reuse_new_ids[k]), character(1)))
  if (length(inst) == 0L || length(ports) == 0L) return(invisible(state))
  
  for (P in ports) {
    key = make_port_key(P, sort_members = sort_members)
    old = if (exists(key, envir = state$elite_seen_env, inherits = FALSE)) {
      get(key, envir = state$elite_seen_env, inherits = FALSE)
    } else integer(0)
    assign(key, unique(c(old, occ)), envir = state$elite_seen_env)
  }
  invisible(state)
}

get_seen_instances_from_port_env = function(port_seen_env) {  
  keys = ls(envir = port_seen_env, all.names = TRUE)
  if (length(keys) == 0L) return(integer(0))
  allv = unlist(lapply(keys, function(k) get(k, envir = port_seen_env, inherits = FALSE)))
  parsed = parse_ir_key(allv)
  base_i = parsed$i
  sort(unique(base_i[is.finite(base_i)]))
}

get_seen_occurrences_from_port_env = function(port_seen_env) {
  keys = ls(envir = port_seen_env, all.names = TRUE)
  if (length(keys) == 0L) {
    return(data.frame(i = integer(0), r = integer(0), u = integer(0), key = character(0), stringsAsFactors = FALSE))
  }
  allv = unique(unlist(lapply(keys, function(k) get(k, envir = port_seen_env, inherits = FALSE))))
  parsed = parse_ir_key(allv)
  keep = is.finite(parsed$i) & is.finite(parsed$r) & is.finite(parsed$u)
  if (!any(keep)) {
    return(data.frame(i = integer(0), r = integer(0), u = integer(0), key = character(0), stringsAsFactors = FALSE))
  }
  parsed = parsed[keep, , drop = FALSE]
  parsed$key = make_ir_key(parsed$i, parsed$r, parsed$u)
  unique(parsed)
}

get_instance_repeat_count = function(state, i) {
  key = as.character(as.integer(i))
  if (is.null(state$inst_repeat_env) ||
      !exists(key, envir = state$inst_repeat_env, inherits = FALSE)) {
    return(-1L)
  }
  as.integer(get(key, envir = state$inst_repeat_env, inherits = FALSE))
}

compute_elite_thresholds = function(state, elite_keys, Tnew = 1L) {
  elite_keys = as.character(elite_keys)
  if (length(elite_keys) == 0L) return(integer(0))
  
  e_vec = vapply(elite_keys, function(k) get_port_e(state, k), integer(1))
  thr  = e_vec + as.integer(Tnew)
  names(thr) = elite_keys
  thr
}

make_instance_queue = function(state,
                               race_id,
                               Tnew = 1L,
                               seed = 1L,
                               training_instances = NULL,
                               min_required = NULL) {
  Tnew = as.integer(Tnew)
  stopifnot(Tnew >= 0L)
  if (!is.null(min_required)) min_required = as.integer(min_required)
  
  if (race_id <= 1L) {
    return(list(instances = integer(0), repeats = integer(0), reuse_new_ids = integer(0), kinds = character(0)))
  }

  build_occurrences = function(base_ids, rep_ids, reuse_ids, kind) {
    n = length(base_ids)
    list(
      instances = as.integer(base_ids),
      repeats = as.integer(rep_ids),
      reuse_new_ids = as.integer(reuse_ids),
      kinds = rep.int(as.character(kind), n)
    )
  }
  combine_occurrences = function(parts) {
    if (length(parts) == 0L) {
      return(list(instances = integer(0), repeats = integer(0), reuse_new_ids = integer(0), kinds = character(0)))
    }
    list(
      instances = as.integer(unlist(lapply(parts, `[[`, "instances"), use.names = FALSE)),
      repeats = as.integer(unlist(lapply(parts, `[[`, "repeats"), use.names = FALSE)),
      reuse_new_ids = as.integer(unlist(lapply(parts, `[[`, "reuse_new_ids"), use.names = FALSE)),
      kinds = as.character(unlist(lapply(parts, `[[`, "kinds"), use.names = FALSE))
    )
  }
  alloc_reuse_new = function(base_pool, n_need, seed_shift = 0L) {
    if (n_need <= 0L || length(base_pool) == 0L) {
      return(build_occurrences(integer(0), integer(0), integer(0), "reuse_new"))
    }
    counter = setNames(vapply(base_pool, function(i) get_instance_repeat_count(state, i), integer(1)), as.character(base_pool))
    set.seed(as.integer(seed + seed_shift))
    chosen = sample(base_pool, size = n_need, replace = TRUE)
    reps = integer(n_need)
    reuse_ids = integer(n_need)
    for (idx in seq_along(chosen)) {
      key = as.character(chosen[idx])
      counter[[key]] = as.integer(counter[[key]] + 1L)
      reps[idx] = as.integer(counter[[key]])
      reuse_ids[idx] = next_reuse_new_id(state)
    }
    build_occurrences(chosen, reps, reuse_ids, "reuse_new")
  }

  seen_old = get_seen_instances_from_port_env(state$port_seen_env)
  seen_occ = get_seen_occurrences_from_port_env(state$port_seen_env)
  if (is.null(training_instances)) {
    max_old = if (length(seen_old) == 0L) 0L else max(seen_old)
    new_inst = if (Tnew > 0L) seq.int(max_old + 1L, max_old + Tnew) else integer(0)
    set.seed(seed)
    old_shuffled = if (length(seen_occ$i) > 0L) sample(seq_len(nrow(seen_occ)), size = nrow(seen_occ), replace = FALSE) else integer(0)
    parts = list(
      build_occurrences(new_inst, integer(length(new_inst)), integer(length(new_inst)), "new_base")
    )
    if (length(old_shuffled) > 0L) {
      parts[[length(parts) + 1L]] = build_occurrences(seen_occ$i[old_shuffled], seen_occ$r[old_shuffled], seen_occ$u[old_shuffled], "history")
    }
    out = combine_occurrences(parts)
    if (!is.null(min_required) && length(out$instances) < min_required) {
      pad = min_required - length(out$instances)
      pad_occ = build_occurrences(seq.int(max_old + length(new_inst) + 1L, max_old + length(new_inst) + pad), integer(pad), integer(pad), "new_base")
      out = combine_occurrences(list(out, pad_occ))
    }
    return(out)
  }

  base_pool = unique(as.integer(training_instances))
  seen_old = intersect(seen_old, base_pool)
  unseen_old = setdiff(base_pool, seen_old)
  n_new_base = min(Tnew, length(unseen_old))
  new_inst = if (n_new_base > 0L) unseen_old[seq_len(n_new_base)] else integer(0)
  reuse_new = alloc_reuse_new(base_pool, n_need = max(0L, Tnew - n_new_base), seed_shift = 101L)
  seen_occ = seen_occ[seen_occ$i %in% base_pool, , drop = FALSE]
  set.seed(seed)
  old_shuffled = if (nrow(seen_occ) > 0L) sample(seq_len(nrow(seen_occ)), size = nrow(seen_occ), replace = FALSE) else integer(0)
  parts = list(
    build_occurrences(new_inst, integer(length(new_inst)), integer(length(new_inst)), "new_base"),
    reuse_new
  )
  if (length(old_shuffled) > 0L) {
    parts[[length(parts) + 1L]] = build_occurrences(seen_occ$i[old_shuffled], seen_occ$r[old_shuffled], seen_occ$u[old_shuffled], "history")
  }
  out = combine_occurrences(parts)
  if (!is.null(min_required) && length(out$instances) < min_required) {
    out = combine_occurrences(list(out, alloc_reuse_new(base_pool, n_need = min_required - length(out$instances), seed_shift = 1001L)))
  }
  out
}

apply_elite_protection_per_elite = function(out_test_alive,
                                            port_keys_alive_order,
                                            elite_thresholds,
                                            state) {
  if (length(elite_thresholds) == 0L) return(out_test_alive)
  
  for (j in seq_along(out_test_alive)) {
    key = port_keys_alive_order[j]
    if (key %in% names(elite_thresholds)) {
      e_cur = get_port_e(state, key)  
      if (e_cur < elite_thresholds[[key]]) {
        out_test_alive[j] = TRUE     
      }
    }
  }
  out_test_alive
}


get_cfg_sd = function(state, tick, default_sd) {
  k = as.character(as.integer(tick))
  if (exists(k, envir = state$sd_env, inherits = FALSE)) {
    get(k, envir = state$sd_env, inherits = FALSE)
  } else default_sd
}

set_cfg_sd = function(state, tick, sd) {
  k = as.character(as.integer(tick))
  assign(k, as.numeric(sd), envir = state$sd_env)
  invisible(state)
}



#####################DEBUGGING##########
DBG = isTRUE(getOption("DBG", FALSE))

dbg_vec = function(x, name, max_show = 10) {
  if (!DBG) return(invisible(NULL))
  x0 = x
  x = as.numeric(x)
  cat("\n----", name, "----\n")
  cat("len:", length(x0), "\n")
  cat("head:", paste(head(x0, max_show), collapse = ", "), "\n")
  cat("tail:", paste(tail(x0, max_show), collapse = ", "), "\n")
  cat("NA:", sum(is.na(x)), " NaN:", sum(is.nan(x)),
      " Inf:", sum(is.infinite(x)), "\n")
  cat("min/max:", suppressWarnings(min(x, na.rm=TRUE)), "/",
      suppressWarnings(max(x, na.rm=TRUE)), "\n")
  cat("sum:", suppressWarnings(sum(x, na.rm=TRUE)), "\n")
  invisible(NULL)
}

dbg_probs = function(x, prob, tag) {
  if (!DBG) return(invisible(NULL))
  cat("\n==================== PROB CHECK:", tag, "====================\n")
  cat("len(x)   =", length(x), "\n")
  cat("len(prob)=", length(prob), "\n")
  dbg_vec(x,   paste0(tag, " :: x"))
  dbg_vec(prob, paste0(tag, " :: prob"))
  if (length(x) != length(prob)) {
    cat(">>> LENGTH MISMATCH!\n")
  }
  if (any(is.na(prob)) || any(!is.finite(prob))) {
    cat(">>> prob contains NA/NaN/Inf!\n")
  }
  if (length(prob) > 0 && sum(prob, na.rm=TRUE) <= 0) {
    cat(">>> prob sum <= 0!\n")
  }
  cat("=============================================================\n")
  invisible(NULL)
}
