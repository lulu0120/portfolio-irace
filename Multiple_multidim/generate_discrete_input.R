

generate_configs_discrete = function(param, n, seed) { #this seed is for global, no need to dependent on instance etc
  set.seed(seed)
  vals = param$values
  if (n > length(vals)) {
    stop("cannot sample without replacement")
  }
  sample(vals, size = n, replace = FALSE, prob = param$weights)
}

weights_uniform = function(values) {
  rep(1 / length(values), length(values))
}
# ---- counter (Dirichlet-like) state: inherited across races ----
init_counter = function(values, eps = 1) {
  values_chr = as.character(as.integer(values))
  setNames(rep(eps, length(values_chr)), values_chr)
}

update_counter_by_portfolio_count = function(cnt, alive_portfolios) {
  if (length(alive_portfolios) == 0) return(cnt)
  for (P in alive_portfolios) {
    for (c in P) {
      key = as.character(c)
      if (key %in% names(cnt)) cnt[key] = cnt[key] + 1
    }
  }
  cnt
}



update_counter_by_best_count = function(cnt, best_config_per_instance) { #先不用这个了
  if (length(best_config_per_instance) == 0) return(cnt)
  for (c in best_config_per_instance) {
    key = as.character(c)
    if (!is.na(key) && key %in% names(cnt)) cnt[key] = cnt[key] + 1
  }
  cnt
}

counter_to_weights = function(cnt, values) {
  values_chr = as.character(as.integer(values))
  w = as.numeric(cnt[values_chr])
  if (any(!is.finite(w)) || sum(w) <= 0) {
    w = rep(1, length(values_chr))
  }
  w / sum(w)
}

# helper: apply domain + inherited-weight update
apply_domain_and_weight_update = function(param_cur,
                                          next_values,
                                          weight_type = c("uniform", "portfolio_count", "best_count", "keep"),
                                          cnt_state,
                                          alive_portfolios = NULL,
                                          best_config_per_instance = NULL) {
  weight_type = match.arg(weight_type)
  
  next_values = sort(unique(as.integer(next_values)))
  if (length(next_values) == 0L) stop("apply_domain_and_weight_update: empty next_values")
  
  if (weight_type == "uniform") {
    param_cur$values  = next_values
    param_cur$weights = rep(1 / length(next_values), length(next_values))
    return(list(param = param_cur, cnt_state = cnt_state))
  }
  
  if (weight_type == "keep") {
    old_values = as.integer(param_cur$values)
    idx = match(next_values, old_values)
    idx = idx[!is.na(idx)]
    if (length(idx) == 0L) stop("keep: empty intersection")
    param_cur$values  = old_values[idx]
    w = param_cur$weights[idx]
    param_cur$weights = w / sum(w)
    return(list(param = param_cur, cnt_state = cnt_state))
  }
  
  if (weight_type == "portfolio_count") {
    if (is.null(alive_portfolios)) stop("portfolio_count needs alive_portfolios")
    cnt_state = update_counter_by_portfolio_count(cnt_state, alive_portfolios)
  } else if (weight_type == "best_count") {
    if (is.null(best_config_per_instance)) stop("best_count needs best_config_per_instance")
    cnt_state = update_counter_by_best_count(cnt_state, best_config_per_instance)
  }
  
  param_cur$values  = next_values
  param_cur$weights = counter_to_weights(cnt_state, next_values)
  
  list(param = param_cur, cnt_state = cnt_state)
}



# generator: discrete weighted sampling (no replacement)
generator = function(param, n, seed) {
  generate_configs_discrete(param, n = n, seed = seed)
}