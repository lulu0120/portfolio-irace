make_param_discrete = function(name="1D", values=1:100, weights = NULL) { #in case we have more dim,先不管这个了因为先不考虑categorical
  values = as.integer(values)
  if (is.null(weights)) {
    weights = rep(1 / length(values), length(values))
  }
  stopifnot(length(weights) == length(values))
  weights = weights / sum(weights)
  
  list(
    name = name,
    type = "discrete",
    values = values,
    weights = weights
  )
}


cost_config = function(c, i,
                       w = 1, p = 1,
                       sigma2 = 50,
                       instance_random = 0,
                       eps = NULL) {
  if (is.null(eps)) eps = rnorm(1, mean = 0, sd = sigma2)
  c_int = as.integer(round(c))
  i_int = as.integer(i)
  penalty = if ((c_int %% 2L) == (i_int %% 2L)) 0 else p
  instance_random + w * abs(c - 50) + penalty + eps
}

scores_configs_on_instance = function(configs, i, base_seed,
                                      w = 1, sigma1 = 1, sigma2 = 0.1,
                                      p = 10) {
  configs = as.numeric(configs)
  cfg_u = sort(unique(configs))
  set.seed(as.integer(base_seed))
  instance_random = rnorm(1, mean = 0, sd = sigma1)
  eps_u = rnorm(length(cfg_u), mean = 0, sd = sigma2)
  scores_u = vapply(seq_along(cfg_u), function(k) {
    cost_config(
      c = cfg_u[k],
      i = i,
      w = w,
      p = p,
      sigma2 = sigma2,
      instance_random = instance_random,
      eps = eps_u[k]
    )
  }, numeric(1))
  scores_u[match(configs, cfg_u)]
}
