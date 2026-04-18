rowRanks_base = function(y, cols) {
  r = t(apply(y[, cols, drop = FALSE], 1, rank, ties.method = "average"))
  storage.mode(r) = "double"
  r
}

colSums2_base = function(x) {
  as.numeric(colSums(x))
}

aux2_friedman_irace = function(y, I, alive, conf.level = 0.95) { #if we have more than 2 candidate then conover
  
  #y: result matrix, instance*candidate 
  #I: column index, 1:k
  #alive: a placeholder now, may change 
  
  dropped.any = FALSE
  n = nrow(y)
  k = length(I)
  n_num = as.numeric(n)
  k_num = as.numeric(k)
  
  r = rowRanks_base(y, cols = I) # n × k, rank mat
  R = colSums2_base(r) #  sum of rank
  o = order(R)
  best = I[o[1L]] #current best
  
  TIES = tapply(c(r), row(r), table)
  
  tie_corr = sum(unlist(lapply(TIES, function(u) { u^3 - u })))
  denom = n_num * k_num * (k_num + 1) - (tie_corr / max(1, (k_num - 1)))
  numer = 12 * sum((R - n_num * (k_num + 1) / 2)^2)
  STATISTIC = numer / denom
  PARAMETER = k - 1L
  #copied from irace package 
  
  if (!is.finite(STATISTIC) || !is.finite(denom) || denom <= 0 || PARAMETER <= 0) {
    PVAL = NA_real_
  } else {
    PVAL = pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  }
  
  alpha = 1 - conf.level
  if (is.finite(PVAL) && !is.na(PVAL) && PVAL < alpha) {
    A = sum(as.vector(r)^2)
    df_now = (n_num - 1) * (k_num - 1)
    tcrit = qt(1 - alpha / 2, df = df_now) *
      (2 * (n_num * A - sum(R^2)) / df_now)^(1 / 2)
    
    J = best
    for (j in 2L:k) {
      if (abs(R[o[j]] - R[o[1L]]) > tcrit) { #if difference is big enough then drop
        break
      } else {
        J = c(J, I[o[j]])
      }
    }
    alive[-J] = FALSE
    dropped.any = TRUE
  }
  
  list(ranks = R, alive = alive, dropped.any = dropped.any, p.value = PVAL) #alive can indicate what is kept. But this vector is only applying on the current iteration, and the result matrix is not yet globally saved 
}

aux_friedman_irace = function(results, conf.level) {
  k = ncol(results)
  alive = rep(TRUE, k)
  which_alive = seq_len(k)
  if (k > 2L) {
    return(aux2_friedman_irace(results, which_alive, alive, conf.level = conf.level))
  }
  
  V1 = results[, 1L]
  V2 = results[, 2L]
  diffs = V1 - V2
  
  dropped.any = TRUE
  PVAL = 0
  
  if (all(diffs <= 0)) {
    ranks = c(1L, 2L)
  } else if (all(diffs >= 0)) {
    ranks = c(2L, 1L)
  } else {
    ZEROES = any(diffs == 0)
    if (ZEROES) diffs = diffs[diffs != 0]
    
    r = rank(abs(diffs))
    TIES = length(r) != length(unique(r))
    
    diffs2 = outer(diffs, diffs, "+")
    diffs2 = sort(diffs2[!lower.tri(diffs2)]) / 2
    
    PVAL = wilcox.test(V1, V2, paired = TRUE,
                       exact = if (ZEROES || TIES) FALSE else NULL)$p.value
    if (!is.finite(PVAL) || is.na(PVAL) || PVAL >= 1 - conf.level) dropped.any = FALSE
    
    ranks = if (median(diffs2) <= 0) c(1L, 2L) else c(2L, 1L)
  }
  
  if (dropped.any) alive[ranks[2L]] = FALSE
  
  list(ranks = ranks, alive = alive, dropped.any = dropped.any, p.value = PVAL)
}
