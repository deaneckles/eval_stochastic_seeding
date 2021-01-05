# lim_exogenous = function(v, seeds, intercept, b_spill) {
#   adj = as_adjacency_matrix(v$g)
#   z = rep(0, v$n_nodes)
#   z[seeds] = 1
#   
#   # fraction treated neighbors
#   z_frac_trt = as.vector(.row_norm(adj) %*% z)
#   
#   # ticket redemption: seeds can't do it
#   y = ifelse(z == 1, 0, intercept + beta * z_frac_trt + rnorm(v$n_nodes))
#   ifelse(y > 0, 1, 0) %>% mean
# }

.row_norm = function(adj) {
  N = nrow(adj)
  degrees = apply(adj, 1, sum)
  inv.degrees = ifelse(degrees == 0, 0, 1 / degrees)
  adj * matrix(rep(inv.degrees, N), nrow=N)
}

# linear-in-means
lim_endogenous = function(v, seeds, intercept, b_degree, b_spill, n_steps) {
  adj = as_adjacency_matrix(v$g)
  
  # initial seeds
  z = rep(0, v$n_nodes)
  z[seeds] = 1
  
  # seeds are the only initial adopters
  y = z
  
  norm_adj = .row_norm(adj)
  sum_degrees = sum(v$degree[seeds]) # sum of degrees in seed set
  
  for (t in 1:n_steps) {
    avg_nbr_y = as.vector(norm_adj %*% y)
    y_star = intercept + b_degree * sum_degrees + b_spill * avg_nbr_y + rnorm(v$n_nodes)
    y = ifelse(y == 1, y, 1 * (y_star > 0)) # adopters can't revert
  }
  
  return(mean(y))
}

# SIR with deterministic I->R transition
independent_cascade = function(v, seeds, p_infect, p_spontaneous = 0, n_steps = Inf) {
  adj = as_adjacency_matrix(v$g)
  if (is.na(p_spontaneous)) p_spontaneous = 0
  if (is.na(n_steps)) n_steps = Inf
  
  # initial seeds
  z = rep(0, v$n_nodes)
  z[seeds] = 1
  
  # seeds are the only initial adopters
  y = z
  y_new = y
  
  norm_adj = .row_norm(adj)
  sum_degrees = sum(v$degree[seeds]) # sum of degrees in seed set

  t = 0
  while (TRUE) {
    nbr_y_new = as.vector(adj %*% y_new)
    prob_adopt = ifelse(y == 1, 0, 1 - (1 - p_spontaneous) * (1 - p_infect)^nbr_y_new)
    y_new = rbinom(length(y), 1, prob_adopt)
    
    y = pmax(y, y_new)
    
    if (sum(y_new) == 0) break
    t = t + 1
    if (t >= n_steps) break
  }
  
  return(mean(y))

}
