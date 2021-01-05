# Village object.

compute_p_nom_one_seed = function(g, valid_seeds = NULL) {
  if (is(g, 'sparseMatrix')) {
    adj = g
  } else {
    adj = as_adj(g)
  }
  
  if (!is.null(valid_seeds)) {
    # remove edges to invalid seeds
    adj[, -valid_seeds] = 0
  }
  out_degree = Matrix::rowSums(adj)
  
  valid_start = out_degree > 0
  adj = adj[valid_start, , drop = FALSE]
  out_degree = out_degree[valid_start]

  one_seed_p = Matrix::colSums(adj * (1 / out_degree)) / nrow(adj)
  one_seed_p
}

compute_p_nom = function(g, 
                         seed_sets, 
                         complete = TRUE, 
                         approx_pi = NULL, 
                         pi_dp = TRUE,
                         valid_seeds = NULL, 
                         return_pi = FALSE,
                         R=10000) {
  one_seed_p = compute_p_nom_one_seed(g, valid_seeds)
  names(one_seed_p) <- c()
  # The probability of a set, with replacement, is the probability 
  # of a single unordered set. To get the probability of the unordered set
  # requires multiplying by k!.
  p_nom = apply(seed_sets, 1, function(ss) { prod(one_seed_p[ss]) })
  k = ncol(seed_sets)
  if (complete) {
    pi = sum(p_nom)
    p_nom = p_nom / pi
    # put k! back in for use elsewhere
    pi = factorial(k) * pi
  } else if (pi_dp) {
    pi = pi_dynamicprogram(g=g, k=k)
  } else if (!is.null(approx_pi)) {
    pi = approx_pi(g=g, k=k, R=R)
    p_nom = p_nom / pi
  } else {
    stop("No normalization selected.")
  }
  if (sum(p_nom) < 1e-10) stop("Rounding error")
  if (return_pi) {
    return(list(
      pi = pi,
      p_nom = p_nom
    ))
  }
  p_nom
}

# bounds on variance of pi
pi_variance_bound_graph = function(village, k, R) {
  m = village$n_nodes
  
  cvec = sort(compute_p_nom_one_seed(village$g))
  cmink = prod(cvec[1:k])
  cmaxk = prod(cvec[(m-k+1):m])

  bound = (
    (1/(4 * R)) * (factorial(k) * choose(m, k))^2
    *(cmaxk - cmink)^2
    )
  bound
}

pi_variance_bound_degree = function(village, k, R) {
  m = village$n_nodes
  
  dinmax = max(degree(village$g,mode='in'))
  doutmax = max(degree(village$g,mode='out'))
  dinmin = 1
  doutmin = 1

  bound = (
    (1/(4 * R)) * (factorial(k) * choose(m, k) / m^k)^2
    * ((dinmax / doutmin)^k - (dinmin / doutmax)^k)^2
  )
  bound
}

# estimate probability that seed set drawn with replacement has no repeats
approx_pi = function(g, k, R = 1000) {
  one_seed_p = compute_p_nom_one_seed(g)
  names(one_seed_p) <- c()
  seed_sets = matrix(0, R, k)
  for (i in 1:R) {
    seed_sets[i,] = sample(one_seed_p, size = k, replace = FALSE)
    }
  prob_with_replacement = apply(seed_sets, 1, prod)
  
  pi_hat = factorial(k)*(choose(vcount(g), k)) * mean(prob_with_replacement)
  pi_hat
}

# compute probability that seed set drawn with replacement has no repeats via DP
pi_dynamicprogram = function(g, k) {
  # g - graph
  # k - size of seed set 
  fvec = compute_p_nom_one_seed(g)
#  fvec = sample(fvec) # Sanity check to make sure value was perm invariant: yes
  ni = length(fvec)
  names(fvec) <- c()

  S = matrix(0,ni+1,k+1)
  S[1,1] = 0 # Doesn't matter, but to stay sane
  for (j in 1:k) {
    S[j+1,j+1] = prod(fvec[1:j])
  }
  for (j in 1:(ni+1)) {
    S[j,1] = 1
  }
  for (j in 2:ni) { #2:ni
      for (ell in 1:min((j-1),k)) { #1:(j-1)
#        print(c(j,ell))
        S[j+1,ell+1] = S[(j-1)+1,(ell-1)+1]*fvec[j] + S[(j-1)+1,ell+1] 
      }
    }
  factorial(k)*S[ni+1,k+1]  # pi
}


approx_pi_dependent = function(g, k, R = 1000, truncate = TRUE) {
  one_seed_p = compute_p_nom_one_seed(g)
  names(one_seed_p) <- c()
  sets_per = floor(length(one_seed_p) / k)
  reps = ceiling(R / sets_per)

  seed_sets = matrix(0, sets_per * reps, k)
  for (i in 1:reps) {
    row = (i - 1) * sets_per + 1
    seed_sets[row:(row + sets_per - 1), ] = matrix(
      sample(one_seed_p, size = k * sets_per, replace = FALSE),
      ncol = k
      )
  }
  prob_with_replacement = apply(seed_sets, 1, prod)

  if (truncate)
    prob_with_replacement = prob_with_replacement[1:R]
  
  pi_hat = factorial(k)*(choose(vcount(g), k)) * mean(prob_with_replacement)
  pi_hat
}



#' @param g graph
#' @param name village name
#' @param ss_size num individuals in each seed set
#' @param approx_pi whether to use DP approximation
#' @param R DP reps
#' @param n_ss_samples number of seed sets to sample; if NULL, get all combinations
Village = function(g, name, ss_size, approx_pi = FALSE, R = 1000, n_ss_samples = NULL) {
  n_nodes = vcount(g)
  n_ss = choose(n_nodes, ss_size)
  if (is.null(n_ss_samples)) {
    seed_sets = t(combn(n_nodes, ss_size))
    p_rand = rep(1 / n_ss, n_ss)
    # combn is faster than gtools::combinations() but returns transpose,
    # so had to make some other small changes
  } else {
    # samples seed sets independently, then drops duplicates.
    seed_sets = sapply(1:n_ss_samples, function(i) sort(sample(n_nodes, ss_size, replace=FALSE)))
    if (ss_size == 1) {
      seed_sets = unique(as.matrix(seed_sets))
    } else {
      seed_sets = unique(t(seed_sets))
    }
    p_rand = rep(1 / nrow(seed_sets), nrow(seed_sets))
  }
  p_nom_results = compute_p_nom(g, 
                                seed_sets, 
                                pi_dp = TRUE, 
                                approx_pi = approx_pi,
                                return_pi = TRUE,
                                R = R)

  p_nom = p_nom_results$p_nom
  pi = p_nom_results$pi
  p_opt = abs(p_rand - p_nom)
  p_opt = p_opt / sum(p_opt)
  structure(
    list(
      ss_size = ss_size,
      g = g,
      degree = degree(g, mode = "in"),
      name = name,
      n_nodes = n_nodes,
      n_ss = n_ss,
      seed_sets = seed_sets,
      p_rand = p_rand,
      p_nom = p_nom,
      p_opt = p_opt,
      p_nom_pi = pi
    ),
    class = 'village'
  )
}

create_village_cai = function(data, village_name, ss_size=2, approx_pi = FALSE, R = 1000, n_ss_samples = NULL) {
  edgelist_within = data %>% 
    dplyr::filter(
      address == village_name,
      network_address == village_name,
      !is.na(network_id),
      id != 99,
      network_id != 99
    )
  if(nrow(edgelist_within) > 0) {
    g <- edgelist_within %>% 
      dplyr::select(id, network_id) %>%
      as.matrix %>% 
      apply(2, as.character) %>%
      matrix(ncol = 2) %>% 
      graph_from_edgelist(directed=TRUE)
    
    return(Village(g, village_name, ss_size = ss_size, approx_pi = approx_pi, R = R, n_ss_samples = n_ss_samples))
  } else {
    return(NA)
  }
}
