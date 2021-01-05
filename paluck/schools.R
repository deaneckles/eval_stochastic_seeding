# functions for creating school objects, computing probs

create_school = function(school_i) {
  in_i = which(pd$schid == school_i)
  pd_i = pd[in_i, ]
  adj_i = adj[in_i, in_i]

  #browser()
  g = igraph::graph_from_adjacency_matrix(adj_i) %>%
    igraph::set_vertex_attr("referent", value = pd_i$referent) %>%
    igraph::set_vertex_attr("treat", value = pd_i$treat) %>%
    igraph::set_vertex_attr("student_block", value = pd_i$student_block)

  # summarise n, n_treated for each block
  block_stats = pd_i %>%
    filter(student_block > 0) %>%
    group_by(student_block) %>%
    summarise(
      n = n(),
      n_treated = sum(treat == 1)
    ) %>%
    arrange(student_block)

  # for each block, compute all seed sets
  # and probabilities.
  by_block = foreach(b = 1:4) %do% {
    p_rand_block = choose(
      block_stats$n[b],
      block_stats$n_treated[b]
    )^-1

    # generate all seed sets in block
    valid_seeds = which(pd_i$student_block == b)
    seed_sets_block = combn(
      valid_seeds,
      block_stats$n_treated[b]
    ) %>% t()

    # compute p_nom for each possible seed set
    p_nom_block = compute_p_nom(
      g = adj_i,
      seed_sets = seed_sets_block,
      complete = TRUE,
      valid_seeds = valid_seeds,
      return_pi = TRUE
    )
    
    list(seed_sets = seed_sets_block,
         p_rand = p_rand_block,
         p_nom = p_nom_block$p_nom,
         p_nom_pi = p_nom_block$pi
         )
  }

  # combine probs with product (blocks are independent)
  p_rand = prod(sapply(by_block, function(block) { block$p_rand }))

  school = list(
    schid = school_i,
    adj = adj_i,
    g = g,
    pd = pd_i,
    p_rand = p_rand,
    by_block = by_block
  )

  school['p_nom'] = compute_p_nom_paluck(school)
  school
}

# compute probability of given seed set
# uses precomputed possible seed sets by block
compute_p_nom_paluck = function(school, seed_set = NULL) {
  if (is.null(seed_set)) {
    seed_set = which(school$pd$treat == 1)
  }
  
  seed_set = sort(seed_set)
  
  p_nom_block = foreach(b = 1:4, .combine = c) %do% {
    # subset of seeds in this block
    ss_block = seed_set[seed_set %in% which(school$pd$student_block == b)]

    # find index for this seed set
    obs_idx = apply(
      school$by_block[[b]]$seed_sets,
      1,
      function(x) identical(x, ss_block)
    ) %>% which()
    
    school$by_block[[b]]$p_nom[obs_idx]
  }
  # combine probs with product (blocks are independent)
  prod(p_nom_block)
}

# draw seed sets according to the design
draw_seed_sets = function(school, R = 10) {
  Z = school$pd$treat == 1
  Z[is.na(Z)] = 0
  treat_perm = ri::genperms(
    Z = Z,
    blockvar = school$pd$student_block,
    maxiter = R
  )
  p_nom = apply(treat_perm, 2, function(zp) {
    compute_p_nom_paluck(school = school, seed_set = which(zp == 1))
  })
  list(treat_perm = treat_perm, p_nom = p_nom)
}
