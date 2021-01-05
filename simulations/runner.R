# Runs experiment

source('outcomes.R')

library(foreach)

draw_experiment_data = function(n, villages,
                                model = "LIM",
                                intercept = 0, b_degree = 0,
                                b_spill = 0, n_steps = 0,
                                p_spontaneous = NA, p_infect = NA,
                                prob_assign_nom = 0.5,
                                with_replacement = FALSE, 
                                return_villages = TRUE) {

  if (class(villages) == "function") {
    villages = villages(n)
  } else {
    if (n > length(villages)) stop("n larger than number of available villages.")
    villages = sample(villages, n, replace = with_replacement)
  }
  treatments = rbinom(n, 1, prob_assign_nom) # Bernoulli design

  draw_one_village = function(i) {
    v = villages[[i]]
    
    w = treatments[i]
    
    nom_idx = sample.int(v$n_ss, size=1, prob=v$p_nom)
    rand_idx = sample.int(v$n_ss, size=1, prob=v$p_rand)
    idx = if (w == 1) nom_idx else rand_idx

    result = list(
      p_rand = v$p_rand[idx],
      p_nom = v$p_nom[idx],
      w = w
    ) %>% as.data.frame() %>% mutate(
      p_design = (1 - prob_assign_nom) * p_rand + prob_assign_nom * p_nom
    )

    if (!is.na(model)) {
      nom_seeds = v$seed_sets[nom_idx, ]
      rand_seeds = v$seed_sets[rand_idx, ]

      if (model == "LIM") {
        outcome_function = function(v, seeds) {
          lim_endogenous(v, seeds, intercept, b_degree, b_spill, n_steps)
        }
      } else if (model == "IC") {
        outcome_function = function(v, seeds) {
          independent_cascade(v, seeds, p_infect, p_spontaneous, n_steps)
        }
      }
      
      result = result %>% mutate(
        y_nom = outcome_function(v, nom_seeds),
        y_rand = outcome_function(v, rand_seeds),
        y = y_rand * (w == 0) + y_nom * (w == 1)
        )
    }

    return(result)
  }

  #debug(draw_one_village)
  df = foreach(i = 1:n, .combine = bind_rows) %do% draw_one_village(i)

  return(list(
    sampled_villages = if (return_villages) villages else NA,
    df = df
  ))
}

# Comparison of estimators for given response model
run_experiment = function(villages, rep, model,
                          intercept = NA, b_degree = NA, b_spill = NA, n_steps = NA,
                          p_spontaneous = NA, p_infect = NA,
                          prob_assign_nom = 0.5,
                          with_replacement = FALSE, 
                          n_villages = length(villages)
                          ) {

  data = draw_experiment_data(
    n = n_villages,
    villages = villages,
    model = model,
    intercept = intercept, b_degree = b_degree,
    b_spill = b_spill, n_steps = n_steps,
    p_spontaneous = p_spontaneous,
    p_infect = p_infect,
    with_replacement = with_replacement,
    prob_assign_nom = prob_assign_nom
  )

  df = data$df
  
  mu_rand = mean(df$y_rand)
  mu_nom = mean(df$y_nom)

  df_op = df %>%
    filter(w == 0) %>%
    mutate(p_design = p_rand)

  ## hb = df %>% hajek_bootstrap_variance_estimate
  ## hb_op = df_op %>% hajek_bootstrap_variance_estimate

  ## hperm = hajek_permutation(
  ##   df,
  ##   data$sampled_villages,
  ##   prob_assign_nom = 0.5,
  ##   R = 1000
  ## )
  ## hperm_op = hajek_permutation(
  ##   df_op,
  ##   data$sampled_villages[which(df$w == 0)],
  ##   prob_assign_nom = 0,
  ##   R = 1000
  ## )
  
  return(data.frame(
    rep = rep,
    prob_assign_nom = prob_assign_nom,
    n_villages = n_villages,
    with_replacement = with_replacement,
    model = model,
    intercept = intercept,
    b_degree = b_degree,
    b_spill = b_spill,
    n_steps = n_steps,
    p_infect = p_infect,
    mu_rand = mu_rand,
    mu_nom = mu_nom,
    tau = mu_nom - mu_rand,
    rand_p_rand_sd = sd(df_op$p_rand),
    rand_p_nom_sd = sd(df_op$p_nom),
    rand_p_nom_max = max(df_op$p_nom),
    dm = df %>% difference_in_means,
    ht = df %>% horvitz_thompson,
    hajek = df %>% hajek,
    dm_var_est = df %>% dm_variance_estimate,
    ht_var_est = df %>% ht_variance_estimate,
    hajek_var_est = df %>% hajek_variance_estimate,
    ## ht_off_policy = df_op %>% horvitz_thompson,
    ## ht_off_policy_var_est = df_op %>% ht_variance_estimate,
    ## hajek_off_policy = df_op %>% hajek,
    ## hajek_off_policy_var_est = df_op %>% hajek_variance_estimate,
    ## hajek_boot_var_est = hb$variance,
    ## hajek_boot_ci_low = hb$low,
    ## hajek_boot_ci_high = hb$high,
    hajek_sandwich_var_est = df %>% hajek_sandwich_variance_estimate
    ## hajek_off_policy_boot_var_est = hb_op$variance,
    ## hajek_off_policy_boot_ci_low = hb_op$low,
    ## hajek_off_policy_boot_ci_high = hb_op$high,
    ## hajek_off_policy_sandwich_var_est = df_op %>% hajek_sandwich_variance_estimate#,
    ## hajek_perm_var_est = hperm$variance,
    ## hajek_perm_null_low = hperm$null_low,
    ## hajek_perm_null_high = hperm$null_high,
    ## hajek_perm_p_value = hperm$p_value,
    ## hajek_off_policy_perm_var_est = hperm_op$variance,
    ## hajek_off_policy_perm_null_low = hperm_op$null_low,
    ## hajek_off_policy_perm_null_high = hperm_op$null_high,
    ## hajek_off_policy_perm_p_value = hperm_op$p_value
  ))
}

# Optimal design under a null response.
run_optimal_design = function(villages, rep) {
  n = length(villages)
  
  treatments = rbinom(n, 1, 0.5) # Bernoulli design
  n_nodes = sapply(villages, function(v) v$n_nodes)
  ys = rbinom(n, n_nodes, 0.5) / n_nodes
  
  df_half = foreach(i = 1:n, .combine=rbind) %do% {
    v = villages[[i]]
    w = treatments[i]
    p = if (w == 1) v$p_nom else v$p_rand
    ss_index = sample.int(v$n_ss, size=1, prob=p)
    seeds = v$seed_sets[ss_index, ]
    
    #y = mean(rbinom(v$n_nodes, 1, 0.5)) # pure noise
    y = ys[i]
    
    # probabilities of the observed seed set
    p_rand_ = v$p_rand[ss_index]
    p_nom_ = v$p_nom[ss_index]
    p_half = 0.5 * (p_rand_ + p_nom_)
    
    result = c(
      y = y,
      p_rand = p_rand_,
      p_nom = p_nom_,
      p_half = p_half,
      w = w
    )
    
    return(result)
  } %>% data.frame

  df_opt = foreach(i = 1:n, .combine=rbind) %do% {
    v = villages[[i]]
    ss_index = sample.int(v$n_ss, size=1, prob=v$p_opt)
    seeds = v$seed_sets[ss_index, ]
    
    #y = mean(rbinom(v$n_nodes, 1, 0.5)) # pure noise
    y = ys[i]
    
    # probabilities of the observed seed set
    p_rand_ = v$p_rand[ss_index]
    p_nom_ = v$p_nom[ss_index]
    p_opt_ = v$p_opt[ss_index]
    
    result = c(
      y = y,
      p_rand = p_rand_,
      p_nom = p_nom_,
      p_opt = p_opt_,
      w = w
    )
    
    return(result)
  } %>% data.frame
  
  return(c(
    rep=rep,
    dm = df_half %>% difference_in_means,
    ht_half = df_half %>% mutate(p_design=p_half) %>% horvitz_thompson,
    hajek_half = df_half %>% mutate(p_design=p_half) %>% hajek,
    ht_opt = df_opt %>% mutate(p_design=p_opt) %>% horvitz_thompson,
    hajek_opt = df_opt %>% mutate(p_design=p_opt) %>% hajek,
    p_opt = sum(df_opt$p_opt),
    p_half = sum(df_half$p_half)
  ))
}
