difference_in_means = function(df) {
  with(df, mean(y[w==1]) - mean(y[w==0]))
}

horvitz_thompson = function(df) {
  n = nrow(df)
  df = df %>% filter(p_design > 0)
  with(df, 1 / n * sum(y * (p_nom - p_rand) / p_design))
}

hajek = function(df) {
  df = df %>% filter(p_design > 0)
  mu_nom = with(df, weighted.mean(y, p_nom / p_design))
  mu_rand = with(df, weighted.mean(y, p_rand / p_design))
  mu_nom - mu_rand
}

dm_variance_estimate = function(df) {
  with(df, var(y[w==1]) / sum(w==1) + var(y[w==0]) / sum(w==0))
}

ht_variance_estimate = function(df) {
  ht = horvitz_thompson(df)
  n = nrow(df)
  wt = with(df, (p_nom - p_rand) / p_design)
  1 / n * with(df, mean((wt * y - ht)^2))
}

hajek_variance_estimate = function(df) {
  mu_nom = with(df, weighted.mean(y, p_nom / p_design))
  mu_rand = with(df, weighted.mean(y, p_rand / p_design))
  n = nrow(df)
  1 / n * with(df, mean(((mu_nom * p_nom - mu_rand * p_rand - y * (p_nom - p_rand)) / p_design)^2))
}

# simple nonparametric bootstrap
loadNamespace("rsample")
loadNamespace("purrr")
hajek_bootstrap_variance_estimate = function(df, R = 1000, alpha = .1) {
  mu_nom = with(df, weighted.mean(y, p_nom / p_design))
  mu_rand = with(df, weighted.mean(y, p_rand / p_design))
  tau_hat = hajek(df)

  boot_df = df %>%
    rsample::bootstraps(times = R)
    
  ests = purrr::map(boot_df$splits, function(s) hajek(rsample::analysis(s))) %>% unlist()

  # bias corrected percentile bootstrap intervals
  z_0 = qnorm(mean(ests <= tau_hat))

  p = c(alpha / 2, 1 - alpha / 2)
  p = pnorm(2 * z_0 + qnorm(p))

  data.frame(
    variance = var(ests),
    low = quantile(ests, p[1]),
    high = quantile(ests, p[2])
    )
}

loadNamespace("estimatr")
hajek_sandwich_variance_estimate = function(df) {
  df2 <- rbind(df, df) %>%
    mutate(
      w = rep(c(0, 1), each = nrow(df)),
      id = rep(1:nrow(df), 2)
    )
  fit <- estimatr::lm_robust(
    y ~ w,
    weights = ifelse(w==1, p_nom, p_rand)/p_design,
    data = df2,
    cluster = id
  )
  fit$vcov[2, 2]
}


test_beta = function(df, p_nom_perms, beta) {
  df = df %>% mutate(
    score = (p_nom - p_rand) / p_design
  )
  approx_tau = function(beta, p_nom_perm) {
    s = df %>% mutate(
      score = (p_nom - p_rand) / p_design,
      score_perm = (p_nom_perm - p_rand) / p_design,
      y_perm = y + (score_perm - score) * beta
    ) %>%
      ungroup() %>%
      summarise(
        mu_nom = weighted.mean(y_perm, p_nom_perm),
        mu_rand = weighted.mean(y_perm, p_rand),
        tau = mu_nom - mu_rand
      )
    s$tau
  }

  results = apply(p_nom_perms, 2, function(p_nom_perm) {
    approx_tau(betas[1], p_nom_perm)
  })
  results
}

find_acceptance_region = function(est, est_var, est_perm, est_var_perm, alpha) {
  browser()
  get_p_value <- function(tau_0) {
    #t_obs = (est - tau_0) / sqrt(est_var)
    t_obs = (est) / sqrt(est_var)
    t_perm = (est_perm - tau_0) / sqrt(est_var_perm)
    p_value = 2 * min(mean(t_obs > t_perm), mean(t_obs < t_perm))
    return(p_value)
  }
  objective = function(tau_0, side = 1) {
    p = get_p_value(tau_0, est, est_var, est_perm, est_var_perm)
    ifelse(
      sign(tau_0 - est) == side & p - alpha > 0,
      p - alpha,
      2 + (p - alpha)
    )
  }
  result_upper = optim(
    par = est,
    fn = objective,
    method = "BFGS",
    side = 1
  )
  result_lower = optim(
    par = est,
    fn = objective,
    method = "BFGS",
    side = -1
  )
  c(ci_low = result_lower$par, ci_high = result_upper$par)
}


hajek_permutation = function(df, villages, prob_assign, R = 100, alpha = .1) {
  tau_hat = hajek(df)
  tau_hat_var = hajek_sandwich_variance_estimate(df)
  t_obs = (tau_hat) / sqrt(tau_hat_var)

  perm = foreach(r = 1:R, .combine = rbind) %do% {
    df_perm = draw_experiment_data(
      villages,
      0, 0, 0, 0,
      with_replacement = FALSE,
      prob_assign = prob_assign,
      return_villages = FALSE
    )$df
    df_perm$y = df$y
    
    data.frame(
      est = hajek(df_perm),
      var_est = hajek_sandwich_variance_estimate(df_perm)
    )
  } %>%
    mutate(
      t = (est) / sqrt(var_est)
    )

  p = c(alpha / 2, 1 - alpha / 2)

  #find_acceptance_region(tau_hat, tau_hat_var, perm$est, perm$var_est, alpha)
  
  data.frame(
    variance = var(perm$t * sqrt(tau_hat_var)),
    null_low = quantile(perm$t, p[1]) * sqrt(tau_hat_var),
    null_high = quantile(perm$t, p[2]) * sqrt(tau_hat_var),
    p_value = 2 * min(mean(t_obs > perm$t), mean(t_obs < perm$t))
    )
}
