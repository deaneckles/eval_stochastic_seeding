# permutation test logic for Paluck et al. design
source("schools.R")

generate_permutations_paluck = function(df, schools, R) {
  p_nom_perm =  foreach(i = df$schid, .combine = rbind) %dopar% {
      # ri::genperms prints a lot
      invisible(capture.output(
        p <- draw_seed_sets(schools[[as.character(i)]], R = R)$p_nom
      ))
      p
  }
  p_nom_perm
}

hajek_permutation_paluck = function(df, schools, y_var = "pconf", R = 1000, alpha = 0.05,
                                    p_nom_perm = NULL, keep_perms = FALSE,
                                    keep_perm_results = FALSE) {
  df$y = df[[y_var]]
  est = hajek(df)
  est_var = hajek_variance_estimate(df)
  t_obs = est / sqrt(est_var)

  # compute permutations, unless provided
  if (is.null(p_nom_perm)) {
    p_nom_perm = generate_permutations_paluck(df, schools, R)
  }
  # compute stats for each permutation
  perm = foreach(r = 1:ncol(p_nom_perm), .combine = rbind) %dopar% { 
    df$p_nom = p_nom_perm[, r]
    data.frame(
      est = hajek(df),
      est_var = hajek_variance_estimate(df)
      )
  } %>%
    mutate(
      t = est / sqrt(est_var)
    )

  p = c(alpha / 2, 1 - alpha / 2)
  
  result = list(
    est = est,
    t_naive = t_obs,
    variance = var(perm$t * sqrt(est_var)),
    null_low = quantile(perm$t, p[1]),
    null_high = quantile(perm$t, p[2]),
    p_value = 2 * min(mean(t_obs > perm$t), mean(t_obs < perm$t))
  )
  if (keep_perms) {
    result$p_nom_perm = p_nom_perm
  }
  if (keep_perm_results) {
    result$perm_df = perm
  }
  return(result)
}


# naive permutation, just for comparison
dumb_perm = function(df, R = 1000, alpha = 0.05) {
  est = hajek(df)
  est_var = hajek_variance_estimate(df)
  t_obs = est / sqrt(est_var)
  
  perm = foreach(r = 1:R, .combine = rbind) %dopar% {
    df$y = sample(df$y)
    data.frame(
      est = hajek(df),
      est_var = hajek_variance_estimate(df)
      )
  } %>%
    mutate(
      t = est / sqrt(est_var)
    )

  p = c(alpha / 2, 1 - alpha / 2)
  
  data.frame(
    t_naive = t_obs,
    variance = var(perm$t * sqrt(est_var)),
    null_low = quantile(perm$est, p[1]),
    null_high = quantile(perm$est, p[2]),
    p_value = 2 * min(mean(t_obs > perm$t), mean(t_obs < perm$t))
    )
}
