# plotting and coverage

library(ggplot2)

strategy_name_map = c(
  'nom' = 'one-hop',
  'rand' = 'random'
)
estimator_name_map = c(
  'dm' = 'difference-in-means',
  'hajek' = 'Hajek',
  'ht' = 'Horvitz-Thompson',
  'truth' = 'truth'
)
palette1 = c(
  "difference-in-means" = "#E69F00",
  "Horvitz-Thompson" = "#56B4E9",
  "Hajek" = "#009E73",
  "truth" = 'black',
  'random' = 'black',
  'one-hop' = 'red',
  "Bernoulli(1/2)" = "black",
  "optimized" = "#018571",
  "random (off-policy)" = "#a6611a",
  "on-policy" = NA
)


shapes1 = c(
  "difference-in-means" = 22, "Horvitz-Thompson" = 21, "Hajek" = 24,
  "truth" = 16,
  'random' = 17,
  'one-hop' = 19,
  "Bernoulli(1/2)" = 16,
  "optimized" = 17,
  "random (off-policy)" = 18
)

shapes2 = c(
  'random' = 2,
  'one-hop' = 1
)

shapes_off_policy = c(
  "difference-in-means (0.5)" = 22, "Horvitz-Thompson (0.5)" = 21, "Hajek (0.5)" = 24,
  #"difference-in-means (0)" = 15,
  "Horvitz-Thompson (0)" = 19, "Hajek (0)" = 17
)

palette_off_policy = rep(palette1, 2)
names(palette_off_policy) = sprintf("%s (%s)", names(palette_off_policy), rep(c("0", "0.5"), each = length(palette1)))
palette_off_policy = palette_off_policy[names(shapes_off_policy)]

palette2 = c("#0072B2", "#D55E00", "#CC79A7")

dataset_name_map = c(
  "cai" = "Cai et al.",
  "paluck" = "Paluck et al.",
  "addhealth" = "AddHealth",
  "microfinance" = "Banerjee et al.",
  "chami" = "Chami et al."
  )


type_map = c(
  'dm' = 'est',
  'ht' = 'est',
  'hajek' = 'est',
  'dm_var_est' = 'var_est',
  'ht_var_est' = 'var_est',
  'hajek_var_est' = 'var_est',
  ## 'ht_off_policy' = 'est',
  ## 'ht_off_policy_var_est' = 'var_est',
  ## 'hajek_off_policy' = 'est',
  ## 'hajek_off_policy_var_est' = 'var_est',
  'hajek_boot_var_est' = 'var_est_bootstrap',
  'hajek_boot_ci_low' = 'boot_ci_low',
  'hajek_boot_ci_high' = 'boot_ci_high',
  'hajek_sandwich_var_est' = 'var_est_sandwich',
  ## 'hajek_off_policy_boot_var_est' = 'var_est_bootstrap',
  ## 'hajek_off_policy_boot_ci_low' = 'boot_ci_low',
  ## 'hajek_off_policy_boot_ci_high' = 'boot_ci_high',
  ## 'hajek_off_policy_sandwich_var_est' = 'var_est_sandwich'#,
  'hajek_perm_var_est' = 'var_est_perm',
  'hajek_perm_null_low' = 'perm_null_low',
  'hajek_perm_null_high' = 'perm_null_high',
  'hajek_perm_p_value' = 'perm_p_value'#,
  ## 'hajek_off_policy_perm_var_est' = 'var_est_perm',
  ## 'hajek_off_policy_perm_null_low' = 'perm_null_low',
  ## 'hajek_off_policy_perm_null_high' = 'perm_null_high',
  ## 'hajek_off_policy_perm_p_value' = 'perm_p_value'
  )

rmse = function(x, y) { sqrt(mean((x-y)^2)) }

plot_CI = function(results, alpha) {
  # create summary df
  means = results %>% dplyr::select(rep, dm, ht, hajek) %>% gather(est, mean, -rep)
  vars = results %>% dplyr::select(rep, dm_var_est, ht_var_est, hajek_var_est) %>% gather(est, var, -rep) %>%
    mutate(est = gsub('dm_var_est', 'dm', est), est = gsub('ht_var_est', 'ht', est), est = gsub('hajek_var_est', 'hajek', est))
  df = means %>% left_join(vars) %>% mutate(
    min=mean-qnorm(1 - alpha/2)*sqrt(var), 
    max=mean+qnorm(1 - alpha/2)*sqrt(var),
    covers=(min < 0) & (max > 0)
  ) %>% group_by(est) %>% arrange(mean) %>% mutate(sorted_rep=1:max(rep)) %>% ungroup
  
  # plot
  df %>% ggplot(aes(x=sorted_rep, y=mean, ymin=min, ymax=max)) + 
    geom_errorbar(aes(colour=covers), alpha=0.2) + 
    geom_hline(aes(yintercept=0)) + facet_grid(. ~ est) +
    theme_bw()
}

.coverage = function(est, var_est, true_mean, alpha=0.2, df = Inf) {
  accept = ifelse(
    var_est == 0,
    est == true_mean,
    pt(abs((est - true_mean) / sqrt(var_est)), df) < 1 - alpha/2
  )
  mean(accept)
  #z_stat = (est - true_mean) / sqrt(var_est)
  #mean(pnorm(abs(z_stat)) < 1 - alpha/2)
}

coverage = function(results, n, alpha) {
  results %>% summarise(
    dm_coverage = .coverage(dm, dm_var_est, mean(dm), alpha),
    ht_coverage = .coverage(ht, ht_var_est, mean(ht), alpha),
    hajek_coverage = .coverage(hajek, hajek_var_est, mean(hajek), alpha)
  )
}
