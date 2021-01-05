# Analyze simulations with independent cascade model

library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(reshape2)
library(gridExtra)
library(grid)
library(scales)
library(readr)
library(arrow)

source('analyze.R')


# read in estimates, either from csv or parquet

results = arrow::read_parquet(
  "results/main_results.parquet",
  as_data_frame = TRUE
)

results = results %>% filter(model == "IC")


# construct data frame for true means and tau
df_tau = results %>% 
  group_by(model, n_steps, p_infect) %>% 
  summarise(rand=mean(mu_rand), nom=mean(mu_nom), tau=mean(tau))

# plot true means
p_truth = df_tau %>%
  gather(strategy, mu, rand, nom) %>%
  mutate(
    strategy = strategy_name_map[strategy]
    ) %>%
  ggplot(
    aes(
      x = p_infect, y = mu,
      group = strategy,
      colour = strategy,
      shape = strategy
    )) +
  geom_point() +
  geom_line() +
  scale_colour_manual(values = palette1) +
  scale_shape_manual(values = shapes1) +
  scale_x_continuous(breaks = unique(df_tau$p_infect)) +
  guides(col = guide_legend(nrow = 1)) +
  xlab(expression('cascade probability '*p)) +
  ylab('mean adoption') +
  #ggtitle('True simulated mean adoption rates') + theme_bw() + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_truth

# construct data frame for estimates

results_long = results %>%
  melt(id.vars = names(results)[c(1:5, 9:16)]) %>%
  rename(estimator = variable) %>%
  mutate(estimator = as.character(estimator)) %>%
  mutate(
    type = type_map[estimator],
    estimator = gsub("_(var_est|boot_var_est|sandwich_var_est|boot_ci_high|boot_ci_low|perm_var_est|perm_null_low|perm_null_high|perm_p_value)", "", estimator)
  ) %>%
  filter(!is.na(type)) %>%
  dcast(
    rep + prob_assign_nom + n_villages + with_replacement + model +
      n_steps + p_infect + 
      mu_rand + mu_nom + tau +
      rand_p_rand_sd + rand_p_nom_sd + rand_p_nom_max +
      estimator ~ type
  ) %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    n_steps,
    p_infect
  ) %>%
  mutate(
    tau_mean = mean(tau),
    tau_se = sd(tau) / sqrt(n())
    )
  
df = results_long %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    n_steps,
    p_infect,
    estimator
  ) %>%
  summarise(
    tau = mean(tau),
    tau_hat = mean(est),
    sd = sd(est),
    sd_est = mean(sqrt(var_est)),
    #sd_est_boot = mean(sqrt(var_est_bootstrap)),
    sd_est_sandwich = mean(sqrt(var_est_sandwich)),
    rmse = rmse(tau, est)
  ) %>%
  mutate(
    bias = tau_hat - tau
    )

df_with_truth = df %>%
  dplyr::select(
    prob_assign_nom, n_villages, with_replacement, model,
    n_steps,
    p_infect, value = tau_hat, estimator
  ) %>%
  rbind(
    df_tau %>% mutate(estimator='truth', off_policy = F) %>%
      dplyr::select(
        model,
        n_steps,
        p_infect, value = tau, estimator
      )
  )

# plot estimates + true tau
tmp_truth = df_with_truth %>% filter(estimator == "truth") %>%
  ungroup() %>% select(estimator, p_infect, value)

p_tau = df_with_truth %>%
  mutate(
    estimator = estimator_name_map[estimator]
  ) %>%
  filter(prob_assign_nom == 0.5, estimator != "truth") %>%
  ggplot(aes(
    x = p_infect, y = value,
    group = estimator, colour = estimator,
    shape = estimator
  )) +
  geom_hline(yintercept=0, linetype='dashed') +
  geom_line(data = tmp_truth) + geom_point(data = tmp_truth) +
  geom_point() + geom_line() +
  facet_grid(n_villages ~ ., labeller=label_bquote(
    #cols=rho == .(prob_assign_nom),
    rows=N == .(n_villages)
  )) +
  theme_bw() +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_x_continuous(breaks = unique(df_tau$p_infect)) +
  guides(col = guide_legend(nrow = 2)) +
  xlab(expression('cascade probability '*p)) +
  ylab('mean treatment effect estimate') +
  #ggtitle('Treatment effects and estimates') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_tau

p_true_tau = grid.arrange(p_truth, p_tau, nrow=1, ncol = 2, widths = c(3, 5))
p_true_tau
ggsave('../figures/sim-IC-true-tau.pdf', p_true_tau, width = 7.5, height = 4.5)


####
# Error: bias, SE, RMSE

# plot bias
p_bias = df %>%
  mutate(
    estimator = estimator_name_map[estimator],
    est_design = sprintf("%s (%s)", estimator, prob_assign_nom)
  ) %>%
  filter(!is.na(tau_hat)) %>%
  ggplot(aes(p_infect, bias, group=est_design, colour=est_design, shape=est_design)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(n_villages ~ ., labeller=label_bquote(
    #cols=rho == .(prob_assign_nom),
    rows=N == .(n_villages)
  )) +
  #scale_y_continuous(limits=c(-0.02, 0.02)) + 
  theme_bw() +
  scale_colour_manual(name = expression(paste("estimator (", rho, ")")), values=palette_off_policy) +
  scale_shape_manual(name = expression(paste("estimator (", rho, ")")), values=shapes_off_policy) +
  scale_x_continuous(breaks = unique(df_tau$p_infect)) +
  xlab(expression('cascade probability '*p)) +
  ylab('bias') +
  #ggtitle('Bias') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 3))
p_bias

ggsave('../figures/sim-IC-bias.pdf', p_bias, width=6, height=6)

# plot standard error
p_se = df %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(p_infect, sd, group=estimator, colour=estimator, shape=estimator)) +
  geom_point() + geom_line() +
  facet_grid(n_villages ~ prob_assign_nom, labeller=label_bquote(
    cols=rho == .(prob_assign_nom),
    rows=N == .(n_villages)
  )) +
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed') + 
  #scale_y_continuous(limits=c(0, 0.1)) + 
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  xlab(expression('cascade probability '*p)) +
  ylab('standard error') +
  ggtitle('Standard error') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_se

# plot RMSE
p_rmse = df %>%
  mutate(
    estimator = estimator_name_map[estimator],
    est_design = sprintf("%s (%s)", estimator, prob_assign_nom)
  ) %>%
  filter(!is.na(tau_hat)) %>%
  ggplot(aes(p_infect, rmse * sqrt(n_villages), group=est_design, colour=est_design, shape=est_design)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(n_villages ~ ., labeller=label_bquote(
    #cols=rho == .(prob_assign_nom),
    rows=N == .(n_villages)
  )) +
  #scale_y_continuous(limits=c(-0.02, 0.02)) + 
  theme_bw() +
  scale_colour_manual(name = expression(paste("estimator (", rho, ")")), values=palette_off_policy) +
  scale_shape_manual(name = expression(paste("estimator (", rho, ")")), values=shapes_off_policy) +
  scale_x_continuous(breaks = unique(df_tau$p_infect)) +
  xlab(expression('cascade probability '*p)) +
  ylab(expression(RMSE %*% sqrt(N))) +
  #ggtitle('Bias') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 3))
p_rmse

ggsave('../figures/sim-IC-rmse.pdf', p_rmse, width=6, height=6)


###
# Coverage rates
nominal = 0.9
nominal.null.95 = nominal - 1.96 * sqrt(nominal * (1 - nominal) / max(results$rep))
alpha = 1 - nominal

df_coverage = results_long %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    n_steps, p_infect, estimator
  ) %>%
  summarise(
    tau_se = sd(tau) / sqrt(n()),
    tau = mean(tau),
    clt = .coverage(est, var_est, tau, alpha = alpha),
    sandwich = .coverage(est, var_est_sandwich, tau, alpha = alpha)#,
    #boot = .coverage(est, var_est_bootstrap, tau, alpha = alpha),
    #boot_percentile = mean(tau <= boot_ci_high & tau >= boot_ci_low),
    #perm = .coverage(est, var_est_perm, tau, alpha = alpha),
    #perm_alt_power = mean(perm_p_value < alpha)
  )

df_coverage %>% filter(sandwich <= .88) %>% print(n = 100)

df_coverage_null = results_long %>%
  filter(p_infect == 0) %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    n_steps, p_infect, estimator
  ) %>%
  summarise(
    tau = mean(tau),
    clt = .coverage(est, var_est, 0, alpha = 1 - nominal),
    sandwich = .coverage(est, var_est_sandwich, 0, alpha = 1 - nominal)#,
    #boot = .coverage(est, var_est_bootstrap, 0, alpha=1 - nominal),
    #boot_percentile = mean(0 <= boot_ci_high & 0 >= boot_ci_low)
  )

p_cov = df_coverage %>%
  mutate(
    estimator = estimator_name_map[estimator],
    est_design = sprintf("%s (%s)", estimator, prob_assign_nom)
    ) %>% filter(!is.na(clt)) %>%
  ggplot(aes(p_infect, clt, group=est_design, colour=est_design, shape=est_design)) +
  geom_point() + geom_line() +
  facet_grid(n_villages ~ ., labeller=label_bquote(
    #cols=rho == .(prob_assign_nom),
    rows=N == .(n_villages)
  )) +
  geom_hline(yintercept=nominal, linetype='dashed') +
  geom_rect(xmin = -1, xmax = 10, ymin=nominal.null.95, ymax = 1, inherit.aes = FALSE, alpha = .005) + 
  scale_colour_manual(name = expression(paste("estimator (", rho, ")")), values=palette_off_policy) +
  scale_shape_manual(name = expression(paste("estimator (", rho, ")")), values=shapes_off_policy) +
  #scale_y_continuous(limits=c(0.5, 1)) + 
  xlab(expression('cascade probability '*p)) +
  ylab('coverage rate') +
  #ggtitle('Coverage rates for 90% nominal confidence interval') + theme_bw() +
  guides(col = guide_legend(nrow = 1)) +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 3))
p_cov

ggsave('../figures/sim-IC-coverage-only.pdf', p_cov, width=6, height=6)

###
# power
df_power = results %>% 
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    n_steps, p_infect
  ) %>%
  filter(p_infect > 0) %>%
  summarise(
    tau = mean(tau),
    dm = .coverage(dm, dm_var_est, 0, alpha = 1 - nominal),
    ht = .coverage(ht, ht_var_est, 0, alpha = 1 - nominal),
    hajek = .coverage(hajek, hajek_var_est, 0, alpha = 1 - nominal)
  ) %>%
  gather(estimator, coverage, dm, ht, hajek) %>%
  mutate(power = 1 - coverage) %>% select(-coverage) %>%
  inner_join(
    df_coverage %>% select(coverage = clt)
  )

p_power_basic = df_power %>%
  #filter(coverage >= .85) %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(
    aes(x = p_infect, y = power,
        group = estimator, colour = estimator, shape = estimator
        )) +
  geom_point() + geom_line() +
  facet_grid(n_villages ~ prob_assign_nom, labeller=label_bquote(
    cols=rho == .(prob_assign_nom),
    rows=N == .(n_villages)
  )) +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_y_continuous(breaks = seq(0, 1, by = .2)) +
  guides(col = guide_legend(nrow = 1)) +
  xlab(expression('spillover effect '*beta)) + ylab('power') +
  #ggtitle('Power of estimators') +
  theme_bw() + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))
p_power_basic

ggsave('../figures/sim-IC-power-basic.pdf', p_power_basic, width=6, height=5)
