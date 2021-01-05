# Analyze simulations with linear-in-means utility model

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

# read in estimates

results = arrow::read_parquet(
  "results/main_results.parquet",
  as_data_frame = TRUE
)

results = results %>% filter(model == "LIM")

# SUBSET
#results = results %>% filter(intercept %in% c(-3, -2, -1, 0), b_degree != 0.2)


# construct data frame for true means and tau
df_tau = results %>% 
  group_by(model, intercept, b_degree, b_spill, n_steps) %>% 
  summarise(rand=mean(mu_rand), nom=mean(mu_nom), tau=mean(tau))

df_tau_long = df_tau %>%
  gather(strategy, mu, rand, nom) %>%
  mutate(
    strategy = strategy_name_map[strategy]
  )

# plot true means
p_truth = df_tau_long %>%
  filter(b_degree == 0) %>%
  ggplot(
    aes(
      x = b_spill, y = mu,
      group = strategy,
      colour = strategy,
      shape = strategy
    )) +
  geom_point() +
  geom_line() +
  facet_grid( ~ intercept, labeller=label_bquote(
    cols=alpha == .(intercept)
  )) +
  scale_colour_manual(values = palette1) +
  scale_shape_manual(values = shapes1) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  guides(col = guide_legend(nrow = 2)) +
  xlab(expression('spillover effect '*beta)) + ylab('mean adoption') +
  #ggtitle('True simulated mean adoption rates') + theme_bw() + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_truth

p_truth_all = p_truth %+% df_tau_long +
  facet_grid(b_degree ~ intercept, labeller=label_bquote(
    cols=alpha == .(intercept),
    rows=gamma == .(b_degree)
  ))
p_truth_all

# construct data frame for estimates

results_long = results %>%
  melt(id.vars = names(results)[c(1:16)]) %>%
  rename(estimator = variable) %>%
  mutate(estimator = as.character(estimator)) %>%
  mutate(
    type = type_map[estimator],
    estimator = gsub("_(var_est|boot_var_est|sandwich_var_est|boot_ci_high|boot_ci_low|perm_var_est|perm_null_low|perm_null_high|perm_p_value)", "", estimator)
  ) %>%
  dcast(
    rep + prob_assign_nom + n_villages + with_replacement + model +
      intercept + b_degree + b_spill + n_steps +
      mu_rand + mu_nom + tau +
      rand_p_rand_sd + rand_p_nom_sd + rand_p_nom_max +
      estimator ~ type
  ) %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps
  ) %>%
  mutate(
    tau_mean = mean(tau),
    tau_se = sd(tau) / sqrt(n())
    )
  
df = results_long %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps,
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

df_se_comparison = df %>%
  filter(estimator %in% c("dm", "hajek")) %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps
  ) %>%
  mutate(
    sd_dm_ratio = sd / sd[estimator == "dm"]
  ) %>%
  filter(estimator == "hajek")

summary(1-df_se_comparison$sd_dm_ratio)
summary(1/df_se_comparison$sd_dm_ratio^2 - 1)

df_with_truth = df %>%
  dplyr::select(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps,
    value=tau_hat, estimator
  ) %>%
  rbind(
    df_tau %>% mutate(estimator='truth', off_policy = F) %>%
      dplyr::select(
        model,
        intercept, b_degree, b_spill, n_steps,
        value=tau, estimator
      )
  )

# plot estimates + true tau

tmp_truth = df_with_truth %>%
  filter(estimator == "truth") %>%
  ungroup() %>% select(estimator, intercept, b_degree, b_spill, value)

p_tau = df_with_truth %>%
  filter(prob_assign_nom == 0.5, b_degree == 0, estimator != "truth") %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(
    x = b_spill, y = value,
    group = estimator, colour = estimator,
    shape = estimator
  )) +
  geom_line(data = tmp_truth %>% filter(b_degree == 0)) +
  geom_point(data = tmp_truth %>% filter(b_degree == 0)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(n_villages ~ intercept, labeller=label_bquote(
    rows = N == .(n_villages),
    cols = alpha == .(intercept)
  )) +
  theme_bw() +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  guides(col = guide_legend(nrow = 2)) +
  xlab(expression('spillover effect '*beta)) + ylab('treatment effect estimate') +
  #ggtitle('Treatment effects and estimates') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_tau

p_tau_less = df_with_truth %>%
  filter(prob_assign_nom == 0.5, n_villages == 150,
         b_degree == 0, estimator != "truth") %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(
    x = b_spill, y = value,
    group = estimator, colour = estimator,
    shape = estimator
  )) +
  geom_line(data = tmp_truth %>% filter(b_degree == 0)) +
  geom_point(data = tmp_truth %>% filter(b_degree == 0)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(n_villages ~ intercept, labeller=label_bquote(
    rows = N == .(n_villages),
    cols = alpha == .(intercept)
  )) +
  theme_bw() +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  guides(col = guide_legend(nrow = 2)) +
  xlab(expression('spillover effect '*beta)) + ylab('treatment effect estimate') +
  #ggtitle('Treatment effects and estimates') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_tau_less

p_tau_all = df_with_truth %>%
  filter(estimator != "truth") %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(
    x = b_spill, y = value,
    group = estimator, colour = estimator,
    shape = estimator
  )) +
  geom_line(data = tmp_truth) +
  geom_point(data = tmp_truth) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(
    n_villages + b_degree ~ prob_assign_nom + intercept,
    labeller = label_bquote(
      rows = list(N == .(n_villages), gamma == .(b_degree)),
      cols = list(rho == .(prob_assign_nom), alpha == .(intercept))
      )
  ) +
  theme_bw() +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  guides(col = guide_legend(nrow = 2)) +
  xlab(expression('spillover effect '*beta)) + ylab('treatment effect estimate') +
  #ggtitle('Treatment effects and estimates') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_tau_all

p_true_tau_less = grid.arrange(p_truth, p_tau_less, nrow=1)
p_true_tau_less
ggsave('../figures/sim-true-tau-simple.pdf', p_true_tau_less, width = 9, height = 4)

p_true_tau = grid.arrange(p_truth_all, p_tau_all, nrow=1, widths = c(1, 2))
p_true_tau
ggsave('../figures/sim-true-tau.pdf', p_true_tau, width = 12, height = 8)


####
# Error: bias, SE, RMSE

df_basic = df %>%
  filter(prob_assign_nom == 0.5, n_villages == 150, b_degree == 0)

# plot bias
p_bias = df_basic %>%
  mutate(
    estimator = estimator_name_map[estimator]
  ) %>%
  ggplot(aes(b_spill, bias, group=estimator, colour=estimator, shape=estimator)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(n_villages + prob_assign_nom ~ intercept, labeller=label_bquote(
    cols=alpha == .(intercept),
    rows = list(N == .(n_villages), rho == .(prob_assign_nom))
  )) +
  #scale_y_continuous(limits=c(-0.02, 0.02)) + 
  theme_bw() +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  xlab(expression('spillover effect '*beta)) + ylab('bias') +
  #ggtitle('Bias') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_bias

# plot standard error
p_se = df_basic %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(b_spill, sd, group=estimator, colour=estimator, shape=estimator)) +
  geom_point() + geom_line() +
  facet_grid(n_villages + prob_assign_nom ~ intercept, labeller=label_bquote(
    cols = alpha == .(intercept),
    rows = list(N == .(n_villages), rho == .(prob_assign_nom))
  )) +
  theme_bw() +
  geom_hline(yintercept=0, linetype='dashed') + 
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  xlab(expression('spillover effect '*beta)) + ylab('true standard error') +
  #ggtitle('Standard error') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_se

# plot RMSE
p_rmse = df %>%
  filter(prob_assign_nom == 0.5, b_degree == 0) %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(b_spill, rmse, group=estimator, colour=estimator, shape=estimator)) +
  geom_point() + geom_line() +
  facet_grid(n_villages + prob_assign_nom ~ intercept, labeller=label_bquote(
    cols = alpha == .(intercept),
    rows = list(N == .(n_villages), rho == .(prob_assign_nom))
  )) +
  theme_bw() +
  coord_cartesian(ylim=c(0, 0.08)) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  xlab(expression('spillover effect '*beta)) + ylab('RMSE') +
  #ggtitle('Estimation error') + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_rmse


p_bias_se = grid.arrange(p_bias, p_se, nrow=1)
p_bias_se
ggsave('../figures/sim-bias-se.pdf', p_bias_se, width=12, height=4)

p_bias_rmse = grid.arrange(p_bias, p_rmse, nrow=1)
p_bias_rmse
ggsave('../figures/sim-bias-rmse.pdf', p_bias_rmse, width=12, height=4)



p_bias_all = p_bias %+% (df %>% mutate(
    estimator = estimator_name_map[estimator]
  )) +
  facet_grid(
    n_villages + b_degree ~ prob_assign_nom + intercept,
    labeller = label_bquote(
      rows = list(N == .(n_villages), gamma == .(b_degree)),
      cols = list(rho == .(prob_assign_nom), alpha == .(intercept))
    )
  )
p_bias_all

p_bias_all = df %>%
  mutate(
    estimator = estimator_name_map[estimator],
    est_design = sprintf("%s (%s)", estimator, prob_assign_nom)
  ) %>%
  filter(!is.na(tau_hat)) %>%
  ggplot(aes(b_spill, bias, colour=est_design, shape=est_design,
             group=est_design)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=0, linetype='dashed') + 
  facet_grid(
    n_villages + b_degree ~ intercept,
    labeller = label_bquote(
      rows = list(N == .(n_villages), gamma == .(b_degree)),
      cols = list(alpha == .(intercept))
    )
  ) +
  theme_bw() +
  scale_colour_manual(name = expression(paste("estimator (", rho, ")")), values=palette_off_policy) +
  scale_shape_manual(name = expression(paste("estimator (", rho, ")")), values=shapes_off_policy) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  xlab(expression('spillover effect '*beta)) + ylab('bias') +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 3))

p_bias_all
ggsave('../figures/sim-bias-all.pdf', p_bias_all, width=8, height=8)

p_rmse_all = df %>%
  mutate(
    estimator = estimator_name_map[estimator],
    est_design = sprintf("%s (%s)", estimator, prob_assign_nom)
  ) %>%
  filter(!is.na(tau_hat)) %>%
  ggplot(aes(b_spill, rmse * sqrt(n_villages), colour=est_design, shape=est_design,
             group=est_design)) +
  geom_point() + geom_line() +
  facet_grid(
    n_villages + b_degree ~ intercept,
    labeller = label_bquote(
      rows = list(N == .(n_villages), gamma == .(b_degree)),
      cols = list(alpha == .(intercept))
    )
  ) +
  theme_bw() +
  scale_colour_manual(name = expression(paste("estimator (", rho, ")")), values=palette_off_policy) +
  scale_shape_manual(name = expression(paste("estimator (", rho, ")")), values=shapes_off_policy) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  coord_cartesian(ylim=c(0, .7)) +
  xlab(expression('spillover effect '*beta)) + ylab(expression(RMSE %*% sqrt(N))) +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 3))

p_rmse_all
ggsave('../figures/sim-rmse-all.pdf', p_rmse_all, width=8, height=8)

p_rmse_op = p_rmse %+% (
  df %>%
    filter(prob_assign_nom == 0, b_degree == 0, estimator != "dm") %>%
    mutate(estimator = estimator_name_map[estimator])
) +
  geom_point(aes(fill = estimator)) +
  scale_fill_manual(values=palette1) +
  coord_cartesian(ylim=c(0, 0.17))

p_rmse_op
ggsave('../figures/sim-rmse-op.pdf', p_rmse_op, width=12, height=4)

###
# Coverage rates
nominal = 0.9
nominal.null.95 = nominal - 1.96 * sqrt(nominal * (1 - nominal) / max(results$rep))
alpha = 1 - nominal

df_coverage = results_long %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps,
    estimator
  ) %>%
  summarise(
    tau_se = sd(tau) / sqrt(n()),
    tau = mean(tau),
    var_est_finite = mean(is.finite(var_est)),
    se_est_mean = mean(sqrt(var_est)),
    clt = .coverage(est, var_est, tau, alpha = alpha),
    sandwich = .coverage(est, var_est_sandwich, tau, alpha = alpha)#,
    #boot = .coverage(est, var_est_bootstrap, tau, alpha = alpha),
    #boot_percentile = mean(tau <= boot_ci_high & tau >= boot_ci_low),
    #perm = .coverage(est, var_est_perm, tau, alpha = alpha),
    #perm_alt_power = mean(perm_p_value < alpha)
  )

df_coverage %>% filter(sandwich <= .88, estimator == "hajek_off_policy") %>% print(n = 100)

df_coverage_null = results_long %>%
  filter(b_degree == 0, b_spill == 0) %>%
  group_by(intercept, b_degree, b_spill, estimator) %>%
  summarise(
    tau = mean(tau),
    clt = .coverage(est, var_est, 0, alpha=1 - nominal),
    sandwich = .coverage(est, var_est_sandwich, 0, alpha=1 - nominal)#,
    #boot = .coverage(est, var_est_bootstrap, 0, alpha=1 - nominal),
    #boot_percentile = mean(0 <= boot_ci_high & 0 >= boot_ci_low)
  )

p_cov = df_coverage %>%
  filter(b_degree == 0, prob_assign_nom == 0.5) %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(aes(b_spill, clt, group=estimator, colour=estimator, shape=estimator)) +
  geom_point() + geom_line() +
  facet_grid(n_villages + prob_assign_nom ~ intercept, labeller=label_bquote(
    cols=alpha == .(intercept),
    rows = list(N == .(n_villages), rho == .(prob_assign_nom))
  )) +
  geom_hline(yintercept=nominal, linetype='dashed') +
  geom_rect(xmin = -1, xmax = 11, ymin=nominal.null.95, ymax = 1, inherit.aes = FALSE, alpha = .01) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  scale_y_continuous() + 
  xlab(expression('spillover effect '*beta)) + ylab('coverage rate') +
  #ggtitle('Coverage rates for 90% nominal confidence interval') + theme_bw() +
  guides(col = guide_legend(nrow = 1)) +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank())
p_cov

ggsave('../figures/sim-coverage-only.pdf', p_cov, width=7, height=5)


p_cov_op = p_cov %+% (
  df_coverage %>%
    filter(prob_assign_nom == 0, b_degree == 0, estimator != "dm") %>%
    mutate(estimator = estimator_name_map[estimator])
) +
  geom_point(aes(fill = estimator)) +
  scale_fill_manual(values=palette1)
p_cov_op

ggsave('../figures/sim-op-coverage-only.pdf', p_cov_op, width=7, height=5)


p_rmse_cov_no_legend = grid.arrange(
    p_rmse + theme(legend.position='none'),
    p_cov + theme(legend.position='none'),
    ncol = 2
    )
g =  ggplotGrob(p_cov+ theme(legend.position="bottom"))$grobs
leg = g[[which(sapply(g, function(x) x$name) == "guide-box")]]
p_rmse_cov = grid.arrange(
  p_rmse_cov_no_legend,
  leg,
  nrow = 2, heights = c(9, 1))
p_rmse_cov

ggsave('../figures/sim-rmse-coverage.pdf', p_rmse_cov, width = 9, height = 5)

p_rmse_cov_op_no_legend = grid.arrange(
    p_rmse_op + theme(legend.position='none'),
    p_cov_op + theme(legend.position='none'),
    ncol = 2
    )
g =  ggplotGrob(p_cov_op + theme(legend.position="bottom"))$grobs
leg = g[[which(sapply(g, function(x) x$name) == "guide-box")]]
p_rmse_cov_op = grid.arrange(
  p_rmse_cov_op_no_legend,
  leg,
  nrow = 2, heights = c(9, 1))
p_rmse_cov_op

ggsave('../figures/sim-op-rmse-coverage.pdf', p_rmse_cov_op, width = 9, height = 5)


###
# power
df_power = results %>% 
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps
  ) %>%
  filter(b_degree > 0 | b_spill > 0) %>%
  summarise(
    tau = mean(tau),
    dm = .coverage(dm, dm_var_est, 0, alpha=1 - nominal),
    ht = .coverage(ht, ht_var_est, 0, alpha=1 - nominal),
    hajek = .coverage(hajek, hajek_var_est, 0, alpha=1 - nominal)
  ) %>%
  gather(estimator, coverage, dm, ht, hajek) %>%
  mutate(power = 1 - coverage) %>% select(-coverage) %>%
  inner_join(
    df_coverage %>% select(coverage = clt)
  )


p_power_basic = df_power %>%
  #filter(coverage >= .85) %>%
  filter(b_degree == 0) %>%
  mutate(
    estimator = estimator_name_map[estimator]
    ) %>%
  ggplot(
    aes(x = b_spill, y = power,
        group = estimator, colour = estimator, shape = estimator
        )) +
  geom_point() + geom_line() +
  facet_grid(n_villages + prob_assign_nom ~ intercept, labeller=label_bquote(
    cols=alpha == .(intercept),
    rows = list(N == .(n_villages), rho == .(prob_assign_nom))
  )) +
  geom_hline(yintercept=1 - nominal, linetype='dashed') +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes1) +
  guides(col = guide_legend(nrow = 1)) +
  xlab(expression('spillover effect '*beta)) + ylab('power') +
  #ggtitle('Power of estimators') +
  theme_bw() + 
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5))
p_power_basic

ggsave('../figures/sim-power-basic.pdf', p_power_basic, width=8, height=8)


p_power_all = df_power %>%
  mutate(
    estimator = estimator_name_map[estimator],
    est_design = sprintf("%s (%s)", estimator, prob_assign_nom)
  ) %>%
  filter(!is.na(power)) %>%
  ggplot(aes(b_spill, power, colour=est_design, shape=est_design,
             group=est_design)) +
  geom_point() + geom_line() +
  facet_grid(
    n_villages + b_degree ~ intercept,
    labeller = label_bquote(
      rows = list(N == .(n_villages), gamma == .(b_degree)),
      cols = list(alpha == .(intercept))
    )
  ) +
  theme_bw() +
  scale_colour_manual(name = expression(paste("estimator (", rho, ")")), values=palette_off_policy) +
  scale_shape_manual(name = expression(paste("estimator (", rho, ")")), values=shapes_off_policy) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  coord_cartesian(ylim=c(0, 1)) +
  xlab(expression('spillover effect '*beta)) + ylab('power') +
  theme(legend.position='bottom', plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(nrow = 3))

p_power_all

ggsave('../figures/sim-power-all.pdf', p_power_all, width=8, height=8)


p_power_ecdf = df_power %>%
  filter(prob_assign_nom == 0.5) %>%
  mutate(
    estimator = estimator_name_map[estimator]
  ) %>%
  ggplot(aes(x = power, group=estimator, colour=estimator)) + 
  stat_ecdf() +
  scale_colour_manual(values=palette1) +
  theme_bw() + theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)) + 
  #scale_x_log10(breaks=10^((-4):(-1)), labels=comma) + 
  geom_hline(yintercept=0, linetype='dashed') + geom_hline(yintercept=1, linetype='dashed') + 
  geom_hline(yintercept=1 - nominal, linetype='dotted') + 
  #scale_y_continuous(limits=c(-0.1, 1.1)) + 
  xlab("power") + ylab("CDF") + ggtitle('Power of estimators')
p_power_ecdf

df_power_diff = df_power %>%
  group_by(
    prob_assign_nom, n_villages, with_replacement, model,
    intercept, b_degree, b_spill, n_steps, tau
  ) %>%
  mutate(
    power_diff = power - power[estimator == 'dm']
  )

ecdf(df_power_diff$power_diff)(.1)

p_power_diff_ecdf = df_power_diff %>%
  filter(prob_assign_nom == 0.5, estimator == "hajek") %>%
  mutate(
    estimator = estimator_name_map[estimator]
  ) %>%
  ggplot(aes(x = power_diff, group=estimator, colour=estimator)) + 
  stat_ecdf() +
  scale_colour_manual(values=palette1) +
  theme_bw() + theme(legend.position='bottom', plot.title = element_text(hjust = 0.5)) + 
  xlab("difference in power\n(vs. difference-in-means)") + ylab("CDF")
p_power_diff_ecdf

p_power = grid.arrange(p_power_all, p_power_ecdf, nrow=1)
p_power

ggsave('../figures/sim-power.pdf', p_power, width=12, height=5)


p_power_diff = grid.arrange(p_power_basic, p_power_diff_ecdf,
                            layout_matrix = matrix(c(1, 1, 2), nrow = 1))
p_power_diff

ggsave('../figures/sim-power-diff.pdf', p_power_diff, width=9, height=6)


p_power_basic_combined_no_legend = grid.arrange(
  p_power_basic + theme(legend.position='none'),
  p_power_all + theme(legend.position='none'),
  p_power_diff_ecdf + theme(legend.position='none'),
  layout_matrix = rbind(c(1, 2), c(1, 3)),
  ncol = 2,
  widths = c(5, 3)
)

g =  ggplotGrob(p_power_basic + theme(legend.position="bottom"))$grobs
leg = g[[which(sapply(g, function(x) x$name) == "guide-box")]]
p_power_basic_combined = grid.arrange(
  p_power_basic_combined_no_legend,
  leg,
  nrow = 2, heights = c(8, 1))
p_power_basic_combined

ggsave('../figures/sim-power-combined.pdf', p_power_basic_combined, width=9, height=6)
