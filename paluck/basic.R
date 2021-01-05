# basic analysis of seed sets in Paluck et al.

library(foreach)
library(igraph)
require(sandwich)
require(lmtest)
library(ri)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(ggrepel)
library(xtable)
library(doMC)
registerDoMC(cores = 38)
require(estimatr)

source("../simulations/village.R")
source("../simulations/estimators.R")
source("../simulations/analyze.R")
source("schools.R")
source("perm.R")
source("data_paths_default.R")
source("data_paths.R")
path_to_data = path_to_data_extract

# load school-level data, main Paluck data, and adjacency matrix
sld = readRDS(paste0(path_to_data, "school_level_data_some_cols.RDS"))
pd = readRDS(paste0(path_to_data, "student_level_data_some_cols.RDS"))
adj = readRDS(paste0(path_to_data, "adjacency_matrix.RDS"))


# create school objects, which include probabilities
schools = foreach(
  i = sld$schid,
  .final = function(x) setNames(x, sld$schid)
) %dopar% {
  create_school(i)
}

# extract probabilities for observed seed sets, merge with data
schools_p_obs = foreach(s = schools, .combine = rbind) %do% {
  data.frame(
    schid = s$schid,
    p_nom = s$p_nom,
    p_rand = s$p_rand,
    n_edges = ecount(s$g),
    mean_indegree=mean(degree(s$g, mode = "in")),
    std_indegree=sd(degree(s$g, mode = "in"))
  )
}

sldp = inner_join(sld, schools_p_obs)

# subset to treated and control units, choose outcome
sldp <- sldp %>% mutate(
    y = pconf,
    p_design = p_rand
  )

sldp_treated = sldp %>%
  filter(treated == 1) %>%
  arrange(school_block) %>%
  mutate(
    foldid = school_block + c(0, 14)
  )

sldp_treated <- sldp_treated %>%
  mutate(
    w_nom = n() * (p_nom / p_design) / sum(p_nom / p_design),
    w_rand = n() * (p_rand / p_design) / sum(p_rand / p_design),
    w_diff = w_nom - w_rand
  )

sldp_control = sldp %>%
  filter(treated == 0) %>%
  arrange(school_block) %>%
  mutate(
    foldid = school_block + c(0, 14)
    )

# descriptives
desc_summary = sldp_treated %>% 
  mutate(
    percent_treated = 100 * n_treated / n,
    `$\\log_{10} {m_i \\choose k_i}$` = log(choose(n, n_treated)),
    pconf_100 = 100 * pconf
  ) %>%
  select(
    Edges = n_edges,
    `Nodes $m_i$` = n,
    `Treated count $k_i$` = n_treated,
    `Treated \\%` = percent_treated,
    `Peer conflict rate ($\\times 100$)` = pconf_100,
    `$\\log_{10} {m_i \\choose k_i}$`,
    `Mean in-degree` = mean_indegree,
    `St. dev. in-degree` = std_indegree,
  ) %>%
  gather %>%
  group_by(key) %>% 
  summarise_all(funs(mean, sd, min, max))
desc_summary = desc_summary[c(2, 4, 3, 6, 7, 8, 1, 5), ]
desc_summary

desc_summary %>%
  xtable(digits = 1) %>%
  print.xtable(
    file = "paluck_descriptives_table.tex",
    include.rownames = FALSE,
    sanitize.text.function = function(x) x)

table(sldp_treated$p_nom == 0)

# probs for each school's seed sets
p_probs_ratio = ggplot(
  aes(
    x = p_rand, y = p_nom / p_rand,
    size = n_treated
  ),
  data = sldp_treated
) +
  geom_point(alpha=0.3) +
  scale_y_log10(breaks = 10^(-2:1)) +
  scale_x_log10() +
  xlab(expression(p[B]* "  (random)")) +
  ylab(expression(p[A] / p[B] * "  (one-hop to random ratio)")) +
  labs(size = 'seed set\n size (k)') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position= c(0.15, 0.87),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=10),
        legend.background = element_rect(fill = "transparent")
        ) + 
  scale_size_area(breaks = c(20, 32), limits = c(1, 32)) +
  geom_abline(slope=0, intercept = log10(1), linetype='dashed')
p_probs_ratio

ggsave(
  "../figures/paluck_probs_ratio_scatterplot.pdf",
  p_probs_ratio,
  width = 3, height = 4.5
  )

# visualize analysis from Paluck et al of composition of seed set
p_pconf_referent = ggplot(
  aes(x = 100 * frac_treat_referent, y = pconf,
      label = schid),
  data = sldp_treated
) +
  geom_point(alpha = .3) +
  geom_text_repel(size = 4, point.padding = .001, segment.alpha = 0) +
  geom_smooth(method = "lm", se = F) +
  xlab("% of seeds who are social referents") +
  ylab("Peer conflict reports per student")
p_pconf_referent

ggsave(
  "../figures/paluck_school_pconf_by_referent.pdf",
  p_pconf_referent,
  width = 4.5, height = 4
  )

# outcomes plotted by p_nom/p_rand ratio
ggplot(
  aes(x = p_nom/p_rand, y = pconf,
      label = schid),
  data = sldp_treated
) +
  scale_x_log10(lim = c(5e-4, 10), na.value = .2) +
  geom_point(alpha = .3) +
  geom_text_repel(size = 4, point.padding = .001, segment.alpha = 0) +
  xlab(expression(p[A] / p[B] * "  (one-hop to random ratio)")) +
  ylab("Peer conflict reports per student")

ggplot(
  aes(x = 100 * frac_treat_referent,
      y = w_nom - w_rand,
      label = schid,
      size = n_treated
      ),
  data = sldp_treated
) +
  ylab(expression(w[A] - w[B])) +
  xlab("% of seeds who are social referents") +
  geom_point(alpha = .3) +
  geom_text_repel(size = 4, point.padding = .001, segment.alpha = 0)

sldp_treated %>% filter(schid == 35) %>% select(frac_treat_referent, w_diff) %>% round(2)

# compare weights with social referent fraction
p_w_diff_referent = ggplot(
  aes(x = 100 * frac_treat_referent,
      y = w_diff,
      label = schid,
      color = pconf,
      size = n_treated
      ),
  data = sldp_treated %>% filter(schid != 35)
) +
  geom_hline(yintercept = 0, alpha = .7, lty = 2) +
  geom_point(alpha = .5) +
  ylab(expression(tilde(w)[i]^A - tilde(w)[i]^B)) +
  xlab("% of seeds who are social referents") +
  scale_color_distiller(palette = "Spectral", breaks = c(0, .2, .4)) +
  scale_y_continuous(breaks = seq(-1, 1, by = .25)) +
  labs(color = "Peer conflict rate  ") +
  theme_bw() +
  guides(color = guide_colorbar(barwidth = 5, barheight = .4, alpha = .5)) +
  theme(legend.position= "top",        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_text(size=10),
        legend.background = element_rect(fill = "transparent")
        ) + 
  scale_size_area(guide = FALSE, limits = c(1, 32))
p_w_diff_referent

ggsave(
  "../figures/paluck_school_weight_diff_by_referent.pdf",
  p_w_diff_referent,
  width = 3, height = 4
)


estimand_names = c(w_nom = "one-hop", w_rand = "random")

estimates = sldp_treated %>%
  select(
    y = pconf,
    w_nom, w_rand
  ) %>%
  gather(
    key = "estimand", value = "weight",
    w_nom, w_rand
  ) %>%
  mutate(
    estimand = estimand_names[estimand]
    ) %>%
  group_by(estimand) %>%
  summarise(
    estimate = weighted.mean(y, weight)
  )

p_outcome_by_weights = sldp_treated %>%
  gather(estimand, weight, w_nom, w_rand) %>%
  mutate(
    estimand = estimand_names[estimand]
  ) %>%
  ggplot(
    aes(pconf, weight,# * 28,
        group = estimand,
        colour = estimand,
        shape = estimand)
  ) +
  geom_point(size = 2.5, alpha = .8) +
  geom_vline(
    data = estimates,
    aes(xintercept = estimate, colour = estimand),
    linetype='dashed'
  ) +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes2) +
  scale_x_continuous(breaks = c(0, .25, .5, .75)) +
  scale_y_log10(
    breaks = c(10, 1, .1, .01, .001)
  ) +
  xlab("Peer conflict rate") +
  #ylab(expression(weight %*% 28)) +
  theme(
    legend.position=c(.8, .2),
    panel.grid.minor = element_blank()
  )
p_outcome_by_weights

ggsave(
  '../figures/paluck_weights_by_pconf_hajek.pdf',
  p_outcome_by_weights,
  width = 3.5, height = 4)

# flipped
p_outcome_by_weights = sldp_treated %>%
  gather(estimand, weight, w_nom, w_rand) %>%
  mutate(
    estimand = estimand_names[estimand]
  ) %>%
  ggplot(
    aes(y = pconf, x = weight,# * 28,
        group = estimand,
        colour = estimand,
        shape = estimand)
  ) +
  geom_point(size = 2.5, alpha = .8) +
  geom_hline(
    data = estimates,
    aes(yintercept = estimate, colour = estimand),
    linetype='dashed'
  ) +
  scale_colour_manual(values=palette1) +
  scale_shape_manual(values=shapes2) +
  scale_y_continuous(breaks = c(0, .25, .5, .75)) +
  scale_x_log10(
    breaks = c(10, 1, .1, .01, .001)
  ) +
  ylab("peer conflict rate") +
  xlab("normalized weight") +
  theme(
    legend.position=c(.2, .8),
    panel.grid.minor = element_blank()
  )
p_outcome_by_weights

ggsave(
  '../figures/paluck_weights_by_pconf_hajek_flipped.pdf',
  p_outcome_by_weights,
  width = 3.5, height = 4)

# illustrate differences via example networks

schools[[27]]$schid
schools[[31]]$schid

s29 = schools[[27]]
s31 = schools[[31]]

sldp_treated %>% filter(schid %in% c(29, 31)) %>% select(schid, frac_treat_referent, p_nom, w_nom, w_diff)


png("../figures/paluck_network_29_with_seeds.png", width = 2000, height = 1500, pointsize = 24)
par(mai=c(0, 0, 0, 0))
V(s29$g)$color = "#888888"
V(s29$g)$color[V(s29$g)$treat == 1] = "red"
V(s29$g)$color[V(s29$g)$treat == 2] = "black"
plot.igraph(
  s29$g,
  vertex.size = sqrt(compute_p_nom_one_seed(s29$g)) * 30, # size sets diameter
  vertex.label = NA,
  vertex.frame.color = "white",
  vertex.shape = ifelse(V(s29$g)$referent, "square", "circle"),
  edge.arrow.width = .5, edge.arrow.size = .5,
  edge.alpha = .5,
  layout = layout_with_graphopt(s29$g, charge = 0.005,
                                spring.length = 1, spring.constant = 2),
  asp = .5
)
dev.off()

png("../figures/paluck_network_31_with_seeds.png", width = 2000, height = 1500, pointsize = 24)
par(mai=c(0, 0, 0, 0))
V(s31$g)$color = "#888888"
V(s31$g)$color[V(s31$g)$treat == 1] = "red"
V(s31$g)$color[V(s31$g)$treat == 2] = "black"
plot.igraph(
  s31$g,
  vertex.size = sqrt(compute_p_nom_one_seed(s31$g)) * 30,
  vertex.label = NA,
  vertex.frame.color = "white",
  vertex.shape = ifelse(V(s31$g)$referent, "square", "circle"),
  edge.arrow.width = .5, edge.arrow.size = .3,
  edge.line.alpha = .5,
  layout = layout_with_graphopt(s31$g, charge = 0.005,
                                spring.length = 1, spring.constant = 2),
  asp = .7
)
dev.off()


###
# inference about nomination strategy

get_results = function(df, y_var = 'pconf') {
  df$y = df[[y_var]]
  
  set.seed(2818)
  boot_1 = hajek_bootstrap_variance_estimate(
    df,
    alpha = .05, R = 1000
  )
  boot_1
  
  results = data.frame(
    estimator = "hajek",
    est = hajek(df),
    var = hajek_variance_estimate(df),
    var_boot = boot_1$variance,
  ci_lower_boot = boot_1$low,
  ci_upper_boot = boot_1$high
  ) %>%
    mutate(
    se = sqrt(var),
    se_boot = sqrt(var_boot), 
    p_value = 2 * (1 - pnorm(abs(est/se))),
    ci_lower = est - 1.96 * se,
    ci_upper = est + 1.96 * se
    )
  results
}

results = get_results(sldp_treated)
results


###
# Fisherian randomization inference

# generate permutations: this can take a while for large R
p_nom_perms_1 = try(
{
  readRDS("p_nom_perms.RDS")
},
error = function(x) {
  perms = generate_permutations_paluck(sldp_treated, schools, R = 1e4)
  saveRDS(perms, "p_nom_perms.RDS")
  return(perms)
}

p1 = hajek_permutation_paluck(
  sldp_treated, schools,
  p_nom_perm = p_nom_perms_1,
  keep_perms = FALSE,
  keep_perm_results = TRUE
)
p1$p_value

write_results_table = function(results, perm_results, file) {
  results = results %>%
    mutate(
      p_value_perm = perm_results$p_value
    )

  results_table = results %>%
    mutate(
      `95\\% CI (analytic)` = sprintf("[%0.4f, %0.4f]", ci_lower, ci_upper),
      `95\\% CI (bootstrap)` = sprintf("[%0.4f, %0.4f]", ci_lower_boot, ci_upper_boot),
    ) %>%
    select(
      `estimate  (one-hop $-$ rand)` = est,
      `SE (analytic)` = se,
      `SE (bootstrap)` = se_boot,
      `95\\% CI (analytic)`,
      `95\\% CI (bootstrap)`,
      `p-value (analytic)` = p_value,
      `p-value (Fisherian)` = p_value_perm
    ) %>%
    format(
      digits=2, width = 4, justify = "right",
      nsmall = 4,
      zero.print = TRUE
    ) %>%
    t() %>%
    xtable(digits=8)
  
  results_table
  
  results_table %>% 
    print.xtable(
      file = file,
      sanitize.text.function = function(x) x
    )
}

write_results_table(results, p1, file = "paluck_results_table.tex")


p_perm_t_hist = ggplot(
  aes(x = t, y = ..density..),
  data = p1$perm_df
) +
  geom_histogram(bins = 50) +
  geom_vline(
    xintercept = p1$t_naive,
    color = "red"
  ) +
  xlab("Studentized estimate") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
p_perm_t_hist

ggsave(
  "../figures/paluck_perm_t_hist.pdf",
  p_perm_t_hist,
  width = 4, height = 2.5
)


# fraction of students at out-degree cap
sum(sapply(schools, function(v) sum(degree(v$g, mode = "out") == 10))) / sum(sapply(schools, function(v) vcount(v$g)))
