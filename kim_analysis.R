library(igraph)
library(dplyr)
library(tidyr)
library(foreach)
library(RColorBrewer)
library(ggplot2)
theme_set(theme_bw())

####
# Table S3

vs = read.csv("kim_data/Kim Table S3.csv")
names(vs) = tolower(names(vs))
vs = vs %>%
  mutate_at(
    vars(indegree_rand, indegree_nom, indegree_top),
    gsub,
    pattern = "(\\(|\\))", replacement = ""
    ) %>%
  separate(
    indegree_rand,
    sep = " ",
    c("indegree_rand_mean", "indegree_rand_sd"),
    extra = "drop", fill = "right", convert = TRUE
  ) %>%
  separate(
    indegree_nom,
    sep = " ",
    c("indegree_nom_mean", "indegree_nom_sd"),
    extra = "drop", fill = "right", convert = TRUE
  ) %>%
  separate(
    indegree_top,
    sep = " ",
    c("indegree_top_mean", "indegree_top_sd"),
    extra = "drop", fill = "right", convert = TRUE
  )

p_obs_seed_sets_rand_vs_nom = ggplot(
  aes(x = indegree_rand_mean,
      y = indegree_nom_mean,
      size = n_seeds),
  data = vs %>%
    filter(!is.na(indegree_rand_mean), !is.na(indegree_nom_mean))
) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  geom_point(alpha = .5) +
  scale_size_area(breaks = c(1, 5, 9)) +
  coord_fixed(xlim = c(1, 10), ylim = c(1, 10)) +
  xlab("random seed set (mean in-degree)") +
  ylab("one-hop seed set (mean in-degree)") +
  labs(size = 'seed set\n size (k)') +
  theme(
    legend.position= c(0.88, 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_text(size=10),
    legend.background = element_rect(fill = "transparent")
  )
p_obs_seed_sets_rand_vs_nom

ggsave(
  "figures/kim_mean_indegree_rand_vs_nom.pdf",
  p_obs_seed_sets_rand_vs_nom,
  width = 4, height = 4
  )


vs_long = vs %>%
  select(
    -indegree_rand_sd, -indegree_nom_sd, -indegree_top_sd
    ) %>%
  gather(
    key = "strategy", value = "indegree_mean",
    indegree_rand_mean, indegree_nom_mean, indegree_top_mean
  )

ggplot(
  aes(x = reorder(factor(village), n_nodes),
      y = indegree_mean,
      color = strategy,
      shape = strategy,
      group = village
      ),
  data = vs_long %>% filter(!is.na(indegree_mean))
) +
  geom_line(color = "grey") +
  geom_point(size = 2) +
  scale_x_discrete() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_text(size=10),
    legend.background = element_rect(fill = "transparent")
  )
