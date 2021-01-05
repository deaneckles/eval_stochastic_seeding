library(igraph)
library(dplyr)
library(tidyr)
library(foreach)
library(RColorBrewer)
library(ggplot2)
theme_set(theme_bw())

el = read.csv("kim_data/village1_edgelist.csv", stringsAsFactors = F)[, 1:3]

el[el[,1] == el[,2],]

g <- graph_from_edgelist(as.matrix(el[, 1:2]), directed = TRUE)
E(g)$relationship_type = el$edge_type

edge_color = c(
  friend = "black",
  sibling = "#880000",
  spouse = brewer.pal(4, "Set1")[4]
)

nd = read.csv("kim_data/village1_nodes.csv", stringsAsFactors = F)[, 1:2]
names(nd) <- c("node_id", "treatment")

V(g)$treatment = nd$treatment[match(V(g)$name, nd$node_id)]

node_color = c(
  multivitamins = brewer.pal(3, "Set1")[1],
  chlorine = brewer.pal(3, "Set1")[3],
  both = brewer.pal(3, "Set1")[2],
  none = "grey"
)

pdf("figures/kim_network_1.pdf", width = 13, height = 13)
par(mai=c(0, 0, 0, 0))
plot.igraph(
  g,
  #vertex.size = 4,
  vertex.size = 4 + 3 * (V(g)$treatment != "none"),
  vertex.label = NA,
  vertex.frame.color = "black",
  vertex.color = node_color[V(g)$treatment],
  edge.arrow.width = 1, edge.arrow.size = .4,
  edge.line.alpha = .5,
  edge.color = edge_color[E(g)$relationship_type],
  edge.curved = FALSE,
  layout = layout_with_kk(g, kkconst = 500),
  margin = rep(0, 4)
)
dev.off()



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


###
# TODO: Compare villages with Table S3

v_in_fig_1 = c(1, 6, 25, 26)

vs[v_in_fig_1, ]

make_kim_village = function(village_num) {
  el = read.csv(
    sprintf("kim_data/village%s_edgelist.csv", village_num),
    stringsAsFactors = F
  )[, 1:3]

  g <- graph_from_edgelist(as.matrix(el[, 1:2]), directed = TRUE)
  E(g)$relationship_type = el$edge_type

  nd = read.csv(
    sprintf("kim_data/village%s_nodes.csv", village_num),
    stringsAsFactors = F
  )[, 1:2]
  names(nd) <- c("node_id", "treatment")
  nd = nd %>%
    mutate(
      multivitamins_seed = treatment %in% c("multivitamins", "both"),
      chlorine_seed = treatment %in% c("chlorine", "both")
    ) %>%
    filter(node_id != "")

  m = match(nd$node_id, V(g)$name)

  nd$indegree = degree(g, mode = "in")[m]

  list(
    village_num = village_num,
    g = g,
    el = el,
    df = nd
    )
}

v1 = make_kim_village(1)

villages = foreach(vi = v_in_fig_1) %do% make_kim_village(vi)
names(villages) = as.character(v_in_fig_1)


seed_summary = foreach(v = villages, .combine = rbind) %do% {
  multivitamins_seed_summary = v$df %>%
    filter(multivitamins_seed) %>%
    summarise(
      n_seeds = n(),
      indegree_mean = mean(indegree),
      indegree_sd = sd(indegree)
    ) %>%
    mutate(set = "multivitamins")
  chlorine_seed_summary = v$df %>%
    filter(chlorine_seed) %>%
    summarise(
      n_seeds = n(),
      indegree_mean = mean(indegree),
      indegree_sd = sd(indegree)
    ) %>%
    mutate(set = "chlorine")
  result = rbind(
    multivitamins_seed_summary,
    chlorine_seed_summary
  ) %>%
    mutate(
      village_num = v$village_num
      )
}


table_to_match = vs[v_in_fig_1, ] %>%
  mutate(
    multivitamin = tolower(multivitamin),
    chlorine = tolower(chlorine)
  )

tmp = vs[v_in_fig_1, ] %>%
  gather(variable, value, contains("indegree")) %>%
  mutate(
    strategy = ifelse(
      grepl("rand", variable), "rand",
      ifelse(
        grepl("top", variable), "indeg", "nom"
      )),
    statistic = ifelse(
      grepl("mean", variable), "indegree_mean", "indegree_sd"
    )
  ) %>%
  select(-variable)
  
  
table_to_match_2 = table_to_match %>%
  select(village, multivitamins = multivitamin, chlorine) %>%
  gather(
    behavior, strategy, multivitamins, chlorine
  ) %>%
  inner_join(
    tmp %>% select(village, strategy, statistic, value)
  ) %>%
  group_by(village, strategy) %>%
  spread(
    statistic, value
  )

seed_summary_2 = seed_summary %>%
  select(
    village = village_num,
    behavior = set,
    indegree_mean, indegree_sd, n_seeds
  ) %>%
  arrange(village, behavior)

output = table_to_match_2 %>%
  ungroup() %>%
  select(-strategy) %>%
  inner_join(
    seed_summary_2,
    by = c("village", "behavior"),
    suffix = c("_table_s3", "_ours")
  )

output


write.csv(output, file = "kim_table_to_match.csv", row.names = FALSE)


output %>%
  mutate(
    indegree_total_table_s3 = indegree_mean_table_s3 * n_seeds,
    indegree_total_ours = indegree_mean_ours * n_seeds
  ) %>%
  select(village, behavior, indegree_total_table_s3, indegree_total_ours)
