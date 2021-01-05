
library(foreign)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(igraph)
library(xtable)
library(Matrix)
library(RColorBrewer)
library(ggplot2)
theme_set(theme_bw())
theme_update(
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

source('village.R')
source('runner.R')
source('analyze.R')
source('estimators.R')


cai_all = read.dta('../cai-data/data/0422allinforawnet.dta')
village_names = unique(cai_all$address)
village_names = village_names[village_names != '']

n_edges = foreach(v = village_names) %dopar% {
  edges = cai_all %>% filter(
    address == v, network_address == v,
    !is.na(network_id),
    id != 99, network_id != 99
  )
  nrow(edges)
}
village_names = village_names[n_edges > 25]
n = length(village_names)

villages = foreach(i = 1:n) %do% create_village_cai(cai_all, village_names[i])

#vi = 80
vi = 115
village = villages[[vi]]
g1 = simplify(village$g)


l1 = layout_with_kk(g1, kkconst = 100)

pdf(sprintf("../figures/example_village_%s_network.pdf", vi), width = 5, height = 7)

V(g1)$color = "black"
par(mar=c(0,0,0,0))
plot(
    g1,
    vertex.size = 7,
    vertex.label = NA,
    vertex.frame.color = NA,
    edge.width = 1.4,
    edge.arrow.size = .6,
    edge.arrow.width = .8,
    layout = l1,
    asp = 7/5
)

dev.off()

tmp <- as.matrix(as_adj(g1))
tmpo <- order(-degree(g1, mode = "in"))
tmp <- tmp[tmpo, tmpo]
tmp <- tmp[, ncol(tmp):1]

squash::savemat(
  x = tmp,
  filename = sprintf("../figures/example_village_%s_network_adj.png", vi),
  dev = 'png',
  map = list(
    breaks = c(0, .5, 2, 10),
    colors = c("black", "white", "black", "black"),
    right = FALSE, include.lowest = FALSE, col.na = NA, base = NA)
)

#set.seed(4217)
set.seed(3001)

#example_ss_i_top = which(rank(-village$p_nom) < 4)
#example_ss_i = c(example_ss_i_top, sample.int(nrow(village$seed_sets), 24))
#example_ss_i = order(-village$p_nom)
example_ss_i = order(-village$p_nom)
example_ss = village$seed_sets[example_ss_i, ]
example_ss_p_nom = village$p_nom[example_ss_i]
seed_str = apply(example_ss, 1, function(x) paste(sort(x), collapse = ", "))

#highlight_i = c(1, length(example_ss_i_top) + 1:2)
highlight_i = c(1, 50, 300)


c1 = rep("#dddddd", length(seed_str))
c1[highlight_i] = brewer.pal(3, "Set1")
names(c1) = seed_str

pdf(sprintf("../figures/example_village_%s_network_colored_seeds.pdf", vi), width = 5, height = 7)

V(g1)$color = "black"
V(g1)[example_ss[highlight_i[1], ]]$color = c1[seed_str[highlight_i[1]]]
V(g1)[example_ss[highlight_i[2], ]]$color = c1[seed_str[highlight_i[2]]]
V(g1)[example_ss[highlight_i[3], ]]$color = c1[seed_str[highlight_i[3]]]

par(mar=c(0,0,0,0))
plot(
    g1,
    vertex.size = 7,
    vertex.label = NA,
    vertex.frame.color = NA,
    edge.width = 1.4,
    edge.arrow.size = .6,
    edge.arrow.width = .8,
    layout = l1,
    asp = 7/5
)

dev.off()


example_ss_df = data.frame(
  seed_str = seed_str,
  p_nom = example_ss_p_nom
) %>%
  filter(!duplicated(seed_str))

nrow(example_ss_df) - 6 - 2

#par(mar=c(50.1, 4.1, 4.1, 2.1))
par(mar=c(0,0,0,0))
p1 = ggplot(
  data = example_ss_df, 
  aes(
    x = reorder(seed_str, p_nom),
    y = 1e3 * p_nom, fill = seed_str
  )) +
  coord_flip() +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 1000 * village$p_rand, linetype = 2, alpha = .4) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +#seq(0, 6, 2)) +
  scale_fill_manual(values = c1) +
  xlab(NULL) + ylab(expression(1000 %*% p[A])) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.y=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks.y=element_blank()
        )
p1

ggsave(sprintf('../figures/example_village_%s_seeds_p_nom.pdf', vi), p1, width=2.2, height=5)

