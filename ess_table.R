library(foreach)
library(igraph)
library(network)
library(intergraph)
library(foreign)
library(dplyr)
library(doParallel)
library(xtable)
library(ggplot2)
library(tidyr)

source('simulations/village.R')
source('simulations/ess.R')
source('simulations/analyze.R')


N_cores=8
registerDoParallel(cores=N_cores)

# Some of the parallelization requires adjustments of system settings: 
options(expressions=500000)
# If having issues with C stack limit, check Cstack_info()
# On MacOS the following system-level adjustments may be necessary:
# https://stackoverflow.com/questions/43181989/how-to-set-the-size-of-the-c-stack-in-r

# Compute population ESS for a bunch of networks.

villages = list()

## CAI ET AL. 2015
cai_all = read.dta('cai-data/data/0422allinforawnet.dta')
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

villages$cai = foreach(i = 1:n) %do% create_village_cai(cai_all, village_names[i])


## INDIAN MICROFINANCE VILLAGES
# Takes all 12 relations and flattens them to a single, simple graph.

create_microfinance_village = function(i) {
  data_dir = 'indian-microfinance/datav4.0/Data/1. Network Data/Adjacency Matrices/'
  adjs = foreach(relation_name = relation_names) %do% {
    path = paste(data_dir, 'adj_', relation_name, '_HH_vilno_', i, '.csv', sep='')
    adj = read.csv(path, header=FALSE)
  }
  g = Reduce('+', adjs) %>%
    as.matrix %>%
    graph_from_adjacency_matrix %>% 
    simplify
  Village(g, i, ss_size=2)
}

relation_names = c(
  'borrowmoney',
  'giveadvice',
  'helpdecision',
  'keroricecome',
  'keroricego',
  'lendmoney',
  'medic',
  'nonrel',
  'rel',
  'templecompany',
  'visitcome',
  'visitgo'
)

village_ids = c(1:12, 14:21, 23:77)
n = length(village_ids)
villages$microfinance = foreach(i = village_ids) %dopar% create_microfinance_village(i)


## ADD HEALTH NETWORKS
# Read in the full (grade 9-12 or grade 7-12, where applicable) networks.

create_addhealth_village = function(i) {
  path = paste('add_health_nocontract/structure_nocontract/comm', i, '.paj', sep='')
  g = read.paj(path) %>% (function(l) {l$networks[[paste('comm', i, '.net', sep='')]]}) %>% asIgraph
  Village(g, i, ss_size=2)
}

village_ids = 1:85
start = proc.time()
# takes a while
villages$addhealth = foreach(i = village_ids) %dopar% create_addhealth_village(i)
print(proc.time() - start)

## CHAMI NETWORKS
# Read in the `friendship' networks. Not currently using the `health-advice' networks.

create_chami_village = function(i) {
  path = paste('chami-data/friendship/', i, '', sep='')
  chami_edgelist_i = read.csv(path, header=FALSE)
  g = graph_from_edgelist(as.matrix(chami_edgelist_i,directed=TRUE))
  Village(g, i, ss_size=2)
}
village_ids = 1:17
villages$chami = foreach(i = village_ids) %dopar% create_chami_village(i)

## PALUCK SCHOOLS

source("paluck/data_paths.R")
path_to_data = path_to_data_extract
if (TRUE) {

  # load school-level data, main Paluck data, and adjacency matrix
  sld = readRDS(paste0(path_to_data, "school_level_data_some_cols.RDS"))
  pd = readRDS(paste0(path_to_data, "student_level_data_some_cols.RDS"))
  adj = readRDS(paste0(path_to_data, "adjacency_matrix.RDS"))

  create_paluck_village = function(i) {
    in_i = pd$schid == i
    pd_i = pd[in_i, ]
    adj_i = adj[in_i, in_i]
    g = graph_from_adjacency_matrix(adj_i, mode='directed')
    Village(g, i, ss_size=2)
  }
  village_ids = sld$schid
  villages$paluck = foreach(i = village_ids) %dopar% create_paluck_village(i)
}

# Compute ESS

designs = c('off-policy', 'bernoulli', 'optimal')
ess_table = foreach (dataset = names(villages), .combine=rbind) %do% {
  ess = foreach (design = designs, .combine=c) %do% {
    compute_population_ess(villages[[dataset]], design)
  }
  data.frame(dataset=dataset, 
             design=designs, 
             n = length(villages[[dataset]]), 
             n_dm = 1/(1/(length(villages[[dataset]])/2)+1/(length(villages[[dataset]])/2)), 
             n_eff=ess)
}

ess_table

saveRDS(ess_table, "ess_table.RDS")

# not as in paper
ess_table %>% xtable %>% print.xtable(include.rownames=FALSE)

# make table for paper
ess_table_wide = ess_table %>%
  mutate(
    ess_mult = n / n_dm,
    n_eff_scaled = ess_mult * n_eff
    ) %>%
  pivot_wider(dataset:N, values_from = n_eff, names_from = design) %>%
  select(dataset, n, bernoulli, optimal, `off-policy`) %>%
  arrange(-n)

ess_table_wide$dataset <- dataset_name_map[ess_table_wide$dataset]

ess_table_wide %>% xtable(digits=1) %>% print.xtable(include.rownames=FALSE)
  
sum(ess_table_wide$n)

#### ESS as a function of k, for Cai

cai_all = read.dta('cai-data/data/0422allinforawnet.dta')
village_names = unique(cai_all$address)
village_names = village_names[village_names != '']
n_edges = foreach(v = village_names) %dopar% {
  edges = cai_all %>% filter(address == v, network_address == v)
  nrow(edges)
}
village_names = village_names[n_edges > 25]
n = length(village_names)


villages = list()

start = proc.time()
for (i in 1:15) {
  villages[[paste0('cai', i)]] = foreach(j = 1:n) %dopar%
  create_village_cai(
    cai_all, village_names[j], ss_size=i,
    approx_pi = TRUE,
    R = 1000, n_ss_samples = 1000
  )
}
print(proc.time() - start)

villages_to_ess = names(villages)


designs = c('off-policy', 'bernoulli', 'optimal')
ess_table_k = foreach (dataset = villages_to_ess, .combine=rbind) %do% {
  ess = foreach (design = designs, .combine=c) %do% {
    compute_population_ess(villages[[dataset]], design)
  }
  data.frame(
    dataset=dataset, 
    design=designs, 
    n = length(villages[[dataset]]), 
    n_dm = 1/(1/(length(villages[[dataset]])/2)+1/(length(villages[[dataset]])/2)), 
    n_eff=ess
  )
}


ess_table_k = ess_table_k %>% mutate(ratio = n_eff / 150, deff = 1/ratio, k=rep(1:15, each=3))

design_names = c(
  'off-policy' = 'random (off-policy)',
  'bernoulli' = 'Bernoulli(1/2)',
  'optimal' = 'optimized'
  )

p1 = ess_table_k %>%
  mutate(
    design = design_names[as.character(design)],
    design = forcats::fct_reorder(design, -ratio)
    ) %>%
  ggplot(aes(k, ratio, color = design, shape = design)) +
  geom_point() + geom_line() +
  geom_hline(yintercept=1, linetype='dashed') +
  scale_colour_manual(values = palette1) +
  scale_shape_manual(values = shapes1) +
  ylim(0, 11) +
  theme_bw() + theme(legend.position = c(.75, .8), plot.title = element_text(hjust = 0.5)) +
  xlab('seed set size, k') + ylab('relative efficiency')
p1

p2 = ess_table_k %>% filter(design == 'off-policy') %>%
  ggplot(aes(k, ratio, color=design)) + geom_point() + geom_line() + geom_hline(yintercept=1, linetype='dashed') +
  scale_colour_manual(values = palette1) +
  theme_bw() + theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5)) +
  xlab('seed set size, k') + ylab('relative efficiency') + ggtitle('Relative efficiency for off-policy evaluation')
p2

ggsave('figures/ess_k_ate.pdf', p1, width=4.5, height=3.5)
ggsave('figures/ess_k_off_policy.pdf', p2, width=4, height=4)
