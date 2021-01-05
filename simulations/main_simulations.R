# Simulations using Cai dataset, with both independent cascade and linear-in-means utility model

library(foreign)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
library(igraph)
library(xtable)
library(Matrix)

source('village.R')
source('runner.R')
source('analyze.R')
source('estimators.R')


# load data and construct village objects
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

# parameters governing simulation
n_reps = 5000

model_params = bind_rows(
  expand.grid(
    model = "LIM",
    intercept = seq(-3, 0, by = 1),
    b_degree = c(0, .1),
    b_spill = seq(0, 10, by = 2),
    n_steps = 3
  ),
  expand.grid(
    model = "IC",
    n_steps = Inf,
    p_spontaneous = 0,
    p_infect = seq(0, 1, by = .1)
  )
)
    
params = expand.grid(
  n_villages = c(50, length(villages)),
  model = c("LIM", "IC"),
  prob_assign_nom = c(0, 0.5)
) %>%
  inner_join(model_params)


# do simulations and write files
start = proc.time()
registerDoParallel(cores = 38)
foreach(i = 1:nrow(params),
        .combine = rbind) %do% {
  param = params[i,]
  print(param)
  results = foreach(rep = 1:n_reps, .combine = rbind) %dopar% {
    with(param, run_experiment(
      n_villages = n_villages,
      villages = villages,
      rep = rep,
      model = model,
      intercept = intercept, b_degree = b_degree,
      b_spill = b_spill,
      n_steps = n_steps,
      p_spontaneous = p_spontaneous,
      p_infect = p_infect,
      with_replacement = FALSE,
      prob_assign_nom = prob_assign_nom
    ))
  } %>% data.frame
  print(proc.time() - start)
  write.table(
    results, file='results/main_results.csv',
    append=TRUE, col.names=FALSE, row.names=FALSE, sep=','
  )
}

write.table(
  names(results),
  file='results/main_results_names.csv',
  col.names=FALSE, row.names=FALSE, sep=','
)

results = read.csv('results/main_results.csv', header=FALSE)
names(results) <- read.csv('results/main_results_names.csv', header = F)[,1]

arrow::write_parquet(
  results,
  "results/main_results.parquet",
  compression = "gzip"
)
