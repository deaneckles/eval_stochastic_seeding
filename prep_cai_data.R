# this reads the Cai data and does things like estimate probabilities pi

options(warn=-1)
library(foreign)
library(dplyr)
library(igraph)
library(doParallel)
library(ggplot2)
library(broom)
library(tidyr)
library(xtable)
library(gridExtra)
library(Matrix)
registerDoParallel(cores = 8)
source("simulations/village.R")
source("simulations/estimators.R")
source("cai_analysis_functions.R")



registerDoParallel(cores=32)

cai_all = read.dta('cai-data/data/0422allinforawnet.dta')
cai_survey = read.dta('cai-data/data/0422survey.dta')
cai_price = read.dta('cai-data/data/0422price.dta')
cai_structure = read.dta('cai-data/data/0422structure_all.dta')

cai_survey = cai_survey %>% 
    select(id, address, hh_size=agpop, rice_area=ricearea_2010,
           takeup_survey, insurance_prediction = insurance_buy, understanding,
           delay, intensive, default) %>% 
    mutate(id = as.character(id)) %>%
    filter(!is.na(hh_size), !is.na(rice_area))

get_strata = function(cai_df) {
    strata_df = cai_df %>% group_by(address) %>% summarise(
        median_hh_size = median(hh_size),
        median_rice_area = median(rice_area)
    )
    with(
        cai_df %>% left_join(strata_df, by='address'),
        2 * (hh_size <= median_hh_size) + 1 * (rice_area <= median_rice_area)
    )
}
cai_survey$stratum = get_strata(cai_survey)

table(cai_survey$stratum)

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

n

villages = foreach(i = 1:n) %dopar% create_village_for_analysis(village_names[i])

required_pi_hat_R = foreach(v = villages, .combine = c) %dopar% {
    se = Inf
    R = 1e4
    while (se > .1) {
        R = R * 10
        se = sqrt(pi_variance_bound_graph(v, v$ss_size, R))
    }
    R
}

pi_hat_with_empirical_se = foreach(v = villages, .combine = rbind) %dopar% {
  v$pi_hat_R = 1e7
  R_each = v$pi_hat_R / 10
  pi_hat_many = replicate(10, approx_pi_dependent(
    g = v$g,
    k = v$ss_size,
    R = R_each
  ))
  
  bound = sqrt(pi_variance_bound_graph(
    v,
    k = v$ss_size,
    R = R_each
  ))
  
  data.frame(
    name = v$name,
    se_bound = bound,
    se = sd(pi_hat_many),
    inv_se = sd(1 / pi_hat_many),
    pi_hat = mean(pi_hat_many)
    )
}

saveRDS(pi_hat_with_empirical_se, "cai_pi_hat_with_empirical_se.RDS")
# pi_hat_with_empirical_se <- readRDS("cai_pi_hat_with_empirical_se.RDS")

pi_hat_with_empirical_se = pi_hat_with_empirical_se %>%
  mutate(
    inv_pi_hat = 1 / pi_hat,
    inv_se_rel = inv_se / inv_pi_hat
  )

registerDoParallel(cores = 39)
pi_hat_sequential = foreach(v = villages, .combine = rbind) %dopar% {
    inv_se_rel = Inf
    R_each = 1e6
    pi_hat_many = c()
    R = 0

    while (R < 5e7 & (is.na(inv_se_rel) | inv_se_rel > .1)) {
      pi_hat_many = c(
        pi_hat_many,
        approx_pi_dependent(
          g = v$g,
          k = v$ss_size,
          R = R_each
        ))
      R = R + R_each
      inv_se = sd(1 / pi_hat_many)
      pi_hat = mean(pi_hat_many)
      inv_se_rel = inv_se / pi_hat
    }
    
    data.frame(
      name = v$name,
      R = R_each * length(pi_hat_many),
      se = sd(pi_hat_many),
      inv_se = sd(1 / pi_hat_many),
      pi_hat = mean(pi_hat_many)
    )
}


saveRDS(pi_hat_with_empirical_se, "cai_pi_hat_with_empirical_se_seq.RDS")
# pi_hat_with_empirical_se <- readRDS("cai_pi_hat_with_empirical_se_seq.RDS")

summary(pi_hat_with_empirical_se)

villages[[1]]$seed_set

villages = foreach(v = villages) %dopar% {
  tmp = pi_hat_with_empirical_se[pi_hat_with_empirical_se$name == v$name, ]
  v$pi_hat = tmp$pi_hat
  v$pi_hat_se = tmp$se
  v$pi_hat_inv_se = tmp$inv_se
  v$p_nom = v$p_nom_unnormalized / v$pi_hat
  v
}

saveRDS(villages, "cai_villages_with_pi_hat.RDS")
