# Computes population ESS for single village.  To compute ESS for all villages take the sum.
.compute_population_ess_for_single_village = function(p_rand, p_nom, design) {
  if (design == 'off-policy') {
    # Condition on positive nomination seeds -- is this the right thing to do?
    p_rand = p_rand[p_nom > 0]
    p_nom = p_nom[p_nom > 0]
    ess = 1 / sum(p_nom^2 / p_rand) 
  } else if (design == 'bernoulli') {
    p_design = 0.5 * (p_nom + p_rand)
    ess = 1 / sum((p_nom - p_rand)^2 / p_design)
  } else if (design == 'optimal') {
    # Avoids problem with 1/sum( ... + 0^2*(1/0) + ) by 
    # first computing inverse 1/0 and then zeroing out infinites
    inv_p_design = abs(p_nom - p_rand) %>% (function(p) sum(p)/p)
    inv_p_design[is.infinite(inv_p_design)]=0
    ess = 1 / sum((p_nom - p_rand)^2 * inv_p_design)
  } else stop() # should not reach here
  
  ess
}

compute_population_ess = function(villages, design) {
  
  if (!(design %in% c('off-policy', 'bernoulli', 'optimal'))) stop('invalid design')
  
  ess = foreach(v = villages, .combine=c) %do% {
    .compute_population_ess_for_single_village(v$p_rand, v$p_nom, design=design)
  }
  if (design == 'off-policy') {
    ess_all = length(villages) / (mean(1 / ess) - 1)
  } else {
    ess_all = length(villages) / mean(1 / ess) 
  }
  ess_all / (0.5 * (1 - 0.5))
}


compute_sample_ess = function(p_rand, p_nom, p_design) {
  w = (p_rand - p_nom) / p_design

  1 / (mean(w^2) / mean(abs(w))^2)

}
