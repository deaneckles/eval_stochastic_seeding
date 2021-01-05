.compute_p_rand = function(g, ss) {
    1 / choose(vcount(g), length(ss))
}

.compute_p_nom_unnormalized = function(g, ss) {
  adj <- as_adj(g)
  out_degree <- degree(g, mode = "out")
  valid_start <- out_degree > 0
  adj <- adj[valid_start, , drop = FALSE]
  out_degree <- out_degree[valid_start]
  one_seed_p = colSums(adj * (1 / out_degree)) / nrow(adj)
  p_nom = prod(one_seed_p[as.character(ss)])
  factorial(length(ss)) * p_nom
}

.compute_p_nom = function(g, ss, pi) {
    .compute_p_nom_unnormalized(g, ss) / pi
}


create_village_for_analysis = function(village_name) {
    nodes = cai_survey %>% 
        filter(address == village_name) %>% 
        mutate(
            trt = delay == 0 & intensive == 1,            
        ) %>%
        select(id, delay, intensive, trt, takeup_survey, stratum)
    edges = cai_all %>% 
        filter(
            address == village_name, network_address == village_name,
            id %in% nodes$id, network_id %in% nodes$id
        )
    nodes = nodes %>% filter(id %in% edges$id | id %in% edges$network_id)
    graph = edges %>% 
        select(id, network_id) %>% 
        as.matrix %>%
        apply(2, as.character) %>% 
        graph_from_edgelist(directed=TRUE)
    
    seed_set = with(nodes, id[trt == 1])
    structure(
        list(
            name=village_name,
            nodes=nodes,
            edges=edges,
            graph=graph,
            n_nodes=vcount(graph),
            n_edges=ecount(graph),
            seed_set=seed_set,
            ss_size=length(seed_set),
            y=mean(nodes$takeup_survey),
            p_rand=.compute_p_rand(graph, seed_set),
            p_nom_unnormalized=.compute_p_nom_unnormalized(graph, seed_set)
        ),
        class = 'village'
    )
}
