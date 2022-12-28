#############################################################################################################################
## Extract 87 summary statistics from a tree and a current trait distribution
## Creation - jeremy.andreoletti@ens.fr - 10/2019
##
## Used with kind permission of JÃ©rÃ©my AndrÃ©ole
##
## This is not actually implemented included in the MESS output yet
## iao 12/2023
#############################################################################################################################

extract_stats <- function(tree, traits) {
  # gsp_extant = result from a simulation (all/gsp_fossil/gsp_extant)
  # Contains tree and tips
  
  nb_stats <- 87
  if (length(traits) < 3){    # ERROR : Tree too small or extinct
    return (rep(NA, nb_stats))
  }
  
  traits <- traits[match(tree$tip.label, names(traits))]  # Reorder traits according to the tree
  
  ## raw trait distribution
  traits_min <- min(traits)
  traits_max <- max(traits)
  traits_mean <- mean(traits)
  traits_median <- median(traits)
  traits_sd <- sd(traits)
  traits_skewness <- moments::skewness(traits)
  traits_kurtosis <- moments::kurtosis(traits)
  traits_mean_median <- mean(traits)-median(traits)
  
  ## independent contrasts
  IC <- ape::pic(traits,tree)        # Phylogenetically Independent Contrasts
  
  tree_traits_IC_min <- min(IC)
  tree_traits_IC_max <- max(IC)
  tree_traits_IC_mean <- -mean(IC)
  tree_traits_IC_median <- median(IC)
  tree_traits_IC_sd <- sd(IC)
  tree_traits_IC_skewness <- moments::skewness(IC)
  tree_traits_IC_kurtosis <- moments::kurtosis(IC)
  tree_traits_IC_mean_median <- mean(IC)-median(IC)
  
  ## nearest-neighbor distance
  traits_NN_dist <- apply(as.matrix(dist(as.numeric(traits))), 2,
                          function(x){sort(x)[2]})   # Take the closest in each column of the pairwise distance matrix
  
  traits_NN_dist_min <- min(traits_NN_dist)
  traits_NN_dist_max <- max(traits_NN_dist)
  traits_NN_dist_mean <- mean(traits_NN_dist)
  traits_NN_dist_median <- median(traits_NN_dist)
  traits_NN_dist_sd <- sd(traits_NN_dist)
  traits_NN_dist_range <- diff(range(traits_NN_dist))
  traits_NN_dist_skewness <- moments::skewness(traits_NN_dist)
  traits_NN_dist_kurtosis <- moments::kurtosis(traits_NN_dist)
  
  ## phylogenetic signal
  cophenetic_tree <- ape::cophenetic.phylo(tree)    # Cophenetic tree
  weights <- 1/cophenetic_tree
  diag(weights) <- 0                                # Diag : Inf -> 0
  
  tree_Moran_I <- ape::Moran.I(as.numeric(traits), weights, scaled = T)$observed
  tree_Blomberg_K <- phylosig(tree, traits)
  
  ## trait distance matrix
  traits_pair_dist <- as.numeric(dist(traits)) # pairwise distance matrix for traits
  
  traits_pair_dist_min <- min(traits_pair_dist)
  traits_pair_dist_max <- max(traits_pair_dist)
  traits_pair_dist_mean <- mean(traits_pair_dist)
  traits_pair_dist_median <- median(traits_pair_dist)
  traits_pair_dist_sd <- sd(traits_pair_dist)
  traits_pair_dist_range <- diff(range(traits_pair_dist))
  
  traits_pair_dist_eig <- eigen(dist(traits))  # eigen values of the pairwise distance matrix
  
  traits_pair_dist_eig_min <- min(traits_pair_dist_eig$values)
  traits_pair_dist_eig_max <- max(traits_pair_dist_eig$values)
  traits_pair_dist_eig_mean <- mean(traits_pair_dist_eig$values)
  traits_pair_dist_eig_median <- median(traits_pair_dist_eig$values)
  traits_pair_dist_eig_sd <- sd(traits_pair_dist_eig$values)
  traits_pair_dist_eig_skewness <- moments::skewness(traits_pair_dist_eig$values)
  traits_pair_dist_eig_kurtosis  <- moments::kurtosis(traits_pair_dist_eig$values)

  ## dist. between ordered trait values
  traits_ordered_gaps <- diff(sort(as.numeric(traits)))   # gaps between ordered trait values
  
  traits_ordered_gaps_min <- min(traits_ordered_gaps)
  traits_ordered_gaps_max <- max(traits_ordered_gaps)
  traits_ordered_gaps_mean <- mean(traits_ordered_gaps)
  traits_ordered_gaps_median <- median(traits_ordered_gaps)
  traits_ordered_gaps_sd <- sd(traits_ordered_gaps)
  traits_ordered_gaps_range <- diff(range(traits_ordered_gaps))
  
  ## spectral density profile (SDP)
  tree_SDP <- spectR(tree)                             # spectral density profile : tree

  tree_SDP_max <- tree_SDP$principal_eigenvalue
  tree_SDP_asymmetry<- tree_SDP$asymmetry
  tree_SDP_peakedness <- tree_SDP$peakedness
  tree_SDP_eigengap <- tree_SDP$eigengap
  
  tree_SDP_normal <- spectR(tree, meth = "normal")     # normalized SDP : tree

  tree_SDP_normal_max <- tree_SDP_normal$principal_eigenvalue
  tree_SDP_normal_asymmetry<- tree_SDP_normal$asymmetry
  tree_SDP_normal_peakedness <- tree_SDP_normal$peakedness
  tree_SDP_normal_eigengap <- tree_SDP_normal$eigengap
  
  tree_traits_SDP <- spectR_t(tree, traits)            # spectral density profile : tree + traits
  
  tree_traits_SDP_splitter <- tree_traits_SDP$splitter
  tree_traits_SDP_tracer <- tree_traits_SDP$tracer
  tree_traits_SDP_fragmenter <- tree_traits_SDP$fragmenter
  
  ## regression of trait distances ~ sqrt (phylogenetic distances)
  cophenetic_tree_lower <- cophenetic_tree[lower.tri(cophenetic_tree)]   # Phylogenetic distances
  
  reg_trait_phylo_dist <- lm(traits_pair_dist~sqrt(cophenetic_tree_lower))
  
  reg_trait_phylo_dist_intercept <- reg_trait_phylo_dist$coefficients["(Intercept)"]
  reg_trait_phylo_dist_beta <- reg_trait_phylo_dist$coefficients[2]  # Beta

  ## tree metrics
  tree_gamma <- ape::gammaStat(tree)                        # Gamma-Statistic of Pybus and Harvey
  tree_beta <- apTreeshape::maxlik.betasplit(tree)$max_lik  # Beta-splitting model
  
  # tree topology
  tree_size <- length(tree$tip.label)       # number of tips

  topo_summaries <- function(tree){
    leafs <- 1:Ntip(tree)
    df_tree_topo_int <- data.frame(external = leafs, depth = NA)
    
    for (l in leafs){
      # depth (number of ancestors) by leaf
      df_tree_topo_int[l,][,"depth"] <- length(phangorn::Ancestors(tree, l))
      }
    
    internal.nodes = (Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
    df_tree_topo_ext <- data.frame(internal = internal.nodes, children_asym = NA, depth = NA, width = NA)
    
    for (n in internal.nodes){
      children_nodes <- phangorn::Children(tree, n) 
      
      # absolute difference between descendants from the right and the left children nodes
      df_tree_topo_ext[df_tree_topo_ext$internal==n,][,"children_asym"] <- abs(length(phangorn::Descendants(tree, children_nodes[1])) -
                                                                         length(phangorn::Descendants(tree, children_nodes[2])))
      
      # depth = number of Ancestors
      depths_table <- table(df_tree_topo_ext$depth)
      df_tree_topo_ext[df_tree_topo_ext$internal==n,][,"depth"] <- length(phangorn::Ancestors(tree, n))
      }
    
    # width = number of nodes at the same depth
    df_tree_topo_ext[,"width"] <- depths_table[paste(df_tree_topo_ext$depth)]
    
    return (list(df_tree_topo_int, df_tree_topo_ext, depths_table))
    }
  
  L = topo_summaries(tree)
  df_tree_topo_int <- data.frame(L[1])
  df_tree_topo_ext <- data.frame(L[2])
  depths_table <- table(df_tree_topo_ext$depth)
  
  tree_leaf_depth_mean <- mean(df_tree_topo_int$depth)                       # mean depth (nb of ancestors) by leaf
  tree_nodes_children_asym_mean <- mean(df_tree_topo_ext$children_asym)      # mean absolute difference between descendants from the right and the left children nodes
  tree_nodes_depth_mean <- mean(df_tree_topo_ext$depth)                      # mean depth by internal node
  tree_nodes_width_mean <- mean(df_tree_topo_ext$width)                      # mean width (nb of nodes at the same depth)
  tree_WD_ratio <- max(df_tree_topo_ext$width)/max(df_tree_topo_ext$depth)   # max width / max depth
  tree_delta_width_max <- max(diff(depths_table))                            # maximal difference in width between 2 consecutive depths
  tree_asym_nodes_ratio <- mean(df_tree_topo_ext$children_asym>0)            # proportion of imbalanced internal nodes

  ## branch lengths
  tree_branch_length_min <- min(tree$edge.length)                            # all branches
  tree_branch_length_max <- max(tree$edge.length)
  tree_branch_length_mean <- mean(tree$edge.length)
  tree_branch_length_median <- median(tree$edge.length)
  tree_branch_length_sd <- sd(tree$edge.length)
  tree_branch_length_skewness <- moments::skewness(tree$edge.length)
  tree_branch_length_kurtosis  <- moments::kurtosis(tree$edge.length)
  tree_branch_length_sum <- sum(tree$edge.length)
  
  tree_int_branch_length <- tree$edge.length[tree$edge[,2] > Ntip(tree)]     # internal branches
  tree_int_branch_length_mean <- mean(tree_int_branch_length)
  tree_int_branch_length_median <- median(tree_int_branch_length)
  tree_int_branch_length_sd <- sd(tree_int_branch_length)

  tree_ext_branch_length <- tree$edge.length[tree$edge[,2] <= Ntip(tree)]    # external branches
  tree_ext_branch_length_mean <- mean(tree_ext_branch_length)
  tree_ext_branch_length_median <- median(tree_ext_branch_length)
  tree_ext_branch_length_sd <- sd(tree_ext_branch_length)
  
  tree_int_ext_branch_length_ratio_mean <- tree_int_branch_length_mean/tree_ext_branch_length_mean     # ratio between internal and external branches
  tree_int_ext_branch_length_ratio_median <- tree_int_branch_length_mean/tree_ext_branch_length_median
  tree_int_ext_branch_length_ratio_sd <- tree_int_branch_length_mean/tree_ext_branch_length_sd
  
  ## LTT = Lineage Through Time plot
  tree_LTT <- phytools::ltt(tree, plot = F)                # number of lineages through time
  tree_LTT_time_gaps <- diff(tree_LTT$times)[-1]           # times between 2 consecutive up steps
  
  tree_LTT_time_gaps_min <- min(tree_LTT_time_gaps)
  tree_LTT_time_gaps_max <- max(tree_LTT_time_gaps)
  tree_LTT_time_gaps_mean <- mean(tree_LTT_time_gaps)
  tree_LTT_time_gaps_median <- median(tree_LTT_time_gaps)
  tree_LTT_time_gaps_sd <- sd(tree_LTT_time_gaps)
  
  ## pairwise phylogenetic distances
  tree_distances_min <- min(cophenetic_tree_lower)
  tree_distances_mean <- mean(cophenetic_tree_lower)
  tree_distances_median <- median(cophenetic_tree_lower)
  tree_distances_sd <- sd(cophenetic_tree_lower)
  tree_distances_skewness <- moments::skewness(cophenetic_tree_lower)
  tree_distances_kurtosis  <- moments::kurtosis(cophenetic_tree_lower)

  summaries <- c(traits_min,                     ## raw trait distribution
                 traits_max,
                 traits_mean,
                 traits_median,
                 traits_sd,
                 traits_skewness,
                 traits_kurtosis,
                 traits_mean_median,
                 
                 tree_traits_IC_min,             ## independent contrasts
                 tree_traits_IC_max,
                 tree_traits_IC_mean,
                 tree_traits_IC_median,
                 tree_traits_IC_sd,
                 tree_traits_IC_skewness,
                 tree_traits_IC_kurtosis,
                 tree_traits_IC_mean_median,
                 
                 traits_NN_dist_min,             ## nearest-neighbor distance
                 traits_NN_dist_max,
                 traits_NN_dist_mean,
                 traits_NN_dist_median,
                 traits_NN_dist_sd,
                 traits_NN_dist_range,
                 traits_NN_dist_skewness,
                 traits_NN_dist_kurtosis,
                 
                 tree_Moran_I,                   ## phylogenetic signal
                 tree_Blomberg_K,
                 
                 #traits_pair_dist_min, (already with NN)
                 traits_pair_dist_max,           ## trait distance matrix
                 #traits_pair_dist_mean, (>0.999 correlation with traits_sd)
                 #traits_pair_dist_median, (>0.99 correlation with traits_sd)
                 #traits_pair_dist_sd, (>0.99 correlation with traits_sd)
                 
                 traits_pair_dist_eig_min,      ## eigen values of the trait distance matrix
                 traits_pair_dist_eig_max,
                 traits_pair_dist_eig_mean,
                 traits_pair_dist_eig_median,
                 traits_pair_dist_eig_sd,
                 traits_pair_dist_eig_skewness,
                 #traits_pair_dist_eig_kurtosis, (>0.999 correlation with tree_size)
                 
                 #traits_ordered_gaps_min, (already with NN)
                 traits_ordered_gaps_max,        ## distances between ordered trait values
                 traits_ordered_gaps_mean,
                 traits_ordered_gaps_median,
                 traits_ordered_gaps_sd,
                 traits_ordered_gaps_range,
                 
                 #tree_SDP_max, (>0.999 correlation with tree_size)
                 tree_SDP_asymmetry,             ## spectral density profile : tree
                 tree_SDP_peakedness,
                 tree_SDP_eigengap,
                 
                 tree_SDP_normal_max,            ## normalized spectral density profile : tree
                 tree_SDP_normal_asymmetry,
                 tree_SDP_normal_peakedness,
                 tree_SDP_normal_eigengap,
                 
                 tree_traits_SDP_splitter,       ## spectral density profile : tree + traits
                 tree_traits_SDP_tracer,
                 tree_traits_SDP_fragmenter,
                 
                 reg_trait_phylo_dist_intercept, ## regression trait distances ~ sqrt(phylogenetic distances)
                 reg_trait_phylo_dist_beta,
                 
                 tree_gamma,                     ## tree metrics
                 tree_beta,
                 
                 tree_size,                      ## tree topology
                 tree_leaf_depth_mean,
                 tree_nodes_children_asym_mean,
                 tree_nodes_depth_mean,
                 tree_nodes_width_mean,
                 tree_WD_ratio,
                 tree_delta_width_max,
                 tree_asym_nodes_ratio,
                 
                 tree_branch_length_min,         ## branch lengths : all branches
                 tree_branch_length_max,
                 tree_branch_length_mean,
                 tree_branch_length_sd,
                 tree_branch_length_skewness,
                 tree_branch_length_kurtosis,
                 tree_branch_length_sum,
                 
                 tree_int_branch_length_mean,    ## branch lengths : internal branches
                 tree_int_branch_length_median,
                 tree_int_branch_length_sd,
                 
                 tree_ext_branch_length_mean,    ## branch lengths : internal branches
                 tree_ext_branch_length_median,
                 tree_ext_branch_length_sd,
                 
                 tree_int_ext_branch_length_ratio_mean,      ## branch lengths : int-ext ratio
                 tree_int_ext_branch_length_ratio_median,
                 tree_int_ext_branch_length_ratio_sd,
                 
                 tree_LTT_time_gaps_min,         ## LTT = Lineage Through Time plot
                 tree_LTT_time_gaps_max,
                 tree_LTT_time_gaps_mean,
                 tree_LTT_time_gaps_median,
                 tree_LTT_time_gaps_sd,
                 
                 tree_distances_min,             ## pairwise phylogenetic distances
                 #tree_distances_max, (always 2*t)
                 tree_distances_mean,
                 tree_distances_median,
                 tree_distances_sd,
                 tree_distances_skewness,
                 tree_distances_kurtosis)

  
  names(summaries) <- c("traits_min",                     ## raw trait distribution
                        "traits_max",
                        "traits_mean",
                        "traits_median",
                        "traits_sd",
                        "traits_skewness",
                        "traits_kurtosis",
                        "traits_mean_median",
                        
                        "tree_traits_IC_min",             ## independent contrasts
                        "tree_traits_IC_max",
                        "tree_traits_IC_mean",
                        "tree_traits_IC_median",
                        "tree_traits_IC_sd",
                        "tree_traits_IC_skewness",
                        "tree_traits_IC_kurtosis",
                        "tree_traits_IC_mean_median",
                        
                        "traits_NN_dist_min",             ## nearest-neighbor distance
                        "traits_NN_dist_max",
                        "traits_NN_dist_mean",
                        "traits_NN_dist_median",
                        "traits_NN_dist_sd",
                        "traits_NN_dist_range",
                        "traits_NN_dist_skewness",
                        "traits_NN_dist_kurtosis",
                        
                        "tree_Moran_I",                   ## phylogenetic signal
                        "tree_Blomberg_K",
                        
                        #"traits_pair_dist_min",
                        "traits_pair_dist_max",           ## trait distance matrix
                        #"traits_pair_dist_mean",
                        #"traits_pair_dist_median",
                        #"traits_pair_dist_sd",
                        
                        "traits_pair_dist_eig_min",         ## eigen values of the trait distance matrix
                        "traits_pair_dist_eig_max",
                        "traits_pair_dist_eig_mean",
                        "traits_pair_dist_eig_median",
                        "traits_pair_dist_eig_sd",
                        "traits_pair_dist_eig_skewness",
                        #"traits_pair_dist_eig_kurtosis",
                        
                        #"traits_ordered_gaps_min",
                        "traits_ordered_gaps_max",        ## distances between ordered trait values
                        "traits_ordered_gaps_mean",
                        "traits_ordered_gaps_median",
                        "traits_ordered_gaps_sd",
                        "traits_ordered_gaps_range",
                        
                        #"tree_SDP_max",
                        "tree_SDP_asymmetry",             ## spectral density profile : tree
                        "tree_SDP_peakedness",
                        "tree_SDP_eigengap",
                        
                        "tree_SDP_normal_max",            ## normalized spectral density profile : tree
                        "tree_SDP_normal_asymmetry",
                        "tree_SDP_normal_peakedness",
                        "tree_SDP_normal_eigengap",
                        
                        "tree_traits_SDP_splitter",       ## spectral density profile : tree + traits
                        "tree_traits_SDP_tracer",
                        "tree_traits_SDP_fragmenter",
                        
                        "reg_trait_phylo_dist_intercept", ## regression trait distances ~ sqrt(phylogenetic distances)
                        "reg_trait_phylo_dist_beta",
                        
                        "tree_gamma",                     ## tree metrics
                        "tree_beta",
                        
                        "tree_size",                      ## tree topology
                        "tree_leaf_depth_mean",
                        "tree_nodes_children_asym_mean",
                        "tree_nodes_depth_mean",
                        "tree_nodes_width_mean",
                        "tree_WD_ratio",
                        "tree_delta_width_max",
                        "tree_asym_nodes_ratio",
                        
                        "tree_branch_length_min",         ## branch lengths : all branches
                        "tree_branch_length_max",
                        "tree_branch_length_mean",
                        "tree_branch_length_sd",
                        "tree_branch_length_skewness",
                        "tree_branch_length_kurtosis",
                        "tree_branch_length_sum",
                        
                        "tree_int_branch_length_mean",    ## branch lengths : internal branches
                        "tree_int_branch_length_median",
                        "tree_int_branch_length_sd",
                        
                        "tree_ext_branch_length_mean",    ## branch lengths : internal branches
                        "tree_ext_branch_length_median",
                        "tree_ext_branch_length_sd",
                        
                        "tree_int_ext_branch_length_ratio_mean",      ## branch lengths : int-ext ratio
                        "tree_int_ext_branch_length_ratio_median",
                        "tree_int_ext_branch_length_ratio_sd",
                        
                        "tree_LTT_time_gaps_min",         ## LTT = Lineage Through Time plot
                        "tree_LTT_time_gaps_max",
                        "tree_LTT_time_gaps_mean",
                        "tree_LTT_time_gaps_median",
                        "tree_LTT_time_gaps_sd",
                        
                        "tree_distances_min",             ## pairwise phylogenetic distances
                        #"tree_distances_max",
                        "tree_distances_mean",
                        "tree_distances_median",
                        "tree_distances_sd",
                        "tree_distances_skewness",
                        "tree_distances_kurtosis")
  return(summaries)
}
