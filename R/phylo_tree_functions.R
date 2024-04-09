#' Construct phylogenetic trees from consensus gene loadings across iNMF K sweep
#'

# build phylogenetic trees based on activity program gene loadings ####
build_phylo_tree <- function(nmf_results){
  # add in dummy program as outgroup to root tree (for depth-first search)
  W_list = lapply(nmf_results$consensus.results, '[[', "W")
  W_list = lapply(W_list, t)
  W_list[["root"]] <- as.matrix(data.frame(root = rowMedians(W_list$R2)))
  W_list = W_list[c("root", names(x$consensus.results))]

  # build phylogenetic tree
  cor_mat = cor(list.cbind(W_list))
  phylo_tree = fastme.bal(1-cor_mat)
  phylo_tree = root(phylo_tree, outgroup = "root", resolve.root = T)
  phylo_tree = drop.tip(phylo_tree, "root")
  # convert negative branches to zero and filter
  phylo_tree_pruned = phylo_tree
  phylo_tree_pruned$edge.length[phylo_tree_pruned$edge.length<0]=0
  dist_mat = cophenetic.phylo(phylo_tree_pruned)

  return(list(phylo_tree = phylo_tree, phylo_tree_pruned = phylo_tree_pruned, distance_matrix = dist_mat))
}

# suggest distance threshold for phylogenetic trees ####
suggest_dist_thresh <- function(phylo_trees){
  res_list = lapply(phylo_trees, function(x){
    mat = x$distance_matrix
    diag(mat) = NA
    res = list(min = quantile(colMins(mat, na.rm = T), 0.95), max = max(colMins(mat, na.rm = T)))
  })
  common_min_thresh = max(unlist(lapply(res_list, `[[`, "min")))
  common_max_thresh = min(unlist(lapply(res_list, `[[`, "max")))
  res = c(common_min_thresh, common_max_thresh)
  names(res) = c("min", "max")
  return(res)
}

# Identify "outlier" activity programs representing rare contaminating cells ####
identify_outlier_programs <- function(nmf_results, ncells = 50){
  test = lapply(nmf_results$consensus.results, `[[`, "H")
  res = lapply(test, function(y) {colMaxs(y)/apply(y, 2, function(z) {mean(sort(z, decreasing = T)[(1:ncells)+1])})})
  return(res)
}

#Visualization function to see programs filtered as outliers
plot_outlier_hists <- function(outlier_scores_by_cell, cell_type){
  outlier_sub <- outlier_scores_by_cell[[cell_type]]
  plots <- lapply(outlier_sub, function(X){
  data <- as.data.frame(X)
  rownames(data) <- 1:length(X)
  ggplot(data, aes(y = X, x = as.numeric(rownames(data)))) + geom_point(aes(y = X, color = X>5)) +
    xlab('Programs') + ylab('Outlier Influence') +
    theme_classic() +
    geom_hline(yintercept = 5, linetype = 2) + ylim(0, 20)  +
    theme(legend.position = 'none') +
    ggtitle(paste0('B cells ', regmatches(names(x), regexpr('^R[(0-9)]*', names(x)))))
  })
  marrangeGrob(grobs = outlier_plots, ncol = 3, nrow = 4, as.table = FALSE)
}

# Partition phylogenetic trees ####
partition_phylo_tree <- function(x, y, dist.thresh = NULL, outlier.thresh = 5){
  if (is.null(dist.thresh)){
    stop("run suggest_dist_thresh")
  }
  tree = x$phylo_tree_pruned
  nodes = 1:tree$Nnode+length(tree$tip.label)
  names(nodes) = nodes
  dist_mat = x$distance_matrix
  res = lapply(nodes, function(x){
    tiplist = tips(tree, x)
    dist = median(dist_mat[tiplist,tiplist][upper.tri(dist_mat[tiplist,tiplist], diag = F)])
  })
  dist = unlist(res)
  ntips<-Ntip(tree)
  cnum <- 0 ## cluster number
  assign <- rep(0,ntips+tree$Nnode) ## cluster assignment
  igraph.tree <- graph.edgelist(tree$edge) ## tree in igraph form
  dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',
                   order=TRUE,dist=TRUE)
  distvec = dist
  ## transverse the tree in depth first order
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    ## skip leaves
    if(node < ntips+1){ next }
    ## If the node's subtree is below the threshold, mark it and
    ## its subtree as members of a new cluster
    if(distvec[node-ntips]<=dist.thresh && assign[node]<=0){
      cnum <- cnum+1
      subtree <- graph.dfs(igraph.tree,node,
                           neimode='out',unreachable=FALSE)$order
      subtree <- subtree[! is.na(subtree)]
      assign[subtree] <- cnum
    }}
  ans <- as.character(assign)
  ans <- ans[1:ntips]
  names(ans) <- tree$tip.label

  # identify outlier activity programs
  outliers = unlist(lapply(y, function(z) names(which(z>outlier.thresh))))
  ans[outliers] = "outlier"

  # set minimum subtree size to at least 2
  ans = plyr::mapvalues(ans, names(which(table(ans)<=2)), rep("0", length(names(which(table(ans)<=2)))))
  return(ans)

}
