#' Cluster Stability Testing/Quantification
#' Subset and cluster multiple times (store cluster assignments in table)
#' Then 1) find a way to assign universal names for clusters
#' and 2) quantify cluster stability for each cluster, and
#' also 3) quantify often each data point (cell) is assigned to its primary cluster
#' @import igraph
#' @param cluster_iter_list a list of lists.  First level list is just all the iterations.  Second
#'  level (inner) list is a list of each cluster, each one containing the indices of all the
#'  original datapoints assigned to that cluster
#' @param cluster_iter_table a matrix with each column as separate iteration of cluster assignments
#'  (same length as full or initially sampled dataset) --note that there should be NAs in this
#'  table because we want clustering iterations from samples of the dataset (to see how stable
#'  clusters are)
#' @export

consensus_and_stability <- function(clust_iter_list=NULL, clust_iter_table=NULL) {
  # There are 2 options for input, clust_iter_list or clust_iter_table:
  # 1) clust_iter_table is a matrix with each column as separate iteration of cluster assignments (same length as full or initially sampled dataset)
  # --note that there should be NAs in this table because we want clustering iterations from samples of the dataset (to see how stable clusters are)
  # 2) clust_iter_list is a list of lists.  First level list is just all the iterations.  Second level (inner) list is a list of each cluster, each
  # one containing the indices of all the original datapoints assigned to that cluster
  # ------------------- See example code above for examples of how the input to this function should be formatted


  # if clusters are only provided in table format, generate the list too
  if(is.null(clust_iter_list)) {
    clust_iter_list <- apply(clust_iter_table, 2, function(x){
      iter_groups <- list()
      for (c in sort(unique(x))) {
        iter_groups[[c]] <- which(x==c)
      }
      return(iter_groups)
    })
  }

  # if clusters are only provided in list format, generate the table too
  if(is.null(clust_iter_table)) {
    clust_iter_table <- matrix(data=NA, nrow=max(unlist(clust_iter_list)), ncol=length(clust_iter_list))
    for (i in 1:length(clust_iter_list)){
      for (j in 1:length(clust_iter_list[[i]])){
        clust_iter_table[clust_iter_list[[i]][[j]],i] <- j
      }
    }
  }

  #------use "clustering of clusters" to identify universal cluster names
  # first, flatten the list of lists into just a single list (clusters from all iterations, each item in the list containing the indices of cells assigned to that cluster)
  all_clusters <- unlist(clust_iter_list, recursive=FALSE)
  # now calculate the jaccard distance between all clusters in flattened list of clusters.  this is an all-by-all comparison, so the runtime increases exponentially and
  # it can take a long time if the number of iterations is high.  With a dataset that gives ~12 clusters, 100 iteration input = 60 seconds for this step, while 20 iterations = 0.5 seconds
  jac_dists <- matrix(data=unlist(lapply(all_clusters, function(y) unlist(lapply(all_clusters, function(x) length(which(x %in% y))/length(unique(c(x, y))))))),
                      nrow=length(all_clusters), ncol=length(all_clusters))
  diag(jac_dists) <- 0 #don't want to count the self-self distances, so make them zero

  # prepare vectors for graph construction.  jac_el = edge list, and jac_ed = edge distances
  jac_el <- as.vector(t(which(jac_dists !=0, arr.ind = T)))
  jac_ed <- as.vector(jac_dists)[which(jac_dists !=0)]

  # make graph that includes clusters from all iterations, connected by their jaccard distances
  # note: jaccard distance is already set up for more positive=closer together, so we don't need to transform from distance to weight here (it's already good to use as edge weight)
  jac_gr <- make_empty_graph(n=length(all_clusters), directed=FALSE)
  jac_gr <- add_edges(jac_gr, edges=jac_el, attr=list(weight=jac_ed))

  # now run clustering on the clusters, and get the global cluster assignments that we'll output as consensus clusters assignments
  cl_cl <- cluster_louvain(jac_gr)
  cl_cl_assignments <- cl_cl$membership

  # transform cluster numbering from local (each iteration) to global (consensus numbering)
  new_cluster_iter_table <- clust_iter_table #start with old table (local, non-consensus numbering) - this will be written over with consensus numbering
  counter <- 0 #this keeps track of the break points
  for (i in 1:length(clust_iter_list)) {
    keys <- cl_cl_assignments[(counter+1):(counter+length(cluster_groups[[i]]))] #take the subset for this iteration, based on counter and number of clusters in current iteration
    new_cluster_iter_table[,i] <- keys[new_cluster_iter_table[,i]] #write over with new consensus numbering
    counter <- counter + length(cluster_groups[[i]]) #update counter
  }
  # calculate the cluster occupancy.  For every cell, how many times was it assigned to each cluster?
  cl_occ <- t(apply(new_cluster_iter_table, 1, function(x) tabulate(x, nbins=max(cl_cl_assignments))))

  # go through each cell one by one, determine which cluster it was assigned to most frequently, then store this as its consensus cluster assignment
  final_assigns <- apply(cl_occ, 1, which.max)
  # for each cell, calculate what percentage of the time it was assigned to its consensus cluster, then store this as stability
  final_stability <- apply(cl_occ, 1, function(x) x[which.max(x)]/sum(x) )
  # return a two column matrix, containing final consensus cluster assignments and cluster stability values for every cell.
  return(cbind(consensus_clusters=final_assigns, stability=final_stability))
}
