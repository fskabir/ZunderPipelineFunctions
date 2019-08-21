#' Cluster Stability Testing/Quantification
#' Subset and cluster multiple times (store cluster assignments in table)
#' Then 1) find a way to assign universal names for clusters
#' and 2) quantify cluster stability for each cluster, and
#' also 3) quantify often each data point (cell) is assigned to its primary cluster

#' Get dbscan clusters across multiple iterations
#' @param input_dataset matrix to be clustered on
#' @param subsample_size number of subsamples to take in each clustering iteration
#' @param num_iterations number of rounds of clustering to test before stability analysis
#' @param minPts hdbscan parameter
#' @importFrom dbscan hdbscan
#' @export
iterate_dbscan_clusters <- function(input_dataset=NULL, subsample_size=NULL, num_iterations=50,
                                    minPts=50) {
  cluster_assigns <- c()
  cluster_probs <- c()
  outlier_scores <- c()

  for (i in 1:num_iterations) {
    #further subset for interations of clustering
    to_sample_for_hdbscan <- sample(1:nrow(input_dataset), size=subsample_size, replace=FALSE)
    to_hdbscan <- input_dataset[to_sample_for_hdbscan,]

    #cluster based on knn graph
    hdbscan_output <- hdbscan(to_hdbscan, minPts=minPts)

    #convert from subsample back to "regular" sample so numbers will all line up in the end
    full_cluster_assigns = full_cluster_probs = full_outlier_scores = rep(NA, nrow(input_dataset))

    full_cluster_assigns[to_sample_for_hdbscan] <- hdbscan_output$cluster
    full_cluster_probs[to_sample_for_hdbscan] <- hdbscan_output$membership_prob
    full_outlier_scores[to_sample_for_hdbscan] <- hdbscan_output$outlier_scores

    #record cluster assignments
    cluster_assigns <- cbind(cluster_assigns, full_cluster_assigns)
    cluster_probs <- cbind(cluster_probs, full_cluster_probs)
    outlier_scores <- cbind(outlier_scores, full_outlier_scores)
    cat(paste("iteration ", i, "\n", sep=""))
  }
  colnames(cluster_assigns) <- colnames(cluster_probs) <- colnames(outlier_scores) <-
    1:ncol(cluster_assigns)

  return(list(cluster_assigns=cluster_assigns,
              cluster_probs=rowMeans(cluster_probs, na.rm=TRUE),
              outlier_scores=rowMeans(outlier_scores, na.rm=TRUE)))
}

#' Iterate through multiple rounds of louvain clustering
#' @param input_dataset matrix to be clustered on
#' @param subsample_size number of subsamples to take in each clustering iteration
#' @param num_iterations number of rounds of clustering to test before stability analysis
#' @param graph_knn number of neighbors on knn graph for clustering
#' @importFrom FNN get.knn
#' @importFrom dplyr groups
#' @export
iterate_louvain_clusters <- function(input_dataset=NULL, subsample_size=NULL, num_iterations=50,
                                     graph_knn=5) {
  cluster_assigns <- c()
  #cluster_groups <- vector(mode="list", length=num_iterations)
  cluster_groups <- list()

  for (i in 1:num_iterations) {
    #further subset for interations of clustering
    to_sample_for_knn <- sample(1:nrow(input_dataset), size=subsample_size, replace=FALSE)
    to_knn <- input_dataset[to_sample_for_knn,]

    #get nearest neighbors and process/prep for graph making
    knn_output <- get.knn(to_knn, k=graph_knn)
    el_i <- as.vector(sapply(1:nrow(knn_output$nn.index), function(x) rep(x,graph_knn)))
    el_j <- as.vector(t(knn_output$nn.index))
    el_d <- as.vector(t(knn_output$nn.dist))

    #convert distances to weight
    #el_d <- -el_d
    #el_d <- (el_d-min(el_d))/(max(el_d)-min(el_d))
    #el_d <- 1/el_d #could use other transforms here, see above lines that are commented out, etc.

    el_d <- el_d^(-0.01)

    #make knn graph with weights
    el <- cbind(el_i, el_j)
    gr <- make_empty_graph(n=subsample_size, directed=FALSE)
    gr <- add_edges(gr, edges=as.vector(t(el)), attr=list(weight=el_d))

    #cluster based on knn graph
    cl <- cluster_louvain(gr, weights=el_d)

    #convert from subsample back to "regular" sample so numbers will all line up in the end
    full_clust_assigns = rep(NA, nrow(input_dataset))
    full_clust_assigns[to_sample_for_knn] <- cl$membership

    #record cluster assignments
    cluster_assigns <- cbind(cluster_assigns, full_clust_assigns)

    #convert from subsample back to "regular" sample so numbers will all line up in the end
    full_groups <- lapply(groups(cl), function(x) to_sample_for_knn[x])

    #record cluster groups
    cluster_groups[[i]] <- full_groups
  }
  return(list(cluster_assigns=cluster_assigns, cluster_groups=cluster_groups,
              ig_graph=gr))
}

#' Get consensus clusters and measures of cluster stability
#' @param clust_iter_list a list of lists.  First level list is just all the iterations.  Second
#'  level (inner) list is a list of each cluster, each one containing the indices of all the
#'  original datapoints assigned to that cluster
#' @param clust_iter_table a matrix with each column as separate iteration of cluster assignments
#'  (same length as full or initially sampled dataset) --note that there should be NAs in this
#'  table because we want clustering iterations from samples of the dataset (to see how stable
#'  clusters are)         
#' @importFrom igraph make_empty_graph
#' @importFrom igraph add_edges
#' @importFrom igraph cluster_louvain
#' @export
consensus_and_stability <- function(clust_iter_list=NULL, clust_iter_table=NULL) {
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
    clust_iter_table <- matrix(data=NA, nrow=max(unlist(clust_iter_list)),
                               ncol=length(clust_iter_list))
    for (i in 1:length(clust_iter_list)){
      for (j in 1:length(clust_iter_list[[i]])){
        clust_iter_table[clust_iter_list[[i]][[j]],i] <- j
      }
    }
  }

  #------use "clustering of clusters" to identify universal cluster names
  # first, flatten the list of lists into just a single list (clusters from
  # all iterations, each item in the list containing the indices of cells assigned to that cluster)
  all_clusters <- unlist(clust_iter_list, recursive=FALSE)
  # now calculate the jaccard distance between all clusters in flattened list of clusters.  this is
  # an all-by-all comparison, so the runtime increases exponentially and
  # it can take a long time if the number of iterations is high.  With a dataset that gives ~12
  # clusters, 100 iteration input = 60 seconds for this step, while 20 iterations = 0.5 seconds
  jac_dists <- matrix(data=unlist(lapply(all_clusters, function(y)
    unlist(lapply(all_clusters, function(x) length(which(x %in% y))/length(unique(c(x, y))))))),
    nrow=length(all_clusters), ncol=length(all_clusters))
  diag(jac_dists) <- 0 #don't want to count the self-self distances, so make them zero

  # prepare vectors for graph construction.  jac_el = edge list, and jac_ed = edge distances
  jac_el <- as.vector(t(which(jac_dists !=0, arr.ind = T)))
  jac_ed <- as.vector(jac_dists)[which(jac_dists !=0)]

  # make graph that includes clusters from all iterations, connected by their jaccard distances
  # note: jaccard distance is already set up for more positive=closer together, so we don't need to
  # transform from distance to weight here (it's already good to use as edge weight)
  jac_gr <- make_empty_graph(n=length(all_clusters), directed=FALSE)
  jac_gr <- add_edges(jac_gr, edges=jac_el, attr=list(weight=jac_ed))

  # now run clustering on the clusters, and get the global cluster assignments that we'll output as
  # consensus clusters assignments
  cl_cl <- cluster_louvain(jac_gr)
  cl_cl_assignments <- cl_cl$membership
        
  # transform cluster numbering from local (each iteration) to global (consensus numbering)
  #start with old table (local, non-consensus numbering) - this will be written over with
  #consensus numbering
  new_cluster_iter_table <- clust_iter_table
  counter <- 0 #this keeps track of the break points
  for (i in 1:length(clust_iter_list)) {
    #take the subset for this iteration, based on counter and number of clusters in
    #current iteration
    keys <- cl_cl_assignments[(counter+1):(counter+length(clust_iter_list[[i]]))]
    #write over with new consensus numbering
    new_cluster_iter_table[,i] <- keys[new_cluster_iter_table[,i]]
    counter <- counter + length(clust_iter_list[[i]]) #update counter
  }
  # calculate the cluster occupancy.  For every cell, how many times was it assigned to
  # each cluster?
  cl_occ <- t(apply(new_cluster_iter_table, 1, function(x) tabulate(x,
                                                                    nbins=max(cl_cl_assignments))))

  # go through each cell one by one, determine which cluster it was assigned to most frequently,
  # then store this as its consensus cluster assignment
  final_assigns <- apply(cl_occ, 1, which.max)
  # for each cell, calculate what percentage of the time it was assigned to its consensus cluster,
  # then store this as stability
  final_stability <- apply(cl_occ, 1, function(x) x[which.max(x)]/sum(x) )
  # return a two column matrix, containing final consensus cluster assignments and cluster
  # stability values for every cell.
  return(cbind(consensus_clusters=final_assigns, stability=final_stability))
}

#' Get UMAP layout for clustering of clusters
#' @param clust_iter_list a list of lists.  First level list is just all the iterations.  Second
#'  level (inner) list is a list of each cluster, each one containing the indices of all the
#'  original datapoints assigned to that cluster
#' @param clust_iter_table a matrix with each column as separate iteration of cluster assignments
#'  (same length as full or initially sampled dataset) --note that there should be NAs in this
#'  table because we want clustering iterations from samples of the dataset (to see how stable
#'  clusters are)          
#' @importFrom igraph make_empty_graph
#' @importFrom igraph add_edges
#' @importFrom igraph cluster_louvain
#' @export                           
consensus_cluster_umap <- function(clust_iter_list=NULL, clust_iter_table=NULL) {
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
    clust_iter_table <- matrix(data=NA, nrow=max(unlist(clust_iter_list)),
                               ncol=length(clust_iter_list))
    for (i in 1:length(clust_iter_list)){
      for (j in 1:length(clust_iter_list[[i]])){
        clust_iter_table[clust_iter_list[[i]][[j]],i] <- j
      }
    }
  }

  #------use "clustering of clusters" to identify universal cluster names
  # first, flatten the list of lists into just a single list (clusters from
  # all iterations, each item in the list containing the indices of cells assigned to that cluster)
  all_clusters <- unlist(clust_iter_list, recursive=FALSE)
  # now calculate the jaccard distance between all clusters in flattened list of clusters.  this is
  # an all-by-all comparison, so the runtime increases exponentially and
  # it can take a long time if the number of iterations is high.  With a dataset that gives ~12
  # clusters, 100 iteration input = 60 seconds for this step, while 20 iterations = 0.5 seconds
  jac_dists <- matrix(data=unlist(lapply(all_clusters, function(y)
    unlist(lapply(all_clusters, function(x) length(which(x %in% y))/length(unique(c(x, y))))))),
    nrow=length(all_clusters), ncol=length(all_clusters))
  diag(jac_dists) <- 0 #don't want to count the self-self distances, so make them zero

  # prepare vectors for graph construction.  jac_el = edge list, and jac_ed = edge distances
  jac_el <- as.vector(t(which(jac_dists !=0, arr.ind = T)))
  jac_ed <- as.vector(jac_dists)[which(jac_dists !=0)]

  # make graph that includes clusters from all iterations, connected by their jaccard distances
  # note: jaccard distance is already set up for more positive=closer together, so we don't need to
  # transform from distance to weight here (it's already good to use as edge weight)
  jac_gr <- make_empty_graph(n=length(all_clusters), directed=FALSE)
  jac_gr <- add_edges(jac_gr, edges=jac_el, attr=list(weight=jac_ed))

  # now run clustering on the clusters, and get the global cluster assignments that we'll output as
  # consensus clusters assignments
  cl_cl <- cluster_louvain(jac_gr)
  cl_cl_assignments <- cl_cl$membership
  
  #------use distance matrix to generate umap layout of clusters
  # currently "dist" are actually weights (this is what is reqd for igraph)
  # set the self-self distances back to almost 1, then do 1.0001-weights to get distances (so nothing ends up being zero)
  diag(jac_dists) <- 1.0 
  jac_dists <- 1.0001-jac_dists
  # run umap with default setting except for specifying that input is distance matrix         
  umap.settings <- umap.defaults
  umap.settings$input <- 'dist'
  clusters_umap <- umap(jac_dists, config = umap.settings)
  return(clusters_umap)
}
                           
#' Take a dataset along with its cluster assignments and cluster certainty/stability
#' (cutoff_values), and calculate statistics related to how "good" the clusters are.
#'
#' Outputs
#'
#' In the outputted list, "cluster_stats"
#'
#'  "min_dist_btwn_means"
#'  Take the centroid of every cluster, and then for each cluster, look at how far away its closest
#'  neighbor is in n-dimensional space.  If a cluster is far away from its closest neighbor, then
#'  it's well separated more likely to be a "good cluster."  Distances here are calculated by city
#'  block/manhattan.
#'
#'  "min_dist_btwn_means_single_param"
#'  Same as above, but calculate for each individual dimension, and then output the biggest
#'  distance to a closest neighbor in 1D space.  The purpose of this is to identify if there's one
#'  cluster that's very far away from all the clusters (well separated) by just one parameter -
#'  this is another way of identifying a "good" cluster
#'
#'  "sum_dists_btwn_means"
#'  This is the same as "min_dist_btwn_means", but instead of taking just the closest neighbor, it
#'  looks at the n-closest neighbors, and sums/averages the distance to all of them.  This
#'  statistic is intended to identify "good" clusters that are very close to each other, but far
#'  from everything else.  Uses a diminishing scale, so 1st neighbor is counted more than 2nd
#'  neighbor, which is counted more than 3rd neighbor, etc.  This statistic will be useful for
#'  automated parameter detection (where we're trying to see if we get "good" clustering with a
#'  given set of parameters) but will not be as useful as a general readout of individual cluster
#'  goodness.
#'
#'  "combined_stdevs"
#'  Standard deviation for each cluster.  This is calculated for each parameter/dimension
#'  individually, and then summed/averaged.  Clusters with lower standard deviation are "better"
#'  (if everything else is equal).
#'
#'  "combined_dist_from_mean_by_cluster"
#'  This is another way of looking at variance within a cluster in addition to standard deviation.
#'  Find the centroid of each cluster, and then calculate how far every point in the cluster is
#'  away from its centroid.  Take the average of these distances and that tells you a bit about how
#'  spread out the cluster is in n-dimensional space.  Clusters with low variance here are "better"
#'  (when everything else is equal).
#'
#'  "c_dist_to_nearest_point"
#'  This is similar to "min_dist_btwn_means" in that it wants to know the distances between
#'  clusters, but it calculates it in a different way.  Instead of looking at means, this looks at
#'  two clusters and finds the shortest distance between a single point in each cluster.  So what
#'  this does is for each cluster, calculate the distance to the nearest point outside that
#'  cluster.  If this distance is far, then the cluster is "well separated" and this makes it a
#'  better cluster.
#'
#'  "c_sum_dists_to_nearest_points"
#'  Like above, but this should sum/average the distance not just to the nearest point/cluster, but
#'  find the n-nearest clusters (by individual points) and then sum average over all these
#'  distances.  Uses a diminishing scale like "sum_dists_btwn_means."  Rationale for looking over
#'  multiple closest clusters instead of just the single closest is the same as for
#'  "sum_dists_btwn_means".
#'
#'  "c_dist_to_nearest_point_trim"
#'  Same as "c_dist_to_nearest_point", but using clusters that are trimmed by stability here.  The
#'  idea for this is that you may have two clusters that are very far apart (well separated),
#'  except that they're connected by a few low stability points that fall in between, so it seems
#'  like they're actually very close to each other.
#'
#'  "c_sum_dists_to_nearest_points_trim"
#'  Like "c_dist_to_nearest_point_trim", but summed/averaged over the n-nearest points, rather than
#'  looking just at the single nearest cluster.  See other parameters for more explanation.
#'
#'  "stability_by_cluster"
#'  See "stability_by_cell" and average this over the whole cluster.  1) Iterate clustering over
#'  many different subsamples. 2) Find the universal cluster/groups and re-assign.  3) For each
#'  cell, what percentage of the time does it fall into it's universal cluster?  That's the cluster
#'  stability by cell.  4) Average over the whole cluster.  Higher stability means better/more
#'  reproducible clusters.
#'
#'  "c_dist_to_nearest_point_single_param"
#'  Same as above - but looking one parameter at a time.  To take into account the scenario where a
#'  cluster is VERY different from its neighbors in one parameter, but this gets drowned out by
#'  being similar in everything else.
#'
#'  "c_dist_to_nearest_point_trim_single_param"
#'  Same as above, but trimmed to remove low stability cells.
#'
#'  "combined_dist_from_mean_by_cell"
#'  Same as above ("combined_dist_from_mean_by_cluster") but looking at each cell individually
#'  rather than by cluster
#'
#'  "stability_by_cell"
#'  Same as described above in "stability_by_cluster", but by cell.  This is what we calculated
#'  and output in the older version of this script.
#'
#' In the outputted list, "param_names"
#'
#' These are to identify which marker is driving the biggest difference between clusters.  This
#' should give an idea about which marker is most important for defining each cell type.
#'
#'  "c_dist_to_nearest_point_single_param_name"
#'  The parameter that drives "c_dist_to_nearest_point_single_param" from above, for each cluster
#'
#'  "c_dist_to_nearest_point_trim_single_param_name"
#'  The parameter that drives "c_dist_to_nearest_point_trim_single_param" from above, for each
#'  cluster
#'
#'  "min_dist_btwn_means_single_param_name"
#'  The parameter that drives "min_dist_btwn_means_single_param" from above, for each cluster
#'
#' In the outputted list, "mst_gap_params"
#'
#' This section is mainly to be used for selecting optimal clustering parameters, and is perhaps
#' not as useful as a general tool for cluster description.  The idea is to find the biggest "gap"
#' in the dataset, with the idea being that if there is a big gap between sections of the dataset,
#' then the parameters must be identifying some useful modality/variation.  How it works is like
#' this: identify the cluster means, and then treat these as individual points and connect them
#' with a minimum spanning tree (mst).  After this is done, the longest edge in the mst should
#' identify the biggest gap between regions of cell identity in the dataset.
#'
#'  "params_ranked"
#'  Once the biggest "gap" is identified, look at the two vertices (cell clusters) in the mst that
#'  form that edge, and see what parameters they're most different for.  I just want to know what
#'  parameters are driving the biggest separation between neighboring clusters.
#'
#'  "v1_stats" and "v2_stats"
#'  One concern about using the MST-gap to identify differences between clusters is that it might
#'  be dominated by a "garbage" cluster that is made up of say, debris and not real cells.  This
#'  could be very different from everything else, but not useful for our analysis.  To protect
#'  against this, check for each of the connected cell clusters 1) how big they are (as a
#'  percentage of total cells), and 2) how stable they are.  Presumably, spurious/debris clusters
#'  would be small and/or have low cluster stability. . .
#'
#' @param clust_dataset original dataset which was clustered
#' @param clust_assigns cluster assignments from consensus clustering
#' @param cutoff_values cluster stability confidence values
#' @param cutoff_trim_factor minimum cluster stability confidence value to keep cell
#' @param knn_to_sum_factor fraction of clusters to be compared to each clusetr
#' @param knn_to_sum_min minimum number of clusters to be compared to each cluster
#' @importFrom stats sd
#' @importFrom stats dist
#' @importFrom FNN get.knnx
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph E
#' @importFrom igraph E<-
#' @importFrom igraph mst
#' @importFrom igraph ends
#'
#' @export
cluster_goodness <- function(clust_dataset=NULL, clust_assigns=NULL, cutoff_values=NULL,
                             cutoff_trim_factor=0.8, knn_to_sum_factor=0.3, knn_to_sum_min=5) {

  # If there's less than 2 clusters present, don't perform analysis of cluster goodness
  if(length(sort(unique(clust_assigns[which(clust_assigns!=0)]))) < 2) {
    return(list(cluster_stats=NA, global_stats=NA, param_names=NA, mst_gap_params=NA))
  }

  # Identify all clusters present, to loop through
  final_clusters <- sort(unique(clust_assigns))
  # This accounts for hdbscan output, which uses zero to denote outliers.  Don't want to use these.
  final_clusters <- final_clusters[which(final_clusters > 0)]

  # First initialize all the things that will be filled in during the for loop
  c_means <- c()
  c_stdevs <- c()
  c_dist_from_mean_by_cell <- matrix(data=NA, nrow=nrow(clust_dataset), ncol=ncol(clust_dataset)) #NAs will be written over
  colnames(c_dist_from_mean_by_cell) <- colnames(clust_dataset)
  c_dist_from_mean_by_cluster <- c()
  c_dist_to_nearest_point <- c()
  c_sum_dists_to_nearest_points <- c()
  c_dist_to_nearest_point_single_param <- c()
  c_dist_to_nearest_point_single_param_name <- c()
  c_dist_to_nearest_point_trim <- c()
  c_sum_dists_to_nearest_points_trim <- c()
  c_dist_to_nearest_point_trim_single_param <- c()
  c_dist_to_nearest_point_trim_single_param_name <- c()
  cutoff_vals_by_cluster <- c()

  # Go through clusters 1 by 1 and calculate cluster "goodness" statistics - filling in the empty initialized containers from above
  for (c in final_clusters) {
    c_inds <- which(clust_assigns==c) #get all indices for this cluster
    c_subset <- clust_dataset[c_inds, ] #isolate this cluster only, will be reusing this later

    # Get cluster means and stdevs for every measurement parameter
    # Note: these are by individual measurement parameters still, so will want to take row averages later
    c_means <- rbind(c_means, colMeans(c_subset))
    c_stdevs <- rbind(c_stdevs, apply(c_subset, 2, sd))

    # Get the distance between cluster means - note that this is by individual cells now
    # Note: these are by individual measurement parameters still, so will want to take row averages later
    tmean <- t(c_dist_from_mean_by_cell) #need to do this transpozing because of how R works when changing multiple matrix rows at once.
    tmean[,c_inds] <- t(abs(t(t(c_subset)-colMeans(c_subset)))) #add in distance from mean, only to the right indices from c_inds (the rest stay NA until it's their turn to be overwritten)
    c_dist_from_mean_by_cell <- t(tmean) #transpose back.

    # Go from individual cells numbers to avergage over all cells in that cluster
    # Note- this is still by parameter, so will need to take row averages later
    c_dist_from_mean_by_cluster <- rbind(c_dist_from_mean_by_cluster,
                                         colMeans(c_dist_from_mean_by_cell[c_inds,]))

    # This is a faster way to find the single nearest point than below, but we'll want to
    # sum up nearest point from all other clusters, so use a for loop for each cluster instead
    #non_outlier_others <- which(clust_assigns !=c & clust_assigns !=0) #non-zero is important for dbscan clustering - don't wont to include outliers as potential neighbors
    #nearest_point_knn <- get.knnx(data=clust_dataset[non_outlier_others,], query=c_subset, k=1)$nn.dist

    # Find the distance to the nearest point for each other cluster, then find the min of all these
    nearests_points_by_cluster <- list()
    for (non_c in final_clusters[which(final_clusters != c)]) {
      non_c_inds <- which(clust_assigns==non_c & clust_assigns!=0) #non-zero is important for dbscan clustering - don't wont to include outliers as potential neighbors
      nearest_point_by_cluster_knn <- get.knnx(data=clust_dataset[non_c_inds, , drop=FALSE],
                                               query=clust_dataset[c_inds, , drop=FALSE], k=1)$nn.dist # Note: use drop=FALSE here for the edge case to keep R from converting one row matrices into vectors!
      nearests_points_by_cluster[[as.character(non_c)]] <- min(nearest_point_by_cluster_knn)
    }
    c_dist_to_nearest_point <- c(c_dist_to_nearest_point, min(unlist(nearests_points_by_cluster)))

    # Calculate the number of clusters to average the distance to from currernt cluster c
    num_clusters_to_compare <- ceiling(length(final_clusters)*knn_to_sum_factor)
    if (num_clusters_to_compare < knn_to_sum_min) { num_clusters_to_compare <- knn_to_sum_min}
    if (num_clusters_to_compare > (length(final_clusters)-1)) { num_clusters_to_compare <- length(final_clusters)-1}
    # Calculate the average distance to the n closest cluster (by nearest point)
    c_sum_dists_to_nearest_points <- c(c_sum_dists_to_nearest_points,
                                       mean(sort(unlist(nearests_points_by_cluster))[1:num_clusters_to_compare]))

    # Repeat the calculations from above (closest cluster by point), but look for the dimension that contributes the most
    nearest_points_by_cluster_and_p <- list()
    # Look (and record) the min distance for each parameter one at a time
    for (p in colnames(clust_dataset)) {
      for (non_c in final_clusters[which(final_clusters != c)]) {
        non_c_inds <- which(clust_assigns == non_c & clust_assigns != 0) #non-zero is important for dbscan clustering - don't wont to include outliers as potential neighbors

        #note: use drop=FALSE here to keep R from converting one row matrices into vectors!
        nearest_point_by_cluster_and_p_knn <- get.knnx(data=clust_dataset[non_c_inds, p, drop=FALSE],
                                                       query=clust_dataset[c_inds, p, drop=FALSE], k=1)$nn.dist
        nearest_points_by_cluster_and_p[[p]][[as.character(non_c)]] <- min(nearest_point_by_cluster_and_p_knn)
      }
    }
    # Out of all the min point distances, choose the biggest one, and record the parameter name
    c_dist_to_nearest_point_single_param <- c(c_dist_to_nearest_point_single_param,
                                              max(sapply(nearest_points_by_cluster_and_p,min)))
    c_dist_to_nearest_point_single_param_name <- c(c_dist_to_nearest_point_single_param_name,
                                                   names(sapply(nearest_points_by_cluster_and_p,min))[which.max(sapply(nearest_points_by_cluster_and_p,min))])

    # Repeat the distances by point calculations from above, this time with trimmed clusters
    # Trimming low stability/confidence points helps give "truer distances". . .?
    # This makes it more similar to distance between means, but not all the way (depending on cutoff choice)

    # First do some tests to make sure there's anything left after trimming
    num_trimmed_cells <- length(which(cutoff_values > cutoff_trim_factor))
    if(num_trimmed_cells > 1) { #want to have at least two cells!
      trimmed_clust_names <- sort(unique(clust_assigns[which(cutoff_values > cutoff_trim_factor)]))
      trimmed_clust_names <- trimmed_clust_names[which(trimmed_clust_names != 0)] # Shouldn't need this test here, because zero-labeled hdbscan outliers shouldn't make it past cutoff)
      num_trimmed_clust_names <- length(trimmed_clust_names)
      if(num_trimmed_clust_names > 1) { #check that there are at least two clusters that make it through the cutoff

        # Trim the dataset, cluster assignments, and cutoff values (by the same index set)
        trim_clust_dataset <- clust_dataset[which(cutoff_values > cutoff_trim_factor),]
        trim_clust_assigns <- clust_assigns[which(cutoff_values > cutoff_trim_factor)]
        trim_cutoff_values <- cutoff_values[which(cutoff_values > cutoff_trim_factor)]

        # Go through the same calculations as above, now on trimmed clusters

        trim_c_inds <- which(trim_clust_assigns==c & trim_clust_assigns!=0)
        if(length(trim_c_inds) > 0) { #this test is important for edge case where current cluster c has no cells after trimming

          # Find the distance to the nearest point for each other cluster, then find the min of all these
          nearests_points_by_cluster_trim <- list()
          all_non_c <- trim_clust_assigns[which(trim_clust_assigns != c)]
          for (non_c in all_non_c) {
            trim_non_c_inds <- which(trim_clust_assigns==non_c & trim_clust_assigns!=0) #non-zero is important for dbscan clustering - don't wont to include outliers as potential neighbors
            #note: use drop=FALSE here to keep R from converting one row matrices into vectors!
            nearest_point_by_cluster_knn_trim <- get.knnx(data=trim_clust_dataset[trim_non_c_inds, , drop=FALSE],
                                                          query=trim_clust_dataset[trim_c_inds, , drop=FALSE], k=1)$nn.dist
            nearests_points_by_cluster_trim[[as.character(non_c)]] <- min(nearest_point_by_cluster_knn_trim)
          }
          c_dist_to_nearest_point_trim <- c(c_dist_to_nearest_point_trim,
                                            min(unlist(nearests_points_by_cluster_trim)))

          # Calculate the number of clusters, to average the distance to from currernt cluster c
          num_clusters_to_compare <- ceiling(num_trimmed_clust_names*knn_to_sum_factor)
          if (num_clusters_to_compare < knn_to_sum_min) {
            num_clusters_to_compare <- knn_to_sum_min
          }
          if (num_clusters_to_compare > (num_trimmed_clust_names-1)) {
            num_clusters_to_compare <- num_trimmed_clust_names-1
          }
          c_sum_dists_to_nearest_points_trim <- c(c_sum_dists_to_nearest_points_trim,
                                                  mean(sort(unlist(nearests_points_by_cluster_trim))[1:num_clusters_to_compare]))

          # Repeat the calculations from above (closest cluster by point), but look for the dimension that contributes the most
          nearests_points_trim_by_cluster_and_p <- list()
          # Look (and record) the min distance for each parameter one at a time
          for (p in colnames(clust_dataset)) {
            for (non_c in trimmed_clust_names[which(trimmed_clust_names != c)]) {
              trim_non_c_inds <- which(trim_clust_assigns==non_c & trim_clust_assigns!=0) #non-zero is important for dbscan clustering - don't wont to include outliers as potential neighbors
              #note: use drop=FALSE here to keep R from converting one row matrices into vectors!
              nearest_point_trim_by_cluster_and_p_knn <- get.knnx(data=trim_clust_dataset[trim_non_c_inds, p, drop=FALSE],
                                                                  query=trim_clust_dataset[trim_c_inds, p, drop=FALSE], k=1)$nn.dist
              nearests_points_trim_by_cluster_and_p[[p]][[as.character(non_c)]] <- min(nearest_point_trim_by_cluster_and_p_knn)
            }
          }
          # Out of all the min point distances, choose the biggest one, and record the parameter name
          c_dist_to_nearest_point_trim_single_param <- c(c_dist_to_nearest_point_trim_single_param,
                                                         max(sapply(nearests_points_trim_by_cluster_and_p,min)))
          c_dist_to_nearest_point_trim_single_param_name <- c(c_dist_to_nearest_point_trim_single_param_name,
                                                              names(sapply(nearests_points_trim_by_cluster_and_p,min))[which.max(sapply(nearests_points_trim_by_cluster_and_p,min))])
        } else {
          # This is what happens if the cluster c is completely trimmed away and no cells are left
          c_dist_to_nearest_point_trim <- c(c_dist_to_nearest_point_trim,0)
          c_sum_dists_to_nearest_points_trim <- c(c_sum_dists_to_nearest_points_trim,0)
          c_dist_to_nearest_point_trim_single_param <- c(c_dist_to_nearest_point_trim_single_param,0)
          c_dist_to_nearest_point_trim_single_param_name <- c(c_dist_to_nearest_point_trim_single_param_name,NULL)
        }
      } else {
        # Put in zeros if not enough clusters pass the stability trim cutoff
        c_dist_to_nearest_point <- c(c_dist_to_nearest_point, 0)
        c_sum_dists_to_nearest_points <- c(c_sum_dists_to_nearest_points, 0)
        c_dist_to_nearest_point_single_param <- c(c_dist_to_nearest_point_single_param, 0)
        c_dist_to_nearest_point_single_param_name <- c(c_dist_to_nearest_point_single_param_name,NULL)
        c_dist_to_nearest_point_trim <- c(c_dist_to_nearest_point_trim, 0)
        c_sum_dists_to_nearest_points_trim <- c(c_sum_dists_to_nearest_points_trim, 0)
        c_dist_to_nearest_point_trim_single_param <- c(c_dist_to_nearest_point_trim_single_param, 0)
        c_dist_to_nearest_point_trim_single_param_name <- c(c_dist_to_nearest_point_trim_single_param_name,NULL)
      }
    } else {
      # Put in zeros if not enough cells pass the stability trim cutoff
      c_dist_to_nearest_point <- c(c_dist_to_nearest_point, 0)
      c_sum_dists_to_nearest_points <- c(c_sum_dists_to_nearest_points, 0)
      c_dist_to_nearest_point_single_param <- c(c_dist_to_nearest_point_single_param, 0)
      c_dist_to_nearest_point_single_param_name <- c(c_dist_to_nearest_point_single_param_name,NULL)
      c_dist_to_nearest_point_trim <- c(c_dist_to_nearest_point_trim, 0)
      c_sum_dists_to_nearest_points_trim <- c(c_sum_dists_to_nearest_points_trim, 0)
      c_dist_to_nearest_point_trim_single_param <- c(c_dist_to_nearest_point_trim_single_param, 0)
      c_dist_to_nearest_point_trim_single_param_name <- c(c_dist_to_nearest_point_trim_single_param_name,NULL)
    }
    # Finally, output the average values for stability/confidence/cutoff for each cluster
    cutoff_vals_by_cluster <- c(cutoff_vals_by_cluster, mean(cutoff_values[c_inds]))
  }

  # Find distance to nearest neighbor by cluster mean
  dists_btwn_means <- get.knn(data=c_means, k=1)$nn.dist
  min_dist_btwn_means <- dists_btwn_means[,1]
  sum_dists_btwn_means <- rowMeans(dists_btwn_means)

  # Find (and record) distance to nearest neighbor by cluster mean, one parameter at a time
  nearest_mean_1d_list <- list()
  for (p in colnames(clust_dataset)) {
    nearest_mean_1d_list[[p]] <- get.knn(data=c_means[,p], k=1)$nn.dist
  }
  nearest_mean_1d_table <- sapply(nearest_mean_1d_list,function(x) x) #this is how I get it into a matrix, so the max and param names can be found in the next step
  # Get the parameter which has the largest minimum distance between cluster means
  min_dist_btwn_means_single_param <- apply(nearest_mean_1d_table, 1, max)
  min_dist_btwn_means_single_param_name <- colnames(nearest_mean_1d_table)[apply(nearest_mean_1d_table, 1, which.max)]

  # This next section is set up to find the biggest "gap" in a dataset.
  # Looking for a big gap that defines separation from two clusters/regions

  # Make a graph with distances between cluster means
  means_graph <- graph_from_adjacency_matrix(adjmatrix=as.matrix(dist(c_means, upper=TRUE, diag=TRUE)),
                                             mode="undirected", weighted=TRUE, diag=FALSE)
  # Add in weights as inverse of distance
  E(means_graph)$weight <- 1/E(means_graph)$weight
  # Calculate the minimum spanning tree (MST)
  means_mst <- mst(graph=means_graph, weights=E(means_graph)$weight)
  # Find the largest distance between two means in the MST (lowest weight edge)
  lowest_mst_edge_weight_index <- which.min(E(means_mst)$weight)

  # Calculate and record stats about the largest gap in the MST
  largest_mst_gap_distance <- 1/E(means_mst)$weight[lowest_mst_edge_weight_index] #transform back
  largest_mst_gap_vertices <- as.numeric(ends(means_mst, E(means_mst)[lowest_mst_edge_weight_index]))
  largest_mst_gap_vertex1_colmeans <- colMeans(clust_dataset[which(clust_assigns==largest_mst_gap_vertices[1]),])
  largest_mst_gap_vertex2_colmeans <- colMeans(clust_dataset[which(clust_assigns==largest_mst_gap_vertices[2]),])
  largest_mst_gap_param_diffs <- abs(largest_mst_gap_vertex1_colmeans-largest_mst_gap_vertex2_colmeans)
  # This will give info on what's driving the biggest gap in the dataset
  largest_mst_gap_params_ranked <- largest_mst_gap_param_diffs[order(largest_mst_gap_param_diffs, decreasing=TRUE)]
  # This will give an idea about how trustworthy thee MST gap info is, based on cluster size/quality
  # -----------Might want to add the ability to trim by cluster quality in the future
  largest_mst_gap_v1_stats <- c(size=length(which(clust_assigns==largest_mst_gap_vertices[1]))/nrow(clust_dataset)*100,
                                confidence=cutoff_vals_by_cluster[largest_mst_gap_vertices[1]])
  largest_mst_gap_v2_stats <- c(size=length(which(clust_assigns==largest_mst_gap_vertices[2]))/nrow(clust_dataset)*100,
                                confidence=cutoff_vals_by_cluster[largest_mst_gap_vertices[2]])
  # This is the list that will be output by cluster_goodness() function
  output_mst_gap <- list(params_ranked=largest_mst_gap_params_ranked,
                         v1_stats=largest_mst_gap_v1_stats,
                         v2_stats=largest_mst_gap_v2_stats)
  # This will also be output for visualization
  mst_gap_by_cluster <- rep(0, length.out=length(final_clusters))
  mst_gap_by_cluster[largest_mst_gap_vertices] <- 1 #high values for the MST gap clusters

  # Go from individual marker values to combined values from all markers

  # Average stdev across all parameters
  combined_stdevs <- rowMeans(c_stdevs)
  # Average dists from mean - note: these are by cell, not by cluster
  combined_dist_from_mean_by_cell <- rowMeans(c_dist_from_mean_by_cell)
  # Same as above for mean, but averaged over clusters
  combined_dist_from_mean_by_cluster <- rowMeans(c_dist_from_mean_by_cluster)

  # Take all the parameters we want to visualize that are by cluster, and expand to all cells (for visualization)
  expand_from_cluster_to_cells_params <- c("min_dist_btwn_means", "sum_dists_btwn_means", "min_dist_btwn_means_single_param",
                                           "mst_gap_by_cluster", "c_dist_to_nearest_point", "c_dist_to_nearest_point_trim",
                                           "c_sum_dists_to_nearest_points", "c_sum_dists_to_nearest_points_trim",
                                           "c_dist_to_nearest_point_single_param", "c_dist_to_nearest_point_trim_single_param",
                                           "combined_dist_from_mean_by_cluster", "combined_stdevs", "cutoff_vals_by_cluster")

  expand_from_cluster_to_cells_orig <- cbind(min_dist_btwn_means, sum_dists_btwn_means, min_dist_btwn_means_single_param,
                                             mst_gap_by_cluster, c_dist_to_nearest_point, c_dist_to_nearest_point_trim,
                                             c_sum_dists_to_nearest_points, c_sum_dists_to_nearest_points_trim,
                                             c_dist_to_nearest_point_single_param, c_dist_to_nearest_point_trim_single_param,
                                             combined_dist_from_mean_by_cluster, combined_stdevs, cutoff_vals_by_cluster)

  # Initialize matrix for expanding from clusters to cells
  expanded_cluster_stats <- matrix(data=NA, nrow=nrow(clust_dataset), ncol=length(expand_from_cluster_to_cells_params))
  colnames(expanded_cluster_stats) <- expand_from_cluster_to_cells_params
  # Add in everything one cluster at a time, using transposed matric because of how R adds by column (not by row) when replacing multiple values
  t_expanded_cluster_stats <- t(expanded_cluster_stats)
  for (i in 1:length(final_clusters)) {
    t_expanded_cluster_stats[,which(clust_assigns==i)] <- expand_from_cluster_to_cells_orig[i,] #change each row (really transposed column)
  }
  expanded_cluster_stats <- t(t_expanded_cluster_stats) #transpose back

  # Add in the others that didn't need to be expanded
  output_cluster_stats <- cbind(expanded_cluster_stats, combined_dist_from_mean_by_cell,
                                cutoff_values, clust_assigns)
  colnames(output_cluster_stats) <- c(expand_from_cluster_to_cells_params, "combined_dist_from_mean_by_cell",
                                      "cutoff_vals_by_cell", "consensus_clusters")

  # Set up global stats - this will be used for parameter selection workflow
  output_global_stats <- c("max_min_dist_btwn_means"=max(min_dist_btwn_means),
                           "max_sum_dists_btwn_means"=max(sum_dists_btwn_means),
                           "max_min_dist_btwn_means_single_param"=max(min_dist_btwn_means_single_param),
                           "max_c_dist_to_nearest_point"=max(c_dist_to_nearest_point),
                           "max_c_dist_to_nearest_point_trim"=max(c_dist_to_nearest_point_trim),
                           "max_c_sum_dists_to_nearest_points"=max(c_sum_dists_to_nearest_points),
                           "max_c_sum_dists_to_nearest_points_trim"=max(c_sum_dists_to_nearest_points_trim),
                           "max_c_dist_to_nearest_point_single_param"=max(c_dist_to_nearest_point_single_param),
                           "max_c_dist_to_nearest_point_trim_single_param"=max(c_dist_to_nearest_point_trim_single_param),
                           "mean_combined_dist_from_mean_by_cluster"=mean(combined_dist_from_mean_by_cluster),
                           "mean_combined_dist_from_mean_by_cell"=mean(combined_dist_from_mean_by_cell[which(!is.na(combined_dist_from_mean_by_cell))]),
                           "mean_combined_stdevs"=mean(combined_stdevs),
                           "mean_cutoff_vals_by_cluster"=mean(cutoff_vals_by_cluster),
                           "mean_cutoff_vals_by_cell"=mean(cutoff_values),
                           "pct_trimmed_by_cutoff"=length(which(cutoff_values < cutoff_trim_factor))/length(cutoff_values),
                           "mst_gap_distance"=largest_mst_gap_distance,
                           "num_clusters"=length(final_clusters))
  # Just get this ready
  output_param_names <- cbind(c_dist_to_nearest_point_single_param_name,
                              c_dist_to_nearest_point_trim_single_param_name,
                              min_dist_btwn_means_single_param_name)
  # Final wrapup for output
  output_list <- list(cluster_stats=output_cluster_stats,
                      global_stats=output_global_stats,
                      param_names=output_param_names,
                      mst_gap_params=output_mst_gap)

  return(output_list)
}
