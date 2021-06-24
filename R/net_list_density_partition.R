#' Calculate the mean subgraph density matrix after partitioning
#'
#' @param net A network (adjacency matrix)
#' @param minsize Minimum size of a network for determining K and tau. Default is the size of the network net. 
#' @param configuration Collection of network statistics (topological features)
#' @param directed Default is FALSE (undirected network).
#' @param tau Minimum size of a subgraph
#' @param K The number of subgraphs
#' @param seed Seed number. Default is 1234.
#' 
#' @export
#'
net_list_density_partition = function(netlist, directed = FALSE, configuration, tau = NULL, K = NULL, seed = 1234)
{
  # Density matrix
  out_list = lapply(netlist, function(x){
    net_stat_partition(net = x, directed = directed, configuration = configuration, 
                       tau = tau, K = K, seed = seed)})
  
  M = data.frame(do.call(rbind, out_list)) # N x p matrix
  
  colnames(M) = configuration
  rownames(M) = names(netlist)
  
  return(M)
}


