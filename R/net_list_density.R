#' Calculate the mean subgraph density matrix
#'
#' @param netlist A list of networks in a form or adjacency matrix
#' @param directed Default is FALSE (undirected network).
#' @param configuration Collection of network statistics (topological features)
#' 
#' @export
#'

net_list_density = function(netlist, directed = FALSE, configuration)
{
  # Density matrix
  out_list = lapply(netlist, function(x){net_stat(net = x, directed = directed, configuration = configuration)})
  
  M = data.frame(do.call(rbind, out_list)) # N x p matrix
 
  colnames(M) = configuration
  rownames(M) = names(netlist)
  
  return(M)
}
