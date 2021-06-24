#' Partition a network and calculate configuration
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

net_stat_partition = function(net, directed = FALSE, configuration, tau = NULL, K = NULL, seed = 1234)
{
  netsize = nrow(net)
  
  # Partition size
  subgraphsize = rep(tau, K) # Minimum size
  
  if(netsize-tau*K > 0)
  {
    set.seed(seed)
    subgraphsize = subgraphsize + rmultinom(1, size = netsize-tau*K, prob = rep(1/K, K))
  }
  
  # Partitioning
  nodenum = 1:netsize
  set.seed(seed)
  subgraphlabel = sample(rep(1:K, subgraphsize))
  subgraphlist = vector("list", K)
  for(k in 1:K){
    seq.k = which(subgraphlabel == k)
    subgraphlist[[k]] = net[seq.k, seq.k]
  }
  
  # Calculate density
  out_list = lapply(subgraphlist, function(x){net_stat(net = x, directed = directed, configuration = configuration)})
  
  out0 = data.frame(do.call(rbind, out_list)) # K x p matrix
  
  # Average over subnetworks: p x 1 matrix
  out = matrix(apply(out0, 2, mean), nrow = 1)
  colnames(out) = colnames(out0)
  
  return(out)
}

