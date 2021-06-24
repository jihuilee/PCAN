#' Calculate network statistics
#'
#' @param net Input network (adjacency matrix)
#' @param stat Statistics
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph distances
#' @importFrom igraph transitivity
#' @importFrom network network
#' @importFrom ergm summary_formula
#'
#' @export

net_stat = function(net, directed = FALSE,
                    configuration = c("edges", "isolate", 
                                  "kstar(2)", "kstar(3)", "kstar(4)", "kstar(5)", 
                                  "cycle(4)", "cycle(5)", "cycle(6)"))
{
  if(directed) {stop("We have not implemented the fuction of calculating network statistics for directed network.")}

#  net0 = matrix(1, nrow = nrow(net), ncol = ncol(net)); diag(net0) = 0
  
  if(!(directed)){net[lower.tri(net)] = 0}
  
  N = nrow(net)
  
  # edges
  edges = NULL
  if("edges" %in% configuration){
    edgeout =  sum(net) / choose(N, 2)
    names(edgeout) = "edges"
  }
  
  # isolate
  isolateout = NULL
  if("isolate" %in% configuration){
    isolateout = sum(apply(net, 1, sum) + apply(net, 2, sum) == 0) / N
    names(isolateout) = "isolate"
  }

  # kstar
  kstarout = NULL
  if(sum(grepl("kstar", configuration)) > 0) {
    kstarstat = configuration[grepl("kstar", configuration)]
    kstarnum = as.numeric(gsub(")", "", gsub("kstar\\(", "", kstarstat)))
    
    # ERGM network statistics
    net.ergm = network(net, directed = directed)
    obj = formula(paste("net.ergm", paste(kstarstat, collapse=" + "), sep=" ~ "))

    kstarout = summary_formula(obj) / (choose(N, kstarnum + 1) * (kstarnum + 1))
    names(kstarout) = kstarstat
  }

  # triangle
  triangleout = NULL
  if("triangle" %in% configuration){
    # ERGM network statistics
    net.ergm = network(net, directed = directed)
    obj = formula("net.ergm ~ triangle")

    triangleout = summary_formula(obj) / choose(N, 3)
    names(triangleout) = "triangle"
  }
  
  # cycles
  cycleout = NULL
  if(sum(grepl("cycle", configuration)) > 0) {
    cyclestat = configuration[grepl("cycle", configuration)]
    cyclenum = as.numeric(gsub(")", "", gsub("cycle\\(", "", cyclestat)))
    
    # ERGM network statistics
    net.ergm = network(net, directed = directed)
    obj = formula(paste("net.ergm", paste(cyclestat, collapse=" + "), sep=" ~ "))
    
    fullcycle = choose(N, cyclenum) * factorial(cyclenum - 1); if(!directed){fullcycle = fullcycle/2}
    cycleout = summary_formula(obj) / fullcycle
    names(cycleout) = cyclestat
  }
  
  out = c(edgeout, isolateout, kstarout, triangleout, cycleout)
  return(out)
}
