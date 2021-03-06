% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/net_stat_partition.R
\name{net_stat_partition}
\alias{net_stat_partition}
\title{Partition a network and calculate configuration}
\usage{
net_stat_partition(net, directed = FALSE, configuration, tau = NULL,
  K = NULL, seed = 1234)
}
\arguments{
\item{net}{A network (adjacency matrix)}

\item{directed}{Default is FALSE (undirected network).}

\item{configuration}{Collection of network statistics (topological features)}

\item{tau}{Minimum size of a subgraph}

\item{K}{The number of subgraphs}

\item{seed}{Seed number. Default is 1234.}

\item{minsize}{Minimum size of a network for determining K and tau. Default is the size of the network net.}
}
\description{
Partition a network and calculate configuration
}
