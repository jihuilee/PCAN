#' PCA on a sample of networks and calculate PCA contribution of variables to explaining variabilities in PCs
#'
#' @param netlist A list of networks in a form or adjacency matrix
#' @param directed Default is FALSE (undirected network).
#' @param configuration Collection of subgraphs (topological features)
#' @param K The number of subgraphs
#' @param tau Minimum size of a subgraph
#' @param seed Seed number. Default is 1234.
#' @param subgroup Vector of subgroup membership. Default is NULL (i.e. no subgroup).
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom factoextra get_pca_var
#' @importFrom factoextra fviz_contrib
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 ggplot_build
#' @importFrom gridExtra arrangeGrob
#'
#' @export
#'

PCAN = function(netlist, directed = FALSE, configuration, subgroup = NULL)
{
  numdim = length(configuration)

  start = Sys.time()
  # Configuration density matrix
  M0 = net_list_density(netlist = netlist, directed = directed, configuration = configuration)

  # Standardize
  M = apply(M0, 2, function(x){(x-mean(x))/sd(x)})
  rownames(M0) = rownames(M) = names(netlist)

  # PCA
  PCA = prcomp(M, scale. = FALSE)
  end = Sys.time()

  # PCA variability
  variability = (100*PCA$sdev^2/sum(PCA$sdev^2))
  numdim = max(which(variability >= 1)) # >= 1% of variability
  variability = variability[1:numdim]

  # PCA contribution
  contrib0 = get_pca_var(PCA)$contrib[,1:min(numdim, length(PCA$sdev))]
  contribution = expand.grid(x = rownames(contrib0), y = colnames(contrib0))
  contribution$value = c(contrib0)

  # K-means if subgroup is NULL
  if(is.null(subgroup)){subgroup = kmeans(x = PCA$x[, 1:numdim], centers = 2)$cluster}

  # Plot1: PCA summary
  PLIST1 = PCAN_plot(PCA, subgroup)
  Plot1 = grid.arrange(grobs = PLIST1, nrow = 2)

  # Plot2 & Plot3: Plot each PC separately
  plotdat = data.frame(PCA$x[,1:numdim], Subgroup = subgroup)
  PLIST2 = PLIST3 = vector("list", numdim)
  for(p in 1:numdim)
  {
    PLIST2[[p]] = fviz_contrib(PCA, choice = "var", axes = p) +
      labs(x = "", y = "", title = "") +
      theme(plot.title = element_text(hjust = 0.5, size = 15), legend.position = "none",
            axis.text.x = element_text(size = 15, angle = 90), axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank())

    PLIST3[[p]] = ggplot(data = plotdat, aes(x = !!as.name(paste0("PC", p)), fill = Subgroup, color = Subgroup)) +
      geom_density(alpha = 0.3) + theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 15),
            legend.position = "bottom", legend.margin = margin(-25, 0, 0, 0), legend.text = element_text(size = 15),
            axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "", y = "", color = "", fill = "", title = "") +
      scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")
  }

  Plot2 = grid.arrange(grobs = PLIST2, nrow = 1)

  g_legend = function(a.gplot) {
    tmp = ggplot_gtable(ggplot_build(a.gplot))
    leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend = tmp$grobs[[leg]]
    return(legend)}

  mylegend = g_legend(PLIST3[[1]])

  PLIST3_2 = vector("list", numdim)
  for (i in 1:numdim) {PLIST3_2[[i]] = PLIST3[[i]] + theme(legend.position = "none")}

  Plot3 = grid.arrange(do.call("arrangeGrob", c(PLIST3_2, nrow = 1)), mylegend, heights = c(9/10, 1/10))

  return(list(M0 = M0, M = M, PCA = PCA,
              contribution = contribution, variability = variability, computingtime = end - start,
              PLIST1 = PLIST1, PLIST2 = PLIST2, PLIST3 = PLIST3,
              Plot1 = Plot1, Plot2 = Plot2, Plot3 = Plot3))
}

