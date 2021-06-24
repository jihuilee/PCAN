#' Plot PCA result
#'
#' @param pca PCA result
#'
#' @importFrom factoextra fviz_eig
#' @importFrom factoextra fviz_contrib
#' @importFrom factoextra fviz_pca_var
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_color_brewer
#' @importFrom ggplot2 expansion
#'
#' @export
#'

PCAN_plot = function(PCA, subgroup, geom = c("point", "text"))
{
  PLIST = vector("list", 4)

  PLIST[[1]] = fviz_eig(PCA, addlabels = TRUE) +
                labs(x = "Principal Components", title = "") +
                scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
                theme(plot.title = element_text(size = 15), legend.position = "none",
                      axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
                      axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  PLIST[[2]] = fviz_contrib(PCA, choice = "var", axes = 1:2) +
                labs(x = "", y = "Contribution (%)", title = "Contribution of variables to the first 2 PCs") +
                scale_y_continuous(expand = expansion(mult = c(0, 0.07))) +
                theme(plot.title = element_text(hjust = 0.5, size = 15), legend.position = "none",
                      axis.text.x = element_text(size = 15, angle = 90), axis.text.y = element_text(size = 15),
                      axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  PLIST[[3]] = fviz_pca_var(PCA, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)  +
                labs(color = "Contribution (%)", title = "") +
                theme(plot.title = element_text(size = 15, hjust = 0.5), legend.position = "bottom",
                      axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
                      axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  PLIST[[4]] = fviz_pca_ind(PCA, axes = c(1, 2),
                            geom = geom,
                            col.ind = subgroup,
                            repel = TRUE) +
                labs(color = "Quality", title = "") + scale_color_brewer(palette = "Dark2") +
                theme(plot.title = element_text(size = 15), legend.position = "none",
                      axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
                      axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                      panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(PLIST)
}
