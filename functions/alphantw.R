#modified by I. Bartomeus from Will Petry public shiny app (https://github.com/wpetry/StructuralCoexistence)

library(igraph)
mk_graph_3sp <- function(alphamat, rs, title = "Species interaction network", 
                         unname_ = FALSE, superunname = FALSE, line_col = "black"){
    if(unname_){
        alphamat <- unname(alphamat)
    }
    g <- igraph::graph_from_adjacency_matrix(alphamat != 0)
    E(g)$weight <- as.numeric(alphamat)
    widths <- E(g)$weight * 5
    #widths[widths > 1] <- sqrt(widths)
    if(superunname){
        plot(g,
             main = title,
             margin = c(0, -0.15, -0.3, -0.15),
             xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
             vertex.label = " ",
             vertex.size = 50 * rs,
             vertex.color = "grey80",
             vertex.frame.color = "transparent",
             edge.curved = TRUE,
             edge.width = abs(widths),
             edge.arrow.size = 0.7,
             edge.arrow.mode = c(0, 2, 2,
                                 2, 0, 2,
                                 2, 2, 0),
             edge.color = line_col,
             edge.lty = ifelse(widths > 0, 1, 3),
             edge.loop.angle = 0.75,
             layout = matrix(c(4, 0, 0, 0, 2, sqrt(3)/2), ncol = 2,
                             byrow = TRUE))
        }  else {  
    plot(g,
         main = title,
         margin = c(0, -0.15, -0.3, -0.15),
         xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
         vertex.label.cex = 2,
         vertex.label.color = "black",
         vertex.size = 50 * rs,
         vertex.color = "grey80",
         vertex.frame.color = "transparent",
         edge.curved = TRUE,
         edge.width = abs(widths),
         edge.arrow.size = 0.7,
         edge.arrow.mode = c(0, 2, 2,
                             2, 0, 2,
                             2, 2, 0),
         edge.color = line_col,
         edge.lty = ifelse(widths > 0, 1, 3),
         edge.loop.angle = 0.75,
         layout = matrix(c(4, 0, 0, 0, 2, sqrt(3)/2), ncol = 2,
                         byrow = TRUE))
        }
}