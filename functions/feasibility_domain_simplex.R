#Function to calculate and plot a simplex (AUthor S. Saavedra)

if(!require(ggtern)) {install.packages("ggtern"); library(ggtern)}

# function to sample from simplex
simplex_sampling <- function(m, n) {
  r <- list()
  for (j in 1:m) {
    dist <- c(sort(runif(n-1, 0, 1)), 1)
    r[[j]] <- c(dist[1], diff(dist))
  }
  return(r)
}

simplex_plot <- function(alpha, r, color = "red", line_col = "black"){
  # sample 3 vectors (these are the column vectors of A, which is transposed here)
  A_vectors <- t(alpha) #simplex_sampling(3, 3)  alpha vs alpha_link (and alpha_no_link)
  A <- as.data.frame(matrix(unlist(A_vectors), nrow = 3, ncol = 3, byrow = TRUE))
  #A <- as.data.frame(alpha)
  R <- as.data.frame(t(r/(sum(r)))) # r vs r_link
  colnames(R) <- c('V1','V2','V3')
  # plot feasibility domain on simplex with r vector
  ggtern() + 
    geom_polygon_closed(data = A, aes(x = V2, y = V3, z = V1), 
                        color = line_col, fill = "gray50", alpha = 0.5) +
    geom_point(data = R, aes(x = V2, y = V3, z = V1), size = 2.5,
               color = color) +
    labs(x = "Raphanus", y = "Vicea", z = "Solanum") +
    theme_classic() + 
    theme_nolabels() +
    theme_noticks() +
    theme(axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),)
}
  

