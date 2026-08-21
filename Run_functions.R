library(MASS)
library(mvtnorm)
source("mgpd_simulation_mixture_logistic.R")
source("mgpd_simulation_mixture_HR.R")
source("Helper_functions.R")
source("angular_measure_simulation.R")

d<-3
r<-3
A<-rbind(c(1,0,0), c(1/2, 1/2, 0), c(1/3, 1/3, 1/3))

#Sample of size 100 from the multivariate generalized Pareto distribution associated to a mixture logistic model with matrix A and alpha=(0.5, 0.5, 0.5)
alpha<-0.5 
sample_mixture_logistic<-function(d,r,alpha,A,N){
  final<-replicate(N,mgpd_simulation_mixture_logistic(d,r,rep(alpha, r),A))
  return(final)
}
N<-100
set.seed(79)
Y_mix_log<-sample_mixture_logistic(d,r,alpha,A,N)
Z_mix_log<-4*(exp(Y_mix_log/4)-1)

new<-rep(0,100)
new[which(Z_mix_log[1, ] > -4 | Z_mix_log[2, ] > -4 | Z_mix_log[3, ] > -4)]<-1
new[which(Z_mix_log[1, ] == -4 & Z_mix_log[2, ] > -4 & Z_mix_log[3,]>-4)]<-2
new[which(Z_mix_log[1, ] == -4 & Z_mix_log[2, ] == -4 & Z_mix_log[3,]>-4)]<-3

Z_mix_log<-as.data.frame(t(Z_mix_log))
Z_mix_log$new<-new

# Define the colors for the points
my_cols <- c("#00AFBB", "#E7B800", "#FC4E07")

# Set up the layout with an extra space for the legend
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1))

# Adjust the margin and distance between the axis and labels
par(mgp = c(4, 1, 0))  # Adjust distance: mgp[1] is axis title, mgp[2] is axis labels, mgp[3] is axis line

# Plot 1: Z2 vs Z1
plot(Z_mix_log[, 2]~Z_mix_log[, 1], 
     pch = c(1,2,3)[Z_mix_log$new], 
     col = my_cols[Z_mix_log$new], 
     cex = 1.5, 
     xlim = c(-4, 7), 
     ylim = c(-4, 7), 
     ylab = "", 
     xlab = "")
title(xlab = expression(Z[1]), ylab = expression(Z[2]), line = 2.3, cex.lab = 1.8)

# Plot 2: Z3 vs Z1
plot(Z_mix_log[, 3]~Z_mix_log[, 1], 
     pch = c(1,2,3)[Z_mix_log$new], 
     col = my_cols[Z_mix_log$new], 
     cex = 1.5, 
     xlim = c(-4, 7), 
     ylim = c(-4, 7), 
     ylab = "", 
     xlab = "")
title(xlab = expression(Z[1]), ylab = expression(Z[3]), line = 2.3, cex.lab = 1.8)

# Plot 3: Z3 vs Z2
plot(Z_mix_log[, 3]~Z_mix_log[, 2], 
     pch = c(1,2,3)[Z_mix_log$new], 
     col = my_cols[Z_mix_log$new], 
     cex = 1.5, 
     xlim = c(-4, 7), 
     ylim = c(-4, 7), 
     ylab = "", 
     xlab = "")
title(xlab = expression(Z[2]), ylab = expression(Z[3]), line = 2.3, cex.lab = 1.8)

# Create a new plot space for the legend
plot.new()

# Draw the legend in the center of the plotting area with a box around it
legend("center", 
       pch = c(1,2,3), 
       col = my_cols, 
       legend = c(expression(A["{1,2,3}"]^{-4}), expression(A["{2,3}"]^{-4}), expression(A["{3}"]^{-4})), 
       cex = 2, 
       box.col = "black", 
       box.lty = "solid")


#Sample of size 100 from the multivariate generalized Pareto distribution associated to a mixture Hüsler-Reiss model with matrix A and variogram matrix Sigma on each column

Sigma<-  rbind(c(1.6, 10 / 11, 10 / 11), c(10 / 11, 1.6, 10 / 11), c(10 / 11, 10 / 11, 1.6) )
Sigma <- list (Sigma, Sigma, Sigma)


sample_mixture_HR<-function(d,r,Sigma,A,N){
  final<-replicate(N,mgpd_simulation_mixture_HR(d,r,Sigma,A))
  return(final)
}
set.seed(7)
Y_mix_HR<-sample_mixture_HR(d,r,Sigma,A,N)
Z_mix_HR<-4*(exp(Y_mix_HR/4)-1)


new<-rep(0,100)
new[which(Z_mix_HR[1, ] > -4 | Z_mix_HR[2, ] > -4 | Z_mix_HR[3, ] > -4)]<-1
new[which(Z_mix_HR[1, ] == -4 & Z_mix_HR[2, ] > -4 & Z_mix_HR[3,]>-4)]<-2
new[which(Z_mix_HR[1, ] == -4 & Z_mix_HR[2, ] == -4 & Z_mix_HR[3,]>-4)]<-3


Z_mix_HR<-as.data.frame(t(Z_mix_HR))
Z_mix_HR$new<-new


# Define the colors for the points
my_cols <- c("#00AFBB", "#E7B800", "#FC4E07")

# Set up the layout with an extra space for the legend
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), widths = c(1, 1), heights = c(1, 1))

# Adjust the margin and distance between the axis and labels
par(mgp = c(4, 1, 0))  # Adjust distance: mgp[1] is axis title, mgp[2] is axis labels, mgp[3] is axis line

# Plot 1: Z2 vs Z1
plot(Z_mix_HR[, 2]~Z_mix_HR[, 1], 
     pch = c(1,2,3)[Z_mix_HR$new], 
     col = my_cols[Z_mix_HR$new], 
     cex = 1.5, 
     xlim = c(-4, 7), 
     ylim = c(-4, 7), 
     ylab = "", 
     xlab = "")
title(xlab = expression(Z[1]), ylab = expression(Z[2]), line = 2.3, cex.lab = 1.8)

# Plot 2: Z3 vs Z1
plot(Z_mix_HR[, 3]~Z_mix_HR[, 1], 
     pch = c(1,2,3)[Z_mix_HR$new], 
     col = my_cols[Z_mix_HR$new], 
     cex = 1.5, 
     xlim = c(-4, 7), 
     ylim = c(-4, 7), 
     ylab = "", 
     xlab = "")
title(xlab = expression(Z[1]), ylab = expression(Z[3]), line = 2.3, cex.lab = 1.8)

# Plot 3: Z3 vs Z2
plot(Z_mix_HR[, 3]~Z_mix_HR[, 2], 
     pch = c(1,2,3)[Z_mix_HR$new], 
     col = my_cols[Z_mix_HR$new], 
     cex = 1.5, 
     xlim = c(-4, 7), 
     ylim = c(-4, 7), 
     ylab = "", 
     xlab = "")
title(xlab = expression(Z[2]), ylab = expression(Z[3]), line = 2.3, cex.lab = 1.8)

# Create a new plot space for the legend
plot.new()

# Draw the legend in the center of the plotting area with a box around it
legend("center", 
       pch = c(1,2,3), 
       col = my_cols, 
       legend = c(expression(A["{1,2,3}"]^{-4}), expression(A["{2,3}"]^{-4}), expression(A["{3}"]^{-4})), 
       cex = 2, 
       box.col = "black", 
       box.lty = "solid")


########Simulation from the L1-angular measure
A <- rbind(c(1/3 , 0 , 1/3 , 1/3),
           c(1/2 , 0 , 1/2 , 0), 
           c(1/2 , 1/2 , 0 , 0 ),
           c(1/2 , 1/2 , 0 , 0 ) )
Sigma<- rho <- 0.5 
d <- 4
Sigma <- matrix(rho, nrow = d, ncol = d)
diag(Sigma) <- 1

sample_angular_measure_mixture_HR<-function(d,r,Sigma,A,N){
  final<-replicate(N,sample_W(d , r , Sigma , alpha = NULL ,  A , model = "HR"))
  return(final)
}

###Simulation on the clique C1
A_C1 = A[1:3 ,  ]   
A_C1 <- A_C1[, colSums(A_C1 != 0) > 0, drop = FALSE]
Sigma_C1 <- Sigma[1:3 , 1:3]
list_Sigma_C1 <- list (Sigma_C1, Sigma_C1, Sigma_C1 , Sigma_C1)
d <- nrow(A_C1)
r <- ncol(A_C1)


set.seed(7)
N <- 100
W1 <- t(sample_angular_measure_mixture_HR(d,r,list_Sigma_C1,A_C1,N)) 
###Simulation on the clique C2
A_C2 = A[2:4 ,  ]   
A_C2 <- A_C2[, colSums(A_C2 != 0) > 0, drop = FALSE]
Sigma_C2 <- Sigma[2:4 , 2:4]
list_Sigma_C2 <- list (Sigma_C2, Sigma_C2, Sigma_C2 , Sigma_C2)
d <- nrow(A_C2)
r <- ncol(A_C2)



set.seed(7)
N <- 100
W2 <- t(sample_angular_measure_mixture_HR(d,r,list_Sigma_C2,A_C2,N)) 


#########Plot the angular measure in a simplex


###Vertical plot
plot_vertical_simplex <- function(
    W1,
    W2,
    point_size = 4.0,
    point_color = "#003366"
) {
  # Check the inputs
  if (!is.matrix(W1) || ncol(W1) != 3) {
    stop("Input 'W1' must be an n x 3 matrix.")
  }
  
  if (!is.matrix(W2) || ncol(W2) != 3) {
    stop("Input 'W2' must be an n x 3 matrix.")
  }
  
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop(
      "Package 'ggtext' is required. Install it with ",
      "install.packages('ggtext')."
    )
  }
  
  # Height of an equilateral triangle
  h <- sqrt(3) / 2
  
  # Vertical position of the upper triangle
  upper_base <- 1.25
  
  # ============================================================
  # Upper triangle: components 1, 2, 3
  #
  # Component 1: top
  # Component 2: bottom-left
  # Component 3: bottom-right
  #
  # Columns of W1 are ordered as 1, 2, 3
  # ============================================================
  
  pts_upper <- data.frame(
    x1 = W1[, 1],
    x2 = W1[, 2],
    x3 = W1[, 3]
  )
  
  # Barycentric transformation
  pts_upper$x_2d <- 0.5 * pts_upper$x1 + pts_upper$x3
  pts_upper$y_2d <- upper_base + h * pts_upper$x1
  
  # Upper triangle boundary
  triangle_upper <- data.frame(
    x = c(0.5, 0, 1, 0.5),
    y = c(
      upper_base + h,
      upper_base,
      upper_base,
      upper_base + h
    )
  )
  
  # Interior point for the upper triangle
  interior_upper_x <- 0.5
  interior_upper_y <- upper_base + h / 3
  
  # ============================================================
  # Lower triangle: components 2, 3, 4
  #
  # Component 2: top-left
  # Component 3: top-right
  # Component 4: bottom
  #
  # Columns of W2 are ordered as 2, 3, 4
  # ============================================================
  
  pts_lower <- data.frame(
    x2 = W2[, 1],
    x3 = W2[, 2],
    x4 = W2[, 3]
  )
  
  # Barycentric transformation:
  #
  # x2 -> (0, h)
  # x3 -> (1, h)
  # x4 -> (0.5, 0)
  pts_lower$x_2d <- pts_lower$x3 + 0.5 * pts_lower$x4
  pts_lower$y_2d <- h * (pts_lower$x2 + pts_lower$x3)
  
  # Lower triangle boundary
  triangle_lower <- data.frame(
    x = c(0, 1, 0.5, 0),
    y = c(h, h, 0, h)
  )
  
  # Interior point for the lower triangle
  interior_lower_x <- 0.5
  interior_lower_y <- 2 * h / 3
  
  # ============================================================
  # Double-struck A corresponding to LaTeX \mathbb{A}
  # ============================================================
  
  bb_A <- intToUtf8(0x1D538)
  
  # Label sizes
  vertex_label_size <- 4.7
  edge_label_size <- 4.2
  interior_label_size <- 4.7
  
  # ============================================================
  # Labels for the upper triangle
  # ============================================================
  
  labels_upper <- data.frame(
    x = c(
      0.50,
      -0.04,
      1.04,
      0.17,
      0.79,
      0.50,
      0.94
    ),
    y = c(
      upper_base + h + 0.04,
      upper_base - 0.03,
      upper_base - 0.03,
      upper_base + h / 2,
      upper_base + h / 2,
      upper_base - 0.04,
      upper_base + 0.55
    ),
    label = paste0(
      "<span style='font-size:14pt;'>",
      bb_A,
      "</span>",
      c(
        "<sub>{1}</sub>",
        "<sub>{2}</sub>",
        "<sub>{3}</sub>",
        "<sub>{1,2}</sub>",
        "<sub>{1,3}</sub>",
        "<sub>{2,3}</sub>",
        "<sub>{1,2,3}</sub>"
      )
    ),
    hjust = c(
      0.5,
      1,
      0,
      1,
      0,
      0.5,
      0
    ),
    vjust = c(
      0,
      1,
      1,
      0.5,
      0.5,
      1,
      0.5
    ),
    size = c(
      vertex_label_size,
      vertex_label_size,
      vertex_label_size,
      edge_label_size,
      edge_label_size,
      edge_label_size,
      interior_label_size
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================
  # Labels for the lower triangle
  #
  # Component 2 is at the top-left
  # Component 3 is at the top-right
  # Component 4 is at the bottom
  # ============================================================
  
  labels_lower <- data.frame(
    x = c(
      -0.04,
      1.04,
      0.50,
      0.17,
      0.83,
      0.50,
      0.94
    ),
    y = c(
      h + 0.03,
      h + 0.03,
      -0.05,
      h / 2,
      h / 2,
      h + 0.035,
      0.31
    ),
    label = paste0(
      "<span style='font-size:14pt;'>",
      bb_A,
      "</span>",
      c(
        "<sub>{2}</sub>",
        "<sub>{3}</sub>",
        "<sub>{4}</sub>",
        "<sub>{2,4}</sub>",
        "<sub>{3,4}</sub>",
        "<sub>{2,3}</sub>",
        "<sub>{2,3,4}</sub>"
      )
    ),
    hjust = c(
      1,
      0,
      0.5,
      1,
      0,
      0.5,
      0
    ),
    vjust = c(
      0,
      0,
      1,
      0.5,
      0.5,
      0,
      0.5
    ),
    size = c(
      vertex_label_size,
      vertex_label_size,
      vertex_label_size,
      edge_label_size,
      edge_label_size,
      edge_label_size,
      interior_label_size
    ),
    stringsAsFactors = FALSE
  )
  
  # Combine all labels
  labels <- rbind(
    labels_upper,
    labels_lower
  )
  
  # ============================================================
  # Construct the plot
  # ============================================================
  
  p <- ggplot2::ggplot() +
    
    # Upper triangle
    ggplot2::geom_polygon(
      data = triangle_upper,
      mapping = ggplot2::aes(
        x = x,
        y = y
      ),
      fill = "white",
      color = "black",
      linewidth = 0.6
    ) +
    
    # Simulations W1
    ggplot2::geom_point(
      data = pts_upper,
      mapping = ggplot2::aes(
        x = x_2d,
        y = y_2d
      ),
      color = point_color,
      size = point_size,
      alpha = 0.85
    ) +
    
    # Lower triangle
    ggplot2::geom_polygon(
      data = triangle_lower,
      mapping = ggplot2::aes(
        x = x,
        y = y
      ),
      fill = "white",
      color = "black",
      linewidth = 0.6
    ) +
    
    # Simulations W2
    ggplot2::geom_point(
      data = pts_lower,
      mapping = ggplot2::aes(
        x = x_2d,
        y = y_2d
      ),
      color = point_color,
      size = point_size,
      alpha = 0.85
    ) +
    
    # Upper interior marker
    ggplot2::annotate(
      geom = "point",
      x = interior_upper_x,
      y = interior_upper_y,
      size = 1.5,
      color = "grey30"
    ) +
    
    # Pointer to A_{1,2,3}
    ggplot2::annotate(
      geom = "segment",
      x = interior_upper_x,
      y = interior_upper_y,
      xend = 0.92,
      yend = upper_base + 0.54,
      color = "grey30",
      linewidth = 0.4
    ) +
    
    # Lower interior marker
    ggplot2::annotate(
      geom = "point",
      x = interior_lower_x,
      y = interior_lower_y,
      size = 1.5,
      color = "grey30"
    ) +
    
    # Pointer to A_{2,3,4}
    ggplot2::annotate(
      geom = "segment",
      x = interior_lower_x,
      y = interior_lower_y,
      xend = 0.92,
      yend = 0.32,
      color = "grey30",
      linewidth = 0.4
    ) +
    
    # Double-struck A labels
    ggtext::geom_richtext(
      data = labels,
      mapping = ggplot2::aes(
        x = x,
        y = y,
        label = label,
        hjust = hjust,
        vjust = vjust,
        size = size
      ),
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(0, "pt"),
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    
    ggplot2::scale_size_identity() +
    
    # Preserve triangular geometry
    ggplot2::coord_fixed(
      xlim = c(-0.20, 1.20),
      ylim = c(
        -0.13,
        upper_base + h + 0.14
      ),
      clip = "off"
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  return(p)
}


p <- plot_vertical_simplex(
  W1 = W1,
  W2 = W2
)

print(p)


#########Horizental plot

plot_horizontal_simplex <- function(
    W1,
    W2,
    point_size = 4.0,
    point_color = "#003366"
) {
  # Check the inputs
  if (!is.matrix(W1) || ncol(W1) != 3) {
    stop("Input 'W1' must be an n x 3 matrix.")
  }
  
  if (!is.matrix(W2) || ncol(W2) != 3) {
    stop("Input 'W2' must be an n x 3 matrix.")
  }
  
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop(
      "Package 'ggtext' is required. Install it with ",
      "install.packages('ggtext')."
    )
  }
  
  # Horizontal height of an equilateral triangle
  h <- sqrt(3) / 2
  
  # Horizontal shift of the right triangle
  right_shift <- 1.18
  
  # ============================================================
  # Left triangle: components 1, 2, 3
  #
  # Component 1: left apex
  # Component 2: top-right
  # Component 3: bottom-right
  #
  # Columns of W1 are ordered as 1, 2, 3
  # ============================================================
  
  pts_left <- data.frame(
    x1 = W1[, 1],
    x2 = W1[, 2],
    x3 = W1[, 3]
  )
  
  # Barycentric transformation:
  # x1 -> (0, 0.5)
  # x2 -> (h, 1)
  # x3 -> (h, 0)
  pts_left$x_2d <- h * (pts_left$x2 + pts_left$x3)
  pts_left$y_2d <- 0.5 * pts_left$x1 + pts_left$x2
  
  triangle_left <- data.frame(
    x = c(0, h, h, 0),
    y = c(0.5, 1, 0, 0.5)
  )
  
  interior_left_x <- 2 * h / 3
  interior_left_y <- 0.5
  
  # ============================================================
  # Right triangle: components 2, 3, 4
  #
  # Component 2: top-left
  # Component 3: bottom-left
  # Component 4: right apex
  #
  # Columns of W2 are ordered as 2, 3, 4
  # ============================================================
  
  pts_right <- data.frame(
    x2 = W2[, 1],
    x3 = W2[, 2],
    x4 = W2[, 3]
  )
  
  # Barycentric transformation:
  # x2 -> (right_shift, 1)
  # x3 -> (right_shift, 0)
  # x4 -> (right_shift + h, 0.5)
  pts_right$x_2d <- right_shift + h * pts_right$x4
  pts_right$y_2d <- pts_right$x2 + 0.5 * pts_right$x4
  
  triangle_right <- data.frame(
    x = c(
      right_shift,
      right_shift,
      right_shift + h,
      right_shift
    ),
    y = c(1, 0, 0.5, 1)
  )
  
  interior_right_x <- right_shift + h / 3
  interior_right_y <- 0.5
  
  # ============================================================
  # Double-struck A corresponding to LaTeX \mathbb{A}
  # ============================================================
  
  bb_A <- intToUtf8(0x1D538)
  
  # Label sizes
  vertex_label_size <- 4.7
  edge_label_size <- 4.2
  interior_label_size <- 4.7
  
  # ============================================================
  # Labels for the left triangle
  # ============================================================
  
  labels_left <- data.frame(
    x = c(
      -0.05,
      h + 0.04,
      h + 0.04,
      h / 2 - 0.03,
      h / 2 - 0.03,
      h + 0.04,
      0.43
    ),
    y = c(
      0.50,
      1.03,
      -0.03,
      0.78,
      0.22,
      0.50,
      1.08
    ),
    label = paste0(
      "<span style='font-size:14pt;'>",
      bb_A,
      "</span>",
      c(
        "<sub>{1}</sub>",
        "<sub>{2}</sub>",
        "<sub>{3}</sub>",
        "<sub>{1,2}</sub>",
        "<sub>{1,3}</sub>",
        "<sub>{2,3}</sub>",
        "<sub>{1,2,3}</sub>"
      )
    ),
    hjust = c(
      1,
      0,
      0,
      1,
      1,
      0,
      0.5
    ),
    vjust = c(
      0.5,
      0,
      1,
      0.5,
      0.5,
      0.5,
      0
    ),
    size = c(
      vertex_label_size,
      vertex_label_size,
      vertex_label_size,
      edge_label_size,
      edge_label_size,
      edge_label_size,
      interior_label_size
    ),
    stringsAsFactors = FALSE
  )
  
  # ============================================================
  # Labels for the right triangle
  # ============================================================
  
  labels_right <- data.frame(
    x = c(
      right_shift - 0.04,
      right_shift - 0.04,
      right_shift + h + 0.05,
      right_shift - 0.04,
      right_shift + h / 2 + 0.03,
      right_shift + h / 2 + 0.03,
      right_shift + 0.43
    ),
    y = c(
      1.03,
      -0.03,
      0.50,
      0.50,
      0.78,
      0.22,
      1.08
    ),
    label = paste0(
      "<span style='font-size:14pt;'>",
      bb_A,
      "</span>",
      c(
        "<sub>{2}</sub>",
        "<sub>{3}</sub>",
        "<sub>{4}</sub>",
        "<sub>{2,3}</sub>",
        "<sub>{2,4}</sub>",
        "<sub>{3,4}</sub>",
        "<sub>{2,3,4}</sub>"
      )
    ),
    hjust = c(
      1,
      1,
      0,
      1,
      0,
      0,
      0.5
    ),
    vjust = c(
      0,
      1,
      0.5,
      0.5,
      0.5,
      0.5,
      0
    ),
    size = c(
      vertex_label_size,
      vertex_label_size,
      vertex_label_size,
      edge_label_size,
      edge_label_size,
      edge_label_size,
      interior_label_size
    ),
    stringsAsFactors = FALSE
  )
  
  labels <- rbind(
    labels_left,
    labels_right
  )
  
  # ============================================================
  # Construct the plot
  # ============================================================
  
  p <- ggplot2::ggplot() +
    
    # Left triangle
    ggplot2::geom_polygon(
      data = triangle_left,
      mapping = ggplot2::aes(x = x, y = y),
      fill = "white",
      color = "black",
      linewidth = 0.6
    ) +
    
    # Simulations W1
    ggplot2::geom_point(
      data = pts_left,
      mapping = ggplot2::aes(x = x_2d, y = y_2d),
      color = point_color,
      size = point_size,
      alpha = 0.85
    ) +
    
    # Right triangle
    ggplot2::geom_polygon(
      data = triangle_right,
      mapping = ggplot2::aes(x = x, y = y),
      fill = "white",
      color = "black",
      linewidth = 0.6
    ) +
    
    # Simulations W2
    ggplot2::geom_point(
      data = pts_right,
      mapping = ggplot2::aes(x = x_2d, y = y_2d),
      color = point_color,
      size = point_size,
      alpha = 0.85
    ) +
    
    # Interior marker for the left triangle
    ggplot2::annotate(
      geom = "point",
      x = interior_left_x,
      y = interior_left_y,
      size = 1.5,
      color = "grey30"
    ) +
    
    # Pointer to A_{1,2,3}
    ggplot2::annotate(
      geom = "segment",
      x = interior_left_x,
      y = interior_left_y,
      xend = 0.43,
      yend = 1.04,
      color = "grey30",
      linewidth = 0.4
    ) +
    
    # Interior marker for the right triangle
    ggplot2::annotate(
      geom = "point",
      x = interior_right_x,
      y = interior_right_y,
      size = 1.5,
      color = "grey30"
    ) +
    
    # Pointer to A_{2,3,4}
    ggplot2::annotate(
      geom = "segment",
      x = interior_right_x,
      y = interior_right_y,
      xend = right_shift + 0.43,
      yend = 1.04,
      color = "grey30",
      linewidth = 0.4
    ) +
    
    # Double-struck A labels
    ggtext::geom_richtext(
      data = labels,
      mapping = ggplot2::aes(
        x = x,
        y = y,
        label = label,
        hjust = hjust,
        vjust = vjust,
        size = size
      ),
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(0, "pt"),
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    
    ggplot2::scale_size_identity() +
    
    # Preserve the simplex geometry
    ggplot2::coord_fixed(
      xlim = c(-0.22, right_shift + h + 0.30),
      ylim = c(-0.13, 1.18),
      clip = "off"
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  return(p)
}

p <- plot_horizontal_simplex(
       W1 = W1,
       W2 = W2
   )
print(p)
