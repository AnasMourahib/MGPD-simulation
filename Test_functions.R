source("mgpd_simulation_mixture_logistic.R")
source("mgpd_simulation_mixture_HR.R")


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














