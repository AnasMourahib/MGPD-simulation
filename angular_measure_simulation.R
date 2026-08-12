# ==============================================================================
# Simulation of a random vector W \sim d^{-1} H, where H is the angular measure associated to the L1-norm
# Formula:
#   J ~ Uniform{1, ..., d}
#   Given J = j, W = (Y / ||Y||_1 | Y_j > 1)
#   Result: W ~ d^(-1) * H
# ==============================================================================

#' Sample the random vector W according to the conditional distribution
#'
#' @param n Number of samples to generate
#' @param d Dimension of the vector
#' @param rY Function to generate random vector Y (default: Exponential(1))
#' @return A matrix of size (n x d) where each row is a realization of W
sample_W <- function(d, r , Sigma = NULL , alpha = NULL ,  A , model = c("HR", "logistic")) {
  
  # Matrix to store output vectors W (n samples of length d)
  W <- rep(0 , d)
  
  # --------------------------------------------------------------------------
  # Step 1: Sample J uniformly from {1, 2, ..., d}
  # --------------------------------------------------------------------------
  j <- sample(1:d, size = 1)
  
  # --------------------------------------------------------------------------
  # Step 2: Generate Y conditioned on Y_j > 1 (Acceptance-Rejection)
  # --------------------------------------------------------------------------
  repeat {
    if(model == "HR"){
      Z <- mgpd_simulation_mixture_HR(d,r,Sigma,A) 
    }
    if(model == "logistic"){
      Z <- mgpd_simulation_mixture_logistic(d,r,alpha,A)
    }
    if (Z[j] > 0) {
      Y <- exp(Z)
      break  # Accept Y if the j-th component exceeds 1
    }
  }
  
  # --------------------------------------------------------------------------
  # Step 3: Compute L1 norm and normalize to get W on the simplex
  # W = Y / ||Y||_1
  # --------------------------------------------------------------------------
  L1_norm <- sum(Y)
  W <- Y / L1_norm
  return(W)
}
