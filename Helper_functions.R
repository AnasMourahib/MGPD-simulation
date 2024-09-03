rejection_sampling<- function(x) (return( exp(max(x))/(sum(exp(x))) ) )

transformed_logsitic_generator<-function(i, sign, column, alpha,r){
  T<- -(log(rexp(length(sign),1)))/(1/alpha)+log((column*r)/(gamma(1-(alpha))))
  T[which(sign==i)]<- -(alpha)*(log(rgamma(1,shape=1-alpha,rate=1)))+log((column[which(sign==i)]*r)/(gamma(1-alpha)))
  return(T)
}

matrix_transformation <- function(Sigma) {
  gamma <- (as.matrix(rep(1, nrow(Sigma))) %*% diag(Sigma)) / 2 + (as.matrix(diag(Sigma)) %*% rep(1, nrow(Sigma))) / 2 - Sigma
  diag(gamma) <- diag(Sigma) / 2
  return(gamma)
}

calculus_of_eta <- function(i, mu, gamma) {
  eta <- (gamma[i, -i] / 2)^(1 / 2) + (mu[i] + diag(gamma)[i] - mu[-i] - diag(gamma)[-i]) / (2 * gamma[i, -i])^(1 / 2)
  return(eta)
}

calculus_of_R <- function(i, gamma) {
  R <- (as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1) + t(as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1)) - gamma[-i, -i]) / (2 * (as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1) * t(as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1)))^(1 / 2))
  diag(R) <- 1
  return(R)
}

case_d_2 <- function(mu, Gamma) {
  expec <- exp(mu[1] + diag(Gamma)[1]) * pnorm(((2 * Gamma[1, 2])^(1 / 2)) / (2) - ((mu[2] - mu[1] + (diag(Gamma)[2] - diag(Gamma)[1])) / ((2 * Gamma[1, 2])^(1 / 2))), mean = 0, sd = 1) + exp(mu[2] + diag(Gamma)[2]) * pnorm(((2 * Gamma[1, 2])^(1 / 2)) / (2) - ((mu[1] - mu[2] + (diag(Gamma)[1] - diag(Gamma)[2])) / ((2 * Gamma[1, 2])^(1 / 2))), mean = 0, sd = 1)
  return(expec)
}

calculus_d <- function(d, Sigma, mu) {
  if (d == 1) {
    return(exp(mu + (Sigma / 2)))
  } else {
    expectedvalue <- 0
    gamma <- matrix_transformation(Sigma)
    
    if (d == 2) {
      return(case_d_2(mu, gamma))
    } else {
      for (i in 1:d){
        # Computing eta
        eta <- calculus_of_eta(i, mu, gamma)
        # Computing R
        R <- calculus_of_R(i, gamma)
        expectedvalue <- expectedvalue + exp(mu[i] + gamma[i, i]) * pmvnorm(mean = rep(0, d - 1), corr = R, lower = rep(-Inf, d - 1), upper = eta)
      }
      return(expectedvalue)
    }
  }
}

mass_of_scenario <- function(d, r, list_of_k_matrices, A) # A is a (d x k) matrix and Sigma is a (dxd) matrix
{
  pi<-numeric(r)
  for (k in 1:r) {
    dimension <- length(which(A[, k] > 0))
    mu <- log(r * A[which(A[, k] > 0), k]) - (1 / 2) * diag(list_of_k_matrices[[k]])[which(A[, k] > 0)]
    Sigma <- list_of_k_matrices[[k]][which(A[, k] > 0), which(A[, k] > 0)]
    pi[k] <- calculus_d(dimension, Sigma, mu)
  }
  pi <- pi / sum(pi)
  return(pi)
}
 
HR_generator<-function(i,j, column, sign, covariance_matrix, r){
  T<-mvrnorm(1, mu = -(1 / 2) * diag(covariance_matrix)[sign] + log(r * column) + covariance_matrix[sign, i], Sigma = covariance_matrix[sign, sign])
  return(T)
} 



