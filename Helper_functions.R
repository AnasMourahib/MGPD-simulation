rejection_sampling<- function(x) (return( exp(max(x))/(sum(exp(x))) ) )

transformed_logsitic_generator<-function(i, vec, column, alpha,r){
  T<- -(log(rexp(length(vec),1)))/(1/alpha)+log((column*r)/(gamma(1-(alpha))))
  T[which(vec==i)]<- -(alpha)*(log(rgamma(1,shape=1-alpha,rate=1)))+log((column[which(vec==i)]*r)/(gamma(1-alpha)))
  return(T)
}

mass_of_scenario <- function(d, r, list_of_k_matrices, A) 
{
  mass <- numeric(r)
  
  # the element pi[j] will contain the mass of the j-th scenario
  for (k in 1:r) {
    dimension <- length(which(A[, k] > 0))
    mu <- log(k * A[which(A[, j] > 0), j]) - (1 / 2) * diag(list_of_k_matrices[[j]])[which(A[, j] > 0)]
    Sigma <- list_of_k_matrices[[j]][which(A[, j] > 0), which(A[, j] > 0)]
    pi[j] <- calculus_d(dimension, Sigma, mu)
  }
  
  
  
  pi <- pi / sum(pi)
  return(pi)
}
