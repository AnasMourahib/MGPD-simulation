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
X_mix_log<-sample_mixture_logistic(d,r,alpha,A,N)


#Sample of size 100 from the multivariate generalized Pareto distribution associated to a mixture Hüsler-Reiss model with matrix A and variogram matrix Sigma on each column

Sigma<-  rbind(c(1.6, 10 / 11, 10 / 11), c(10 / 11, 1.6, 10 / 11), c(10 / 11, 10 / 11, 1.6) )
Sigma <- list (Sigma, Sigma, Sigma)



sample_mixture_HR<-function(d,r,Sigma,A,N){
  final<-replicate(N,mgpd_simulation_mixture_HR(d,r,Sigma,A))
  return(final)
}
set.seed(11)
X_mix_HR<-sample_mixture_HR(d,r,Sigma,A,11)
