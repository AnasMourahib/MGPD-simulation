source("mgpd_simulation.R")

Nmain<-function(d,r,alpha,A,N){
  final<-replicate(N,mgpd_simulation(d,r,rep(alpha, r),A))
  return(final)
}
d<-3
r<-3
A<-rbind(c(1,0,0), c(1/2, 1/2, 0), c(1/3, 1/3, 1/3))
alpha<-0.5
N<-10
set.seed(7)
Nmain(d,r,alpha,A,N)



