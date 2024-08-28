source("Helper_functions.R")

mgpd_simulation_mixture_HR<-function(d,r,Sigma,A){
  w<-mass_of_scenario(d, r, Sigma, A)  
  T<-rep(-Inf,d)
  b<-sample(c(1:r), prob=w,size=1)
  sign_column_b<-which(A[,b]>0)
  n_column_b<-(A[sign_column_b,b])/sum(A[sign_column_b,b])  
  accept=FALSE
  while(!(accept)){
    if(length(sign_column_b)==1){a<-sign_column_b}
    else{a<-sample(sign_column_b, prob=n_column_b,size=1)}
    T[sign_column_b]<- HR_generator(a ,b, A[sign_column_b, b], sign_column_b, Sigma[[b]],r)  
    U_0<-runif(1,min=0,max=1)
    if ( (U_0 < rejection_sampling(T)) ) {
      accept = TRUE
      }
  }
  E<-rexp(1,1)
  Y <- exp(T - max(T) + E) - 1
  return(Y)
}


