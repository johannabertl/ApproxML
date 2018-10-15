# The functions in this file can be used to analyse the gradients along the estimation sequence. The input is the matrix lik that is produced by KDKW.FD.comp2.

# 

grad = function(lik, ce, gamma = 1/6, C = 0, lg=T){
  p = ncol(lik)/2
  K = nrow(lik)
  grad = matrix(NA, K, p)
 	ck = ce/(((1:K)+1+C)^gamma)
 	if(lg==T) lik = log(lik)
  for(i in 1:p){
    grad[,i] = (lik[,i*2-1] - lik[,i*2])/(2*ck)
  }
  grad
}