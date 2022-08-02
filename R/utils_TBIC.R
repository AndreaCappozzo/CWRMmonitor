#' Title
#'
#' @param c1
#' @param c2
#' @param k
#' @param p
#'
#' @return
#' @export
#'
#' @examples
pnlt_term<-function(c1,c2,k,p){
  fixed_term<-(k-1) +k*p +k*(p+1)
  #(k-1) mixture weights, k*p means in X,k*(p+1) beta coeff for regression bgx+b0g,
  var_term <- 1 + ( (k*p-1) + k*p*(p-1)/2 ) * (1-1/c1)
  # this is the part relative to modelling X
  # 1+ (k*p-1)  1 free eigenvalue and k*p-1 constrained eigenvalues,
  # k*p*(p-1)/2 rotation matrices for Sigma_g,
  # (1-1/c1) take into account constrained eigenvalues
  var_term <- var_term + 1 + (k-1) * (1-1/c2)
  # this is the part relative to modelling Y|X
  # 1+ (k-1)*(1-1/c2) one free σ^2g and k-1 constrained σ^2g
  return(fixed_term+var_term)
}

#' Title
#'
#' @param a
#' @param Y
#' @param X
#'
#' @return
#' @export
#'
#' @examples
var_dec<-function(a,Y,X){
  # in a$czise we have  sumCols(z_ij)
  rob_mean_Y<-mean(Y*rowSums(a$z_ij))*length(Y)/length(which(rowSums(a$z_ij)==1))
  #  iter$center[k,] = (t(a$z_ij[,k]) %*% Y) / iter$csize[k]
  means_y_g = (t(a$z_ij) %*% Y) / a$csize
  #  means_y_g <- a$b[,1] + a$center * a$b[,2]
  K<-length(a$cw)
  sc_mu_j_da_mu<-as.matrix(means_y_g-rep(rob_mean_Y,K),dim=c(1,K))
  BSS<- sum((sc_mu_j_da_mu)^2 * a$csize)
  #  X1<-cbind(rep(1,length(X)),X)
  #  explained<-X1 %*% t(a$b)- rep(1,length(X)) %*% t(means_y_g)
  #  explained^2
  one=matrix(1,nrow=length(Y),ncol=1)
  X1=cbind(one,X)
  #\mu(x_i;\hat(beta)_g) in formula 30
  # beta_0g+beta_1g'x
  forESS<-one%*%t(means_y_g)-X1%*%t(a$b)
  EWSS<-sum(a$z_ij*(forESS^2))
  #primo modo per calcolare gli errori
  onek<-rep(1,K)
  err=Y%*%t(onek)-X1%*%t(a$b)
  RWSS<-sum(a$z_ij*(err^2))
  # secondo modo per calcolare gli errori
  residual<-Y %*% t(rep(1,K))-X1 %*% t(a$b)
  # sono uguali gli errori=scarti residui?
  # sum(abs(residual-err))
  # cbind(Y,X1 %*% t(a$b),a$z_ij)[1:5,]
  scartirob<-(Y-rob_mean_Y)^2*rowSums(a$z_ij)
  # cbind((Y-mean(Y))^2,rowSums(a$z_ij))
  TSS<-sum(scartirob)
  devtot<-stats::var(Y)*(length(Y)-1)
  # devtot
  # TSS
  # BSS  # the soft between sum of squares
  # # (variability of Y explained by the latent group variable)
  # # BSS can be seen as a separation measure along the Y axis
  # EWSS  # the soft within-group sum of squares explained
  # # by the model (thanks to the covariates)
  # RWSS  # the soft residual within-group sum of squares can
  # # WSS=EWSS+RWSS can be seen as a compactness measure
  # BSS+EWSS+RWSS
  # TSS-(BSS+EWSS+RWSS)
  return(c(BSS,EWSS,RWSS))
}
