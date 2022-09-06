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
  y_bar <- stats::weighted.mean(x = Y,w = rowSums(a$z_ij))
  y_bar_g <- mclust::covw(X = Y,Z = a$z_ij,normalize = F)$mean # vector of 1 x G
  BSS <- sum(a$csize*(y_bar_g-y_bar)^2)
  mu_xi_beta_G <- apply(a$b,1,function(beta) beta%*%t(cbind(1,X)))
  EWSS <- sum(sweep(mu_xi_beta_G,MARGIN = 2,STATS = y_bar_g,FUN = "-")^2*a$z_ij)
  RWSS <- sum((Y-mu_xi_beta_G)^2*a$z_ij)
  # TSS <- sum((Y-y_bar)^2)
  # c(BSS=BSS,EWSS=EWSS,RWSS=RWSS, TSS=TSS)
  c(BSS=BSS,EWSS=EWSS,RWSS=RWSS)
}
