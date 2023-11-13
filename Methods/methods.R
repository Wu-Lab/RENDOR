

# RENDOR
RENDOR <- function(mat, m, eps1, eps2) {
  n <- dim(mat)[1]
  mat <- mat + t(mat)
  
  mat <- mat + eps1 + eps2*diag(n)
  
  P1 <- mat / rowSums(mat)
  
  P2 <- m * ginv((m - 1) * diag(n) + P1) %*% P1
  P2 <- P2 - pmin(apply(P2, 1, min), 0)
  P2 <- P2 / rowSums(P2)
  
  w_out <- diag(rowSums(mat)) %*% P2
  
  m1 <- min(min(w_out))
  m2 <- max(max(w_out))
  w_out <- (w_out-m1)/(m2-m1)
  
  w_out <- w_out + t(w_out)
  w_out <- w_out - diag(diag(w_out))
  
  return(w_out)
}




# partial normalize
partial <- function(mat){
  mat_new <- -diag(1/diag(mat))%*%mat%*%diag(1/diag(mat))
  diag(mat_new) <- 1
  return(mat_new)
}
# ICM
ICM <- function(mat) {
  n <- dim(mat)[1]
  mat <- mat + t(mat)
  
  net_new <- ginv(mat)
  net_new <- net_new + 1e-18
  net_new <- partial(net_new)
  net_new <- abs(net_new)
  
  m1 <- min(min(net_new))
  m2 <- max(max(net_new))
  net_new <- (net_new-m1)/(m2-m1)
  
  net_new <- net_new + t(net_new)
  net_new <- net_new - diag(diag(net_new))
  
  return(net_new)
}



# Silencer
Silencer <- function(mat) {
  n <- dim(mat)[1]
  mat <- mat + t(mat)
  
  temp <- diag(diag((mat-diag(n))%*%mat))
  net_new <- (mat-diag(n)+temp) %*% ginv(mat)
  net_new <- abs(net_new)
  
  m1 <- min(min(net_new))
  m2 <- max(max(net_new))
  net_new <- (net_new-m1)/(m2-m1)
  
  net_new <- net_new + t(net_new)
  net_new <- net_new - diag(diag(net_new))
  
  return(net_new)
}



# network deconvolution
ND <- function(mat) {
  n <- dim(mat)[1]
  mat <- mat + t(mat)
  
  Eigen <- eigen(mat)
  U <- Eigen$vectors
  D <-diag(Eigen$values)
  lam_n <- abs(min(min(diag(D)),0))
  lam_p <- abs(max(max(diag(D)),0))
  
  beta <-  0.99
  m1 <- lam_p*(1-beta)/beta
  m2 <- lam_n*(1+beta)/beta
  m <- max(m1,m2)
  
  for (i in 1:n){
    D[i,i] <- (D[i,i])/(m+D[i,i])
  }
  
  mat_nd <- U%*%D%*%ginv(U)
  
  m1 <- min(min(mat_nd))
  m2 <- max(max(mat_nd))
  mat_nd <- (mat_nd-m1)/(m2-m1)
  
  mat_nd <- mat_nd + t(mat_nd)
  mat_nd <- mat_nd - diag(diag(mat_nd))
  
  return(mat_nd)
}


