


#---------------------------------------------------------------------------------
#                               Methods
#---------------------------------------------------------------------------------
change_into_adj_keep_edges <- function(w, edges){
  n <- dim(w)[1]
  temp <- sort(w,decreasing = T)
  thre <- temp[2*edges+1]
  A <- matrix(0,nrow = n,ncol = n)
  A[w>thre] <- 1
  A <- A+t(A)
  A[A>0] <- 1
  return(A)
}


change_into_adj <- function(w,thre){
  adj_mat <- matrix(0,nrow = n,ncol = n)
  adj_mat[w>thre] <- 1
  return(adj_mat)
}


# compute noisy proportion
compute_error <- function(A_true,A_pre) {
  error_edges <- sum((A_true==0)&(A_pre==1))
  true_edges <- sum((A_true==1)&(A_pre==1))
  
  SNR_inverse <- error_edges/true_edges
  return(SNR_inverse)
  
}


add_noise <- function(A,noise){
  
  n <- dim(A)[1]
  for (i in 1:n){
    for (j in 1:n){
      if(abs(i-j)<8&abs(i-j)>2){
        len <- length(shortest_paths(g1, from=i, to =j)$vpath[[1]])-2
        p <- 0.02+1/((len^2))
        A[i,j] <- sample(c(0,1),1,replace = T,
                         prob=c(1-p,p))
      }
    }
  }
  
  A <- (A+t(A))/2
  A <- (A>0)*1
  return(A)
}



