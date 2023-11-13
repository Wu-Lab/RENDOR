source('./1. Simulated/help_func.R')


#------------------------------------------------------------------
# generate BA noisy graphs
#------------------------------------------------------------------

n <- 50
loops <- 100

noisy_proportion <- seq(0.1,0.5,0.01)  # noise proportion

for (i in 1:length(noisy_proportion)){
  for (j in 1:loops){
    print(paste0('i_',i,'    j',j))
    g1 <- sample_pa(n,m=3,directed = F)
    e1 <- length(E(g1))
    A1 <- matrix(as_adjacency_matrix(g1),nrow = n )
    write.csv(A1,paste0('./1. Simulated/simudata_BA/true_',i,'_loop_',j,'.csv'))
    
    #If the graph has outliers, regenerate the graph
    check <- rowSums(A1)
    if (length(check[check==0])!=0){
      print('error')
      next}
    
    # add noise
    P1 <- A1/rowSums(A1)
    P2 <- P1%*%solve(diag(2,nrow = n)-P1)
    D2 <- as.vector(Null(P2-diag(1,nrow=n)))
    D2 <- D2/sum(D2)
    W2 <- diag(D2)%*%P2
    W2 <- (W2+t(W2))/2
    W2 <- W2-diag(diag(W2))
    A2 <- change_into_adj_keep_edges(W2,e1+noisy_proportion[i]*e1/(1-noisy_proportion[i]))
    
    write.csv(A2,paste0('./1. Simulated/simudata_BA/noisy_',i,'_loop_',j,'.csv'))
  }
}


