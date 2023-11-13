library(igraph)
library(MASS)
library(expm)
library(ggplot2)
source('./1. Simulated/help_func.R')
source('Methods/methods.R')


TP <- rep(0, 4)
FP <- rep(0, 4)
TN <- rep(0, 4)
FN <- rep(0, 4)

# --------- generate original cycle graph ---------------------
# g1 <- make_ring(20)
# n <- length(V(g1))
# e1 <- length(E(g1))
# A1 <- matrix(as_adjacency_matrix(g1),nrow = n)
A1 <- as.matrix(read.csv(
  './1. Simulated/Example/cycle_A1.csv',
  header = 1,
  row.names = 1
))
colnames(A1) <- 1:20
g1 <- graph_from_adjacency_matrix(A1, mode = 'undirected')
e1 <- length(E(g1))
n <- length(V(g1))
E(g1)$color <- 'black'
plot(
  g1,
  edge.width = 2,
  vertex.color = "#56B4E9",
  main = 'original graph',
  vertex.label.cex = 0.8,
  vertex.size = 20,
  layout = layout_in_circle(g1)
)



# -------------- add indirect edges randomly ---------------------
# A2 <- add_noise(A1,e1)
# g2 <- graph_from_adjacency_matrix(A2,mode = 'undirected')
# e2 <- length(E(g2))
A2 <- as.matrix(read.csv(
  './1. Simulated/Example/cycle_A2.csv',
  header = 1,
  row.names = 1
))
colnames(A2) <- 1:20
g2 <- graph_from_adjacency_matrix(A2, mode = 'undirected')
e2 <- length(E(g2))
g2 <- g1 + edge(c(t(get.edgelist(g2 - g1))), color = '#EE2C2C')
plot(
  g2,
  edge.width = 2,
  vertex.color = "#56B4E9",
  main = 'noisy graph',
  vertex.label.cex = 0.8,
  vertex.size = 20,
  layout = layout_in_circle(g2)
)




#  --------------------- RENDOR denoising  ---------------------
w_out <- RENDOR(A2, 2, 1, 1)
A_out <- change_into_adj_keep_edges(w_out, e1)
g3 <- graph_from_adjacency_matrix(A_out, mode = 'undirected')
E(g3)$color <- 'black'
e3 <- length(E(g3))
g3 <- g1 - (g1 - g3) + edge(c(t(get.edgelist(g3 - g1))), color = '#EE2C2C')
plot(
  g3,
  edge.width = 2,
  vertex.color = "#56B4E9",
  main = 'RENDOR denoised graph',
  vertex.label.cex = 0.8,
  vertex.size = 20,
  layout = layout_in_circle(g3)
)

TP[1] <- length(E(g3 - (g3 - g1)))
FP[1] <- length(E(g3 - g1))
FN[1] <- length(E(g1 - g3))
TN[1] <- n * (n - 1) / 2 - TP[1] - FP[1] - FN[1]


#  --------------------- ND denoising  ---------------------
w_out <- ND(A2)
A_out <- change_into_adj_keep_edges(w_out, e1)
g4 <- graph_from_adjacency_matrix(A_out, mode = 'undirected')
E(g4)$color <- 'black'
e4 <- length(E(g4))
g4 <- g1 - (g1 - g4) + edge(c(t(get.edgelist(g4 - g1))), color = '#EE2C2C')
plot(
  g4,
  edge.width = 2,
  vertex.color = "#56B4E9",
  main = 'ND denoised graph',
  vertex.label.cex = 0.8,
  vertex.size = 20,
  layout = layout_in_circle(g4)
)

TP[2] <- length(E(g4 - (g4 - g1)))
FP[2] <- length(E(g4 - g1))
FN[2] <- length(E(g1 - g4))
TN[2] <- n * (n - 1) / 2 - TP[2] - FP[2] - FN[2]



#  --------------------- ICM denoising  ---------------------
w_out <- ICM(A2)
A_out <- change_into_adj_keep_edges(w_out, e1)
g5 <- graph_from_adjacency_matrix(A_out, mode = 'undirected')
E(g5)$color <- 'black'
e5 <- length(E(g5))
g5 <- g1 - (g1 - g5) + edge(c(t(get.edgelist(g5 - g1))), color = '#EE2C2C')
plot(
  g5,
  edge.width = 2,
  vertex.color = "#56B4E9",
  main = 'ICM denoised graph',
  vertex.label.cex = 0.8,
  vertex.size = 20,
  layout = layout_in_circle(g3)
)

TP[3] <- length(E(g5 - (g5 - g1)))
FP[3] <- length(E(g5 - g1))
FN[3] <- length(E(g1 - g5))
TN[3] <- n * (n - 1) / 2 - TP[3] - FP[3] - FN[3]



#  --------------------- Silencer denoising  ---------------------
w_out <- Silencer(A2)
A_out <- change_into_adj_keep_edges(w_out, e1)
g6 <- graph_from_adjacency_matrix(A_out, mode = 'undirected')
E(g6)$color <- 'black'
e6 <- length(E(g6))
g6 <- g1 - (g1 - g6) + edge(c(t(get.edgelist(g6 - g1))), color = '#EE2C2C')
plot(
  g6,
  edge.width = 2,
  vertex.color = "#56B4E9",
  main = 'Silencer denoised graph',
  vertex.label.cex = 0.8,
  vertex.size = 20,
  layout = layout_in_circle(g3)
)
# 400*400

TP[4] <- length(E(g6 - (g6 - g1)))
FP[4] <- length(E(g6 - g1))
FN[4] <- length(E(g1 - g6))
TN[4] <- n * (n - 1) / 2 - TP[4] - FP[4] - FN[4]



# -----------------------------------------------------------
plot_df_P <- data.frame(
  methods = rep(c('RENDOR', 'ND', 'ICM', 'Silencer'), 2),
  count.type = c(rep('TP', 4), rep('FP', 4)),
  counts = c(TP, FP)
)

plot_df_N <- data.frame(
  methods = rep(c('RENDOR', 'ND', 'ICM', 'Silencer'), 2),
  count.type = c(rep('TN', 4), rep('FN', 4)),
  counts = c(TN, FN)
)

plot_df_P$methods <- factor(plot_df_P$methods,
                            c('RENDOR', 'ND', 'ICM', 'Silencer'))
plot_df_N$methods <- factor(plot_df_N$methods,
                            c('RENDOR', 'ND', 'ICM', 'Silencer'))


ggplot(plot_df_P, aes(methods, counts)) +
  geom_bar(stat = "identity", aes(fill = count.type, color = count.type)) +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#4AAD52", "#F4D35E")) +
  scale_color_manual(values = c("#4AAD52", "#F4D35E"))


ggplot(plot_df_N, aes(methods, counts)) +
  geom_bar(stat = "identity", aes(fill = count.type, color = count.type)) +
  theme(legend.position = "top") +
  scale_fill_manual(values = c("#28AFB0", "#EE964B")) +
  scale_color_manual(values = c("#28AFB0", "#EE964B"))
# 400*250
