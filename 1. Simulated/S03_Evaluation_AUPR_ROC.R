library(ggplot2)
library(igraph)
library(MASS)
library(pROC)
library(ROCR)
library(ggplot2)
source('1. Simulated/help_func.R')
source('./Methods/methods.R')



#----------------------------------------------------------------------
#   Compare the result of denoising methods on noisy graph
#            Compute the average result of 100 BA/ER graph
#----------------------------------------------------------------------
data_name <- 'BA'
loops <- 100


# Network settings
if (data_name == 'BA') {
  n <- 50
  m <- 3
  title_AUROC <- paste0('AUROC of BA graphs')
  title_AUPR <- paste0('AUPR of BA graphs')
  
} else if (data_name == 'ER') {
  p <- 0.3
  n <- 50
  title_AUROC <- paste0('AUROC of ER graphs')
  title_AUPR <- paste0('AUPR of ER graphs')
}


# record AUROC and AUPR scores
noisy_proportion <- seq(0.1, 0.5, 0.01) 
auc_ND <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
pr_ND <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
auc_RENDOR2 <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
pr_RENDOR2 <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
auc_RENDOR3 <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
pr_RENDOR3 <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
auc_ICM <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
pr_ICM <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
auc_Silencer <- matrix(0, nrow = loops, ncol = length(noisy_proportion))
pr_Silencer <- matrix(0, nrow = loops, ncol = length(noisy_proportion))



# import BA data
for (i in 1:length(noisy_proportion)) {
  for (j in 1:loops) {
    print(paste0('i_', i, '    j_', j))
    
    if (data_name == 'BA') {
      A1 <- read.csv(paste0('./1. Simulated/simudata_BA/true_', i, '_loop_', j, '.csv'))
      A2 <- read.csv(paste0('./1. Simulated/simudata_BA/noisy_', i, '_loop_', j, '.csv'))
    } else if (data_name == 'ER') {
      A1 <- read.csv(paste0('./1. Simulated/simudata_ER/true_', i, '_loop_', j, '.csv'))
      A2 <- read.csv(paste0('./1. Simulated/simudata_ER/noisy_', i, '_loop_', j, '.csv'))
    }
    
    A1 <- as.matrix(A1[, 2:51])
    A2 <- as.matrix(A2[, 2:51])
    
    
    output_w_RENDOR2 <- RENDOR(A2, 2, 1, 1)
    output_w_RENDOR3 <- RENDOR(A2, 3, 1, 1)
    output_w_ND <- ND(A2)
    output_w_ICM <- ICM(A2)
    output_w_Silencer <- Silencer(A2)
    
    
    test <- data.frame(
      label = as.factor(as.matrix(A1)),
      output_RENDOR2 = as.vector(output_w_RENDOR2),
      output_RENDOR3 = as.vector(output_w_RENDOR3),
      output_ND = as.vector(output_w_ND),
      output_Silencer = as.vector(output_w_Silencer),
      output_ICM = as.vector(output_w_ICM)
    )
    
    pred_ND <- prediction(test$output_ND, test$label)
    pred_RENDOR2 <- prediction(test$output_RENDOR2, test$label)
    pred_RENDOR3 <- prediction(test$output_RENDOR3, test$label)
    pred_ICM <- prediction(test$output_ICM, test$label)
    pred_Silencer <- prediction(test$output_Silencer, test$label)
    
    auc_RENDOR2[j, i] <- unlist(slot(performance(pred_RENDOR2, 'auc'),
                                  "y.values"))
    pr_RENDOR2[j, i] <- unlist(slot(performance(pred_RENDOR2, 'aucpr'),
                                 "y.values"))
    
    auc_RENDOR3[j, i] <- unlist(slot(performance(pred_RENDOR3, 'auc'),
                                  "y.values"))
    pr_RENDOR3[j, i] <- unlist(slot(performance(pred_RENDOR3, 'aucpr'),
                                 "y.values"))
    
    auc_ND[j, i] <- unlist(slot(performance(pred_ND, 'auc'),
                                "y.values"))
    pr_ND[j, i] <- unlist(slot(performance(pred_ND, 'aucpr'),
                               "y.values"))
    
    auc_ICM[j, i] <- unlist(slot(performance(pred_ICM, 'auc'),
                                 "y.values"))
    pr_ICM[j, i] <- unlist(slot(performance(pred_ICM, 'aucpr'),
                                "y.values"))
    
    auc_Silencer[j, i] <-
      unlist(slot(performance(pred_Silencer, 'auc'),
                  "y.values"))
    pr_Silencer[j, i] <-
      unlist(slot(performance(pred_Silencer, 'aucpr'),
                  "y.values"))
    
  }
}


#------------------------------------------------------
data <- data.frame(
  noise_proportion = c(rep(noisy_proportion, 5)),
  group = c(
    rep('ND', length(noisy_proportion)),
    rep('ICM', length(noisy_proportion)),
    rep('Silencer', length(noisy_proportion)),
    rep('RENDOR(m=2)', length(noisy_proportion)),
    rep('RENDOR(m=3)', length(noisy_proportion))
  ),
  score_auc = c(
    apply(auc_ND, 2, mean),
    apply(auc_ICM, 2, mean),
    apply(auc_Silencer, 2, mean),
    apply(auc_RENDOR2, 2, mean),
    apply(auc_RENDOR3, 2, mean)
  ),
  sd_auc = c(
    apply(auc_ND, 2, sd),
    apply(auc_ICM, 2, sd),
    apply(auc_Silencer, 2, sd),
    apply(auc_RENDOR2, 2, sd),
    apply(auc_RENDOR3, 2, sd)
  ),
  score_pr = c(
    apply(pr_ND, 2, mean),
    apply(pr_ICM, 2, mean),
    apply(pr_Silencer, 2, mean),
    apply(pr_RENDOR2, 2, mean),
    apply(pr_RENDOR3, 2, mean)
  ),
  sd_pr = c(
    apply(pr_ND, 2, sd),
    apply(pr_ICM, 2, sd),
    apply(pr_Silencer, 2, sd),
    apply(pr_RENDOR2, 2, sd),
    apply(pr_RENDOR3, 2, sd)
  )
)


data$group <- factor(data$group,
                     c('RENDOR(m=2)', 'RENDOR(m=3)', 'ND', 'ICM', 'Silencer'))


#-----------------------------------------------------------------
# plot performance
ggplot(data,
       aes(
         x = noise_proportion,
         y = score_auc,
         group = group,
         color = group
       )) +
  geom_line(size = 0.8) + geom_point(size = 1.5) +
  geom_ribbon(aes(ymin = score_auc - sd_auc, ymax = score_auc + sd_auc), alpha = 0.05) +
  scale_fill_manual(values = c("#CC2936",
                               "#81282f",
                               "#4BA3C3",
                               "#FFA07A",
                               "#25A161")) +
  scale_color_manual(values = c("#CC2936",
                                "#81282f",
                                "#4BA3C3",
                                "#FFA07A",
                                "#25A161")) +
  theme_bw() +
  labs(title = title_AUROC) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 15)
  )


ggplot(data,
       aes(
         x = noise_proportion,
         y = score_pr,
         group = group,
         color = group
       )) +
  geom_line(size = 0.8) + geom_point(size = 1.5) +
  geom_ribbon(aes(ymin = score_pr - sd_pr, ymax = score_pr + sd_pr), alpha = 0.05) +
  scale_fill_manual(values = c("#CC2936",
                               "#81282f",
                               "#4BA3C3",
                               "#FFA07A",
                               "#25A161")) +
  scale_color_manual(values = c("#CC2936",
                                "#81282f",
                                "#4BA3C3",
                                "#FFA07A",
                                "#25A161")) +
  theme_bw() +
  labs(title = title_AUPR) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 15)
  )


# 600*400
