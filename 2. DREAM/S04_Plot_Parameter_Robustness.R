library(ggplot2)
library(pheatmap)


data_auc <- read.csv('./2. DREAM/Results/para_auc_after_RENDOR.csv', header = F)
colnames(data_auc) <- c('CLR','Relavance','ARACNE','Pearson','Spearman',
                        'GENIE3','TIGRESS','Inferelator','ANOVerence')
rownames(data_auc) <- c(seq(1.2,2,by=0.1), seq(4,10,by=2))

data_pr <- read.csv('./2. DREAM/Results/para_pr_after_RENDOR.csv', header = F)
colnames(data_pr) <- c('CLR','Relavance','ARACNE','Pearson','Spearman',
                        'GENIE3','TIGRESS','Inferelator','ANOVerence')
rownames(data_pr) <- c(seq(1.2,2,by=0.1), seq(4,10,by=2))

data_before_auc <- read.csv('./2. DREAM/Results/para_auc_before.csv', header = F)
data_before_pr <- read.csv('./2. DREAM/Results/para_pr_before.csv', header = F)


pheatmap::pheatmap(data_auc-sapply(data_before_auc,rep,13),
                   color = colorRampPalette(c('#008B8B', 'white', '#FF7F00'), bias=0.35)(100),
                   cluster_rows = F,
                   cluster_cols = F,
                   gaps_col = c(5),
                   angle_col = 45,
                   border_color = "black",
                   # display_numbers = TRUE,
                   display_numbers = matrix(ifelse((data_auc-sapply(data_before_auc,rep,13)) >= 0, "+", "-"),
                                            nrow = nrow(data_auc)),
                   )


pheatmap::pheatmap(data_pr-sapply(data_before_pr,rep,13),
                   color = colorRampPalette(c('#008B8B', 'white', '#FF7F00'), bias=0.6)(100),
                   cluster_rows = F,
                   cluster_cols = F,
                   angle_col = 45,
                   gaps_col = c(5),
                   border_color = "black",
                   # display_numbers = TRUE,
                   display_numbers = matrix(ifelse((data_pr-sapply(data_before_pr,rep,13)) > 0, "+", "-"),
                                            nrow = nrow(data_auc)),
                   )


# 500*400


