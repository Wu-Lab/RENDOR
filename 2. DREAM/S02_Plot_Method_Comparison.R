library(ggplot2)

name_ <- '(in silico)'

AUPR <- read.table("./2. DREAM/Results/AUPR.csv",sep=",")
AUROC <- read.table("./2. DREAM/Results/AUROC.csv",sep=",")
AUPR <- as.matrix(AUPR)
AUROC <- as.matrix(AUROC)

methods_name <- c('CLR','Relavance','ARACNE','Pearson','Spearman',
                  'GENIE3','TIGRESS','Inferelator','ANOVerence')

data <- data.frame(
  algorithm = rep(methods_name, 5),
  methods = c(
    rep('raw', 9),
    rep('RENDOR', 9),
    rep('ND', 9),
    rep('ICM', 9),
    rep('Silencer', 9)
  ),
  AUPR = as.vector(t(AUPR)),
  AUROC = as.vector(t(AUROC))
)


data$algorithm <- factor(data$algorithm,
                         c('CLR','ARACNE','Relavance',
                           'Pearson','Spearman','GENIE3',
                           'TIGRESS','Inferelator','ANOVerence'))

data$methods <-
  factor(data$methods, c('raw', 'RENDOR','ND','ICM','Silencer'))



#-------------------------------------------------------------------------------
# plot result
ggplot(data, mapping = aes(x = algorithm, y = AUPR, fill = methods)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  scale_fill_manual(values = c("#BEBEBE", 
                               "#CC2936",
                               "#4BA3C3",
                               "#FFA07A",
                               "#25A161"
  )) +
  ylab(paste0('AUPR ', name_)) + xlab("") + theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text.x = element_text(size=15,angle = 45,hjust = 0.5,vjust=0.5),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15))


ggplot(data, mapping = aes(x = algorithm, y = AUROC, fill = methods)) +
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.5) +
  scale_fill_manual(values = c("#BEBEBE", 
                               "#CC2936",
                               "#4BA3C3",
                               "#FFA07A",
                               "#25A161"
                               )) +
  ylab(paste0('AUROC ', name_)) + xlab("") + theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text.x = element_text(size=15,angle = 45,hjust = 0.5,vjust=0.5),
        axis.text.y = element_text(size=15),
        axis.title.y = element_text(size=15))

# 1000*400


