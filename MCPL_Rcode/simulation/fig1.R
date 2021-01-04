# load("unbalance_sim_n300p5d3_e3")

n=300

s.est <- sapply(sim.res, function(x) {x$s.est})
length(which(s.est == 2))/500
table(s.est)
theta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$theta.hat})
beta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$beta.hat})
delta.hat <- sapply(sim.res[which(s.est == 2)], function(x) {x$obj$delta.hat})
a.hat <- sapply(sim.res[which(s.est == 2)], function(x){x$obj$a.hat})

cl.hat = sapply(sim.res[which(s.est == 2)], function(x){ x$id.est }, simplify = F )
splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
cl.tr <- sapply(sim.res[which(s.est == 2)], function(x){ splitAt(1:n, x$tr.cp+1) }, simplify = F)


grp.id <- function(x){
  
  a <- NULL
  for (i in 1:length(x)) {
    a[[i]] = rep(i, length(x[[i]]))
  }
  
  grp.id = do.call("c", a)
  
}

group.hat <- sapply(cl.hat, grp.id, simplify = F)
group.tr <- sapply(cl.tr, grp.id, simplify = F)

NMI.fun <- function(clusters, true_labels) {
  
  tbl = table(clusters, true_labels)
  conv_df = as.data.frame.matrix(tbl)
  mutual_information = 0.0
  for (i in 1:nrow(conv_df)) {
    
    for (j in 1:ncol(conv_df)) {
      
      if (conv_df[i,j] > 0.0) {
        
        mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log2((sum(tbl) * 
                                                                                       conv_df[i,j]) / (sum(conv_df[i,]) * sum(conv_df[,j]))))
      }
    }
  }
  entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))
  
  entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log2(sum(x) / sum(tbl))))
  
  NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))
  
  return(NMI)
  
}

boxplot(mapply(NMI.fun, group.hat, group.tr))


dat = mapply(NMI.fun, group.hat, group.tr)
dat <- stack(as.data.frame(dat))

xlab21 = expression(paste("Example 2", " (", Sigma[1], ", ", "n=150", ")" ))
ex2n150p5 <- ggplot(dat) + ggtitle(xlab21)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab22 = expression(paste("Example 2", " (", Sigma[1], ", ", "n=300", ")" ))
ex2n300p5 <- ggplot(dat) + ggtitle(xlab22)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


xlab21s2 = expression(paste("Example 2", " (", Sigma[2], ", ", "n=150", ")" ))
ex2n150p5s2 <- ggplot(dat) + ggtitle(xlab21s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab22s2 = expression(paste("Example 2", " (", Sigma[2], ", ", "n=300", ")" ))
ex2n300p5s2 <- ggplot(dat) + ggtitle(xlab22s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


xlab21s3 = expression(paste("Example 2", " (", Sigma[3], ", ", "n=150", ")" ))
ex2n150p5s3 <- ggplot(dat) + ggtitle(xlab21s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab22s3 = expression(paste("Example 2", " (", Sigma[3], ", ", "n=300", ")" ))
ex2n300p5s3 <- ggplot(dat) + ggtitle(xlab22s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )



xlab41 = expression(paste("Example 4", " (", Sigma[1], ", ", "n=150", ")" ))
ex4n150p5 <- ggplot(dat) + ggtitle(xlab41)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab42 = expression(paste("Example 4", " (", Sigma[1], ", ", "n=300", ")" ))
ex4n300p5 <- ggplot(dat) + ggtitle(xlab42)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


xlab41s2 = expression(paste("Example 4", " (", Sigma[2], ", ", "n=150", ")" ))
ex4n150p5s2 <- ggplot(dat) + ggtitle(xlab41s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab42s2 = expression(paste("Example 4", " (", Sigma[2], ", ", "n=300", ")" ))
ex4n300p5s2 <- ggplot(dat) + ggtitle(xlab42s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


xlab41s3 = expression(paste("Example 4", " (", Sigma[3], ", ", "n=150", ")" ))
ex4n150p5s3 <- ggplot(dat) + ggtitle(xlab41s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab42s3 = expression(paste("Example 4", " (", Sigma[3], ", ", "n=300", ")" ))
ex4n300p5s3 <- ggplot(dat) + ggtitle(xlab42s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}


multiplot( ex2n150p5, ex2n300p5,
           ex2n150p5s2, ex2n300p5s2, 
           ex2n150p5s3, ex2n300p5s3, cols=2)


multiplot( ex4n150p5, ex4n300p5,
           ex4n150p5s2, ex4n300p5s2, 
           ex4n150p5s3, ex4n300p5s3, cols=2)



####
load("single_res_n300p5d2e3")
################# misclassification rate
n=300
id.fun <- function(x) {
  id <- rep(1,n)
  id[x] <- 0
  id
}

##
cl.hat <- apply(res, 2, function(x){ x$fit$nind } )
cl.tr <- apply(res, 2, function(x){ c(1:x$tr.cl)})

group.hat <- sapply(cl.hat, id.fun, simplify = F)
group.tr <- sapply(cl.tr, id.fun, simplify = F)
mean(mapply(function(x,y){mean(abs(x-y))}, group.hat, group.tr ))
sd(mapply(function(x,y){mean(abs(x-y))}, group.hat, group.tr ))




### NMI
NMI.fun <- function(clusters, true_labels) {
  
  tbl = table(clusters, true_labels)
  conv_df = as.data.frame.matrix(tbl)
  mutual_information = 0.0
  for (i in 1:nrow(conv_df)) {
    
    for (j in 1:ncol(conv_df)) {
      
      if (conv_df[i,j] > 0.0) {
        
        mutual_information = mutual_information + ((conv_df[i,j] / sum(tbl)) * log((sum(tbl) * 
                                                                                      conv_df[i,j]) / (sum(conv_df[i,]) * sum(conv_df[,j]))))
      }
    }
  }
  entr_cluster = sum(apply(conv_df, 1, function(x) -(sum(x) / sum(tbl)) * log(sum(x) / sum(tbl))))
  
  entr_class = sum(apply(conv_df, 2, function(x) -(sum(x) / sum(tbl)) * log(sum(x) / sum(tbl))))
  
  NMI = (mutual_information / ((entr_cluster + entr_class) / 2.0))
  
  return(NMI)
  
}

boxplot(mapply(NMI.fun, group.hat, group.tr))

dat = mapply(NMI.fun, group.hat, group.tr)
dat <- stack(as.data.frame(dat))

library(ggplot2)
xlab11 = expression(paste("Example 1", " (", Sigma[1], ", ", "n=150", ")" ))
ex1n150p5 <- ggplot(dat) + ggtitle(xlab11)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab12 = expression(paste("Example 1", " (", Sigma[1], ", ", "n=300", ")" ))
ex1n300p5 <- ggplot(dat) + ggtitle(xlab12)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )


xlab11s2 = expression(paste("Example 1", " (", Sigma[2], ", ", "n=150", ")" ))
ex1n150p5s2 <- ggplot(dat) + ggtitle(xlab11s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab12s2 = expression(paste("Example 1", " (", Sigma[2], ", ", "n=300", ")" ))
ex1n300p5s2 <- ggplot(dat) + ggtitle(xlab12s2)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab11s3 = expression(paste("Example 1", " (", Sigma[3], ", ", "n=150", ")" ))
ex1n150p5s3 <- ggplot(dat) + ggtitle(xlab11s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

xlab12s3 = expression(paste("Example 1", " (", Sigma[3], ", ", "n=300", ")" ))
ex1n300p5s3 <- ggplot(dat) + ggtitle(xlab12s3)+ xlab("") + ylab("NMI")+
  geom_boxplot(aes(x = ind, y = values))+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() )

multiplot( ex1n150p5, ex1n300p5,
           ex1n150p5s2, ex1n300p5s2, 
           ex1n150p5s3, ex1n300p5s3, cols=2)
