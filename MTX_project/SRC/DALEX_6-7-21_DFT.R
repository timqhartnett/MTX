#tim DFT features models MTX compounds Multi-caloric Magnets
#3-4-21
#University of Virginia
#Charlottesville, VA 22903

##### CLEAN WORKSPACE ########
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())

#directroy definitions
DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data'
MODEL_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Models'
PLOT_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Plots'
ANALYSIS_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/analysis'
CP_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots'
BD_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/breakdown'

#feature reduction by material science assumptions
dft.dataset <- read.csv(paste(DATA_DIR,'DFT_combined_4-7-21.csv',sep = '/'))
raw.features <- dft.dataset[,2:31]

#### Predict Functions for DALEX explainer, DALEX requires predict functional form predict(model,data)-> Yhat ####
fn.mean.svr.predict = function(model,data){#SVR bootstrapped prediction function is required, 
  predict.data = NULL
  i=1
  print(nrow(data))
  repeat{
    predict.data[[i]] = predict(model[[i]], data)
    i=i+1
    if(i>length(model)) break ()
  }
  if (nrow(data) == 1){ #previous rowmeans immplementation had a bug for single observation
    predict.data.mean = mean(predict.data)
  }
  else {
    predict.data.df = as.data.frame(predict.data)
    predict.data.mean = rowMeans(as.matrix(predict.data.df))
  }
  return(predict.data.mean)
}

fn.sd.svr.predict = function(model,data){#SVR bootstrapped prediction function is required, 
  predict.data = NULL
  i=1
  print(nrow(data))
  repeat{
    predict.data[[i]] = predict(model[[i]], data)
    i=i+1
    if(i>length(model)) break ()
  }
  if (nrow(data) == 1){ #previous rowmeans immplementation had a bug for single observation
    predict.data.mean = mean(predict.data)
    predict.data.sd = sd(predict.data)
    
  }
  else {
    predict.data.df = as.data.frame(predict.data)
    predict.data.mean = rowMeans(as.matrix(predict.data.df))
    predict.data.sd = rowSds(as.matrix(predict.data.df))
  }
  final = predict.data.sd
  return(final)
}


z = fn.median.svr.predict(best_model,theating.test[1:3,])
#### functions for cp plots #######

#libraries
library(matrixStats) # for rowSds
library(e1071) # the e1071 library contains implementations for a number of statistical learning methods.
library(parallel)
library(randomForest)

### models must have associated bagged samples for interpretibilty so we retrain models using the optimized oob
## in future renditions we need to save the bagged sample as a part of the model

## DALEX
library(gridExtra)
library(DALEX)
library(ggplot2)

reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,28)])
colnames(reduced)[9] = "response"

best_sample <- c("1","12","16","21","23","24","25")
theating.train <- reduced[!(rownames(reduced)) %in% best_sample,]
theating.test <- reduced[best_sample,]
load(file = paste(MODEL_DIR,'best.svr.dft.6-10-21.rda',sep = '/'))

explainer_mean = explain(best_model,data = theating.train[,1:8],y=theating.train[,9],predict_function = fn.mean.svr.predict)

training.features <- theating.train[,1:7]
tt.importance <- feature_importance(explainer_mean)
compositions = na.omit(dft.dataset[1:29])

as.data.frame(strtoi(rownames(theating.train)))

feature.names <- colnames(theating.train)

features.list = colnames(reduced[,-ncol(theating.train)])
list.of.list.of.breakdown.values = NULL #since break down does not give all predictions in the same order, must sort interatively into lists
for (i in 1:length(features.list)){
  print(i)
  list.of.list.of.breakdown.values[[i]] <- c(NA,NA) 
}
for (i in 1:nrow(reduced)){
  print(rownames(reduced)[i])
  break.down <- predict_parts_break_down(explainer = explainer_mean,
                                         new_observation = reduced[i,-(ncol(theating.train))],
                                         verbose = FALSE)
  values <- break.down$contribution[-c(1,length(break.down$variable_value))]
  print(values)
  # the above is an aggressively long line which takes the factor output from DALEX and converts to string only considering features
  #not intercepts or predictions
  features = as.character(break.down$variable_name[-1])
  for (j in 1:length(features.list)){
    list.of.list.of.breakdown.values[[j]] <- c(list.of.list.of.breakdown.values[[j]],values[which(features==features.list[j])])
  }
}
full.breakdown <- na.omit(data.frame(list.of.list.of.breakdown.values))
colnames(full.breakdown) <- features.list
rownames(full.breakdown) <- rownames(reduced)

### PCA of breakdown plots ####
library(factoextra)
pca.pred <- prcomp(full.breakdown,scale. = TRUE)
fviz_eig(pca.pred)
fviz_pca_biplot(pca.pred)

theating.clustering <- scale(full.breakdown)

fviz_nbclust(theating.clustering,kmeans,method = 'wss')
#7

fviz_nbclust(theating.clustering,kmeans, method = 'silhouette')
#9

library(cluster)
set.seed(123)
gap_stat <- clusGap(theating.clustering, FUN = kmeans, nstart = 25,
                    K.max = 20, B = 50)
fviz_gap_stat(gap_stat)
#1 cluster (not conclusive)

#WSS gives 7  the most reasonable to deal with but 9 gives better sorting 
k14 <- kmeans(theating.clustering,9,nstart = 25)
fviz_cluster(k14,theating.clustering)
full.breakdown$cluster <- k14$cluster
reduced$cluster <- k14$cluster
write.csv(full.breakdown, file = paste(BD_DIR,'DFT_Tt_breakdown_silhouette.csv',sep = '/'))
write.csv(reduced, file = paste(BD_DIR,'DFT_Tt_silhouette.csv',sep = '/'))

breakdown.bygroup <- NULL
for (i in 1:9){
  mean.breakdown.group <- NULL
  sd.breakdown.group <- NULL
  median.breakdown.group <- NULL
  upper.quant.breakdown.group <- NULL
  lower.quant.breakdown.group <- NULL
  positive.group <- NULL
  for (j in 1:length(feature.names[-length(feature.names)])){
    mean.breakdown.group[[j]] <- mean(full.breakdown[which(full.breakdown$cluster==i),j])
    median.breakdown.group[[j]] <- median(full.breakdown[which(full.breakdown$cluster==i),j])
    lower.quant.breakdown.group[[j]] <- quantile(full.breakdown[which(full.breakdown$cluster==i),j])[[2]]
    upper.quant.breakdown.group[[j]] <- quantile(full.breakdown[which(full.breakdown$cluster==i),j])[[4]]
    sd.breakdown.group[[j]] <- sd(full.breakdown[which(full.breakdown$cluster==i),j])
  }
  breakdown.bygroup[[i]] <- data.frame(mean=mean.breakdown.group,sd=sd.breakdown.group,median=median.breakdown.group,
                                       upper_quant = upper.quant.breakdown.group, lower_quant = lower.quant.breakdown.group,
                                       row.names = features.list)
  #write.csv(breakdown.bygroup[[i]],file = paste(paste(BD_DIR,paste("DFT_breakdown_cluster",i,sep = '-'),sep = "/"),
  #                                             "csv",sep = "."))
}
values.c9 = reduced[which(reduced$cluster==9),-c(9,10)]
predicted.mean.c9 = fn.mean.svr.predict(best_model,values.c9)
predicted.sd.c9 = fn.sd.svr.predict(best_model,values.c9)

values.c3 = reduced[which(reduced$cluster==3),-c(9,10)]
predicted.mean.c3 = fn.mean.svr.predict(best_model,values.c3)
predicted.sd.c3 = fn.sd.svr.predict(best_model,values.c3)

values.c5 = reduced[which(reduced$cluster==5),-c(9,10)]
predicted.mean.c8 = fn.mean.svr.predict(best_model,values.c5)

values.c4 = reduced[which(reduced$cluster==6),-c(9,10)]
predicted.mean.c8 = fn.mean.svr.predict(best_model,values.c4)
#delta V for compositions in clusters 2 and 9
volume.compositions <- reduced[which(reduced$cluster == 3),-length(reduced)]
delta.v.cp <- seq(min(reduced$Delta_V),max(reduced$Delta_V),(max(reduced$Delta_V)-min(reduced$Delta_V))/100)
for (i in 1:nrow(volume.compositions)){
  predictions <- data.frame(Delta_V=delta.v.cp)
  df <- volume.compositions[i,]
  df <- df[rep(1,length(delta.v.cp)),]
  df$Delta_V <- delta.v.cp
  predictions$mean <- fn.mean.svr.predict(best_model,df)
  predictions$sd <- fn.sd.svr.predict(best_model,df)
  write.csv(predictions,file = paste(paste(CP_DIR,paste('delta_v3_cp-',i,sep = ''),sep = '/'),'csv',sep = '.'),
            row.names=FALSE)
}
 #### T-X profile Dataframe ####
for (i in 1:9){
  print(paste('cluster: ',i,sep = ''))
  tx.compositions <- reduced[which(reduced$cluster == i),-length(reduced)]
  predictions.tx <- cbind(tx.compositions$Pnma.T.X,fn.mean.svr.predict(best_model,tx.compositions))
  write.csv(predictions.tx,file = paste(paste(CP_DIR,paste('tx-',i,'-predictions' ,sep = ''),
                                           sep = '/'),'csv',sep = '.'), row.names = FALSE)
  t.x.cp <- seq(min(reduced$Pnma.T.X),max(reduced$Pnma.T.X),(max(reduced$Pnma.T.X)-min(reduced$Pnma.T.X))/100)
  for (j in 1:nrow(tx.compositions)){
    predictions <- data.frame(Pnma.T.X=t.x.cp)
    df <- tx.compositions[j,]
    df <- df[rep(1,length(t.x.cp)),]
    df$Pnma.T.X <- t.x.cp
    predictions$mean <- fn.mean.svr.predict(best_model,df)
    predictions$sd <- fn.sd.svr.predict(best_model,df)
    composition_name <- as.character(dft.dataset[which(dft.dataset$Theating==tx.compositions$response[j]),1])
    print(composition_name)
    write.csv(predictions,file = paste(paste(CP_DIR,paste('tx',i,'.cp-',composition_name ,sep = ''),
                                             sep = '/'),'csv',sep = '.'), row.names = FALSE)
  }
}

tx.compositions <- reduced[which(reduced$cluster == 9),-length(reduced)]
t.x.cp <- seq(min(reduced$Pnma.T.X),max(reduced$Pnma.T.X),(max(reduced$Pnma.T.X)-min(reduced$Pnma.T.X))/100)
for (i in 1:nrow(tx.compositions)){
  predictions <- data.frame(Pnma.T.X=t.x.cp)
  df <- tx.compositions[i,]
  df <- df[rep(1,length(t.x.cp)),]
  df$Pnma.T.X <- t.x.cp
  predictions$mean <- fn.mean.svr.predict(best_model,df)
  predictions$sd <- fn.sd.svr.predict(best_model,df)
  write.csv(predictions,file = paste(paste(CP_DIR,paste('tx.cp-',i,sep = ''),sep = '/'),'csv',sep = '.'),
              row.names = FALSE)
}

tx.compositions <- reduced[which(reduced$cluster == 5),-length(reduced)]
t.x.cp <- seq(min(reduced$Pnma.T.X),max(reduced$Pnma.T.X),(max(reduced$Pnma.T.X)-min(reduced$Pnma.T.X))/100)
for (i in 1:nrow(tx.compositions)){
  predictions <- data.frame(Pnma.T.X=t.x.cp)
  df <- tx.compositions[i,]
  df <- df[rep(1,length(t.x.cp)),]
  df$Pnma.T.X <- t.x.cp
  predictions$mean <- fn.mean.svr.predict(best_model,df)
  predictions$sd <- fn.sd.svr.predict(best_model,df)
  write.csv(predictions,file = paste(paste(CP_DIR,paste('tx5.cp-',i,sep = ''),sep = '/'),'csv',sep = '.'),
            row.names = FALSE)
}

tx.compositions <- reduced[which(reduced$cluster == 3),-length(reduced)]
t.x.cp <- seq(min(reduced$Pnma.T.X),max(reduced$Pnma.T.X),(max(reduced$Pnma.T.X)-min(reduced$Pnma.T.X))/100)
for (i in 1:nrow(tx.compositions)){
  predictions <- data.frame(Pnma.T.X=t.x.cp)
  df <- tx.compositions[i,]
  df <- df[rep(1,length(t.x.cp)),]
  df$Pnma.T.X <- t.x.cp
  predictions$mean <- fn.mean.svr.predict(best_model,df)
  predictions$sd <- fn.sd.svr.predict(best_model,df)
  write.csv(predictions,file = paste(paste(CP_DIR,paste('tx2.cp-',i,sep = ''),sep = '/'),'csv',sep = '.'),
            row.names = FALSE)
}

tx.compositions <- reduced[which(reduced$cluster == 6),-length(reduced)]
t.x.cp <- seq(min(reduced$Pnma.T.X),max(reduced$Pnma.T.X),(max(reduced$Pnma.T.X)-min(reduced$Pnma.T.X))/100)
for (i in 1:nrow(tx.compositions)){
  predictions <- data.frame(Pnma.T.X=t.x.cp)
  df <- tx.compositions[i,]
  df <- df[rep(1,length(t.x.cp)),]
  df$Pnma.T.X <- t.x.cp
  predictions$mean <- fn.mean.svr.predict(best_model,df)
  predictions$sd <- fn.sd.svr.predict(best_model,df)
  write.csv(predictions,file = paste(paste(CP_DIR,paste('tx4.cp-',i,sep = ''),sep = '/'),'csv',sep = '.'),
            row.names = FALSE)
}

e.compositions <- reduced[which(reduced$cluster == 6),-length(reduced)]
e.cp <- seq(min(reduced$energy),max(reduced$energy),(max(reduced$energy)-min(reduced$energy))/100)
for (i in 1:nrow(tx.compositions)){
  predictions <- data.frame(energy=e.cp)
  df <- e.compositions[i,]
  df <- df[rep(1,length(e.cp)),]
  df$energy <- e.cp
  predictions$mean <- fn.mean.svr.predict(best_model,df)
  predictions$sd <- fn.sd.svr.predict(best_model,df)
  write.csv(predictions,file = paste(paste(CP_DIR,paste('energy.cp-',i,sep = ''),sep = '/'),'csv',sep = '.'),
            row.names = FALSE)
}

#Hysteresis 
load(file = paste(MODEL_DIR,'best.hyst.dft.6-7-21.rda',sep = '/'))
reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,30)])
colnames(reduced)[9] <- "response"
best_sample <- c("2","4","9","11","14","30")
theating.train <- reduced[!(rownames(reduced)) %in% best_sample,]
theating.test <- reduced[best_sample,]

explainer_mean = explain(best_model,data = theating.train[,1:8],y=theating.train[,9],predict_function = fn.mean.svr.predict)
hyst.importance = feature_importance(explainer_mean)

for (i in 1:nrow(dft.dataset)){
  break.down <- predict_parts_break_down(explainer = explainer_mean,new_observation =reduced[i,1:8])
  cp <- predict_profile(explainer = explainer_mean,new_observation = reduced[i,1:8])
  p1 <- plot(break.down)
  p2 <- plot(cp)
  png(paste(paste(BD_DIR,paste('DFT_breakdown_hyst_6-10-21',dft.dataset$Compostion[i], sep = '-'),sep = '/'),'png',sep = '.'),
      width = 1200,height = 822)
  grid.arrange(p1,nrow=1,ncol=1)
  dev.off()
  png(paste(paste(CP_DIR,paste('DFT_CP_hyst_6-10-21',dft.dataset$Compostion[i], sep = '-'),sep = '/'),'png',sep = '.'),
      width = 846,height = 580)
  grid.arrange(p2,nrow=1,ncol=1)
  dev.off()
}

#clustering
library(factoextra)
PCA.for.clustering <- cbind(theating.train,fn.mean.svr.predict(model,theating.train))
pca.pred <- prcomp(PCA.for.clustering,scale. = TRUE)
fviz_eig(pca.pred)
fviz_pca_biplot(pca.pred)
