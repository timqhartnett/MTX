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
#DATA
library(readxl)
library(e1071)
library(matrixStats)
comp.dataset <- read_xlsx(paste(DATA_DIR,'composition_5_27_21.xlsx',sep = '/'))
reduced.dataset <- as.data.frame(na.omit(comp.dataset[,-c(1,12:14)]))
colnames(reduced.dataset)[10] <- 'response'

load(paste(MODEL_DIR,'best.svr.comp.5-27-21.rda',sep = '/'))
best_sample = c("2","3","5","13","24","28","31","32","41","45","50","54","58","62","66","71")
theating.train <- reduced.dataset[!(rownames(reduced.dataset)) %in% best_sample,]
theating.test <- reduced.dataset[best_sample,]

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


### models must have associated bagged samples for interpretibilty so we retrain models using the optimized oob
## in future renditions we need to save the bagged sample as a part of the model

## DALEX
library(gridExtra)
library(DALEX)
library(ggplot2)

explainer_mean = explain(best_model,data = theating.train[,1:9],y=theating.train[,10],predict_function = fn.mean.svr.predict,
                         label = 'Bootstrapped SVR mean prediction')
#explainer_median = explain(best_model,data = theating.train[,1:9],y=theating.train[,10],predict_function = fn.median.svr.predict)
#explainer_upper = explain(best_model,data = theating.train[,1:9],y=theating.train[,10],predict_function = fn.upper.svr.predict)
#explainer_lower = explain(best_model,data = theating.train[,1:9],y=theating.train[,10],predict_function = fn.lower.svr.predict)
svr.importance <- feature_importance(explainer_mean)
as.data.frame(strtoi(rownames(theating.train)))

feature.names <- colnames(theating.train)

for (i in 1:nrow(reduced.dataset)){
  break.down <- predict_parts_break_down(explainer = explainer_mean,new_observation = reduced.dataset[i,])
  cp <- predict_profile(explainer = explainer_mean,new_observation = reduced.dataset[i,])
  p1 <- plot(break.down)
  p2 <- plot(cp)
  png(paste(paste(BD_DIR,paste('comp_breakdown',rownames(reduced.dataset)[i], sep = '-'),sep = '/'),'png',sep = '.'),
      width = 1200,height = 822)
  grid.arrange(p1,nrow=1,ncol=1)
  dev.off()
  png(paste(paste(CP_DIR,paste('comp_CP',rownames(reduced.dataset)[i], sep = '-'),sep = '/'),'png',sep = '.'),
      width = 846,height = 580)
  grid.arrange(p2,nrow=1,ncol=1)
  dev.off()
}

for (i in 1:nrow(reduced.dataset)){
  
  
}

features.list = colnames(reduced.dataset[,-ncol(theating.train)])
list.of.list.of.breakdown.values = NULL #since break down does not give all predictions in the same order, must sort interatively into lists
for (i in 1:length(features.list)){
  print(i)
  list.of.list.of.breakdown.values[[i]] <- c(NA,NA) 
}
for (i in 1:nrow(reduced.dataset)){
  print(rownames(reduced.dataset)[i])
  break.down <- predict_parts_break_down(explainer = explainer_mean,
                                         new_observation = reduced.dataset[i,-(ncol(theating.train))],
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
full.breakdown$intercept <- break.down$contribution[1]
full.breakdown$prediction <- fn.mean.svr.predict(best_model,reduced.dataset)
write.csv(full.breakdown, file = paste(BD_DIR,'breakdown_comosition_8-17-21.csv',sep = '/'))

### PCA of breakdown plots ####
library(factoextra)
pca.pred <- prcomp(full.breakdown,scale. = TRUE)
fviz_eig(pca.pred)
fviz_pca_biplot(pca.pred)

theating.clustering <- scale(full.breakdown)

fviz_nbclust(theating.clustering,kmeans,method = 'wss')
#10
fviz_nbclust(theating.clustering,kmeans, method = 'silhouette')
#9
library(cluster)
set.seed(123)
gap_stat <- clusGap(theating.clustering, FUN = kmeans, nstart = 25,
                    K.max = 20, B = 50)
fviz_gap_stat(gap_stat)

#### try each #####
#14
k14 <- kmeans(theating.clustering,9,nstart = 25)
fviz_cluster(k14,theating.clustering)
full.breakdown$cluster <- k14$cluster
reduced.dataset$cluster <- k14$cluster
write.csv(full.breakdown, file = paste(BD_DIR,'comp_Tt_breakdown_gapstat.csv',sep = '/'))
write.csv(reduced.dataset, file = paste(BD_DIR,'comp_Tt_gapstat.csv',sep = '/'))

breakdown.bygroup <- NULL
for (i in 1:14){
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
  write.csv(breakdown.bygroup[[i]],file = paste(paste(BD_DIR,paste("local_breakdown_silhoutte_cluster",i,sep = '-'),sep = "/"),
                                                "csv",sep = "."))
}

# 10 silhoutte
k=10
k10 <- kmeans(theating.clustering,k,nstart = 25)
fviz_cluster(k10,theating.clustering)
full.breakdown$cluster <- k10$cluster
reduced.dataset$cluster <- k10$cluster
write.csv(full.breakdown, file = paste(BD_DIR,'comp_Tt_breakdown_silhoutte.csv',sep = '/'))
write.csv(reduced.dataset, file = paste(BD_DIR,'comp_Tt_silhoutte.csv',sep = '/'))

breakdown.bygroup <- NULL
for (i in 1:10){
  mean.breakdown.group <- NULL
  sd.breakdown.group <- NULL
  median.breakdown.group <- NULL
  upper.quant.breakdown.group <- NULL
  lower.quant.breakdown.group <- NULL
  positive.group <- NULL
  for (j in 1:length(feature.names[-length(feature.names)])){
    positive.group[[j]] <- as.integer(mean(full.breakdown[which(full.breakdown$cluster==i),j]) > 0)
    mean.breakdown.group[[j]] <- abs(mean(full.breakdown[which(full.breakdown$cluster==i),j]))
    median.breakdown.group[[j]] <- median(full.breakdown[which(full.breakdown$cluster==i),j])
    lower.quant.breakdown.group[[j]] <- quantile(full.breakdown[which(full.breakdown$cluster==i),j])[[2]]
    upper.quant.breakdown.group[[j]] <- quantile(full.breakdown[which(full.breakdown$cluster==i),j])[[4]]
    sd.breakdown.group[[j]] <- sd(full.breakdown[which(full.breakdown$cluster==i),j])
  }
  breakdown.bygroup[[i]] <- data.frame(mean=mean.breakdown.group,sd=sd.breakdown.group,median=median.breakdown.group,
                                       upper_quant = upper.quant.breakdown.group, lower_quant = lower.quant.breakdown.group,
                                       positive=positive.group,row.names = features.list)
  breakdown.bygroup[[i]] <- breakdown.bygroup[[i]][order(-breakdown.bygroup[[i]]$mean),]
  write.csv(breakdown.bygroup[[i]],file = paste(paste(BD_DIR,paste("local_breakdown_silhouette_cluster",i,sep = '-'),sep = "/"),
                                                "csv",sep = "."))
}

values.c2 = reduced.dataset[which(reduced.dataset$cluster==2),-c(10,11)]
predicted.mean.c2 = fn.mean.svr.predict(best_model,values.c2)
predicted.sd.c2 = fn.sd.svr.predict(best_model,values.c2)

values.c5 = reduced.dataset[which(reduced.dataset$cluster==5),-c(10,11)]
predicted.mean.c5 = fn.mean.svr.predict(best_model,values.c5)
predicted.sd.c5 = fn.sd.svr.predict(best_model,values.c5)

values.c6 = reduced.dataset[which(reduced.dataset$cluster==6),-c(10,11)]
predicted.mean.c6 = fn.mean.svr.predict(best_model,values.c6)
predicted.sd.c6 = fn.sd.svr.predict(best_model,values.c6)

values.c7 = reduced.dataset[which(reduced.dataset$cluster==7),-c(10,11)]
predicted.mean.c7 = fn.mean.svr.predict(best_model,values.c7)
predicted.sd.c7 = fn.sd.svr.predict(best_model,values.c7)

#k=9 WSS
k=9
k9 <- kmeans(theating.clustering,k,nstart = 25)
fviz_cluster(k10,theating.clustering)
full.breakdown$cluster <- k9$cluster
reduced.dataset$cluster <- k9$cluster
write.csv(full.breakdown, file = paste(BD_DIR,'comp_Tt_breakdown_WSS.csv',sep = '/'))
write.csv(reduced.dataset, file = paste(BD_DIR,'comp_Tt_WSS.csv',sep = '/'))

breakdown.bygroup <- NULL
for (i in 1:k){
  mean.breakdown.group <- NULL
  sd.breakdown.group <- NULL
  median.breakdown.group <- NULL
  upper.quant.breakdown.group <- NULL
  lower.quant.breakdown.group <- NULL
  positive.group <- NULL
  for (j in 1:length(feature.names[-length(feature.names)])){
    positive.group[[j]] <- as.integer(mean(full.breakdown[which(full.breakdown$cluster==i),j]) > 0)
    mean.breakdown.group[[j]] <- abs(mean(full.breakdown[which(full.breakdown$cluster==i),j]))
    median.breakdown.group[[j]] <- median(full.breakdown[which(full.breakdown$cluster==i),j])
    lower.quant.breakdown.group[[j]] <- quantile(full.breakdown[which(full.breakdown$cluster==i),j])[[2]]
    upper.quant.breakdown.group[[j]] <- quantile(full.breakdown[which(full.breakdown$cluster==i),j])[[4]]
    sd.breakdown.group[[j]] <- sd(full.breakdown[which(full.breakdown$cluster==i),j])
  }
  breakdown.bygroup[[i]] <- data.frame(mean=mean.breakdown.group,sd=sd.breakdown.group,median=median.breakdown.group,
                                       upper_quant = upper.quant.breakdown.group, lower_quant = lower.quant.breakdown.group,
                                       positive=positive.group,row.names = features.list)
  breakdown.bygroup[[i]] <- breakdown.bygroup[[i]][order(-breakdown.bygroup[[i]]$mean),]
  write.csv(breakdown.bygroup[[i]],file = paste(paste(BD_DIR,paste("local_breakdown_WSS_cluster",i,sep = '-'),sep = "/"),
                                                "csv",sep = "."))
}
#write.csv(breakdown.df.bygroup,file = paste(BD_DIR,"breakdown_clustering_6-15-21.csv",sep = "/"))

#### CP plots for X site substition
X = seq(0,33,0.01)
Si =33-X
Mn = rep(21.33,length(X))
Fe = rep(21.33,length(X))
Ni = rep(21.33,length(X))
Co = rep(0,length(X))
Ge = rep(0,length(X))
Sn = rep(0,length(X))
Al = rep(0,length(X))
Ga = rep(0,length(X))

df.cp <- data.frame(Mn,Fe,Co,Ni,Si,Ge,Sn,Al,Ga)
X_sub =  c('Ge','Sn','Al','Ga')
for(i in 6:9){
  CP.predictions <- data.frame(base=X)
  df.temp <- df.cp
  df.temp[,i] <- X
  CP.predictions <- cbind(CP.predictions,fn.mean.svr.predict(best_model,df.temp))
  CP.predictions <-cbind(CP.predictions,fn.sd.svr.predict(best_model,df.temp))
  df.temp[,i] <- 0
  colnames(CP.predictions) <- c("X","mean","sd")
  write.table(CP.predictions,file = paste(paste(CP_DIR,paste(X_sub[[i-5]],"distribution",sep='-'),sep = '/'),"dat",sep = ".")
              ,row.names = FALSE)
  write.csv(CP.predictions,file = paste(paste(CP_DIR,paste(X_sub[[i-5]],"distribution",sep='-'),sep = '/'),"csv",sep = ".")
              ,row.names = FALSE)
}

#### Fe at M vs Fe at M' vs Fe at both
#at M
Fe = seq(5,30,0.01)
Si = rep(33,length(Fe))
Mn = 33-Fe
Ni = Si
Co = rep(0,length(Fe))
Ge = rep(0,length(Fe))
Sn = rep(0,length(Fe))
Al = rep(0,length(Fe))
Ga = rep(0,length(Fe))
df.cp <- data.frame(Mn,Fe,Co,Ni,Si,Ge,Sn,Al,Ga)
Cp.Fe.Mn <- data.frame(base = Fe)
Cp.Fe.Mn <- cbind(Cp.Fe.Mn,fn.mean.svr.predict(best_model,df.cp))
Cp.Fe.Mn <-cbind(Cp.Fe.Mn,fn.sd.svr.predict(best_model,df.cp))
colnames(Cp.Fe.Mn) <- c('X','mean','sd')
write.table(Cp.Fe.Mn,file = paste(CP_DIR,'Fe-Mn.dat',sep = '/'),row.names = FALSE)
write.csv(Cp.Fe.Mn,file = paste(CP_DIR,'Fe-Mn.csv',sep = '/'),row.names = FALSE)
#At T
Mn = Si
Ni = 33-Fe
df.cp <- data.frame(Mn,Fe,Co,Ni,Si,Ge,Sn,Al,Ga)
Cp.Fe.Ni <- data.frame(base = Fe)
Cp.Fe.Ni <- cbind(Cp.Fe.Ni,fn.mean.svr.predict(best_model,df.cp))
Cp.Fe.Ni <-cbind(Cp.Fe.Ni,fn.sd.svr.predict(best_model,df.cp))
colnames(Cp.Fe.Ni) <- c('X','mean','sd')
write.table(Cp.Fe.Ni,file = paste(CP_DIR,'Fe-Ni.dat',sep = '/'),row.names = FALSE)
write.csv(Cp.Fe.Ni,file = paste(CP_DIR,'Fe-Ni.csv',sep = '/'),row.names = FALSE)
#At both
Fe = seq(5,30,0.01)
Mn = 33-Fe/2
Ni = 33-Fe/2
df.cp <- data.frame(Mn,Fe,Co,Ni,Si,Ge,Sn,Al,Ga)
Cp.Fe.Mn.Ni <- data.frame(base = Fe)
Cp.Fe.Mn.Ni <- cbind(Cp.Fe.Mn.Ni,fn.mean.svr.predict(best_model,df.cp))
Cp.Fe.Mn.Ni <-cbind(Cp.Fe.Mn.Ni,fn.sd.svr.predict(best_model,df.cp))
colnames(Cp.Fe.Mn.Ni) <- c('X','mean','sd')
write.table(Cp.Fe.Mn.Ni,file = paste(CP_DIR,'Fe-Mn-Ni.dat',sep = '/'),row.names = FALSE)
write.csv(Cp.Fe.Mn.Ni,file = paste(CP_DIR,'Fe-Mn-Ni.csv',sep = '/'),row.names = FALSE)
#plotting
#for (i in 1:nrow(reduced.dataset)){
#  break.down <- predict_parts_break_down(explainer = explainer_mean,new_observation = observations_to_check[i,])
#  cp <- predict_profile(explainer = explainer_mean,new_observation = observations_to_check[i,])
#  p1 <- plot(break.down)
#  p2 <- plot(cp)
#  png(paste(paste(BD_DIR,paste('comp_breakdown_6-6-21',rownames(observations_to_check)[i], sep = '-'),sep = '/'),'png',sep = '.'),
#      width = 1200,height = 822)
#  grid.arrange(p1,nrow=1,ncol=1)
#  dev.off()
#  png(paste(paste(CP_DIR,paste('comp_CP_6-6-21',rownames(observations_to_check)[i], sep = '-'),sep = '/'),'png',sep = '.'),
#      width = 846,height = 580)
#  grid.arrange(p2,nrow=1,ncol=1)
#  dev.off()
#}
