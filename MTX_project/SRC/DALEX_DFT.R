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

for (i in 1:nrow(dft.dataset)){
  break.down <- predict_parts_break_down(explainer = explainer_mean,new_observation = reduced[i,])
  cp <- predict_profile(explainer = explainer_mean,new_observation = reduced[i,])
  p1 <- plot(break.down)
  p2 <- plot(cp)
  png(paste(paste(BD_DIR,paste('DFT_breakdown_9-20-21',dft.dataset$Compostion[which(dft.dataset$Theating==reduced$response[i])], sep = '-'),sep = '/'),'png',sep = '.'),
      width = 1200,height = 822)
  grid.arrange(p1,nrow=1,ncol=1)
  dev.off()
  png(paste(paste(CP_DIR,paste('DFT_CP_9-20-21',dft.dataset$Compostion[which(dft.dataset$Theating==reduced$response[i])], sep = '-'),sep = '/'),'png',sep = '.'),
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
