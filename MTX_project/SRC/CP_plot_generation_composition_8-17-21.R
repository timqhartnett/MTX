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
CP_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/CP_Plots/homebrew_8-17-21'
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

predictions <- cbind(reduced.dataset[,1:9],fn.mean.svr.predict(best_model,reduced.dataset))
colnames(predictions)[10] <- 'Theating'
write.csv(predictions,file = paste(CP_DIR,'predicted_composition_8-18-21.csv',sep ='/'))


features.list <- colnames(reduced.dataset[,-ncol(reduced.dataset)])
for (i in 1:nrow(reduced.dataset)){
  print(i)
  for (j in 1:length(features.list)){
    print(features.list[j])
    X = seq(min(reduced.dataset[features.list[j]]),max(reduced.dataset[features.list[j]]),0.01)
    Mn = rep(reduced.dataset[i,which(features.list=='Mn')],length(X))
    Fe = rep(reduced.dataset[i,which(features.list=='Fe')],length(X))
    Co = rep(reduced.dataset[i,which(features.list=='Co')],length(X))
    Ni = rep(reduced.dataset[i,which(features.list=='Ni')],length(X))
    Si = rep(reduced.dataset[i,which(features.list=='Si')],length(X))
    Ge = rep(reduced.dataset[i,which(features.list=='Ge')],length(X))
    Sn = rep(reduced.dataset[i,which(features.list=='Sn')],length(X))
    Al = rep(reduced.dataset[i,which(features.list=='Al')],length(X))
    Ga = rep(reduced.dataset[i,which(features.list=='Ga')],length(X))
    virtual.dataset = data.frame(Mn,Fe,Co,Ni,Si,Ge,Sn,Al,Ga)
    virtual.dataset[features.list[j]] <- X
    virtual.dataset.mean.predictions <- fn.mean.svr.predict(best_model,virtual.dataset)
    virtual.dataset.sd.predictions <- fn.sd.svr.predict(best_model,virtual.dataset)
    cp.dataset <- data.frame(X,virtual.dataset.mean.predictions,virtual.dataset.sd.predictions)
    colnames(cp.dataset) <- c(features.list[j],'mean','sd')
    write.csv(cp.dataset,file = paste(paste(CP_DIR,paste('observation',i,features.list[j],sep = '-'),sep ='/'),'csv',sep = '.'))
  }
}
