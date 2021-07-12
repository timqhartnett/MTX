#tim Validation Dataset prediction and comparisons MTX compounds Multi-caloric Magnets
#4-7-21
#University of Virginia
#Charlottesville, VA 22903

##### CLEAN WORKSPACE ########
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014")
# Clean workspace
rm(list=ls())

#### directroy definitions ####
DATA_DIR = paste(getwd(),'DATA',sep='/')
MODEL_DIR = paste(getwd(),'MODELS',sep='/')
PLOT_DIR = paste(getwd(),'PLOTS',sep='/')
ANALYSIS_DIR = paste(getwd(),'ANALYSIS',sep='/')

#### libraries ####
library(matrixStats) #used for rowmeans and sd function
library(readxl)
library(e1071) #contains svr algorithms and prediction functions

#### function definitions ####
fn.composition.to.magpie = function(data,magpie_features){
  magpie.dataframe = data.frame(matrix(NA,nrow = nrow(data), ncol=ncol(magpie_features))) #create blank dataframe for weighted features
  colnames(magpie.dataframe) = colnames(magpie_features)
  mean.features = colnames(magpie_features)
  for (i in 1:length(mean.features)){
    feature.vector = numeric(nrow(data))
    for (j in 1:ncol(data)){
      feature.vector = feature.vector + data[,j]*magpie_features[colnames(data)[j],mean.features[i]]
    }
    magpie.dataframe[,mean.features[i]] <- feature.vector/100
  }
  return(magpie.dataframe)
}

fn.temp.predict = function(data, model){
  predict.data = NULL
  i=1
  repeat{
    predict.data[[i]] = predict(model[[i]], data)
    i=i+1
    if(i>length(model)) break ()
  }
  predict.data.df = as.data.frame(predict.data)
  predict.data.mean = rowMeans(as.matrix(predict.data.df))
  predict.data.sds = rowSds(as.matrix(predict.data.df))
  data.predictions <- as.data.frame(cbind(predict.data.mean, predict.data.sds))
  names(data.predictions) <- c("prediction.mean","prediction.sd")
  return(data.predictions)
}

#### load data ####
composition <- as.data.frame(read_excel(paste(DATA_DIR,'Radhika_comp_5-28-21.xlsx',sep = '/')))
validation.composition <- composition[,-c(1,11:12)]
magpie.elements <- read.csv(paste(DATA_DIR,'element_magpie.csv',sep = '/'),row.names = 'X')
validation.magpie <- fn.composition.to.magpie(validation.composition,magpie.elements)

#### Composition ####
#Heating
load(paste(MODEL_DIR,'best.svr.comp.5-27-21.rda',sep = '/'))
comp.predictions.theating <- fn.temp.predict(validation.composition,best_model)

#hysteresis
load(paste(MODEL_DIR,'hyst.svr.comp.5-27-21.rda',sep = '/'))
comp.predictions.hyst <- fn.temp.predict(validation.composition,best_model)

### Magpie ###
load(paste(MODEL_DIR,'best.svr.magpie.6-1-21.rda',sep = '/'))
magpie.predictions.theating <- fn.temp.predict(validation.magpie,best_model)

load(paste(MODEL_DIR,'best.hyst.magpie.6-1-21.rda',sep = '/'))
magpie.predictions.hysteresis <- fn.temp.predict(validation.magpie,best_model)
