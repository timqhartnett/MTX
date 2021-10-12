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

#directroy definitions
DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data'
MODEL_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Models'
PLOT_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Plots'
ANALYSIS_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/analysis'

fn.composition.to.magpie = function(data,magpie_features,selected_features){#function converts from composition to Magpie features using magpie elements CSV
  magpie.dataframe = data.frame(matrix(NA,nrow = nrow(data), ncol=length(selected_features))) #create blank dataframe for weighted features
  colnames(magpie.dataframe) = selected_features
  mean.features = selected_features
  for (i in 1:length(mean.features)){
    feature.vector = numeric(nrow(data))
    for (j in 1:ncol(data)){
      feature.vector = feature.vector + data[,j]*magpie_features[colnames(data)[j],mean.features[i]]
    }
    magpie.dataframe[,mean.features[i]] <- feature.vector/100
  }
  return(magpie.dataframe)
}

fn.temp.predict = function(data, model){#function for ensemble prediction of property (temp/hyst)
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

#data
library(readxl)
library(e1071)
library(matrixStats)
composition <- as.data.frame(read_excel(paste(DATA_DIR,'Radhika_comp_5-28-21.xlsx',sep = '/')))
validation.composition <- composition[,2:10]
composition.82721 <- as.data.frame(read_excel(paste(DATA_DIR,'VCU_comp_8-27-21.xlsx',sep = '/')))
magpie.82721 <- as.data.frame(read_excel(paste(DATA_DIR,'VCU_magpy_comp_8-27-21.xlsx',sep = '/')))
validation.composition.82721 <- na.omit(composition.82721)[,2:12]

magpie.elements <- read.csv('/home/timothy/magpy_elements.csv',row.names = 'Symbol')
selected_features <- c('CovalentRadius','MeltingT','AtomicWeight','MendeleevNumber','NdValence')
validation.magpie <- na.omit(fn.composition.to.magpie(validation.composition,magpie.elements,selected_features))
validation.82721.magpie <- na.omit(fn.composition.to.magpie(magpie.82721[,-c(1,ncol(magpie.82721),ncol(magpie.82721)-1)],
                                                            magpie.elements,selected_features))
colnames(validation.magpie) <- c('mean_CovalentRadius','mean_MeltingT','mean_AtomicWeight','mean_MendeleevNumber','mean_NdValence')
#validation.82721.magpie <- fn.composition.to.magpie(validation.composition,magpie.elements)
#### Composition ####
#Heating
load(paste(MODEL_DIR,'best.svr.comp.5-27-21.rda',sep = '/'))
comp.predictions.theating <- fn.temp.predict(validation.composition,best_model)
comp.82721.predictions.theating <- fn.temp.predict(validation.composition.82721,best_model)

#hysteresis
load(paste(MODEL_DIR,'hyst.svr.comp.5-27-21.rda',sep = '/'))
comp.predictions.hyst <- fn.temp.predict(validation.composition,best_model)
comp.827.21.predictions.hyst <- fn.temp.predict(validation.composition.82721,best_model)
### Magpie ###
load(paste(MODEL_DIR,'best.svr.magpie.6-1-21.rda',sep = '/'))
magpie.predictions.theating <- fn.temp.predict(validation.magpie,best_model)
validation.82721.magpie.theating <- cbind(magpie.82721[,1],fn.temp.predict(validation.82721.magpie,best_model))
write.csv(validation.82721.magpie.theating,'/home/timothy/VCU_composition_magpie_8-28-21.csv')

load(paste(MODEL_DIR,'best.hyst.magpie.6-1-21.rda',sep = '/'))
magpie.predictions.hysteresis <- fn.temp.predict(validation.magpie,best_model)
validation.82721.magpie.hyst <- cbind(magpie.82721[,1],fn.temp.predict(validation.82721.magpie,best_model))
write.csv(validation.82721.magpie.hyst,'/home/timothy/VCU_magpie_hyst_8-28-21.csv')
