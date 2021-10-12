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

### functions #####
fn.boot.lasso = function(data,seed,B,model_name,model_directory) {
  ptm <- proc.time()
  set.seed(seed)
  n = dim(data)[1]
  lasso.models<- NULL
  boot.sample <- NULL
  print(paste('training models 1-',str(B),sep = ''))
  for(i in 1:B)
  {
    print(i)
    boot.sample = sample(n, replace = TRUE)
    boot.sample = data[boot.sample,]
    y = boot.sample[,length(data)]
    x = boot.sample[,-length(data)]
    x = as.matrix(x)
    lasso_reg <- cv.glmnet(x, y, alpha=1, nfolds=10)
    lambda.best = lasso_reg$lambda.min
    lasso.models[[i]] = glmnet(x, y, alpha = 1, lambda = lambda.best, standardize = TRUE)
  }
  save(lasso.models, file = paste(model_directory,model_name,sep = '/'))
  print('complete')
  print(proc.time() - ptm)
  return(lasso.models)
}

fn.predict.lasso = function(model,data){
  insample.features = as.matrix(data[,-ncol(data)])
  predict.insample = NULL
  i=1
  repeat{
    predict.insample[[i]] = predict(model[[i]], insample.features)
    i=i+1
    if(i>100) break ()
  }
  predict.insample.df = as.data.frame(predict.insample)
  predict.insample.mean.df = rowMeans(as.matrix(predict.insample.df))
  predict.insample.sd.df = rowSds(as.matrix(predict.insample.df))
  data.insample.final = cbind(data, predict.insample.mean.df, predict.insample.sd.df)
  return(data.insample.final)
}

fn.feature.importance = function(model){#returns percentage of models which use a given feature
  usage = integer(length(coef(model[[1]])))
  for (i in 1:length(model)){
    coefficients <- as.array(coef(model[[i]]))
    logic = as.numeric(sapply(coefficients,as.logical))
    usage = usage+logic
  }
  usage = data.frame(usage/length(model)*100,row.names = rownames(coef(model[[1]])))
  return(usage)
}

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

#directroy definitions
DATA_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Data'
MODEL_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Models'
PLOT_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/Plots'
ANALYSIS_DIR = '/home/timothy/Desktop/MnNiSi/TransitionTemp_3-5-21/analysis'

#DATA
#feature reduction by material science assumptions
dft.dataset <- read.csv(paste(DATA_DIR,'DFT_combined_4-7-21.csv',sep = '/'))
raw.features <- dft.dataset[,2:31]
reduced = na.omit(raw.features[,c(5:6,10:12,22,25,27,28)])
colnames(reduced)[9] = "response"

#test train split
set.seed(42)
sample = sample.int(n = nrow(reduced),size = floor(0.8*nrow(reduced)),replace = F)
train = reduced[sample,]
test = reduced[-sample,]

#heating vs cooling, models were constructed and validated as functions of theating not tcooling.
theating.train <- train[,1:ncol(train)]
theating.test <- test[,1:ncol(test)]
colnames(theating.train)[ncol(theating.train)] <- 'response'
colnames(theating.test)[ncol(theating.test)] <- 'response'

#bootstrapped lasso
library(glmnet)
library(matrixStats)
theating.model.lasso = fn.boot.lasso(theating.train,42,1000,'8-10-21.lasso.dft.rda',model_directory = MODEL_DIR)

#performance
load(paste(MODEL_DIR,'8-10-21.lasso.magpie.rda',sep = '/'))

lasso.train.predictions <- fn.predict.lasso(theating.model.lasso,theating.train)[,9:11]
colnames(lasso.train.predictions)[1] <- 'T_heating'
lasso.train.predictions$residuals = (lasso.train.predictions$T_heating-lasso.train.predictions$predict.insample.mean.df)**2
lasso.train.mse = mean(lasso.train.predictions$residuals)
train.rsqr = cor(lasso.train.predictions$T_heating, lasso.train.predictions$predict.insample.mean.df)**2

lasso.test.predictions <- fn.predict.lasso(theating.model.lasso,theating.test)[,9:11]
colnames(lasso.test.predictions)[1] <- 'T_heating'
lasso.test.predictions$residuals = (lasso.test.predictions$T_heating-lasso.test.predictions$predict.insample.mean.df)**2
lasso.test.mse = mean(lasso.test.predictions$residuals)
test.rsqr = cor(lasso.test.predictions$T_heating, lasso.test.predictions$predict.insample.mean.df)**2
lasso.performance = rbind(lasso.train.predictions, lasso.test.predictions)

write.csv(lasso.performance, file = paste(ANALYSIS_DIR,'best_dft_model_lasso_8-10-21.csv',sep = '/'))
