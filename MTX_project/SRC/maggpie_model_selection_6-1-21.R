#tim DFT features models MTX compounds Multi-caloric Magnets
#check against multiple training sets for Magpie features
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

#### Support Vector Machines functions for bootstrap optimized model training ####
#libraries
library(matrixStats) # for rowSds
library(e1071) # the e1071 library contains implementations for a number of statistical learning methods.
library(parallel)

fn.assign.response = function(data,response.name){
  colnumber = which(colnames(data)==response.name)
  colnames(data)[colnumber] <- 'response'
  return(data)
}

fn.cv.svr.rbf = function(data){
  tuneResult.svr.rbf = tune(svm, response ~ ., data = data,
                            kernel = "radial",
                            ranges = list(gamma = c(0.001,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), 
                                          cost = c(0.001, 0.01, 0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100)),
                            tunecontrol = tune.control(cross = 9))
  return(tuneResult.svr.rbf$best.parameters)
}

fn.train = function(data){
  parameter.svr.rbf = fn.cv.svr.rbf(data)
  svr.rbf.ecorr = svm(response ~ ., data = data, kernel = "radial", 
                      gamma =parameter.svr.rbf[1], cost = parameter.svr.rbf[2])
  return(svr.rbf.ecorr)
}

fn.svr.predict = function(model,data){#SVR bootstrapped prediction function is required, 
  predict.data = NULL
  i=1
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

fn.rsqrd = function(model,data){
  predictions <- fn.svr.predict(model,data)
  Rsqrd <- cor(predictions,data[,ncol(data)])**2
  return(Rsqrd)
}

fn.boot.mean.oob = function(data,number_cores) {
  values <- c(15,25,50,75,100,125,150)
  oob.mean.err <- NULL
  oob.sd.err <- NULL
  for (i in 1:length(values)){
    print(values[i])
    oob <- NULL
    data.boot<-NULL
    for (j in 1:values[i]){
      n = dim(data)[1]
      bootsamples <- sample(n, replace = TRUE)
      data.boot[[j]] <- data[bootsamples,]
      oob[[j]] <- data[-bootsamples,]
    }
    oob <<- oob
    data.boot <<- data.boot
    oob.r2 <- NULL
    models <- mclapply(data.boot,fn.train,mc.cores = number_cores)
    models <<- models
    for (j in 1:length(models)){
      oob.r2 <- cor(oob[[j]][,which(colnames(oob[[j]])=='response')],predict(models[[j]],oob[[j]]))**2
    }
    oob.r2<<-oob.r2
    oob.mean.err[[i]] <- mean(oob.r2)
    oob.sd.err[[i]] <- sd(oob.r2)
    oob.mean <<- oob.mean.err
    if (i ==1 ){
      model <- models
    }else if (mean(oob.r2) == max(oob.mean.err)){
      model <- models
    }
  }
  oob.err.total <- as.data.frame(cbind(values,oob.mean.err,oob.sd.err))
  colnames(oob.err.total) <- c('number_bootstraps','mean R^2 error','SD R^2 error')
  oob.err.total <<- oob.err.total
  return(model)
} 

fn.ensemble.predict = function(data, model){
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

fn.k.crossfold = function(data,k){
  new_order = sample.int(nrow(data),nrow(data))
  shuffled.data = data[new_order,]
  k.size = round(nrow(data)/(k+1))
  k.samples = NULL
  training.samples = NULL
  breaks = seq(1,nrow(data),k.size)
  for (i in 1:(length(breaks)-1)){
    samples = breaks[i]:breaks[i+1]
    k.samples[[i]] = data[samples,]
    training.samples[[i]] = data[-samples,]
  }
  cross.fold = c(training.samples,k.samples)
  return(cross.fold)
}

fn.kfold.oob.optimization = function(data,response.name,k,number.cores){
  data <- fn.assign.response(data, response.name)
  kfold.split <<- fn.k.crossfold(data,k)
  performance.k = NULL
  models = NULL
  for (i in 1:k){
    print(i)
    model.data = kfold.split[[i]]
    validation.data = kfold.split[[i+k]]
    models[[i]] <- fn.boot.mean.oob(model.data,number.cores)
    performance.k[[i]] <- (cor(validation.data$response,fn.ensemble.predict(validation.data,models[[i]])[[1]]))**2
  }
  performance.k <<- performance.k
  models <<- models
  best.model = models[which.max(performance.k)]
  info = list(best.model,performance.k)
  return(info)
}

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

#feature selection
library(readxl)
comp.dataset <- as.data.frame(read_xlsx(paste(DATA_DIR,'composition_5_27_21.xlsx',sep = '/')))
magpie.elements <- read.csv('/home/timothy/magpy_elements.csv',row.names = 'Symbol')
magpie.elements <- magpie.elements[,-ncol(magpie.elements)]
selected_features = colnames(magpie.elements)
magpie <- fn.composition.to.magpie(comp.dataset[,-c(1,(ncol(comp.dataset)-3):ncol(comp.dataset))],
                                   magpie.elements,selected_features)
first.to.drop <- c('Number','Row','Column','SpaceGroupNumber')
magpie = magpie[,-which(colnames(magpie) %in% first.to.drop)]
corr.magpie.total = cor(magpie)
na.to.drop = rownames(corr.magpie.total[is.na(corr.magpie.total[,1]),])
magpie.2 <- magpie[,-which(colnames(magpie) %in% na.to.drop)]
corr.magpie.2 <- cor(magpie.2)
highly.correlated = NULL
for (i in 1:ncol(magpie.2)){
  for (j in 1:ncol(magpie.2)){
    if (corr.magpie.2[i,j]>0.8){
      highly.correlated = c(highly.correlated,paste(colnames(magpie.2)[i],colnames(magpie.2)[j],corr.magpie.2[i,j],sep = '-'))
    }
  }
}
next.to.drop <- c('NsValence','NpValence','NdValence',colnames(magpie.2)[10],'sValence',colnames(magpie.2)[11],'NUnfilled','AtomicWeight','GSbandgap')
magpie3 <- magpie.2[,-which(colnames(magpie.2) %in% next.to.drop)]
corr.magpie.3 <- cor(magpie3)
highly.correlated = NULL
for (i in 1:ncol(magpie3)){
  for (j in 1:ncol(magpie3)){
    if (corr.magpie.3[i,j]>0.8){
      highly.correlated = c(highly.correlated,paste(colnames(magpie3)[i],colnames(magpie3)[j],corr.magpie.3[i,j],sep = '-'))
    }
  }
}

fn.scale = function(data){
  feature.means = NULL
  feature.sds = NULL
  for (i in 1:ncol(data)){
    feature.means[[i]] <- mean(data[,i])
    feature.sds[[i]] <- sd(data[,i])
    data[,i] <- (data[,i]-feature.means[[i]])/feature.sds[[i]]
  }
  output <- list(data,feature.means,feature.sds)
  return(output)
}
scaled <- fn.scale(magpie3)
mapgie.3.scaled <- scaled[[1]]
magpie.dataset.heating = as.data.frame(cbind(mapgie.3.scaled,comp.dataset$Theating)[-c(78,79),])
colnames(magpie.dataset.heating)[8] = 'Theating'

magpie.dataset.hyst = na.omit(cbind(mapgie.3.scaled,comp.dataset$hysteresis))
colnames(magpie.dataset.heating)[8] = 'hysteresis'
#theating
k.fold.heating <- fn.kfold.oob.optimization(magpie.dataset.heating,response.name = 'Theating',k=10,number.cores = 8)

#training
sample_performance <- fn.robust.ttsplit(magpie.dataset,number_of_seeds = 20,number_of_cores = 10)

print(which(sample_performance==max(sample_performance)))
rownames(test_sets[[which(sample_performance==max(sample_performance))]])
#"6"  "8"  "9"  "11" "13" "14" "21" "23" "28" "45" "48" "50" "57"
best_model <- models[[which(sample_performance==max(sample_performance))]]
save(best_model,file = paste(MODEL_DIR,'best.svr.magpie.6-1-21.rda',sep = '/'))
best_sample <- c("6","8","9","11","13","14","21","23","28","45","48","50","57")

theating.train <- magpie.dataset[!(rownames(magpie.dataset)) %in% best_sample,]
theating.test <- magpie.dataset[best_sample,]

### model performance ####
load(paste(MODEL_DIR,'best.svr.magpie.6-1-21.rda',sep = '/'))
svr.train.predictions = fn.temp.predict(theating.train,model = best_model)
svr.train.performance = cbind(theating.train$response,svr.train.predictions)
colnames(svr.train.performance) <- c('Theating','prediction.mean','prediction.sd')
svr.train.performance$residuals = (svr.train.performance$Theating - svr.train.performance$prediction.mean)**2
RMSE.train = sqrt(sum(svr.train.performance$residuals)/nrow(svr.train.performance))
Rsqrd = cor(theating.train$response,svr.train.predictions$prediction.mean)^2

svr.test.predictions = fn.temp.predict(theating.test,model = best_model)
svr.test.performance = cbind(theating.test$response,svr.test.predictions)
colnames(svr.test.performance) <- c('Theating','prediction.mean','prediction.sd')
svr.test.performance$residuals = (svr.test.performance$Theating - svr.test.performance$prediction.mean)**2
RMSE.test = sqrt(sum(svr.test.performance$residuals)/nrow(svr.test.performance))
Rsqrd.test = cor(svr.test.predictions$prediction.mean, theating.test$response)
svr.performance.total = rbind(svr.train.performance,svr.test.performance)

#Plotting
library(ggplot2)
labels = c(rep('train',51),rep('test',13))
svr.performance.total$label = labels

p = ggplot(svr.performance.total,aes(x=Theating,y=prediction.mean,col=label))+geom_point()
p_error = p+geom_errorbar(data = svr.performance.total,aes(x=Theating, ymin=prediction.mean-prediction.sd, ymax=prediction.mean+prediction.sd, width=0.25))
p_final = p_error+geom_abline()
p_final

#write for xmgrace plotting
write.csv(svr.performance.total,paste(ANALYSIS_DIR,'bestsvr_Magpie_theating_6-6-21.csv',sep = '/'))

#hysteresis
magpie.elements <- read.csv(paste(DATA_DIR,'element_magpie.csv',sep = '/'),row.names = 'X')
magpie.dataset <- fn.composition.to.magpie(reduced.dataset[,-c(10:12)],magpie.elements)
magpie.dataset$response <- reduced.dataset[,12]

sample_performance <- fn.robust.ttsplit(magpie.dataset,number_of_seeds = 20,number_of_cores = 10)

print(which(sample_performance==max(sample_performance)))
rownames(test_sets[[which(sample_performance==max(sample_performance))]])
#"2"  "13" "17" "24" "26" "29" "34" "35" "38" "45" "46" "54" "62"
best_model <- models[[which(sample_performance==max(sample_performance))]]
save(best_model,file = paste(MODEL_DIR,'best.hyst.magpie.6-1-21.rda',sep = '/'))
best_sample <- c("2","13","17","24","26","29","34","35","38","45","46","54","62")

theating.train <- magpie.dataset[!(rownames(magpie.dataset)) %in% best_sample,]
theating.test <- magpie.dataset[best_sample,]

### model performance ####
load(paste(MODEL_DIR,'best.hyst.magpie.6-1-21.rda',sep = '/'))
svr.train.predictions = fn.temp.predict(theating.train,model = best_model)
svr.train.performance = cbind(theating.train$response,svr.train.predictions)
colnames(svr.train.performance) <- c('Theating','prediction.mean','prediction.sd')
svr.train.performance$residuals = (svr.train.performance$Theating - svr.train.performance$prediction.mean)**2
RMSE.train = sqrt(sum(svr.train.performance$residuals)/nrow(svr.train.performance))
Rsqrd = cor(theating.train$response,svr.train.predictions$prediction.mean)^2

svr.test.predictions = fn.temp.predict(theating.test,model = best_model)
svr.test.performance = cbind(theating.test$response,svr.test.predictions)
colnames(svr.test.performance) <- c('Theating','prediction.mean','prediction.sd')
svr.test.performance$residuals = (svr.test.performance$Theating - svr.test.performance$prediction.mean)**2
RMSE.test = sqrt(sum(svr.test.performance$residuals)/nrow(svr.test.performance))
Rsqrd.test = cor(svr.test.predictions$prediction.mean, theating.test$response)
svr.performance.total = rbind(svr.train.performance,svr.test.performance)
write.csv(svr.performance.total,paste(ANALYSIS_DIR,'bestsvr_Magpie_hsyt_6-6-21.csv',sep = '/'))

#Plotting
library(ggplot2)
labels = c(rep('train',51),rep('test',13))
svr.performance.total$label = labels

p = ggplot(svr.performance.total,aes(x=Theating,y=prediction.mean,col=label))+geom_point()
p_error = p+geom_errorbar(data = svr.performance.total,aes(x=Theating, ymin=prediction.mean-prediction.sd, ymax=prediction.mean+prediction.sd, width=0.25))
p_final = p_error+geom_abline()
p_final


