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


#feature reduction by correlation anaysis (pearson)
#theating
library(readxl)
comp.dataset <- read_xlsx(paste(DATA_DIR,'composition_5_27_21.xlsx',sep = '/'))
reduced.dataset <- as.data.frame(na.omit(comp.dataset[,-c(1,14)]))
magpie.elements <- read.csv(paste(DATA_DIR,'element_magpie.csv',sep = '/'),row.names = 'X')
magpie.dataset <- fn.composition.to.magpie(reduced.dataset[,-c(10:12)],magpie.elements)
magpie.dataset$response <- reduced.dataset[,10]

#test train split
set.seed(42)
sample = sample.int(n = nrow(reduced.dataset),size = floor(0.8*nrow(reduced.dataset)),replace = F)
train = magpie.dataset[sample,]
test = magpie.dataset[-sample,]

#heating vs cooling, models were constructed and validated as functions of theating not tcooling.
theating.train <- train[,1:ncol(train)]
theating.test <- test[,1:ncol(test)]
colnames(theating.train)[ncol(theating.train)] <- 'response'
colnames(theating.test)[ncol(theating.test)] <- 'response'

#RandomForest
#Random Forest training
library(randomForest)
library(pROC)

fn.train.rf = function(data,ntree,mtry){
  oob = NULL
  mtry_values = NULL
  tree_values = NULL
  print("optimizing hyperparamters ntree and mtry")
  i = 1
  for (m in 2:mtry) {## please check if 9 is our dimensionalit
    print(paste('mtree = ', m))
    initial.tree = 100
    for (tree in seq(from = initial.tree, to = ntree, by = 20)){
      print(paste('ntree = ', tree))
      model = randomForest(response ~ ., data = data, ntree = tree, mtry = m, importance = TRUE)
      oob[[i]] = mean(model$mse)
      mtry_values[[i]] = m
      tree_values[[i]] = tree
      i = i +1
    }
  }
  model.performance <<- data.frame(mtry_values,tree_values,oob)
  ordered.model.performance <<- model.performance[order(model.performance$oob,decreasing = FALSE),]
  optimized.model = randomForest(response ~ ., data = data, ntree = ordered.model.performance[1,2],
                                 mtry = ordered.model.performance[1,1])
  return(optimized.model)
}
model <- fn.train.rf(theating.train,ntree=2000,mtr=6)

save(model, file = paste(MODEL_DIR,"rforest_magpie.8-10-21.rda",sep = '/'))
set.seed(42)

##performance
load(paste(MODEL_DIR,"rforest_magpie.8-10-21.rda",sep = '/'))
rf.train.predictions = predict(model,theating.train,type = 'response',predict.all = TRUE)
individual_trees = data.frame(rf.train.predictions[2])
sd = NULL
for (i in 1:nrow(individual_trees)){
  sd[[i]] = sd(individual_trees[i,])
}
rf.train.performance = data.frame(theating.train[,6],rf.train.predictions[1],sd)
colnames(rf.train.performance) <- c('Theating','prediction.mean','prediction.sd')
rf.train.performance$residuals <- (rf.train.performance$Theating-rf.train.performance$prediction.mean)**2
train.rmse = sqrt(mean(rf.train.performance$residuals))
tran.rsqr = cor(rf.train.performance$Theating,rf.train.performance$prediction.mean)**2

rf.test.predictions = predict(model,theating.test,type = 'response',predict.all = TRUE)
individual_trees = data.frame(rf.test.predictions[2])
sd = NULL
for (i in 1:nrow(individual_trees)){
  sd[[i]] = sd(individual_trees[i,])
}
rf.test.performance = data.frame(theating.test[,6],rf.test.predictions[1],sd)
colnames(rf.test.performance) <- c('Theating','prediction.mean','prediction.sd')
rf.test.performance$residuals = (rf.test.performance$Theating-rf.test.performance$prediction.mean)**2
test.rmse = sqrt(mean(rf.test.performance$residuals))
test.rsqr = cor(rf.test.performance$Theating,rf.test.performance$prediction.mean)**2

rf.performance = rbind(rf.train.performance,rf.test.performance)
write.csv(rf.performance,paste(PLOT_DIR,'RF_magpie_8-10-21.csv',sep = '/'))
