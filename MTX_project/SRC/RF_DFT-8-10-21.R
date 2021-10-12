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

save(model, file = paste(MODEL_DIR,"rforest_DFT.8-10-21.rda",sep = '/'))
set.seed(42)

##performance
load(paste(MODEL_DIR,"rforest_DFT.8-10-21.rda",sep = '/'))
rf.train.predictions = predict(model,theating.train,type = 'response',predict.all = TRUE)
individual_trees = data.frame(rf.train.predictions[2])
sd = NULL
for (i in 1:nrow(individual_trees)){
  sd[[i]] = sd(individual_trees[i,])
}
rf.train.performance = data.frame(theating.train[,9],rf.train.predictions[1],sd)
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
rf.test.performance = data.frame(theating.test[,9],rf.test.predictions[1],sd)
colnames(rf.test.performance) <- c('Theating','prediction.mean','prediction.sd')
rf.test.performance$residuals = (rf.test.performance$Theating-rf.test.performance$prediction.mean)**2
test.rmse = sqrt(mean(rf.test.performance$residuals))
test.rsqr = cor(rf.test.performance$Theating,rf.test.performance$prediction.mean)**2

rf.performance = rbind(rf.train.performance,rf.test.performance)
write.csv(rf.performance,paste(PLOT_DIR,'RF_magpie_8-10-21.csv',sep = '/'))