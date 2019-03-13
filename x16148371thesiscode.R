setwd("~/Desktop/NCI/Thesis_Repeat/R")

library(class)
library(caret)
library(randomForest)
library(caTools)
library(corrplot)

#Data Cleaing

#labeling of class variable 
label_bc <- read.csv("labels.csv")
str(label_bc)
label_bc$Class <- factor(label_bc$Class, levels = c('BRCA','KIRC','COAD','LUAD','PRAD'), labels = c(1,2,3,4,5))

#finding correlation between the variables
library(caret)
Breast_Cancer <- read.csv("Breast Cancer.csv",stringsAsFactors = T)
Breast_Cancer$Tumor <- factor(Breast_Cancer$Tumor, levels = c('BRCA','KIRC','COAD','LUAD','PRAD'), labels = c(1,2,3,4,5))
cordata <- Breast_Cancer
sapply(cordata, function(x) sum(is.na(x)))
corm <- cor(cordata[,-1])
cordata_new <- na.omit(corm)
highcor <- findCorrelation(corm,cutoff=0.5)
Breast_Cancer <- Breast_Cancer[,-c(5,17,43,47,49,56,52)] 


write.csv(Breast_Cancer , file= "Breast_Cancer_final1.csv",row.names=FALSE)
################################################################################################################################################################################


dataset <- read.csv("Breast_Cancer_final.csv", stringsAsFactors = T)

#correlation matrix of the dataset
r <- cor(dataset)
corrplot(r, method = "circle", type = "upper")
corrplot.mixed(r)


#checking for missing values.
sapply(dataset,function(x) sum(is.na(x)))

#check for class imbalance
table(dataset$Tumor) 
prop.table(table(dataset$Tumor))
plot(dataset$Tumor)

################################################################################################################################################################################
#factorization of the dependent variable
set.seed(123)
library(caret)
dataset$Tumor <- factor(dataset$Tumor, levels = c(1,2,3,4,5), labels = c('BRCA','KIRC','COAD','LUAD','PRAD'))
model_data <- dataset

################################################################################################################################################################################


#scaling of data
split = sample.split(model_data$Tumor, SplitRatio = 0.8)
train_set = subset(model_data, split == TRUE)
test_set = subset(model_data, split == FALSE)
train_set[-1] = as.data.frame(scale(train_set[-1]))
test_set[-1] = as.data.frame(scale(test_set[-1]))

################################################################################################################################################################################
#Feature selection using Random Forest
forest <- randomForest(x = train_set[-1],y = train_set$Tumor, importance = TRUE,ntree=100)
varImpPlotData <- varImpPlot(forest)  
varImpPlot(forest)
rf <- predict(forest, newdata = test_set, type = "class")
plot(rf, main = "Class distribution", xlab = "Classes", ylab = "Frequency")
confusionMatrix(rf, test_set$Tumor, mode = "prec_recall")


#parameter for tuning are set
tuneParams <- trainControl(method = "cv",number =10,savePredictions = 'final')

#Bagging
#feature slection for bagging
#Kernel SVM 
meanDecGini <- varImpPlotData[, 2]
meanDecGini <- meanDecGini[order(-meanDecGini)]
temp <- c(1:length(meanDecGini))
temp <- temp %% 2 == 1
featuressvm <- names(meanDecGini[temp])
featuressvm

svm_model <- train(train_set[,featuressvm], train_set$Tumor, method="svmRadial", trControl=tuneParams)
pred_svm <- predict(svm_model, newdata = test_set[,featuressvm])
confusionMatrix(pred_svm, test_set$Tumor, mode = "prec_recall")


#kNN
meanDecAcc <- varImpPlotData[, 1]
meanDecAcc <- meanDecAcc[order(-meanDecAcc)]
a <- c(1:length(meanDecAcc))
a <- a %% 2 == 1
featureknn <- names(meanDecAcc[a])
featureknn

knn <- train(train_set[,featureknn],train_set$Tumor,method='knn', trControl=tuneParams)
pred_knn <- predict(knn, newdata = test_set[,featureknn])
confusionMatrix(test_set$Tumor, pred_knn, mode = "prec_recall")

#c5.0
library(C50)
featuresc50 <- names(meanDecAcc[!a])
featuresc50

c50_model <- train(train_set[,featuresc50], train_set$Tumor, method="C5.0", trControl=tuneParams)
pred_c50 <- predict(c50_model, newdata = test_set[,featuresc50])
confusionMatrix(pred_c50, test_set$Tumor,mode = "prec_recall")


#Averaging
test_set$prob_knn <- predict(object = knn, test_set[, featureknn], type = 'prob')
test_set$prob_c50 <- predict(object = c50_model, test_set[,featuresc50], type = 'prob')
test_set$pred_avg <- (test_set$prob_knn+test_set$prob_c50)/2
test_set$pred <- apply(test_set$pred_avg, 1, FUN=function(x) {which.max(x)})
table(test_set$pred)
test_set$pred <- factor(test_set$pred, levels = c(1,2,3,4,5), labels = c('BRCA','KIRC','COAD','LUAD','PRAD'))
confusionMatrix(test_set$Tumor, test_set$pred,mode = "prec_recall")

################################################################################################################################################################################
#Individual model implementation
#training and testing data are created.
model_data <- dataset
split = sample.split(model_data$Tumor, SplitRatio = 0.8)
train_set = subset(model_data, split == TRUE)
test_set = subset(model_data, split == FALSE)

train_set[-1] = as.data.frame(scale(train_set[-1]))
test_set[-1] = as.data.frame(scale(test_set[-1]))
################################################################################################################################################################################
#C5.0
library(C50)
c50_model <- train(train_set[,-1], train_set$Tumor, method="C5.0", trControl=tuneParams)
pred_c50 <- predict(c50_model, newdata = test_set[,-1])
confusionMatrix(pred_c50, test_set$Tumor,mode = "prec_recall")
plot(c50_model)
################################################################################################################################################################################
#radial svm
svm_model <- train(train_set[,-1], train_set$Tumor, method="svmRadial", trControl=tuneParams)
pred_svm <- predict(svm_model, newdata = test_set[,-1])
confusionMatrix(pred_svm, test_set$Tumor,mode = "prec_recall")
plot(svm_model)
################################################################################################################################################################################
#KNN
knn_model <- train(train_set[,-1],train_set$Tumor,method='knn',trControl=tuneParams)
pred_knn <- predict(knn_model, newdata = test_set[,-1])
confusionMatrix(pred_knn, test_set$Tumor,mode = "prec_recall")
plot(knn_model)
################################################################################################################################################################################
#Stacking

train_set$knn_pred<-factor(knn_model$pred$pred[order(knn_model$pred$rowIndex)])
train_set$svm_pred<-factor(svm_model$pred$pred[order(svm_model$pred$rowIndex)])
train_set$c50_pred<-factor(c50_model$pred$pred[order(c50_model$pred$rowIndex)])

predictors <- c("knn_pred", "svm_pred", "c50_pred")
library(gbm)
gbm <- train(train_set[, predictors],train_set$Tumor,method='gbm',trControl=tuneParams)

test_set$knn_pred <- factor(pred_knn)
test_set$svm_pred <- factor(pred_svm)
test_set$c50_pred <- factor(pred_c50)

pred_gbm <- predict(gbm, newdata = test_set[,-1])
confusionMatrix(test_set$Tumor, pred_gbm,mode = "prec_recall")
summary(gbm)
