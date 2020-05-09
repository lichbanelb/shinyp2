#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# if (!requireNamespace("BiocManager", quietly=TRUE))
#    install.packages("BiocManager")
# BiocManager::install("GEOquery")
# library("BiocManager")

library(BiocManager)
options(repos = BiocManager::repositories())

library(devtools)
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery")
library("GEOquery")  ## go to https://github.com/seandavi/GEOquery for installation details
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma)
library(dplyr)
library(caret)
  # Please install caret packages for all classifiers: knn, random forests, naive bayes (klaR), support vector machine and multilayer perceptron (RSNNS)
  # Caret will automatically prompt installation, when the app is run, within the console

clinical_outcome = getGEO("GSE120396")
clinical_outcome = clinical_outcome$GSE120396_series_matrix.txt.gz

rejection_status = clinical_outcome$characteristics_ch1.1
rejection_status = unlist(lapply( strsplit(as.character(rejection_status), ": " ) , `[[` , 2)  )
table(rejection_status)

# Note: please change this dir to point to the folder where your dataset is
datadir = "data/"
fileNames = list.files(datadir)

# for (files in fileNames){  # Run once, then comment out
#   gunzip(file.path(datadir,files))
# }

fileNames = list.files(datadir)
fileNames = fileNames[substr(fileNames, nchar(fileNames)-3, nchar(fileNames)) == ".txt"]

gse = c()
for(i in 1:length(fileNames)){
  temptable = read.delim(file.path(datadir, fileNames[i]), header=TRUE)
  gse = cbind(gse, temptable[,2])
  colnames(gse)[i] = colnames(temptable)[2]
}

rownames(gse) = read.delim(file.path(datadir, fileNames[1]), header=TRUE)[,1]

gse_pca = prcomp(t(gse), scale = TRUE)
pData = gse_pca$x

# Define server 
shinyServer(function(input, output) {
  
  setModel = "knn"  # Default selection is KNN
  setModel = eventReactive(input$modelConfirm, {
    if (input$classifier == "KNN"){
      model = "knn"
    } else if (input$classifier == "Random Forests") {
      model = "rf"
    } else if (input$classifier == "Naive Bayes"){
      model = "naive_bayes"
    } else if (input$classifier == "Multilayer Perceptron"){
      model = "mlp"
    } else if (input$classifier == "Linear Support Vector Machine"){
      model = "svmLinear"
    } else if (input$classifier == "Decision Tree"){
      model = "rpart"
    } 

    return(model)
  })
  output$result <- renderText({
    paste("\nCurrently selected model:", input$classifier)
  })
  output$accPlot <- renderPlot({
    
    nComps = input$nComps
    
    X = as.matrix(pData[,1:nComps])
    y = rejection_status
    
    mlData = data.frame(X)
    mlData$y = y
    
    k1 = 1
    k2 = input$k2
    
    nFolds = input$nFolds
    repeats = input$repeats
    
    train_method = trainControl(method = "repeatedcv", number = nFolds, repeats = repeats)
    
    chosenModel = setModel()
    if (chosenModel == "knn") {
      modelObj <<- caret::train(factor(y) ~ ., data = mlData, method = chosenModel, trControl = train_method, tuneGrid = expand.grid(k = k1:k2))
      cm = modelObj$resampledCM
      cm$acc = (cm$cell1 + cm$cell4)/(cm$cell1 + cm$cell2 + cm$cell3 + cm$cell4)
      cm %>% ggplot(aes(x = k, y = acc, group = k)) + geom_boxplot() + theme_classic() + 
        theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
              plot.title = element_text(size = 20), plot.caption = element_text(size = 16)) + 
        labs(title = "Accuracy for k Values for KNN in Cross Validation", y = "Accuracy", x = "k Value",
             caption = paste("\nAccuracy results by k value for KNN performed on selected principal components.
              Objective is classification of patient rejection status using gene expression data
             
             ", "Estimated Accuracy:", round(max(modelObj$results$Accuracy), 2)))
    } else if (chosenModel == "rf") {
      modelObj <<- caret::train(factor(y) ~ ., data = mlData, method = chosenModel, trControl = train_method)
      cm = modelObj$resampledCM
      cm$acc = (cm$cell1 + cm$cell4)/(cm$cell1 + cm$cell2 + cm$cell3 + cm$cell4)
      cm %>% ggplot(aes(x = mtry, y = acc, group = mtry)) + geom_boxplot() + theme_classic() + 
        theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
              plot.title = element_text(size = 20), plot.caption = element_text(size = 16)) +
        labs(title = "Accuracy for mtry Values\nfor Random Forests in Cross Validation", y = "Accuracy", x = "mtry",
             caption = paste("\nAccuracy results by mtry value for Random Forests performed on selected principal components.
              Objective is classification of patient rejection status using gene expression data
                             
              ", "Estimated Accuracy:", round(max(modelObj$results$Accuracy), 2)))
    } else {
      modelObj <<- caret::train(factor(y) ~ ., data = mlData, method = chosenModel, trControl = train_method)
      cm = modelObj$resampledCM
      cm$acc = (cm$cell1 + cm$cell4)/(cm$cell1 + cm$cell2 + cm$cell3 + cm$cell4)
      cm %>% ggplot(aes(y = acc)) + geom_boxplot() + theme_classic() + 
        theme(axis.text = element_text(size = 15), axis.title = element_text(size = 18), 
              plot.title = element_text(size = 20), plot.caption = element_text(size = 16)) +
        labs(title = "Accuracy Across Cross Validation", y = "Accuracy", x = "",
             caption = paste("\nAccuracy results on a classifier applied to selected principal components.
              Objective is classification of patient rejection status using gene expression data
             
             ", "Estimated Accuracy:", round(max(modelObj$results$Accuracy), 2)))

    }
    #knn
    #output$accuracy = 
    #boxplot(knn$resample$Accuracy, main = "Boxplot of KNN Accuracies", ylab = "Accuracy")
  })
  # output$accuracy <- renderText({
  #   paste("\nEstimated accuracy:",round(max(modelObj$results$Accuracy) ,2))
  # })
})
