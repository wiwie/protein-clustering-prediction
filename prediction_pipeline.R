### Clustering Prediction Pipeline
## Script Version & Date
# 1.0, 2016/09/22

## Script Usage
# This script is intended to be applied to a ClustEval data analysis result of a new data set to predict expected 
# clustering performances of clustering methods when applied to the new data set.
# This script does not require any modifications by you.
# It is free of any size limits; i.e. the limit of 50,000 sequences mentioned on the website and the publication
# is website-specific.

## Requirements
# (1) The script expects the new6.RData file located in the same directory. 
#     This R data file contains the prediction models and should be distributed together with this script.
# (2) This script requires the following R packages:
#     install.packages(c("glmnet","plyr","ggplot2","reshape","grid","pdist"))
# (3) The script requires a ClustEval data analysis result folder with the calculated data statistics of a new data set.
#     The following statistics are expected to be contained in the folder and are necessary to perform predictions:
#       GraphMinCutRDataStatistic
#       NumberOfSamplesDataStatistic
#       GraphDiversityAverageRDataStatistic
#       ClusteringCoefficientRDataStatistic
#       AssortativityDataStatistic
#       MatrixRankRDataStatistic
#       GraphDensityRDataStatistic
#       AssortativityWeightedDataStatistic
#       GraphAdhesionRDataStatistic
#       SimilarityPercentile10DataStatistic
#       SimilarityPercentile20DataStatistic
#       SimilarityPercentile30DataStatistic
#       SimilarityPercentile40DataStatistic
#       SimilarityPercentile50DataStatistic
#       SimilarityPercentile60DataStatistic
#       SimilarityPercentile70DataStatistic
#       SimilarityPercentile80DataStatistic
#       SimilarityPercentile90DataStatistic
#       SimilarityPercentile100DataStatistic
#       AbsoluteZScore01DataStatistic
#       AbsoluteZScore12DataStatistic
#       AbsoluteZScore23DataStatistic
#       AbsoluteZScore3InfDataStatistic
#       MinSimilarityDataStatistic
#       MaxSimilarityDataStatistic
#       MinimalNodeSimilaritySumDataStatistic

## Script Parameters:
# This R script is intended to be invoked from command line with the following parameters:
# (1) Absolute path to the ClustEval result folder that contains all the data analysis statistic results of a new data set.
#     Example: Rscript prediction_pipeline.R /home/username/repository/results/08_13_2016-01_00_34_analysis_new_dataset/

## Script Output:
# (1) The script will write the predicted F1-Scores into a file called "predicted_qualities.txt" in the passed ClustEval result folder.
#     Example: /home/username/repository/results/08_13_2016-01_00_34_analysis_new_dataset/predicted_qualities.txt
# (2) Predicted best parameters will be written to a file called "predicted_parameters.txt" in the passed ClustEval result folder.
#     Example: /home/username/repository/results/08_13_2016-01_00_34_analysis_new_dataset/predicted_parameters.txt
# (3) A plot of the parsed data statistics of the new data set will be stored in a file called "data_statistics.png" in the passed ClustEval result folder.
#     Example: /home/username/repository/results/08_13_2016-01_00_34_analysis_new_dataset/data_statistics.png
# (4) A plot of the data statistics of the most similar data set to the new data set in terms of data statistics will be stored in a file "data_statistics_most_similar.png".
#     Example: /home/username/repository/results/08_13_2016-01_00_34_analysis_new_dataset/data_statistics_most_similar.png

## Authors:
# Christian Wiwie, wiwiec@imada.sdu.dk
# Richard Röttger, roettger@imada.sdu.dk

### SCRIPT BEGIN
library("glmnet")
load("new6.RData")

args <- commandArgs(trailingOnly = TRUE)
resFolder <- args[1]

toolList <-list()
modelTypeList <- list()
modelList <- list()
for (finalModel in cvResultsAllFam3$finalModels) {
  toolList[[length(toolList)+1]] <- finalModel$Program.Config
  modelTypeList[[length(modelTypeList)+1]] <- finalModel$modelType
  modelList[[length(modelList)+1]] <- finalModel$model
}

predictedF1Scores <- parseDataStatisticsAndPredict(
  toolList, modelList, modelTypeList, 
  resFolder)


write.table(predictedF1Scores$f1, file=paste(resFolder,"predicted_qualities.txt",sep="/"), quote = F, row.names = F, sep="\t")



# get best parameters of most similar data set
library(plyr)
statsWithBestParams <- stats[as.character(join(data.frame(Data.Config=c(dataConfigIdsCFamily, dataConfigIdsFoldFamily)),dataConfigToDataSet)$Data.Set),]
statsWithBestParams <- statsWithBestParams[apply(statsWithBestParams, MARGIN=1, function(x) {return (!any(is.na(x)))}),]
predictedParams <- getBestParametersOfSimilarDataset(statsWithBestParams, 
                                                     predictedF1Scores$dataStatistics, 
                                                     bestParameters2.df.avg, resFolder)
write.table(predictedParams, file=paste(resFolder,"predicted_parameters.txt",sep="/"), quote = F, row.names = F, sep="\t")



### FUNCTIONS
parseDataStatisticsAndPredict <- function(toolList, modelList, modelTypeList, dataAnalysisRunResultFolder) {
  dataStatistics <- parseDataStatisticsFromFolder(dataAnalysisRunResultFolder)
  plotDataStatistics(dataAnalysisRunResultFolder, dataStatistics)
  f1 <- predictedToolPerformances(toolList, modelList, modelTypeList, dataStatistics)
  list(dataStatistics=dataStatistics, f1=f1)
}

plotDataStatistics <- function(dataAnalysisRunResultFolder, dataStatistics, isNewDataStatistic=T) {
  library(ggplot2)
  library(reshape)
  stats.melt <- melt(dataStatistics)
  colnames(stats.melt) <- c("Statistic","Value")
  normalizedStats <- gsub("DataStatistic","",c("GraphMinCutRDataStatistic",
                       #"NumberOfSamplesDataStatistic",
                       "GraphDiversityAverageRDataStatistic",
                       "ClusteringCoefficientRDataStatistic",
                       "AssortativityDataStatistic",
                       "MatrixRankRDataStatistic","GraphDensityRDataStatistic",
                       "AssortativityWeightedDataStatistic","GraphAdhesionRDataStatistic",
                       "SimilarityPercentile10DataStatistic","SimilarityPercentile20DataStatistic",
                       "SimilarityPercentile30DataStatistic","SimilarityPercentile40DataStatistic",
                       "SimilarityPercentile50DataStatistic","SimilarityPercentile60DataStatistic",
                       "SimilarityPercentile70DataStatistic","SimilarityPercentile80DataStatistic",
                       "SimilarityPercentile90DataStatistic","SimilarityPercentile100DataStatistic",
                       "AbsoluteZScore01DataStatistic","AbsoluteZScore12DataStatistic",
                       "AbsoluteZScore23DataStatistic",
                       "AbsoluteZScore3InfDataStatistic"#, 
                       #"MinSimilarityDataStatistic",
                       #"MaxSimilarityDataStatistic",
                       #"MinimalNodeSimilaritySumDataStatistic"
  ))
  stats.melt$Normalized[stats.melt$Statistic %in% normalizedStats] <- "Normalized"
  stats.melt$Normalized[!stats.melt$Statistic %in% normalizedStats] <- "Not Normalized"
  
  p <- ggplot(stats.melt, aes(x=Statistic, y = Value, fill=Normalized))
  p <- p + xlab("Data Statistic")
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (isNewDataStatistic)
    p <- p + labs(title = "Data Statistics of Your Dataset")
  else
    p <- p + labs(title = "Data Statistics of the Dataset Most Similar to Yours")
  p <- p + theme(#axis.text=element_text(size=12),
    #axis.title=element_text(size=14),
    plot.margin = unit(c(1, 1, 1, 2), "cm"),
    strip.text.x = element_text(size = 12))
  p <- p + facet_wrap(~ Normalized, scales="free")
  p <- p + geom_boxplot()
  
  # Get the ggplot grob
  gt = ggplotGrob(p)
  # Check for the widths - you need to change the two that are set to 1null
  gt$widths
  # The required widths are 4 and 7
  # Replace the default widths with relative widths:
  library(grid)
  gt$widths[4] = unit(22/25, "null")
  gt$widths[7] = unit(3/25, "null")
  # I think it is better to have some extra space between the two panels
  gt$widths[5] = unit(1, "cm")
  
  if (isNewDataStatistic)
    png(paste(dataAnalysisRunResultFolder, "data_statistics.png", sep="/"), width=14, height=6, units = "in", res = 300)
  else
    png(paste(dataAnalysisRunResultFolder, "data_statistics_most_similar.png", sep="/"), width=14, height=6, units = "in", res = 300)
  #plot(p)
  grid.newpage()
  grid.draw(gt)
  dev.off()
}

predictedToolPerformances <- function(toolList, modelList, modelTypeList, newDataStatisticsDataFrame) {
  f1 <- data.frame(Program.Config=character(), modelType=character(), f1=numeric(), truncated.f1=numeric())
  for (i in 1:length(modelList)) {
    model <- modelList[[i]]
    modelType <- modelTypeList[[i]]
    tool <- toolList[[i]]
    if (modelType %in% c("ridge","lasso")) {
      predictedF1 <- predict(model, as.matrix(newDataStatisticsDataFrame))[1,1]
    } else if (modelType == "ordinary") {
      predictedF1 <- predict(model, newDataStatisticsDataFrame)
    }
    f1Trunc <- max(0.0, min(predictedF1, 1.0))
    f1 <- rbind(f1,
                data.frame(Program.Config=tool, modelType=modelType, f1=predictedF1, truncated.f1=f1Trunc))
  }
  f1
}

parseDataStatisticsFromFolder <- function(folder) {
  statisticsNames <- c("GraphMinCutRDataStatistic",
                       #"ClusteringCoefficientDataStatistic",
                       #"IntraInterOverlapDataStatistic",
                       #"RichClubMaxDataStatistic",
                       "NumberOfSamplesDataStatistic",
                       "GraphDiversityAverageRDataStatistic",
                       #"GraphCohesionRDataStatistic",
                       "ClusteringCoefficientRDataStatistic",
                       "AssortativityDataStatistic",
                       #"HopkinsDataStatistic",
                       "MatrixRankRDataStatistic","GraphDensityRDataStatistic",
                       "AssortativityWeightedDataStatistic","GraphAdhesionRDataStatistic",
                       #"RichClub10DataStatistic",
                       "SimilarityPercentile10DataStatistic","SimilarityPercentile20DataStatistic",
                       "SimilarityPercentile30DataStatistic","SimilarityPercentile40DataStatistic",
                       "SimilarityPercentile50DataStatistic","SimilarityPercentile60DataStatistic",
                       "SimilarityPercentile70DataStatistic","SimilarityPercentile80DataStatistic",
                       "SimilarityPercentile90DataStatistic","SimilarityPercentile100DataStatistic",
                       "AbsoluteZScore01DataStatistic","AbsoluteZScore12DataStatistic",
                       "AbsoluteZScore23DataStatistic",
                       "AbsoluteZScore3InfDataStatistic", "MinSimilarityDataStatistic",
                       "MaxSimilarityDataStatistic",
                       # we use this one to normalize the graph min-cut
                       "MinimalNodeSimilaritySumDataStatistic"
                       )
  
  dataConfigIds <- list.files(paste(folder, "configs", sep="/"), pattern=".dataconfig")
  # only take first dataconfig; we assume only one dataset here
  dataConfigIds <- c(gsub(".dataconfig","",dataConfigIds[1],fixed=T))
  
  dataStatistics <- data.frame(statistic=character(), value=numeric())
  for (absoluteResultFolder in paste(folder,"analyses",sep="/")) {
    for (statisticsName in statisticsNames) {
      for (dataConfig in dataConfigIds) {
        potentialStatFile <- paste(paste(dataConfig, statisticsName,sep="_"), ".txt",sep="")
        absPotentialStatFile <- paste(absoluteResultFolder,potentialStatFile,sep="/")
        if (file.exists(absPotentialStatFile)) {
          if (length(readLines(absPotentialStatFile)) > 0) {
            statFileX <- as.numeric(readLines(absPotentialStatFile))
            dataStatistics <- rbind(dataStatistics, data.frame(statistic  = statisticsName, value=statFileX))
          }
        }
      }
    }
  }
  
  statisticVector <- rbind(dataStatistics$value)
  colnames(statisticVector) <- dataStatistics$statistic
#  print("Original data statistic values:")
#  print(statisticVector)
  statisticVector[is.na(statisticVector)] <- 0
#  print("Data statistic values after filtering:")
#  print(statisticVector)
  origStats <- as.data.frame(statisticVector)
  
  colnames(origStats) <- gsub("value.","",colnames(origStats))
  colnames(origStats) <- gsub("DataStatistic","",colnames(origStats))
  
  # normalize statistics
  if ("Assortativity" %in% colnames(origStats))
    origStats[,"Assortativity"] <- (origStats[,"Assortativity"]+1)/2
  if ("AssortativityWeighted" %in% colnames(origStats))
    origStats[,"AssortativityWeighted"] <- (origStats[,"AssortativityWeighted"]+1)/2
  if ("MatrixRankR" %in% colnames(origStats))
    origStats[,"MatrixRankR"] <- origStats[,"MatrixRankR"]/origStats[,"NumberOfSamples"]
  if ("GraphAdhesionR" %in% colnames(origStats))
    origStats[,"GraphAdhesionR"] <- origStats[,"GraphAdhesionR"]/(origStats[,"NumberOfSamples"]*origStats[,"NumberOfSamples"])
  if ("GraphMinCutR" %in% colnames(origStats)) {
    origStats$GraphMinCutR <- origStats$GraphMinCutR/origStats$MinimalNodeSimilaritySum
    origStats$GraphMinCutR[is.nan(origStats$GraphMinCutR)] <- 0
  }
  
  stats <- origStats
  
  # remove statistics
  if ("GraphCohesionR" %in% colnames(origStats))
    stats <- stats[,-c(which(colnames(stats)=="GraphCohesionR"))]
  if ("MinimalNodeSimilaritySum" %in% colnames(origStats))
    stats <- stats[,-c(which(colnames(stats)=="MinimalNodeSimilaritySum"))]

  stats
}

getBestParametersOfSimilarDataset <- function(stats, newStatistics, bestParameters2.avg, dataAnalysisRunResultFolder) {
  library(pdist)
  featuresForDistCalc <- colnames(newStatistics)
  featuresForDistCalc <- featuresForDistCalc[!featuresForDistCalc %in% c("NumberOfSamples",
                                  "MaxSimilarity",
                                  "MinSimilarity")]
  d <- as.matrix(pdist(as.matrix(stats[, featuresForDistCalc]), as.matrix(newStatistics[,featuresForDistCalc])))
  minDistInd <- which.min(d)
  minDistDataConfig <- as.character(stats$Data.Config[minDistInd])
  plotDataStatistics(dataAnalysisRunResultFolder, stats[stats$Data.Config == minDistDataConfig,-1], isNewDataStatistic = F)
  print(sprintf("Most similar data configuration found is: %s", minDistDataConfig))
  params <- subset(bestParameters2.df.avg, Data.Config==minDistDataConfig)
  params <- params[!is.na(params$valueNum),c("Data.Config", "Program.Config", "parameter", "value")]
  params
}