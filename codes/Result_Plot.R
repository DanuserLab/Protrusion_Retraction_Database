##########################################################################################################################
###### 2022/06/28 PA Analysis
##########################################################################################################################
#### load R packages
#require(depmixS4)    # For hmm estimation
require(ggplot2)     # For graphics   
library(readr)
require(gtools)
library(RColorBrewer)
library(dplyr)       # Database Managing
library(zoo)
library(gridExtra)

########################################################################################################################### 
# Input main working directory & subdirectory
# This had to be setted where the movielist file is located.
# Create the subdirectory for output files and start in the main directory to load moviedata files.
#### Set working directory
mainDir <- "/project/bioinformatics/Danuser_lab/P01biosensor/analysis/Jaewon/Github/Protrusion_Retraction_Database"
subDir  <- "Result"
dataDirectory <- mainDir
dir.create(file.path(mainDir))
setwd(file.path(mainDir)) 
## Set the data directory
options(bitmapType = 'cairo')

########################################################################################################################### 
### Basic inputs for pipeline
AbsoluteVelocityThreshold <- 2  # Absolute threshold value in (nm/s) to filter protrusion/retraction data base based on its average speed.
AbsoluteDurationThreshold <- 3  # Absolute threshold value in duration to filter protrusion/retraction data base

#### Set working directory
mainDir <- getwd()
subDir  <- "Result"

#### Read photo activation indicator information
setwd(file.path(mainDir))
PA_Info <- read.csv("PA_Info.txt", header = T, stringsAsFactors = F, strip.white = T)
maxTime <- max(Velocity.df$frameIndex)
################################
setwd(file.path(mainDir, subDir))
### Get all labels
tempLabel <- unique(Velocity.df$Label)

ProtrusionSummaryData <- list()
RetractionSummaryData <- list()

#############################################################################
for (k in 1:length(tempLabel)){
  temp <- Velocity.df %>% filter(Label == tempLabel[k]) %>% group_by(cellIndex) %>% summarise(winNum   = max(windowIndex, na.rm = T),
                                                                                              frameNum = max(frameIndex, na.rm = T)) %>% ungroup()
  cellwiseWindowNumber <- temp$winNum
  cellWiseFrameNumber  <- temp$frameNum
  
  protrusion <- read.csv(paste0(tempLabel[k], "_Protrusion.txt"), header = T, stringsAsFactors = F, strip.white = T)
  retraction <- read.csv(paste0(tempLabel[k], "_Retraction.txt"), header = T, stringsAsFactors = F, strip.white = T)
  
  ### Additional filters can be applied for duration and speed
  protrusion <- protrusion %>% filter(avgSpeed > AbsoluteVelocityThreshold & Duration > AbsoluteDurationThreshold)
  retraction <- retraction %>% filter(avgSpeed < -AbsoluteVelocityThreshold & Duration > AbsoluteDurationThreshold)
  
  #############Full Combo protrusion
  temp_info <- PA_Info %>% filter(Label == tempLabel[k])
  protrusion <- mutate(protrusion, PA_window = apply(protrusion, 1, function(x){
    tt <- temp_info %>% filter(cellIndex == x[1] & windowIndex == x[2]) %>% dplyr::select(PA_window)
    return(tt$PA_window[1])
  }))
  
  #protrusion %>% group_by(cellIndex, PA_window) %>% summarise(count = n())
  
  dummy <- as.data.frame(matrix(NaN, length(cellwiseWindowNumber) * 2, dim(protrusion)[2]))
  names(dummy) <- names(protrusion)
  dummy$cellIndex <- rep(1:length(cellwiseWindowNumber), each = 2)
  dummy$PA_window <- rep(0:1, length(cellwiseWindowNumber))
  
  temp_data <- rbind(protrusion, dummy) %>% group_by(cellIndex, PA_window) %>% summarise(Count = n(),
                                                                                         Speed = mean(avgSpeed, na.rm = T),
                                                                                         Duration = 5 * mean(Duration, na.rm = T)) %>% ungroup()
  
  winNum   <- PA_Info %>% filter(Label == tempLabel[k]) %>% group_by(cellIndex, PA_window) %>%
    summarise(count = n()) %>% ungroup() %>% dplyr::select(count)
  frameNum <- rep(cellWiseFrameNumber/maxTime, each = 2)
  winNum   <- winNum/frameNum
  
  temp_data <- cbind(temp_data, winNum)
  ProtrusionSummaryData[[k]] <- temp_data %>% mutate(Frequency = Count/count) %>% 
    dplyr::select(cellIndex, PA_window, Speed, Duration, Frequency) %>% mutate(Label = tempLabel[k])
  
  #############Full Combo retraction
  temp_info <- PA_Info %>% filter(Label == tempLabel[k])
  retraction <- mutate(retraction, PA_window = apply(retraction, 1, function(x){
    tt <- temp_info %>% filter(cellIndex == x[1] & windowIndex == x[2]) %>% dplyr::select(PA_window)
    return(tt$PA_window[1])
  }))
  
  #protrusion %>% group_by(cellIndex, PA_window) %>% summarise(count = n())
  
  dummy <- as.data.frame(matrix(NaN, length(cellwiseWindowNumber) * 2, dim(protrusion)[2]))
  names(dummy) <- names(protrusion)
  dummy$cellIndex <- rep(1:length(cellwiseWindowNumber), each = 2)
  dummy$PA_window <- rep(0:1, length(cellwiseWindowNumber))
  
  temp_data <- rbind(retraction, dummy) %>% group_by(cellIndex, PA_window) %>% summarise(Count = n(),
                                                                                         Speed = mean(avgSpeed, na.rm = T),
                                                                                         Duration = 5 * mean(Duration, na.rm = T)) %>% ungroup()
  
  winNum   <- PA_Info %>% filter(Label == tempLabel[k]) %>% group_by(cellIndex, PA_window) %>%
    summarise(count = n()) %>% ungroup() %>% dplyr::select(count)
  frameNum <- rep(cellWiseFrameNumber/maxTime, each = 2)
  winNum   <- winNum/frameNum
  
  temp_data <- cbind(temp_data, winNum)
  RetractionSummaryData[[k]] <- temp_data %>% mutate(Frequency = Count/count) %>% 
    dplyr::select(cellIndex, PA_window, Speed, Duration, Frequency) %>% mutate(Label = tempLabel[k])
}

ProtrusionSummaryData <- do.call(rbind, ProtrusionSummaryData) %>% mutate(Timing = ifelse(grepl("Before", Label, fixed = TRUE), 0, 1),
                                                 Label = ifelse(grepl("Control", Label, fixed = TRUE), 'Control', "Mutant"))
RetractionSummaryData <- do.call(rbind, RetractionSummaryData) %>% mutate(Timing = ifelse(grepl("Before", Label, fixed = TRUE), 0, 1),
                                                                          Label = ifelse(grepl("Control", Label, fixed = TRUE), 'Control', "Mutant"))

write.table(ProtrusionSummaryData, file = paste0("ProtrusionSummary.txt"), sep = ",",
            row.names = F, col.names = T)
write.table(RetractionSummaryData, file = paste0("RetractionSummary.txt"), sep = ",",
            row.names = F, col.names = T)




write.table(temp_data, file=paste0("protrusion_table.csv"), sep=",",
            row.names = F, col.names = T)

#############Retraction Full Combo
#############Full Combo
temp <- Signal.df %>% filter(layerIndex == 1 & frameIndex == 3) %>% dplyr::select(cellIndex, windowIndex, PA_window)

retraction <- mutate(retraction, PA_window = apply(retraction, 1, function(x){
  tt <- temp %>% filter(cellIndex == x[1] & windowIndex == x[2]) %>% dplyr::select(PA_window)
  return(tt$PA_window)
}))

retraction %>% group_by(cellIndex, PA_window) %>% summarise(count = n())

dummy <- as.data.frame(matrix(NaN, cellNumber * 2, dim(retraction)[2]))
names(dummy) <- names(retraction)
dummy$cellIndex <- rep(1:cellNumber, each = 2)
dummy$PA_window <- rep(0:1, cellNumber)

temp_data <- rbind(retraction, dummy) %>% group_by(cellIndex, PA_window) %>% summarise(Count = n(),
                                                                                       Speed = mean(avgSpeed, na.rm = T),
                                                                                       Duration = mean(10 * protDuration, na.rm = T)) %>% ungroup()
winNum <- Signal.df %>% filter(layerIndex == 1 & frameIndex == 3) %>% group_by(cellIndex, PA_window) %>% summarise(count = n()) %>% ungroup() %>% dplyr::select(count)
frameNum <- rep(cellWiseFrameNumber/maxTime, each = 2)
winNum   <- winNum/frameNum

temp_data <- cbind(temp_data, winNum)
temp_data <- temp_data %>% mutate(Frequency = Count/count) %>% dplyr::select(cellIndex, PA_window, Speed, Duration, Frequency)
library(ggpubr)

png(file = paste0("R_Speed.png"), width = 1000, height = 800)
pp <- ggplot(data = temp_data, aes(x = factor(PA_window), y = Speed, group = factor(PA_window))) + geom_boxplot() +
  #geom_jitter(aes(x = layerIndex, y = odds_act_ratio, group = layerIndex, colour = factor(cellIndex)),
  #            size = 5, width = 0.1, height = 0) +
  theme(text = element_text(size = 40), legend.position='none') + stat_compare_means(paired = F, size = 10) +
  labs(x = "PA_window", y = "Speed(nm/s)")
print(pp)
dev.off()

png(file = paste0("R_Duration.png"), width = 1000, height = 800)
pp <- ggplot(data = temp_data, aes(x = factor(PA_window), y = 10 * Duration, group = factor(PA_window))) + geom_boxplot() +
  #geom_jitter(aes(x = layerIndex, y = odds_act_ratio, group = layerIndex, colour = factor(cellIndex)),
  #            size = 5, width = 0.1, height = 0) +
  theme(text = element_text(size = 40), legend.position='none') + stat_compare_means(paired = F, size = 10) +
  labs(x = "PA_window", y = "Time(s)")
print(pp)
dev.off()

png(file = paste0("R_Frequency.png"), width = 1000, height = 800)
pp <- ggplot(data = temp_data, aes(x = factor(PA_window), y = Frequency, group = factor(PA_window))) + geom_boxplot() +
  #geom_jitter(aes(x = layerIndex, y = odds_act_ratio, group = layerIndex, colour = factor(cellIndex)),
  #            size = 5, width = 0.1, height = 0) +
  theme(text = element_text(size = 40), legend.position='none') + stat_compare_means(paired = F, size = 10) +
  labs(x = "PA_window", y = "Frequency")
print(pp)
dev.off()

write.table(temp_data, file=paste0("retraction_table.csv"), sep=",",
            row.names = F, col.names = T)

