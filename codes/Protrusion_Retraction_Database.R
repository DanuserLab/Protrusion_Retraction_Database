##########################################################################################################################
###### 2022/06/27 Protrusion/Retraction Standard Pipeline
#### load packages
library(dplyr)       # Basic Dataframe mainipulation
library(depmixS4)    # For hmm estimation
library(ggplot2)     # For graphics   
library(readr)       # Read files
library(pheatmap)    # Heatmap production
library(gtools)      # Basic TS tools
library(zoo)         # Basic TS tools
#library(RColorBrewer)

#### Functions that compute the HMM modeling for single time series data
markovSingleFunction <- function(timeSeries, nstates) {
  out <- tryCatch(
    {
      hmm      <- depmix(timeSeries ~ 1, family = gaussian(), nstates = nstates, ntimes = length(timeSeries))
      hmmfit   <- fit(hmm, verbose = F)
      return(hmmfit)
    },
    error = function(cond) {
      return(NA)
    },
    finally = {
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      #message(paste("Processed URL:", url))
      #message("Some other message at the end")
    }
  )    
  return(out)
}

### General Input Variables
maxState  <- 8
smoothing <- T

# Input main working directory & subdirectory
# Create the subdirectory for output files and set the main directory to load moviedata files.
#save.image(file = "my_work_space.RData")
#### Set working directory
mainDir <- getwd()
subDir  <- "Result"

## Create directories and set to main directory
dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir))

####################################### Raw data input
## NA information
Velocity.df <- read.csv("GithubExampleData.txt", header = T, stringsAsFactors = F, strip.white = T)
head(Velocity.df)

########################################################################################################################### 
############ HMM Modeling starts
## Go back to main directory
setwd(file.path(mainDir))

## HMM purpose temp data frames
Ch0.hmm.df  <- filter(Velocity.df, Ch0.NAinfo == T)

## List objects to save temporary state data
Ch0.temp.fit  <- list()
tempChan0 <- Ch0.hmm.df$Ch0

##### Run HMM
t1 <- Sys.time()
set.seed(2)
for (i in 2:maxState) {
  tempfitmod             <- markovSingleFunction(tempChan0, i)
  Ch0.temp.fit[[i]]      <- list()
  Ch0.temp.fit[[i]][[1]] <- sample.int(i, length(tempChan0), replace = T)
  Ch0.temp.fit[[i]][[2]] <- tempfitmod
  if(!is.na(tempfitmod)){
    Ch0.temp.fit[[i]][[1]] <- tempfitmod@posterior$state
    print("good")
  }
}
Sys.time() - t1

########## Save data into data frame structures
## Save sorted state information to original data frame for Ch0
temp        <- as.data.frame(do.call(cbind, lapply(Ch0.temp.fit, function(x) return(x[[1]]))))
temp        <- rbind(temp, matrix(NA, nrow = (dim(Velocity.df)[1] - dim(temp)[1]), ncol = dim(temp)[2]))
names(temp) <- paste0("Ch0.S", 2:maxState)

temp.df <- arrange(Velocity.df, desc(Ch0.NAinfo))
temp.df <- cbind(temp.df, temp)
temp.df <- arrange(temp.df, Label, cellIndex, layerIndex, windowIndex, frameIndex)

###### Check
temp.state  <- dplyr::select(temp.df, starts_with("Ch0.S"))
temp.vector <- Velocity.df$Ch0
for (i in 2:maxState) {
  temp.mean <- rep(NaN, i)
  tempstate <- temp.state[, i - 1]
  for (j in 1:i) {
    temp.mean[j] <- mean(temp.vector[tempstate == j], na.rm = T)
  }
  change.vector <- rank(temp.mean)
  temp.state[, i - 1] <- change.vector[temp.state[, i - 1]]
}
temp.df <- dplyr::select(temp.df, -c(paste0("Ch0.S", 2:maxState)))
temp.df <- cbind(temp.df, temp.state)
State.df <- temp.df
### Write Data as text file
write.table(State.df, file=paste0("VelocityState.txt"), sep=",",
            row.names = F, col.names = T)



###################### Part 2 is to extract protrusion/retraction database
setwd(file.path(mainDir, subDir))

### Get all labels
tempLabel <- unique(Velocity.df$Label)
#############################################################################
for (k in 1:length(tempLabel)) {
  temp <- Velocity.df %>% filter(Label == tempLabel[k]) %>% group_by(cellIndex) %>% summarise(winNum   = max(windowIndex, na.rm = T),
                                                                                          frameNum = max(frameIndex, na.rm = T)) %>% ungroup()
  cellwiseWindowNumber <- temp$winNum
  cellWiseFrameNumber  <- temp$frameNum
  
  ###### Build protrusion event data frame
  ## Here we use the 8 state model for edge velocity
  temp.df <- dplyr::select(Velocity.df, Ch0, cellIndex, windowIndex, frameIndex, Ch0.S8, Label)
  temp.df <- filter(temp.df, Label == tempLabel[k])
  temp.df <- mutate(temp.df, binaryInd = 1 * (Ch0.S8 > 4))
  temp.df <- arrange(temp.df, cellIndex, windowIndex, frameIndex)
  
  temp <- temp.df$binaryInd
  tempLength <- length(temp)
  A <- rep(0, tempLength)
  for (i in 2:tempLength) {
    if(((temp[i] == 1) %in% T) & ((temp[i - 1] == 0) %in% T)){
      A[i] <- 1
    }else{
      if(((temp[i] == 1) %in% T) & is.na(temp[i - 1])){
        A[i] <- 1
      }else{
        if(((temp[i] == 0) %in% T) & ((temp[i - 1] == 1) %in% T)){
          A[i] <- 2
        }else{
          if(is.na(temp[i]) & ((temp[i - 1] == 1) %in% T)){
            A[i] <- 2
          }
        }
      }
    }
  }
  temp.df <- cbind(temp.df, A)
  protrusion.df <- filter(temp.df, A == 1)
  if(table(A)[2] > table(A)[3]){
    protrusion.df <- protrusion.df[-nrow(protrusion.df), ]
  }
  terminationTime <- dplyr::select(filter(temp.df, A == 2), frameIndex)
  pointCellIndex  <- dplyr::select(filter(temp.df, A == 2), cellIndex)
  names(terminationTime) <- "termination"
  
  terminationTime[terminationTime == 1] <- cellWiseFrameNumber[pointCellIndex[terminationTime == 1]]
  
  protrusion.df <- cbind(protrusion.df, terminationTime)
  protrusion.df <- dplyr::select(protrusion.df, -c("Ch0", "Ch0.S8", "A", "binaryInd"))
  
  ### Add variables in protrusion.df
  protrusion.df <- protrusion.df %>% mutate(Duration = termination - frameIndex + 1) %>%
    dplyr::select(-"Label")
  
  temp      <- filter(Velocity.df, Label == tempLabel[k])
  Before_protrusion.df <- mutate(protrusion.df, ## Protrusion speed related variables
                                         # average speed
                                         avgSpeed = apply(protrusion.df, 1, function(x){
                                           temp <- filter(temp, cellIndex == x[1] & windowIndex == x[2] &
                                                            frameIndex >= x[3] & frameIndex <= x[4])
                                           return(mean(temp$Ch0, na.rm = T))
                                         }),
                                         maxSpeed = apply(protrusion.df, 1, function(x){
                                           temp <- filter(temp, cellIndex == x[1] & windowIndex == x[2] &
                                                            frameIndex >= x[3] & frameIndex <= x[4])
                                           return(max(temp$Ch0, na.rm = T))
                                         })
  )
  #setwd(file.path(mainDir))
  
  ############## Filters can be manually added
  Before_protrusion.df <- filter(Before_protrusion.df, avgSpeed > 2 & Duration > 3)
  write.table(Before_protrusion.df, file = paste0(tempLabel[k], "_Protrusion.txt"), sep = ",",
              row.names = F, col.names = T)
  
  #############################################################################
  ###### Build retraction event data frame
  ## Here we use the 8 state model for edge velocity
  temp.df <- dplyr::select(Velocity.df, Ch0, cellIndex, windowIndex, frameIndex, Ch0.S8, Label)
  temp.df <- filter(temp.df, Label == tempLabel[k])
  temp.df <- mutate(temp.df, binaryInd = 1 * (Ch0.S8 < 5))
  temp.df <- arrange(temp.df, cellIndex, windowIndex, frameIndex)
  
  temp <- temp.df$binaryInd
  tempLength <- length(temp)
  A <- rep(0, tempLength)
  for (i in 2:tempLength) {
    if(((temp[i] == 1) %in% T) & ((temp[i - 1] == 0) %in% T)){
      A[i] <- 1
    }else{
      if(((temp[i] == 1) %in% T) & is.na(temp[i - 1])){
        A[i] <- 1
      }else{
        if(((temp[i] == 0) %in% T) & ((temp[i - 1] == 1) %in% T)){
          A[i] <- 2
        }else{
          if(is.na(temp[i]) & ((temp[i - 1] == 1) %in% T)){
            A[i] <- 2
          }
        }
      }
    }
  }
  temp.df <- cbind(temp.df, A)
  protrusion.df <- filter(temp.df, A == 1)
  if(table(A)[2] > table(A)[3]){
    protrusion.df <- protrusion.df[-nrow(protrusion.df), ]
  }
  terminationTime <- dplyr::select(filter(temp.df, A == 2), frameIndex)
  pointCellIndex  <- dplyr::select(filter(temp.df, A == 2), cellIndex)
  names(terminationTime) <- "termination"
  
  terminationTime[terminationTime == 1] <- cellWiseFrameNumber[pointCellIndex[terminationTime == 1]]
  
  protrusion.df <- cbind(protrusion.df, terminationTime)
  protrusion.df <- dplyr::select(protrusion.df, -c("Ch0", "Ch0.S8", "A", "binaryInd"))
  
  ### Add variables in protrusion.df
  protrusion.df <- protrusion.df %>% mutate(Duration = termination - frameIndex + 1) %>%
    dplyr::select(-"Label")
  
  temp      <- filter(Velocity.df, Label == tempLabel[k])
  PA_Rac1_Before_protrusion.df <- mutate(protrusion.df, ## Protrusion speed related variables
                                         # average speed
                                         avgSpeed = apply(protrusion.df, 1, function(x){
                                           temp <- filter(temp, cellIndex == x[1] & windowIndex == x[2] &
                                                            frameIndex >= x[3] & frameIndex <= x[4])
                                           return(mean(temp$Ch0, na.rm = T))
                                         }),
                                         maxSpeed = apply(protrusion.df, 1, function(x){
                                           temp <- filter(temp, cellIndex == x[1] & windowIndex == x[2] &
                                                            frameIndex >= x[3] & frameIndex <= x[4])
                                           return(max(temp$Ch0, na.rm = T))
                                         })
  )
  #setwd(file.path(mainDir))
  PA_Rac1_Before_protrusion.df <- filter(PA_Rac1_Before_protrusion.df, avgSpeed < -2 & Duration > 3)
  write.table(PA_Rac1_Before_protrusion.df, file = paste0(tempLabel[k], "_Retraction.txt"), sep = ",",
              row.names = F, col.names = T)
  
}



