# Purpose of this function
# Test if there is bunching at a notch point, i.e. excess mass after and missing mass before the notch point
# The structure of the notch point is infomred by abortion legal forms in Spain in 2011
#
# function follows Kleven and Wassem (2013). Using notches to uncover optimization frictions ...
#
# input
# data: vector of the population of observations whose density need to be tested for bunching
# notchPoint: self-explanatory
# binw: bin width to collapse data
# regFrom - regTo: data points to be included for counterfactual prediction
# excludeFrom - excludeTo: the exclusion regions
# polySize: order of polynomial for regression
# optimalMiss: logical value - whether to optimize the boundary of the missing mass
# graph: 0 if no graph, 1 if line graph, 2 if bar graph

notchEst <- function(data, notchPoint, binw, regFrom, regTo, excludeFrom, excludeTo, polySize, optimalMiss){
  subData <- data[which(data>=regFrom & data <= regTo)]
  
  # data frequency
  dataHist <- graphics::hist(subData, breaks = seq(from=regFrom,to=regTo,by=binw),plot=FALSE)
  
  # extract from dataHist the bins and the mid-points
  regData <- as.data.frame(dataHist$counts)
  regData$mids <- dataHist$mids
  
  # rename the "counts" column for ease of reference
  colnames(regData)[colnames(regData)=="dataHist$counts"] <- "counts"
  
  # add the bin names sequentially for regression
  regData$bin <- seq(from=1,to=length(dataHist$counts),by=1)
  
  # add the exclude column
  # if optimalMiss = true, create a list of possible excludeFrom points
  if (optimalMiss == TRUE){
    search = seq(from=excludeFrom,to=notchPoint,by=binw)
  } else {search = excludeFrom}
  
  excess = rep(NA,length(search))
  total = rep(NA,length(search))
  
  # initialize the optimization for the boundary of the missing mass
  for (i in 1:length(search)){
    # create the indicators for exclude variables
    regData$exclude <- (regData$mids - search[i])/binw + 0.5 # alternatively, just regData$mids - excludeFrom is sufficient; the algebra is purely cosmetic
    regData$exclude[regData$mids < search[i]] <- 0 #if the mid-points < excludeFrom then exclude = 0
    regData$exclude[regData$mids > excludeTo] <- 0 #if the mid-points > excludeTo then exclude = 0
    
    test <- notchReg(regData,polySize,notchPoint,search[i],excludeTo)
    total[i] <- test$total
    excess[i] <- test$excess
  }
  
  #take the absolute value of total to find the minimum deviation from 0
  totalAbs <- abs(total)
  opt <- which(totalAbs == min(totalAbs))
  excludeFromOpt <- search[opt]
  excessOpt <- excess[opt]
  totalOpt <- total[opt]
  regData$cf <- notchReg(regData,polySize,notchPoint,excludeFromOpt,excludeTo)$cf
  
  results <- list("excess" = excessOpt,
                  "excludeFrom" = excludeFromOpt,
                  "totalBunch" = totalOpt,
                  "totalDist" = total,
                  "data" = regData)
  
  return(results)
}