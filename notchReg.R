notchReg <- function(regData,polySize,notchPoint,excludeFrom,excludeTo){
  # create the "new" data frame for prediction
  # add the bin from regData and rename
  predData <- as.data.frame(regData$bin)
  colnames(predData)[colnames(predData)=="regData$bin"] <- "bin"
  
  # add counter-factual exclude values for prediction
  predData$exclude <- 0
  
  # counterfactual values using "new" data frame
  regData$cf <- predict(lm(counts~poly(bin,polySize)+exclude, data=regData), predData)
  
  # difference between actual and counterfactual
  regData$diff <- regData$counts - regData$cf
  # sum over the excess mass region
  excess <- sum(regData$diff[which(regData$mids> notchPoint & regData$mids < excludeTo)])
  # sum over the exclusion region to test whether the missing and excess mass balance
  total <-  sum(regData$diff[which(regData$exclude>0)])
  
  results <- list("excludeFrom" = excludeFrom, "excess" = excess, "total" = total, "cf" = regData$cf)
  return(results)
}
  