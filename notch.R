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
# convergence: convergence criteria if optimalMiss = TRUE
# maxIter: maximum number of interations if optimalMiss = TRUE
# graph: 0 if no graph, 1 if line graph, 2 if bar graph

notch <- function(data, notchPoint, binw, regFrom, regTo, excludeFrom, excludeTo, polySize, 
                  optimalMiss, graph, nboots){
  
  ans <- notchEst(data,notchPoint, binw, regFrom, regTo, excludeFrom, excludeTo, polySize, 
                  optimalMiss)
  regData <- ans$data
  
  #bootstrap
  excessBoot = rep(NA,nboots)
  excludeFromBoot = rep(NA,nboots)
  if (nboots>0){
    for (i in 1:nboots){
      bootData <- data[sample(1:length(data), length(data), replace=TRUE)]
      test <- notchEst(bootData,notchPoint, binw, regFrom, regTo, excludeFrom, excludeTo, polySize, optimalMiss)
      excessBoot[i] <- test$excess
      excludeFromBoot[i] <- test$excludeFrom
    }
  }
  excessSD <- sd(excessBoot)
  excludeFromSD <- sd(excludeFromBoot)
  
  if (graph == 1){
    graphics::plot(regData$mids,regData$counts,col="blue", 
                   xlab = "Abortion Day - 18 bday day", ylab = "counts",type="o")
    graphics::lines(regData$mids,regData$cf,col="red",type="l")
    graphics::abline(v=c(ans$excludeFrom,excludeTo),col="black",lty=2)
    #graphics::abline(v=c(ceiling(notchPoint)),col="green",lty=5,lwd=1)
    legend("topleft",legend=c(paste("ExcessMass (B) = ", format(ans$excess, digits=4), 
                                    " (", format(excessSD, digits=4), ")"), 
                              paste("NotchSize = ", format(ceiling(ans$excludeFrom), digits=4), 
                                    " (", format(excludeFromSD, digits=4), ")")),
           cex=0.7)
  }
  results <- list("excess" = ans$excess,
                  "excludeFrom" = ans$excludeFrom,
                  "totalBunch" = ans$totalBunch,
                  "totalDist" = ans$totalDist,
                  "data" = regData,
                  "excessBoot" = excessBoot,
                  "excessSD" = excessSD)
  
  return(results)
}