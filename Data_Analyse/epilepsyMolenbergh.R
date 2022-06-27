########################################################################################################################
#########################################################  GLM #########################################################
########################################################################################################################
rm(list=ls())
library(RdeGlmLassoGam)
source("RobGLM/glmrob.R")
source("RobGLM/glmrobMqle.R")

library(sas7bdat)
library(dglm)
library(RdeGlmLassoGam)
library(ggplot2)
#######################
# First analysis
#######################
RawData <- read.sas7bdat(file =  "Data_Analyse/epilepsy.sas7bdat")
RawData$id <- as.factor(RawData$id)
RawData$sex <- as.factor(RawData$sex)
RawData$race <- as.factor(RawData$race)
RawData$trt <- as.factor(RawData$trt)
RawData$height <- (RawData$height - median(RawData$height))/mad(RawData$height)
RawData$weight <- (RawData$weight - median(RawData$weight))/mad(RawData$weight)

Data <- data.frame(cbind(RawData[1:2,], NA))
colnames(Data)[ncol(Data)] <- "totseiz"
PatientID <- levels(RawData$id)
count <- 1
for(i in 1:length(PatientID)) {
  tempInd <- which(RawData$id==PatientID[i])
  if (length(tempInd) < 12) { 
    print(PatientID[i])
  } else {
    Data[count, 1:ncol(RawData)] <- RawData[tempInd[1],]
    Data$totseiz[count] <- sum(RawData$nseizw[tempInd[9:12]])
    count <- count + 1 
  }
}
Data <- Data[,-c(9:11)]

response <- Data$trt

formula <- totseiz ~ sex + race + height +  weight + bserate + trt
dformula <- totseiz ~ trt

hist(Data$height)
hist(Data$weight)
hist(Data$bserate)

######################
# CEQL
######################
dglm.est <- dglm(data = Data,
                 formula = formula,
                 dformula = dformula,
                 family = "poisson",
                 method = "reml")
summary.dglm(dglm.est)

round(dglm.est$coefficients, digits = 2)
round(dglm.est$dispersion.fit$coefficients, digits = 2)   #Report with a minus sign as model specification is different from our RDE specification

#######################
# CDE
#######################
# Format inputs
glmCant <- glmrob(formula, family = poisson,  data = Data, control = glmrobMqle.control(maxit = 500))

designX <- model.matrix(glmCant)
designZ <- model.matrix(dglm.est$dispersion.fit)

optionList1 <- list(huberC = 20000,
                    tukeyC = 20000,
                    tol = 1e-4,
                    maxIt = 15)

# Actual fitting
rob.est <- RDE(y = dglm.est$y,
                                X=designX, Z=designZ,
                                weights.on.xz1 = "none", gamma.ini = c(-2, 0.001),
                                weightFunction1 = "Huber",
                                family = "poisson",
                                optionList = optionList1)
mu <- as.numeric(exp(designX %*% rob.est$betaEstimate))
theta <- as.numeric(exp(designZ %*% rob.est$gammaEstimate))
plot(1 / theta)

# Format results
alpha <- c(rob.est$betaEstimate, rob.est$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * rob.est$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

#######################
# HRDE
#######################
# Format inputs
optionList2 <- list(huberC = 1.85,
                   tukeyC = 4.75,
                   tol = 1e-4,
                   maxIt = 50,
                   alpha = 0.75)

# Actual fitting
rob.est <- RDE(y = Data$totseiz,
                                X=designX, Z=designZ,
                                weights.on.xz1 = "covMcd",
                                contX = 4:6,
                                weightFunction1 =  "Huber",
                                family = "poisson",
                                optionList = optionList2)

mu <- as.numeric(exp(designX %*% rob.est$betaEstimate))
theta <- as.numeric(exp(designZ %*% rob.est$gammaEstimate))
plot(1 / theta)

# Format results
alpha <- c(rob.est$betaEstimate, rob.est$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * rob.est$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

#######################
# TRDE
#######################
rob.est <- RDE(y = Data$totseiz,
                                X=designX, Z=designZ,
                                weights.on.xz1 = "covMcd",
                                contX = 4:6, 
                                weightFunction1 =  "Tukey",
                                family = "poisson",
                                optionList = optionList2)

mu <- as.numeric(exp(designX %*% rob.est$betaEstimate))
theta <- as.numeric(exp(designZ %*% rob.est$gammaEstimate))
plot(1 / theta)

# Format results
alpha <- c(rob.est$betaEstimate, rob.est$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * rob.est$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

#######################
# Robust Testing
#######################
# Format inputs
designX <- designX[,c(1,6,2,3,4,5,7)]
designZ <- designZ

# Actual fitting
rob.est_Final <- RDE(y = Data$totseiz,
                                X=designX, Z=designZ,
                                weights.on.xz1 = "covMcd",
                                contX = c(2,5,6), 
                                weightFunction1 = "Tukey",
                                family = "poisson",
                                optionList = optionList2)
rob.est_Final_NULL <- RDE(y = Data$totseiz ,X=designX[,1:2], Z=designZ[,1,drop=FALSE],
                                     weights.on.xz1 = "covMcd",
                                     contX = c(2), 
                                     weightFunction1 = "Tukey",
                                     family = "poisson",
                                     optionList = optionList2)

# Actual testing
DispTest <- DispersionTest(y=Data$totseiz,X=designX, Z=designZ,
                           m = NULL,
                           alphaHAlt = c(rob.est_Final$betaEstimate, rob.est_Final$gammaEstimate),
                           alphaHNull = c(rob.est_Final_NULL$betaEstimate,0,0,0,0,0, rob.est_Final_NULL$gammaEstimate,0),
                           nzeroBeta = 5, #q1
                           nzeroGamma = 1, #q2
                           family="poisson", 
                           weightFunction1 = "Tukey",
                           weights.on.xz1 = rob.est_Final$weights.xz1,
                           weights.on.xz2 = rob.est_Final$weights.xz2,
                           optionList = optionList2,
                           M = rob.est_Final$M, Q = rob.est_Final$Q)
DispTest

# Generate residual plot
PearsonResid <- (Data$totseiz - mu) / sqrt(mu / (theta) )
PlotData <- data.frame(1:nrow(Data), PearsonResid, response )
colnames(PlotData) <- c("index", "residual", "Treatment")
Plot <- ggplot(PlotData) + geom_point(aes(x = index, y = residual, color = Treatment)) 
Plot <- Plot + scale_colour_manual(values = c("blue", "orange"), labels = c("0", "1"))
Plot <- Plot + xlab("Index of observation") + ylab("Pearson residual") +
  ggtitle("Residual plot for robust RDE estimator") + mrfDepth:::mrfDepth_theme() 
Plot
ggsave(plot = Plot, filename = "EpilepsyMol.eps", width = 4 * 1.5, height =  3 * 1.5)

designX <- model.matrix(glmCant)
designZ <- model.matrix(dglm.est$dispersion.fit)

badLev = (rownames(Data[order(PearsonResid)[1:2],])) # Select the patient codes resulting in the two lowest residuals.
Data[order(PearsonResid)[1:2],2] #id
Data[order(PearsonResid)[1:2],8] #bserate
sort(Data[,8]) #sorted bserate
rob.est_Final$weights.xz1[which(rownames(Data)%in%badLev)] # These were indeed indicated with weight 0.

########################################################################################################################
################################################# ADAPTIVE LASSO GLM ###################################################
########################################################################################################################
# lambdaBetaValues = seq(5,10,0.05)
# lambdaGammaValues = seq(1,5,0.05)
# EBICGamMat = matrix(NA, nrow = length(lambdaBetaValues), ncol = length(lambdaGammaValues))
# beta.iniList = NULL
# gamma.iniList = NULL
# beta.iniListStart = NULL
# gamma.iniListStart = NULL
# for(i in rev(1:length(lambdaBetaValues))){
#   lb = lambdaBetaValues[i]
#   for(j in rev(1:length(lambdaGammaValues))){
#     lg = lambdaGammaValues[j]
#     if(j==length(lambdaGammaValues)){
#       res = EBICGam(y= Data$totseiz, X=designX, Z=designZ, family="poisson",
#                     muEstim = 'pen', thetaEstim = 'pen',
#                     beta.ini = beta.iniListStart, gamma.ini = gamma.iniListStart,
#                     lambdaBeta=lb, lambdaGamma=lg,
#                     weights.on.xz1 = "covMcd", contX = 4:6,
#                     weightFunction1 = "Huber",
#                     optionList = optionList2)
#     }else{
#       res = EBICGam(y= Data$totseiz, X=designX, Z=designZ, family="poisson",
#                     muEstim = 'pen', thetaEstim = 'pen',
#                     beta.ini = beta.iniList, gamma.ini = gamma.iniList,
#                     lambdaBeta=lb, lambdaGamma=lg,
#                     weights.on.xz1 = "covMcd", contX = 4:6,
#                     weightFunction1 = "Huber",
#                     optionList = optionList2)
#     }
# 
#     EBICGamMat[i,j] = res$EBICGam
#     beta.iniList = res$fit$betaEstimate
#     gamma.iniList = res$fit$gammaEstimate
# 
#     if(j ==length(lambdaGammaValues)){
#       beta.iniListStart = beta.iniList
#       gamma.iniListStart = gamma.iniList
#     }
#   }
# }
# 
# EBICBetMat = rep(NA, length(lambdaBetaValues))
# 
# beta.iniList = NULL
# gamma.iniList = NULL
# 
# for(i in rev(1:length(lambdaBetaValues))){
#   lb = lambdaBetaValues[i]
#   lg= lambdaGammaValues[which(abs(diff( EBICGamMat[i,]))<10^(-4))[1]]
#   res = EBICBet(y= Data$totseiz, X=designX, Z=designZ, family="poisson",
#                 muEstim = 'pen', thetaEstim = 'pen',
#                 beta.ini = beta.iniList, gamma.ini = gamma.iniList,
#                 lambdaBeta=lb, lambdaGamma=lg,
#                 weights.on.xz1 = "covMcd", contX = 4:6,
#                 weightFunction1 = "Huber",
#                 optionList = optionList2)
# 
#   EBICBetMat[i] = res$EBICBet
#   beta.iniList = res$fit$betaEstimate
#   gamma.iniList = res$fit$gammaEstimate
# }
# 
# lb = lambdaBetaValues[which(abs(diff( EBICBetMat))<10^(-4))[1]] #8.05
# lg = lambdaGammaValues[which(abs(diff( EBICGamMat[which(abs(diff( EBICBetMat))<10^(-4))[1],]))<10^(-4))[1]] #1.7

rob.est_Lasso <- RDE(y = Data$totseiz,
                     X=designX, Z=designZ,
                     muEstim = "pen", thetaEstim = "pen",
                     intercept = TRUE, 
                     weights.on.xz1 = "covMcd",
                     contX = 4:6, 
                     weightFunction1 = "Huber",
                     family = "poisson", lambdaBeta=8.05, lambdaGamma=1.7,
                     optionList = optionList2)
rob.est_Lasso$betaEstimate
rob.est_Lasso$gammaEstimate


# lambdaBetaValues = seq(6,7,0.05)
# lambdaGammaValues = seq(1,5,0.05)
# EBICGamMat = matrix(NA, nrow = length(lambdaBetaValues), ncol = length(lambdaGammaValues))
# beta.iniList = NULL
# gamma.iniList = NULL
# beta.iniListStart = NULL
# gamma.iniListStart = NULL
# for(i in rev(1:length(lambdaBetaValues))){
#   lb = lambdaBetaValues[i]
#   for(j in rev(1:length(lambdaGammaValues))){
#     lg = lambdaGammaValues[j]
#     if(j==length(lambdaGammaValues)){
#       res = EBICGam(y= Data$totseiz, X=designX, Z=designZ, family="poisson",
#                     muEstim = 'pen', thetaEstim = 'pen',
#                     beta.ini = beta.iniListStart, gamma.ini = gamma.iniListStart,
#                     lambdaBeta=lb, lambdaGamma=lg,
#                     weights.on.xz1 = "covMcd", contX = 4:6,
#                     weightFunction1 = "Tukey",
#                     optionList = optionList2)
#     }else{
#       res = EBICGam(y= Data$totseiz, X=designX, Z=designZ, family="poisson",
#                     muEstim = 'pen', thetaEstim = 'pen',
#                     beta.ini = beta.iniList, gamma.ini = gamma.iniList,
#                     lambdaBeta=lb, lambdaGamma=lg,
#                     weights.on.xz1 = "covMcd", contX = 4:6,
#                     weightFunction1 = "Tukey",
#                     optionList = optionList2)
#     }
#     
#     EBICGamMat[i,j] = res$EBICGam
#     beta.iniList = res$fit$betaEstimate
#     gamma.iniList = res$fit$gammaEstimate
#     
#     if(j ==length(lambdaGammaValues)){
#       beta.iniListStart = beta.iniList
#       gamma.iniListStart = gamma.iniList
#     }
#   }
# }
# 
# EBICBetMat = rep(NA, length(lambdaBetaValues))
# 
# beta.iniList = NULL
# gamma.iniList = NULL
# 
# for(i in rev(1:length(lambdaBetaValues))){
#   lb = lambdaBetaValues[i]
#   lg= lambdaGammaValues[which(abs(diff( EBICGamMat[i,]))<10^(-4))[1]]
#   res = EBICBet(y= Data$totseiz, X=designX, Z=designZ, family="poisson",
#                 muEstim = 'pen', thetaEstim = 'pen',
#                 beta.ini = beta.iniList, gamma.ini = gamma.iniList,
#                 lambdaBeta=lb, lambdaGamma=lg,
#                 weights.on.xz1 = "covMcd", contX = 4:6,
#                 weightFunction1 = "Tukey",
#                 optionList = optionList2)
#   
#   EBICBetMat[i] = res$EBICBet
#   beta.iniList = res$fit$betaEstimate
#   gamma.iniList = res$fit$gammaEstimate
# }
# 
# lb = lambdaBetaValues[which(abs(diff( EBICBetMat))<10^(-4))[1]] #6.2
# lg = lambdaGammaValues[which(abs(diff( EBICGamMat[which(abs(diff( EBICBetMat))<10^(-4))[1],]))<10^(-4))[1]] #1.25


rob.est_Lasso <- RDE(y = Data$totseiz,
                     X=designX, Z=designZ,
                     muEstim = "pen", thetaEstim = "pen",
                     intercept = TRUE, 
                     weights.on.xz1 = "covMcd",
                     contX = 4:6, 
                     weightFunction1 = "Tukey",
                     family = "poisson", lambdaBeta=6.2, lambdaGamma=1.25,
                     optionList = optionList2)
rob.est_Lasso$betaEstimate
rob.est_Lasso$gammaEstimate
