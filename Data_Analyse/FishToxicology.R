# Data:
# Author(s): Lisa M. Ganio and Daniel W. Schafer
# Source: Journal of the American Statistical Association, Vol. 87, No. 419 (Sep., 1992), pp. 795-804

#######################
rm(list=ls())
library(RdeGlmLassoGam)
library(ggplot2)
library(dglm)
source("RobGLM/glmrob.R")
source("RobGLM/glmrobMqle.R")

#######################
# Make the data
#######################
a <- c(86,86,88,86,87,90,83,88,90,89,89,87,86,80,89,88,87,88,88,84)
b <- c(3,5,4,2,14,14,9,12,29,31,33,26,44,40,44,43,62,67,59,58)
c <- c(87,86,89,85,86,86,86,88,89,86,90,88,88,89,88,90,86,82,81,89)
d <- c(9,5,2,9,30,41,27,34,54,53,64,55,71,73,65,72,66,75,72,73)
e <- c(0.01,0.01,0.01,0.01,0.025,0.025,0.025,0.025,0.05,0.05,0.05,0.05,0.1,0.1,0.1,.1,.25,.25,.25,.25)
Data <- cbind(c(e,e),c(a,c), c(b,d), c(rep(0,length(e)),rep(1,length(e))))
colnames(Data) <- c("dose", "tot", "cancer", "drug")
Data <- data.frame(Data)
Data$drug <- as.factor(Data$drug)
Data$doseF <- as.factor(Data$dose)
attach(Data)

PlotData <- data.frame(1:40, Data$dose, Data$cancer/Data$tot, Data$drug)
colnames(PlotData) <- c("index", "dosage", "cancerRate", "Treatment")
Plot <- ggplot(PlotData) 
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = Treatment)) 
Plot <- Plot + scale_colour_manual(values = c("blue", "orange"), labels = c("Aflatoxin B1", "Aflatoxicol"))
Plot <- Plot + xlab("dosage (ppm)") + ylab("tumor incidence rate") +
  ggtitle("Fish Toxicology data") + mrfDepth:::mrfDepth_theme() 
Plot
ggsave(filename = "Fish_Data.eps", plot = Plot, width = 4 * 1.5, height =  3 * 1.5)

#######################
# CEQL
#######################
cdglm <- dglm(cbind(cancer, tot - cancer) ~ drug + doseF, ~ drug, family = "binomial", data = Data)
summary.dglm(cdglm)

round(cdglm$coefficients, digits = 2)
round(cdglm$dispersion.fit$coefficients, digits = 2)   #Report with a minus sign as model specification is different from our RDE specification

#######################
# CDE
#######################
# Format inputs
robglm <- glmrob(data = Data, cbind(cancer, tot - cancer) ~ drug + doseF, family = "binomial")

ResponsVariable <- Data$cancer / Data$tot
m <- Data$tot

designX <- model.matrix(robglm)
designZ <- cbind(rep(1, nrow(designX)), designX[,2])
colnames(designZ) <- c("Intercept", "Trt")

optionList1 <- list(huberC = 20000,
                   tukeyC = 20000,
                   tol = 1e-4,
                   maxIt = 100)

# Actual fitting
CDE <- RDE(y = ResponsVariable,
                          X = designX,
                          Z = designZ,
                          mBin = m,
                          family = "binomial",
                          weightFunction1 = "Huber",
                          optionList = optionList1)

linpred <- as.numeric(designX %*% CDE$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
theta <- as.numeric(exp(designZ %*% CDE$gammaEstimate))
1 / theta

# Format results
alpha <- c(CDE$betaEstimate, CDE$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * CDE$ASVar))
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
                    maxIt = 100)

# Actual fitting
HRDE <- RDE(y = ResponsVariable,
                       X = designX,
                       Z = designZ,
                       mBin = m,
                       family = "binomial",
                       weightFunction1 = "Huber",
                       optionList = optionList2 )

linpred <- as.numeric(designX %*% HRDE$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
theta <- as.numeric(exp(designZ %*% HRDE$gammaEstimate))
1 / theta

# Format results
alpha <- c(HRDE$betaEstimate, HRDE$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * HRDE$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

#######################
# TRDE 
#######################
TRDE <- RDE(y = ResponsVariable,
                          X = designX,
                          Z = designZ,
                          mBin = m,
                          family = "binomial",
                          weightFunction1 = "Tukey",
                          optionList = optionList2 )

linpred <- as.numeric(designX %*% TRDE$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
theta <- as.numeric(exp(designZ %*% TRDE$gammaEstimate))
1 / theta

# Format results
alpha <- c(TRDE$betaEstimate, TRDE$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * TRDE$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

############################################################
# ADD SOME OUTLIERS
############################################################
a <- c(86,86,88,86,87,90,83,88,90,89,89,87,86,80,89,88,87,88,88,84)
b <- c(3,5,4,2,14,14,9,12,29,31,33,26,44,40,44,43,62,67,59,58)
c <- c(87,86,89,85,86,86,86,88,89,86,90,88,88,89,88,90,86,82,81,89, 90, 90)
d <- c(9,5,2,9,30,41,27,34,54,53,64,55,71,73,65,72,66,75,72,73, 10, 10)
e <- c(0.01,0.01,0.01,0.01,0.025,0.025,0.025,0.025,0.05,0.05,0.05,0.05,0.1,0.1,0.1,.1,.25,.25,.25,.25)
Data <- cbind(c(e,e,.1,.05),c(a,c), c(b,d), c(rep(0,length(e)),rep(1,length(e)), 0,0))
colnames(Data) <- c("dose", "tot", "cancer", "drug")
Data <- data.frame(Data)
Data$drug <- as.factor(Data$drug)
Data$doseF <- as.factor(Data$dose)

PlotData <- data.frame(1:42, Data$dose, Data$cancer/Data$tot, Data$drug)
colnames(PlotData) <- c("index", "dosage", "cancerRate", "Treatment")
Plot <- ggplot(PlotData[1:40,]) 
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = Treatment))
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = Treatment), data = PlotData[41:42,], shape = 17, size = 2)
Plot <- Plot + scale_colour_manual(values = c("blue", "orange"), labels = c("Aflatoxin B1", "Aflatoxicol"))
Plot <- Plot + xlab("dosage (ppm)") + ylab("tumor incidence rate") +
  ggtitle("Fish Toxicology data") + mrfDepth:::mrfDepth_theme() 
Plot
ggsave(filename = "Fish_Data.eps", plot = Plot, width = 4 * 1.5, height =  3 * 1.5)

#######################
# CEQL 
#######################
cdglm <- dglm(cbind(cancer, tot) ~ drug + doseF, ~ drug, family = "binomial", data = Data)
summary.dglm(cdglm)

round(cdglm$coefficients, digits = 2)
round(cdglm$dispersion.fit$coefficients, digits = 2)   #Report with a minus sign as model specification is different from our RDE specification

#######################
# CDE
#######################
# Format inputs
robglm <- glmrob(data = Data, cbind(cancer, tot - cancer) ~ drug + doseF, family = "binomial")
ResponsVariable <- Data$cancer / Data$tot
m <- Data$tot
designX <- model.matrix(robglm)
designZ <- cbind(rep(1, nrow(designX)), designX[,2])
colnames(designZ) <- c("Intercept", "Trt")

# Actual fitting
CDE <- RDE(y = ResponsVariable,
                      X = designX,
                      Z = designZ,
                      mBin <- m,
                      family = "binomial",
                      weightFunction1 = "Huber",
                      optionList = optionList1 )

# Format results
alpha <- c(CDE$betaEstimate, CDE$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * CDE$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

#######################
# HRDE 
#######################
HRDE <- RDE(y = ResponsVariable,
                          X = designX,
                          Z = designZ,
                          mBin <- m,
                          family = "binomial",
                          weightFunction1 = "Huber",
                          optionList = optionList2 )

# Format results
alpha <- c(HRDE$betaEstimate, HRDE$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * HRDE$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

#######################
# TRDE 
#######################
TRDE <- RDE(y = ResponsVariable,
                       X = designX,
                       Z = designZ,
                       mBin <- m,
                       family = "binomial",
                       weightFunction1 = "Tukey",
                       optionList = optionList2 )

# Format results
alpha <- c(TRDE$betaEstimate, TRDE$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * TRDE$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

############################################################
# Generate residual plots
############################################################
HuberResid <- function(x, c){
  x[x <= -c] <- -c
  x[x >= c] <- c
  return(x)
}

TukeyResid <- function(x, c){
  response <- ((x/c)^2-1)^2*x
  Ind <- which(abs(x)>c)
  if(length(Ind)>0){ response[Ind] <- 0 }
  return( response )
}
#######################
# TRDE
#######################
linpred <- as.numeric(designX %*% TRDE$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
theta <- as.numeric(exp(designZ %*% TRDE$gammaEstimate))
y <- ResponsVariable
PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) / (m * theta))
plot(PearsonResid)
robWeights <- TukeyResid(PearsonResid, optionList2$tukeyC) / PearsonResid
plot(robWeights)

# See what observations are downweighted
PlotData <- data.frame(1:42, Data$dose, Data$cancer/Data$tot, Data$drug, robWeights)
colnames(PlotData) <- c("index", "dosage", "cancerRate", "Treatment", "weight")
Plot <- ggplot(PlotData[1:40,]) 
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = weight)) 
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = weight), data = PlotData[41:42,], shape = 17) 
Plot <- Plot + xlab("dosage (ppm)") + ylab("tumor incidence rate") +
  ggtitle("Robustness weights") + mrfDepth:::mrfDepth_theme() 
Plot
ggsave(filename = "Fish_weights.eps", plot = Plot, width = 4 * 1.5, height =  3 * 1.5)

#######################
# HRDE
#######################
linpred <- as.numeric(designX %*% HRDE$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
theta <- as.numeric(exp(designZ %*% HRDE$gammaEstimate))
y <- ResponsVariable
PearsonResid <- (y - mu) / sqrt(mu * (1 - mu) / (m * theta))
plot(PearsonResid)
robWeights <- TukeyResid(PearsonResid, optionList2$huberC) / PearsonResid
plot(robWeights)

# See what observations are downweighted
PlotData <- data.frame(1:42, Data$dose, Data$cancer/Data$tot, Data$drug, robWeights)
colnames(PlotData) <- c("index", "dosage", "cancerRate", "Treatment", "weight")
Plot <- ggplot(PlotData[1:40,]) 
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = weight)) 
Plot <- Plot + geom_point(aes(x = dosage, y = cancerRate, color = weight), data = PlotData[41:42,], shape = 17) 
Plot <- Plot + xlab("dosage (ppm)") + ylab("tumor incidence rate") +
  ggtitle("Robustness weights") + mrfDepth:::mrfDepth_theme() 
Plot