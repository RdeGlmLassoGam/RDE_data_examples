########################################################################################################################
#########################################################  GLM #########################################################
########################################################################################################################

rm(list=ls())

source("RobGLM/glmrob.R")
source("RobGLM/glmrobMqle.R")

library(sas7bdat)
library(dglm)
library(RdeGlmLassoGam)
library(ggplot2)
#######################
# Initialise Data
#######################
d <- as.data.frame(UCBAdmissions)
d <- tidyr::spread(d, Admit, Freq)
d[order(d$Dept), ]

Data <- d
attach(Data)

formula <- cbind(Admitted, Rejected) ~ Dept + Gender
dformula <- cbind(Admitted, Rejected)  ~  1 
response <- Admitted / (Admitted +  Rejected)

Data$total <- (Admitted +  Rejected)

family <- "binomial"

#######################
# CEQL
#######################
dglm.est <- dglm(data = Data,
                 ind = NULL,
                 formula = formula,
                 dformula = dformula,
                 family = family,
                 method = "reml")
summary.dglm(dglm.est)

round(dglm.est$coefficients, digits = 2)
round(dglm.est$dispersion.fit$coefficients, digits = 2)   #Report with a minus sign as model specification is different from our RDE specification

#######################
# CDE
#######################
# Format inputs
glmCant <- glmrob(formula, family = family,  data = Data)

designX <- model.matrix(glmCant)
designZ <- model.matrix(dglm.est$dispersion.fit)

optionList1 <- list(huberC = 20000,
                    tukeyC = 20000,
                    tol = 1e-4,
                    maxIt = 15)

# Actual fitting
CDE <- RDE(y = response ,X=designX, Z=designZ,
                         weights.on.xz1 = "none", 
                         weights.on.xz2 = "none", 
                         weightFunction1 = "Huber",
                         weightFunction2 = "Huber",
                         mBin = Data$total,
                         family = family,
                         optionList = optionList1)

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
                    maxIt = 15)

# Actual fitting
rob.est <- RDE(y = response ,X=designX, Z=designZ,
                    weights.on.xz1 = "none", weights.on.xz2 = "none", 
                    weightFunction1 = "Huber", weightFunction2 = "Huber",
                    mBin = Data$total,
                    family = family,
                    optionList = optionList2)

linpred <- as.numeric(designX %*% rob.est$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
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
rob.est <- RDE(y = response ,X=designX, Z=designZ,
                                weights.on.xz1 = "none", 
                                weights.on.xz2 = "none", 
                                weightFunction1 = "Tukey",
                                weightFunction2 = "Tukey",
                                mBin = Data$total,
                                family = family,
                                optionList = optionList2)

linpred <- as.numeric(designX %*% rob.est$betaEstimate)
mu <- as.numeric(exp(linpred) / (1 + exp(linpred)))
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

# Generate residual plot
RDEPearson <- (response - mu) / sqrt((1-mu) * mu / ( Data$total * theta))
PlotData <- data.frame(1:12,
                       c("Am", "Bm", "Cm", "Dm", "Em", "Fm",
                         "Af", "Bf", "Cf", "Df", "Ef", "Ff"),
                       RDEPearson, response )
colnames(PlotData) <- c("index", "label", "residual")
Plot <- ggplot(PlotData) + geom_point(aes(x = index, y = residual)) 
Plot <- Plot + geom_text(data = PlotData,
                         aes(x = index, y = residual, label = label),
                         hjust = -0.5, vjust = 0.1)
Plot <- Plot + coord_cartesian(xlim = c(0,13), ylim = c(-2,10))
Plot <- Plot + xlab("Index of observation") + ylab("Pearson residual") +
  ggtitle("Residual plot for RDE estimator") + mrfDepth:::mrfDepth_theme() 
Plot
ggsave(plot = Plot, filename = "UCBAdmissions_residuals.eps", width = 4 * 1.5, height =  3 * 1.5)

#######################
# Robust Testing
#######################
# Format inputs
designX <- designX[,c(1,3,4,5,6,2,7)]
designZ <- designZ

# Actual fitting
rob.est_Final <- RDE(y = response ,X=designX, Z=designZ,
                                weights.on.xz1 = "none", 
                                weights.on.xz2 = "none", 
                                weightFunction1 = "Tukey",
                                weightFunction2 = "Tukey",
                                mBin = Data$total,
                                family = family,
                                optionList = optionList2)
rob.est_Final_NULL <- RDE(y = response ,X=designX[,1:5], Z=designZ,
                                weights.on.xz1 = "none", 
                                weights.on.xz2 = "none", 
                                weightFunction1 = "Tukey",
                                weightFunction2 = "Tukey",
                                mBin = Data$total,
                                family = family,
                                optionList = optionList2)

# Actual testing
DispTest <- DispersionTest(y = response ,X=designX, Z=designZ,
                           m = Data$total,
                           alphaHAlt = c(rob.est_Final$betaEstimate, rob.est_Final$gammaEstimate),
                           alphaHNull = c(rob.est_Final_NULL$betaEstimate,0,0, rob.est_Final_NULL$gammaEstimate),
                           nzeroBeta = 2, #q1
                           nzeroGamma = 0, #q2
                           family=family, 
                           weights.on.xz1 = rob.est_Final$weights.xz1, weights.on.xz2 = rob.est_Final$weights.xz2,
                           weightFunction1 = "Tukey",
                           optionList = optionList2,
                           M = rob.est_Final$M, Q = rob.est_Final$Q)
DispTest

#######################
# CDE without obs 7
#######################
# Format inputs
designX <- model.matrix(glmCant)
designZ <- model.matrix(dglm.est$dispersion.fit)

# Actual fitting
CDENoObs7 <- RDE(y = response[-7] ,X=designX[-7,], Z=designZ[-7,],
                      weights.on.xz1 = "none", 
                      weights.on.xz2 = "none", 
                      weightFunction1 = "Huber",
                      weightFunction2 = "Huber",
                      mBin = Data$total[-7],
                      family = family,
                      optionList = optionList1)
alpha <- c(CDENoObs7$betaEstimate, CDENoObs7$gammaEstimate)
errorbarRDE <- sqrt(diag((1/nrow(Data)) * CDENoObs7$ASVar))
zscore <- alpha / errorbarRDE
cbind(c(colnames(designX), colnames(designZ)),
      alpha,
      errorbarRDE,
      2 * (1 - pnorm(abs(zscore))))

