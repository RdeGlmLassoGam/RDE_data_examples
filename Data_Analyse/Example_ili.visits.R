library(RdeGlmLassoGam)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

############ Real data example
load("ili.visits.RData")
ili.visits <-ili.visits[order(ili.visits$week, ili.visits$season,ili.visits$visits),]
ili.visits$visits=ili.visits$visits

x <- ili.visits$week
## The model is scale sensitive, so SCALE the response to get reasonable models!
y <- ili.visits$visits/1000
n <- length(y)

# lambdaBetaValues = seq(5,10,0.05)
# lambdaGammaValues = seq(5,10,0.05)
# RBICGamMatExtra = matrix(NA, nrow = length(lambdaBetaValues), ncol = length(lambdaGammaValues))
# beta.iniList = NULL
# gamma.iniList = NULL
# beta.iniListStart = NULL
# gamma.iniListStart = NULL
# for(i in rev(1:length(lambdaBetaValues))){
#   print(Sys.time());print(i)
#   lb = lambdaBetaValues[i]
#   for(j in rev(1:length(lambdaGammaValues))){
#     lg = lambdaGammaValues[j]
#     if(j==length(lambdaGammaValues)){
#       res = RBICGam(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
#                     smooth.basis1="cr", smooth.basis2="cr",K1=10, 
#                     beta.ini = beta.iniListStart, gamma.ini = gamma.iniListStart,
#                     weightFunction1 = "Tukey", weightFunction2 = "Tukey",
#                     muEstim = 'gam', thetaEstim = 'gam', sp1 =lb, sp2 =lg)
#     }else{
#       res = RBICGam(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
#                     smooth.basis1="cr", smooth.basis2="cr",K1=10, 
#                     beta.ini = beta.iniList, gamma.ini = gamma.iniList,
#                     weightFunction1 = "Tukey", weightFunction2 = "Tukey",
#                     muEstim = 'gam', thetaEstim = 'gam', sp1 =lb, sp2 =lg)
#     }
#     RBICGamMatExtra[i,j] = res$RBICGam
#     beta.iniList = res$fit$betaEstimate
#     gamma.iniList = res$fit$gammaEstimate
#     save(RBICGamMatExtra, file='RBICGamMatExtra.RData')
#     save(beta.iniList, file='beta.iniList.RData')
#     save(gamma.iniList, file='gamma.iniList.RData')
#     if(j ==length(lambdaGammaValues)){
#       beta.iniListStart = beta.iniList
#       gamma.iniListStart = gamma.iniList
#       save(gamma.iniListStart, file='gamma.iniListStart.RData')
#       save(beta.iniListStart, file='beta.iniList.RData')
#     }
#   }
# }
# 
# RBICBetMat = rep(NA, length(lambdaBetaValues))
# 
# for(i in rev(1:length(lambdaBetaValues))){
#   lb = lambdaBetaValues[i]
#   lg= lambdaGammaValues[which.min( RBICGamMat[i,])]
#   res = RBICBet(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
#                 smooth.basis1="cr", smooth.basis2="cr",K1=10, 
#                 beta.ini = beta.iniList, gamma.ini = gamma.iniList,
#                 weightFunction1 = "Tukey", weightFunction2 = "Tukey",
#                 muEstim = 'gam', thetaEstim = 'gam', sp1 =lb, sp2 =lg)
#   
#   RBICBetMat[i] = res$RBICBet
# }
# 
# lb = lambdaBetaValues[which.min( RBICBetMat)] #9.85
# lg = lambdaGammaValues[which.min( RBICGamMat[which.min( RBICBetMat),])] #6
 
# GLM-GLM model:
modelRef1 = RDE(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
                smooth.basis1="cr", smooth.basis2="cr",K1=10,
                weightFunction1 = "Tukey", weightFunction2 = "Tukey",
                muEstim = 'glm', thetaEstim = 'glm', sp1 =9.85, sp2 =6)

# GAM-GLM model: 
modelRef2 = RDE(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
                smooth.basis1="cr", smooth.basis2="cr",K1=10,
                weightFunction1 = "Tukey", weightFunction2 = "Tukey",
                muEstim = 'gam', thetaEstim = 'glm', sp1 =9.85, sp2 =6)

# GLM-GAM model: 
modelRef2b = RDE(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
                 smooth.basis1="cr", smooth.basis2="cr",K1=10,
                 weightFunction1 = "Tukey", weightFunction2 = "Tukey",
                 muEstim = 'glm', thetaEstim = 'gam', sp1 =9.85, sp2 =6)

# GAM-GAM model: 
modelRef3 = RDE(y, X=cbind(x), Z=cbind(x), family= poisson(link="log"),
                smooth.basis1="cr", smooth.basis2="cr",K1=10, 
                weightFunction1 = "Tukey", weightFunction2 = "Tukey",
                muEstim = 'gam', thetaEstim = 'gam', sp1 =9.85, sp2 =6)

# check the mean models
plot(x,ili.visits$visits, xlab='Week', ylab='ILI visits', pch=19, col=as.factor(ili.visits$season))
legend("topleft", legend = levels(as.factor(ili.visits$season)), pch=19, col= unique(as.factor(ili.visits$season)))
legend("topright", legend = c("gam-gam","gam-glm","glm-gam","glm-glm"), col= c("black","blue","gray","cyan"),lwd=rep(2,4))
lines(x,1000*modelRef1$fitted.mu.values,col='cyan',lwd=2)
lines(x,1000*modelRef2$fitted.mu.values,col='blue',lwd=2)
lines(x,1000*modelRef2b$fitted.mu.values,col='gray',lwd=2)
lines(x,1000*modelRef3$fitted.mu.values,col='black',lwd=2)