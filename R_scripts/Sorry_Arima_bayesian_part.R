# #######################################################
#    Sorry Arima Baysian model part
# #######################################################
install.packages("lubridate")
library(lubridate)
library(bsts)
library(ggplot2)
library(reshape2)

### Fit the model with regressors
data(iclaims)
ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
bsts.reg <- bsts(iclaimsNSA ~ ., state.specification = ss, data =&&
                   initial.claims, niter = 500, ping=0, seed=2016)

### Get the number of burn-ins to discard

burn <- SuggestBurn(0.1, bsts.reg)

### Helper function to get the positive mean of a vector
PositiveMean <- function(b) {
  b <- b[abs(b) > 0]
  if (length(b) > 0) 
    return(mean(b))
  return(0)
}

### Get the average coefficients when variables were selected (non-zero slopes)
coeff <- data.frame(melt(apply(bsts.reg$coefficients[-(1:burn),], 2, PositiveMean)))
coeff$Variable <- as.character(row.names(coeff))
ggplot(data=coeff, aes(x=Variable, y=value)) + 
  geom_bar(stat="identity", position="identity") + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) +
  xlab("") + ylab("") + ggtitle("Average coefficients")



