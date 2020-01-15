
library(bsts)     # load the bsts package
data(iclaims)     # bring the initial.claims data into scope

ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons = 52)
model1 <- bsts(initial.claims$iclaimsNSA,
               state.specification = ss,
               niter = 1000)

plot(model1)
plot(model1, "components")  # plot(model1, "comp") works too!
plot(model1, "help")
plot(model1, "predictors")  # plot(model1, "comp") works too!
PlotBstsCoefficients(model1)


pred1 <- predict(model1, horizon = 12)
plot(pred1, plot.original = 500)


## Example 1:  Time series (ts) data
data(AirPassengers)
y <- log(AirPassengers)
ss <- AddLocalLinearTrend(list(), y)
ss <- AddSeasonal(ss, y, nseasons = 12)
model <- bsts(y, state.specification = ss, niter = 500)
pred <- predict(model, horizon = 12, burn = 100)
par(mfrow = c(1,2))
plot(model)
plot(pred)

## Not run: 

MakePlots <- function(model, ask = TRUE) {
  ## Make all the plots callable by plot.bsts.
  opar <- par(ask = ask)
  on.exit(par(opar))
  plot.types <- c("state", "components", "residuals",
                  "prediction.errors", "forecast.distribution")
  for (plot.type in plot.types) {
    plot(model, plot.type)
  }
  if (model$has.regression) {
    regression.plot.types <- c("coefficients", "predictors", "size")
    for (plot.type in regression.plot.types) {
      plot(model, plot.type)
    }
  }
}
MakePlots(model)




Code perso
```{r, echo = TRUE, message=FALSE, eval=TRUE, fig.width=7, fig.height=5}
burn <- SuggestBurn(0.1, bsts.model)
n_sample <- seq(10,500,20)
l <- length(n_sample)
MAPE_list <- rep(0,l)
for (i in 1:l){
  bsts.model <- bsts(y, state.specification = ss, niter = 3000, ping=0, seed=2016)
  p <- predict.bsts(bsts.model, horizon = 12, burn = burn, quantiles = c(.025, .975))
  ### Actual versus predicted
  d2 <- data.frame(
    # fitted values and predictions
    c(10^as.numeric(-colMeans(bsts.model$one.step.prediction.errors[-(1:burn),])+y),  
      10^as.numeric(p$mean)),
    # actual data and dates 
    as.numeric(AirPassengers),
    as.Date(time(AirPassengers)))
  names(d2) <- c("Fitted", "Actual", "Date")
  ### MAPE
  MAPE_list[i] <- filter(d2, year(Date)>1959) %>% summarise(MAPE=mean(abs(Actual-Fitted)/Actual))
}
plot(n_sample,MAPE_list,type='l')

Test
Les erreur retourné sont calculé selont le filtre de Kalman
```{r, echo = TRUE, message=FALSE, eval=TRUE, fig.width=7, fig.height=5}
# one step ahead edoction errors
errors <- bsts.prediction.errors(bsts.model, burn  = SuggestBurn(0.1, bsts.model))
plot(errors)
PlotDynamicDistribution(errors)

## Compute out of sample prediction errors beyond times 80 and 120.
errors <- bsts.prediction.errors(model1, cutpoints = c(80, 120))
standardized.errors <- bsts.prediction.errors(
  model, cutpoints = c(80, 120), standardize = TRUE)
plot(model, "prediction.errors", cutpoints = c(80, 120))
str(errors)     ## three matrices, with 400 ( = 500 - 100) rows
## and length(y) columns
``````

r_data <- data.frame(
  res = as.numeric(r)
)
ggpl
ggplot(r_data, aes(x=res)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="grey") +
  geom_density(kernel="gaussian",alpha=.2, fill="#FF6666") +
  labs(title="histogramme des résidus")+
  theme_bw()




