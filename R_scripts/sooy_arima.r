##############################################################
#   https://multithreaded.stitchfix.com/blog/2016/04/21/forget-arima/
##############################################################

setwd(
  "/home/pierre/Documents/Courses/project_serie_temporelles/time_series_project/"
)


# library -----------------------------------------------------------------
install.packages("tidyverse")
install.packages("bsts")
install.packages("dplyr")

library(TSA)
# library(dynlm)
library(lmtest)
library(forecast)
library(lubridate)
library(bsts)
library(dplyr)
library(ggplot2)

# data --------------------------------------------------------------------
AirPassengers
str(AirPassengers)

start(AirPassengers)
end(AirPassengers)
plot(AirPassengers)         #le graphe sugg?re un mod?le multiplicatif
frequency(AirPassengers)

Y <- window(AirPassengers,
            start = c(1949, 1),
            end = c(1959, 12))
Y.log = log(Y)

AP.log = log(AirPassengers)  #modele multiplicatif donc on passe au log
tsdisplay(AP.log)

# modele arima 1 -----------------------------------------------------------------------

analyse = Arima(
  AP.log,
  order = c(1, 1, 0),
  seasonal = c(2, 0, 1),
  lambda = 0
)
summary(analyse)
tsdisplay(residuals(analyse))

analyse = Arima(
  Y.log,
  order = c(0, 1, 0),
  lambda = 0,
  seasonal = c(0, 0, 1)
)

#   On peut essayer (1,1,0)(2,0,1)
tsdiag(analyse)
summary(analyse)
tsdisplay(residuals(analyse))

# prediction --------------------------------------------------------------


plot(forecast(analyse, h = 8))



# article -----------------------------------------------------------------


### Fit the ARIMA model
arima <- arima(Y.log,
               order = c(0, 1, 1),
               seasonal=list(order=c(0,1,1), period=12))

summary(arima)

### Actual versus predicted

d1 <- data.frame(c(10 ^ as.numeric(fitted(arima)), # fitted and predicted
                 10 ^ as.numeric(predict(arima, n.ahead = 12)$pred)),
                 as.numeric(AirPassengers),
                 #actual values
                 as.Date(time(AirPassengers)))

head(d1)
names(d1) <- c("Fitted", "Actual", "Date")
plot(d1)


### MAPE (mean absolute percentage error)
MAPE <- filter(d1, Date > 1959-01-01) %>% 


summarise(MAPE = mean(abs(Actual - Fitted) /Actual))

### Plot actual versus predicted
ggplot(data = d1, aes(x = Date)) +
  geom_line(aes(y = Actual, colour = "Actual"), size = 1.2) +
  geom_line(aes(y = Fitted, colour = "Fitted"),
            size = 1.2,
            linetype = 2) +
  theme_bw() + theme(legend.title = element_blank()) +
  ylab("") + xlab("") +
  geom_vline(xintercept = as.numeric(as.Date("1959-12-01")), linetype =
               2) +
  ggtitle(paste0("ARIMA -- Holdout MAPE = ", round(100 * MAPE, 2), "%")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))
