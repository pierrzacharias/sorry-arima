library(forecast)
library(lmtest)
library(TSA)
#install.packages("dynlm")
library(dynlm)



setwd("/home/pierre/Documents/Courses/2A/Modelisation_des_series_temporelles/MODELISATION-SERIES-TEMPORELLES/MODELISATION SERIES TEMPORELLES/Séance 5/")

seriesJ=read.table(file="seriesJ.txt",sep='\t',header=T)

X=ts(seriesJ[,1])
Y=ts(seriesJ[,2])

# cf code TP1:
par(mar=c(3,4,1,1))  
plot(X, col="blue", ylab="",type="o",pch=20)
mtext("Input Gas",side=2,line=2,col="blue")
par(new=T)   # pour superposer le prochain graphe
plot(Y, col="red", type="o", xaxt="n", yaxt="n", xlab="", ylab="", pch=20)
axis(4)  # tracer l'axe ? droite
mtext("Output CO2",side=4,line=2,col="red")



# On approche X par un ARIMA 
tsdisplay(X,main="ACF Input.Gas")
arima_X=Arima(X,order=c(4,0,0))  # AR3 ou AR4 
tsdisplay(residuals(arima_X),main = "residus de Input.Gas pour un AR[4]")
tsdiag(arima_X)
# On approche aussi Y par un ARIMA
arima_Yf=Arima(Y,model=arima_X)
tsdisplay(residuals(arima_Yf),main = "residus de Output.CO2 avec le modèle de X")
# a l'air stationnaire

# On regarde la ccf (cross correlation) entre ces deux séries temporelles
ccf(residuals(arima_X),residuals(arima_Yf),main="corrélogramme croisé")
# Y_t significativement corrélé avec X_{t-3}, X_{t-4}, X_{t-5}, X_{t-6}, X{t-7} (limite)

# Dynamic Linear models
# On cherche un modèle pour expliquer la corrélation entre les variables
# A la lumière de la CCF on constate des pics sur les zones -3, -4, -5, -6 et -7 
out_lm=dynlm(Y ~ L(X,-3)+L(X,-4)+L(X,-5)+L(X,-6)+L(X,-7))
summary(out_lm)
tsdisplay(residuals(out_lm))
# résidus corrélés ~ AR(2)

output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(0,7)),xtransf=X)
summary(output_XY)
tsdiag(output_XY)
tsdisplay(residuals(output_XY))
coeftest(output_XY)
# on retrouve que les coefs 0 çà 2 ne sont pas significatifs
# par curiosité: output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(0,10)),xtransf=X)

output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(0,7)),fixed=c(NA,NA,NA,0,0,0,NA,NA,NA,NA,NA),xtransf=X)
# pour introduire un décalage, on est obligé de jouer avec les coefficients 
# du numérateur en fixant certains à 0, de manière à factoriser B^b, où b est le décalage
summary(output_XY)
# sigma^2 estimated as 0.05823:  log likelihood = -0.82,  aic = 17.64
tsdiag(output_XY)
tsdisplay(residuals(output_XY))
hist(residuals(output_XY), freq = F, col = "grey")
curve(dnorm(x, mean = mean(residuals(output_XY),na.rm=T), sd = sd(residuals(output_XY),na.rm=T)), col = 2, add = TRUE)
qqnorm(residuals(output_XY))    
qqline(residuals(output_XY))  
coeftest(output_XY)



# fonction de transfert à numérateur et dénominateur degré 1:
output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(1,4)),fixed=c(NA,NA,NA,NA,0,0,0,NA,NA),xtransf=X)
summary(output_XY)
tsdiag(output_XY)
tsdisplay(residuals(output_XY))  # des corr?lations
hist(residuals(output_XY), freq = F, col = "grey")
curve(dnorm(x, mean = mean(residuals(output_XY),na.rm=T), sd = sd(residuals(output_XY),na.rm=T)), col = 2, add = TRUE)
qqnorm(residuals(output_XY))    
qqline(residuals(output_XY))
coeftest(output_XY)
# sigma^2 estimated as 0.06113:  log likelihood = -7.82,  aic = 27.64
# Y= 53.351 - (0.4551+0.6320B)/(1-0.6740B) X_{t-3} + 1/(1-1.5177B+0.6267B^2) eps_t 

# fonction de transfert numérateur degré 2 et dénominateur degré 1:
  output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(1,5)),fixed=c(NA,NA,NA,NA,0,0,0,NA,NA,NA),xtransf=X)
  summary(output_XY)
tsdiag(output_XY)

tsdisplay(residuals(output_XY))  # des corr?lations
hist(residuals(output_XY), freq = F, col = "grey")
curve(dnorm(x, mean = mean(residuals(output_XY),na.rm=T), sd = sd(residuals(output_XY),na.rm=T)), col = 2, add = TRUE)
qqnorm(residuals(output_XY))    
qqline(residuals(output_XY))
coeftest(output_XY)
# sigma^2 estimated as 0.0571:  log likelihood = 2.08,  aic = 9.83
# Y= 53.362 - (0.53096+0.38013B+0.51801B^2)/(1-0.54903B) X_{t-3} + 1/(1-1.5272B+0.6288B^2) eps_t 

# fonction de transfer ? num?rateur degr? 1 et d?nominateur degr? 2:
output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(2,4)),fixed=c(NA,NA,NA,NA,NA,0,0,0,NA,NA),xtransf=X)
summary(output_XY)
tsdiag(output_XY)
tsdisplay(residuals(output_XY))  # mieux pour les corr?lations
hist(residuals(output_XY), freq = F, col = "grey")
curve(dnorm(x, mean = mean(residuals(output_XY),na.rm=T), sd = sd(residuals(output_XY),na.rm=T)), col = 2, add = TRUE)
qqnorm(residuals(output_XY))    
qqline(residuals(output_XY))
coeftest(output_XY)
# sigma^2 estimated as 0.05849:  log likelihood = -1.41,  aic = 16.81
# Y= 53.367 - (0.4826+0.3409B)/(1-1.031B-0.2928B^2) X_{t-3} + 1/(1-1.5215B+0.6232B^2) eps_t 
# rem: si num d?gr? 2 et d?nom degr? 2, coef degr? 2 au d?nom non significatif (? tester)

# meilleur modle:
output_XY=arimax(Y,order=c(2,0,0),transfer=list(c(1,5)),fixed=c(NA,NA,NA,NA,0,0,0,NA,NA,NA),xtransf=X)
summary(output_XY)
