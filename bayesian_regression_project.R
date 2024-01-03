library(tidyr)
library(mvtnorm)
library(manipulate)
library(metRology)
library(rstan)
library(bridgesampling)
library(coda)
library(BMS)

setwd("C:/Users/User/Desktop/Szkola/Ekonometria bayesowska/")

rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")

#Kolory uzywane w kodzie
grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)
red_area <- rgb(255, 100, 123, 80, names = NULL, maxColorValue = 255)
red_line <- rgb(200, 0, 30, 160, names = NULL, maxColorValue = 255)

#Wczytywanie danych
fossil_fuels<-read.csv("fossil_fuels.csv")
fossil_fuels<- fossil_fuels[,!names(fossil_fuels) %in% 'ï..Country.Name']

GDPpc<-read.csv("GDPpc.csv")
GDPpc<- GDPpc[,!names(GDPpc) %in% 'ï..Country.Name']

urban<-read.csv("urban_population.csv")
urban<- urban[,!names(urban) %in% 'ï..Country.Name']

emissions<-read.csv("emissions.csv")
emissions<- emissions[,!names(emissions) %in% 'ï..Country.Name']

colnames(fossil_fuels)[1]<-'Code'
colnames(fossil_fuels)[2]<-'2000'
colnames(fossil_fuels)[3]<-'2014'
fossil_fuels<-fossil_fuels[order(fossil_fuels$Code),]

colnames(GDPpc)[1]<-'Code'
colnames(GDPpc)[2]<-'2000'
colnames(GDPpc)[3]<-'2014'
GDPpc<-GDPpc[order(GDPpc$Code),]

colnames(urban)[1]<-'Code'
colnames(urban)[2]<-'2000'
colnames(urban)[3]<-'2014'
urban<-urban[order(urban$Code),]

colnames(emissions)[1]<-'Code'
colnames(emissions)[2]<-'2000'
colnames(emissions)[3]<-'2014'
emissions<-emissions[order(emissions$Code),]

#Literatura, z której zostaly wziete zmienne uzykala wyniki mówiace o jej istotnosci uzywajac metod BMA i LASSO, gdzie zalozono normalnosc 

model.prior<-lm(emissions$'2000'~GDPpc$'2000'+fossil_fuels$'2000'+urban$'2000')
summary(model.prior)

shapiro.test(model.prior$residuals)

#Reszty nie wykazuja normalnosci, nalezy wiec zastosowac inny rozklad.

hist(model.prior$residuals, breaks=10)

qqplot(rcauchy(1000), model.prior$residuals, main="Q-Q plot",
       ylab="Sample Quantiles", xlim=c(-10,10))
abline(0,1)

#Rozklad Cauchy'ego wydaje sie byc dobrym przyblizeniem rozkladu reszt (z wylaczeniem 2 przypadków odstajacych, znajdujacych sie poza wykresem, dla uproszczenia, pozostaja w modelu ale sa ignorowane).
#Zalozone zostanie, ze wszystkie parametry regresji sa rozlozone wg. tego samego rozkladu Cauchy'ego wynikajacego z reszt regresji.

#Parametry rozkladu a priori

#Wartosci oczekiwane parametrów
model.prior$coefficients

beta.p.GDPpc<-0.0001402047
beta.p.fossil_fuels<-0.0726126677
beta.p.urban<-0.0895386803
intercept<--5.6926278646

beta.prior<-c(beta.p.GDPpc, beta.p.fossil_fuels, beta.p.urban, intercept)
rm(beta.p.GDPpc, beta.p.fossil_fuels, beta.p.urban, intercept)

#Estymacja parametru lokalizacji (location) za pomoca mediany
loc<-median(model.prior$residuals)

#Estymacja parametru skali (scale) za pomoca przedzialu miedzykwartylowego
sc<-IQR(model.prior$residuals)

#Lokalizacja zostanie uzyta jako wartosc oczekiwana skladnika losowego
#Skala zostanie uzyta jako "wariancja" skladnika losowego i parametrów beta
#Wartosci oczekiwane a priori beta to wartosci oszacowania OLS dla danych a priori

#Dodatkowe parametry R stan
N<-38
y<-emissions$'2014'
beta<-beta.prior[-4]
inter<-beta.prior[4]
sigma<-sd(y)
GDPpc_X<-GDPpc$'2014'
fossil_fuels_X<-fossil_fuels$'2014'
urban_X<-urban$'2014'

#Wartosci x i y dla których pózniej wykonamy predykcje - dane dla Argentyny
x_tau <- c(12334.798, 87.72241, 91.377)
y_tau <- 4.209096

#R stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

model_fit <- stan(file = "stan.stan",
                     data = c("N", "y", "beta", "sigma","loc", "sc", "GDPpc_X", "fossil_fuels_X", "urban_X", "inter", "x_tau", "y_tau"),
                     iter = 5000,
                     chains = 5)
print(model_fit)
plot(model_fit)

#Rozklady brzegowe zmiennych
model.OLS<-lm(emissions$'2014'~GDPpc$'2014'+fossil_fuels$'2014'+urban$'2014')
beta.OLS<-model.OLS$coefficients
temp<-beta.OLS[1]
beta.OLS<-beta.OLS[-1]
beta.OLS<-c(beta.OLS, temp)
rm(temp)

post_samples <- extract(model_fit)
post_samples<-cbind(post_samples$beta, post_samples$inter)
colnames(post_samples)<-c('GDPpc', 'fossil_fuels', 'urban', 'Intercept')

beta.post<-as.numeric(get_posterior_mean(model_fit)[,6])
beta.post<-beta.post[-6]
beta.post<-beta.post[-5]

DF<-data.frame("beta prior"=round(beta.prior, 5), "beta posterior"=round(beta.post, 5))
rownames(DF)<-c('GDPpc', 'fossil_fuels', 'urban', 'Intercept')
View(DF)

par(mfrow = c(2, 2))

{
plot(density(post_samples[,1]), col=green_line, type='l', main=colnames(post_samples)[1], ylab='gestosc')
abline(v=beta.post[1], col=green_line)
text(beta.post[1], 15000, paste0('E(beta) a posteriori=', round(beta.post[1], 6)), col=green_line)
abline(v=beta.prior[1], col=grey_line)
text(beta.prior[1], 14000, paste0('E(beta) a priori=', round(beta.prior[1], 6)), col=grey_line)
abline(v=beta.OLS[1], col=red_line)
text(beta.OLS[1], 13000, paste0('Oszacowanie OLS=', round(beta.OLS[1], 6)), col=red_line)

for (i in 2:3) {
  plot(density(post_samples[,i]), col=green_line, type='l', main=colnames(post_samples)[i], ylab='gestosc')
  abline(v=beta.post[i], col=green_line)
  text(mean(post_samples[,i]), 9, paste0('E(beta) a posteriori=', round(beta.post[i], 6)), col=green_line)
  abline(v=beta.prior[i], col=grey_line)
  text(beta.prior[i], 8, paste0('E(beta) a priori=', round(beta.prior[i], 6)), col=grey_line)
  abline(v=beta.OLS[i], col=red_line)
  text(beta.OLS[i], 7, paste0('Oszacowanie OLS=', round(beta.OLS[i], 6)), col=red_line)
}

plot(density(post_samples[,4]), col=green_line, type='l', main=colnames(post_samples)[4], ylab='gestosc')
abline(v=beta.post[4], col=green_line)
text(beta.post[4], 0.15, paste0('E(beta) a posteriori=', round(beta.post[4], 6)), col=green_line)
abline(v=beta.prior[4], col=grey_line)
text(beta.prior[4], 0.13, paste0('E(beta) a priori=', round(beta.prior[4], 6)), col=grey_line)
abline(v=beta.OLS[4], col=red_line)
text(beta.OLS[4], 0.12, paste0('Oszacowanie OLS=', round(beta.OLS[4], 6)), col=red_line)
}

#Czynniki Bayesa
#Estymacja modeli bez jednej zmiennej
model_fit2 <- stan(file = "stan2.stan",
                  data = c("N", "y", "beta", "sigma","loc", "sc", "fossil_fuels_X", "urban_X", "inter"),
                  iter = 5000,
                  chains = 5)

model_fit3 <- stan(file = "stan3.stan",
                   data = c("N", "y", "beta", "sigma","loc", "sc", "GDPpc_X", "urban_X", "inter"),
                   iter = 5000,
                   chains = 5)

model_fit4 <- stan(file = "stan4.stan",
                   data = c("N", "y", "beta", "sigma","loc", "sc", "GDPpc_X", "fossil_fuels_X", "inter"),
                   iter = 5000,
                   chains = 5)

log1<-bridge_sampler(model_fit)$logml
log2<-bridge_sampler(model_fit2)$logml
log3<-bridge_sampler(model_fit3)$logml
log4<-bridge_sampler(model_fit4)$logml
log<-cbind(log1, log2, log3, log4)

rm(model_fit2, model_fit3, model_fit4, log1, log2, log3, log4)

BF <- rep(NA, 3)

for (i in 1:3) {
  BF[i] <- exp(as.numeric(log[1])) / exp(as.numeric(log[i+1]))
}
BF<-data.frame('Zmienna'=c('GDPpc', 'fossil_fuels', 'urban'), 'Czynnik Bayesa'=BF)
View(BF)

#Przedzialy HPDI
HPDI1<-round(HPDinterval(mcmc(post_samples[,1])), 6)
HPDI2<-round(HPDinterval(mcmc(post_samples[,2])), 6)
HPDI3<-round(HPDinterval(mcmc(post_samples[,3])), 6)

HPDI<-data.frame('Zmienna'=c('GDPpc', 'fossil_fuels', 'urban'), 'Górna granica 95 proc.'=c(HPDI1[1], HPDI2[1], HPDI3[1]), 'Dolna granica 95 proc.'=c(HPDI1[2], HPDI2[2], HPDI3[2]))
rm(HPDI1, HPDI2, HPDI3)
View(HPDI)

#Rozklad predykcyjny i prognoza punktowa
#Wszystkie prognozy wykonane na danych dla Argentyny z 2014 r.
pred_dist<-extract(model_fit)$y_tau
par(mfrow = c(1, 1))
plot(density(pred_dist), col=green_line, ylab = "Gestosc predykcyjna", main = "Gestosc predykcyjna")
abline(v=as.numeric(get_posterior_mean(model_fit)[,6])[5], col=green_line)
text(x = as.numeric(get_posterior_mean(model_fit)[,6])[5], y = 0.08, 
     paste("prognoza punktowa = ", round(as.numeric(get_posterior_mean(model_fit)[,6])[5], 4)), col = green_line)

#Prognoza przedzialowa - 80% przedzial ufnosci
quantile(density(pred_dist), probs=c(0.1, 0.9), names=T, normalize=T)
