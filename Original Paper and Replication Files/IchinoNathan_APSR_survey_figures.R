############################################################################################
###  Part 3 of 3: 
###  Ichino, Nahomi, and Noah L. Nathan.  2013.  "Crossing the Line: Local Ethnic Geography 
###  and Voting in Ghana," American Political Science Review 107(2): 344-61.
###  http://dx.doi.org/10.1017/S0003055412000664
###
###  IchinoNathan_APSR_survey_figures.R
###
###  Nahomi Ichino and Noah Nathan
###  Department of Government, Harvard University
###  February 2013
###  
###  Note: 
###  This reproduces Figures 4 and 5 of the published article.
###  We used Zelig version 3.5.5 for these plots.
###
##############
###  Calls: 
###  IchinoNathan_APSR_survey_data.csv
###  
###  Data available at:
###  http://dvn.iq.harvard.edu/dvn/dv/nichino
###  http://hdl.handle.net/1902.1/21461  
############################################################################################

#library(foreign)
library(Zelig)
library(car)
library(sandwich)
library(lmtest)
library(calibrate)



#### Figure 4 ####
set.seed(02138)

data<-read.csv("IchinoNathan_APSR_survey_data.csv", header=TRUE) 
data <- data[,c("eth_akan", "eth_ewe", "eth_dagomba", "central", "akan_30km_l_p", "ewe_30km_l_p", "vote_npp_pres", "vote_ndc_pres", "male", "economy_oneyear", "poverty", "dev_factor2", "r4", "fid0725", "popdens5x5")]
data <- na.omit(data)
dim(data)

data <- data[data$popdens5x5<1000,]
dat_akan <- data[data$eth_akan==1,]
dim(dat_akan)
dat_ewe <- data[data$eth_ewe==1,]
dim(dat_ewe)
dat_mole <- data[data$eth_dagomba==1,]
dim(dat_mole)
dat_akan$int<-dat_akan$central*dat_akan$akan_30km_l_p
dat_ewe$int<-dat_ewe$central*dat_ewe$akan_30km_l_p
data$int<-data$central*data$akan_30km_l_p



#################
## AKAN 30 & NPP
#################
### for Akans:

m1 <- glm(vote_npp_pres ~ akan_30km_l_p + int + dev_factor2 + male  +economy_oneyear + poverty + central + r4, data=dat_akan, family=binomial(link="logit"))
w.dist<-cbind(dat_akan$vote_npp_pres, dat_akan$akan_30km_l_p , dat_akan$int, dat_akan$dev_factor2, dat_akan$male, dat_akan$economy_oneyear, dat_akan$poverty,  dat_akan$central, dat_akan$r4, dat_akan$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,ncol(w.dist)])
mat <- estfun(m1)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m1$rank
df <- (length(unique(w.dist[,ncol(w.dist)])) / (length(unique(w.dist[,ncol(w.dist)])) - 1))*((N - 1)/ (N - K))
vcovCL <- df*sandwich(m1, meat=crossprod(u)/N)

a.akan_30km.vals <- seq(quantile(dat_akan$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_akan$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m1<- setx(m1,r4=1,male=1,economy_oneyear=3,central=0,int=0)
x.m1 <- as.matrix(x.m1)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL)

predictions.m1 <- matrix(NA, nrow=1000, ncol=length(a.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(a.akan_30km.vals)){
	x.m1[2] <- a.akan_30km.vals[j]
	predictions.m1[i, j] <- 1 / (1 + exp(-x.m1 %*% betas[i,]))
		}
	}
m1.mean <- apply(X = predictions.m1, MARGIN=2, FUN=mean)
m1.lower <- apply(predictions.m1, 2, quantile, 0.025)
m1.upper <- apply(predictions.m1, 2, quantile, 0.975)



### for Ewes:
m2 <- glm(vote_npp_pres ~ akan_30km_l_p + int  + dev_factor2 + male  +economy_oneyear + poverty + central + r4, data=dat_ewe, family=binomial(link="logit"))
w.dist<-cbind(dat_ewe$vote_npp_pres, dat_ewe$akan_30km_l_p, dat_ewe$int, dat_ewe$dev_factor2, dat_ewe$male, dat_ewe$economy_oneyear, dat_ewe$poverty, dat_ewe$central, dat_ewe$r4, dat_ewe$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,ncol(w.dist)])
mat <- estfun(m2)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m2$rank
df <- (length(unique(w.dist[,ncol(w.dist)])) / (length(unique(w.dist[,ncol(w.dist)])) - 1))*((N - 1)/ (N - K))
vcovCL2 <- df*sandwich(m2, meat=crossprod(u)/N)
final<-coeftest(m2, vcovCL2)

e.akan_30km.vals <- seq(quantile(dat_ewe$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_ewe$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m2<- setx(m2,r4=1,male=1,economy_oneyear=3,central=0)
x.m2 <- as.matrix(x.m2)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL2)
predictions.m2 <- matrix(NA, nrow=1000, ncol=length(e.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(e.akan_30km.vals)){
	x.m2[2] <- e.akan_30km.vals[j]
	predictions.m2[i, j] <- 1 / (1 + exp(-x.m2 %*% betas[i,]))
		}
	}
m2.mean <- apply(X = predictions.m2, MARGIN=2, FUN=mean)
m2.lower <- apply(predictions.m2, 2, quantile, 0.025)
m2.upper <- apply(predictions.m2, 2, quantile, 0.975)

### For Mole/Dagbons:
m3 <- glm(vote_npp_pres ~ akan_30km_l_p  + dev_factor2 + male  +economy_oneyear + poverty + r4, data=dat_mole, family=binomial(link="logit"))
w.dist<-cbind(dat_mole$vote_npp_pres, dat_mole$akan_30km_l_p, dat_mole$dev_factor2, dat_mole$male, dat_mole$economy_oneyear, dat_mole$poverty, dat_mole$r4, dat_mole$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,8])
mat <- estfun(m3)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m3$rank
df <- (length(unique(w.dist[,8])) / (length(unique(w.dist[,8])) - 1))*((N - 1)/ (N - K))
vcovCL3 <- df*sandwich(m3, meat=crossprod(u)/N)
final<-coeftest(m3, vcovCL3)

m.akan_30km.vals <- seq(quantile(dat_mole$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_mole$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m3<- setx(m3,r4=1,male=1,economy_oneyear=3)
x.m3 <- as.matrix(x.m3)
betas <- mvrnorm(n = 1000, mu=summary(m3)$coefficients[,1], vcovCL3)
predictions.m3 <- matrix(NA, nrow=1000, ncol=length(m.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(m.akan_30km.vals)){
	x.m3[2] <- m.akan_30km.vals[j]
	predictions.m3[i, j] <- 1 / (1 + exp(-x.m3 %*% betas[i,]))
		}
	}
m3.mean <- apply(X = predictions.m3, MARGIN=2, FUN=mean)
m3.lower <- apply(predictions.m3, 2, quantile, 0.025)
m3.upper <- apply(predictions.m3, 2, quantile, 0.975)


#################
#### EVERYONE TOGETHER

m4<-glm(vote_npp_pres~ akan_30km_l_p +int + eth_akan + eth_ewe + eth_dagomba + male + economy_oneyear + poverty + dev_factor2 + central + r4, family=binomial(link="logit"), data=data)

### clustered standard errors
w.dist<-cbind(data$vote_npp_pres, data$akan_30km_l_p ,data$int,  data$eth_akan, data$eth_ewe, data$eth_dagomba, data$male, data$economy_oneyear, data$poverty, data$dev_factor2, data$central, data$r4 ,data$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,ncol(w.dist)])
mat <- estfun(m4)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m4$rank
df <- (length(unique(w.dist[,ncol(w.dist)])) / (length(unique(w.dist[,ncol(w.dist)])) - 1))*((N - 1)/ (N - K))
vcovCL4 <- df*sandwich(m4, meat=crossprod(u)/N)
final<-coeftest(m4, vcovCL4)

all.akan_30km.vals <- seq(quantile(data$akan_30km_l_p, probs=c(.1,.9))[1], quantile(data$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
set.seed(02138)
x.m4<- setx(m4,r4=1,male=1,economy_oneyear=3, central=0)
x.m4 <- as.matrix(x.m4)
betas <- mvrnorm(n = 1000, mu=summary(m4)$coefficients[,1], vcovCL4)
predictions.m4 <- matrix(NA, nrow=1000, ncol=length(all.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(all.akan_30km.vals)){
	x.m4[2] <- all.akan_30km.vals[j]
	predictions.m4[i, j] <- 1 / (1 + exp(-x.m4 %*% betas[i,]))
		}
	}
m4.mean <- apply(X = predictions.m4, MARGIN=2, FUN=mean)
m4.lower <- apply(predictions.m4, 2, quantile, 0.025)
m4.upper <- apply(predictions.m4, 2, quantile, 0.975)


#################
## AKAN 30 & NDC
#################
### for Akans:
m1 <- glm(vote_ndc_pres ~ akan_30km_l_p + int+ dev_factor2 + male  +economy_oneyear + poverty + central + r4  , data=dat_akan, family=binomial(link="logit"))
w.dist<-cbind(dat_akan$vote_ndc_pres, dat_akan$akan_30km_l_p , dat_akan$int, dat_akan$dev_factor2, dat_akan$male, dat_akan$economy_oneyear, dat_akan$poverty,  dat_akan$central, dat_akan$r4, dat_akan$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,ncol(w.dist)])
mat <- estfun(m1)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m1$rank
df <- (length(unique(w.dist[,ncol(w.dist)])) / (length(unique(w.dist[,ncol(w.dist)])) - 1))*((N - 1)/ (N - K))
vcovCL <- df*sandwich(m1, meat=crossprod(u)/N)

a.akan_30km.vals <- seq(quantile(dat_akan$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_akan$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m1<- setx(m1,r4=1,male=1,economy_oneyear=3,central=0,int=0)
x.m1 <- as.matrix(x.m1)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL)

predictions.m1 <- matrix(NA, nrow=1000, ncol=length(a.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(a.akan_30km.vals)){
	x.m1[2] <- a.akan_30km.vals[j]
	predictions.m1[i, j] <- 1 / (1 + exp(-x.m1 %*% betas[i,]))
		}
	}
m1b.mean <- apply(X = predictions.m1, MARGIN=2, FUN=mean)
m1b.lower <- apply(predictions.m1, 2, quantile, 0.025)
m1b.upper <- apply(predictions.m1, 2, quantile, 0.975)


### for Ewes:
m2 <- glm(vote_ndc_pres ~ akan_30km_l_p +int + dev_factor2 + male  +economy_oneyear + poverty + central + r4, data=dat_ewe, family=binomial(link="logit"))
w.dist<-cbind(dat_ewe$vote_ndc_pres, dat_ewe$akan_30km_l_p, dat_ewe$int, dat_ewe$dev_factor2, dat_ewe$male, dat_ewe$economy_oneyear, dat_ewe$poverty, dat_ewe$central, dat_ewe$r4, dat_ewe$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,ncol(w.dist)])
mat <- estfun(m2)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m2$rank
df <- (length(unique(w.dist[,ncol(w.dist)])) / (length(unique(w.dist[,ncol(w.dist)])) - 1))*((N - 1)/ (N - K))
vcovCL2 <- df*sandwich(m2, meat=crossprod(u)/N)
final<-coeftest(m2, vcovCL2)

e.akan_30km.vals <- seq(quantile(dat_ewe$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_ewe$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m2<- setx(m2,r4=1,male=1,economy_oneyear=3,central=0,int=0)
x.m2 <- as.matrix(x.m2)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL2)
predictions.m2 <- matrix(NA, nrow=1000, ncol=length(e.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(e.akan_30km.vals)){
	x.m2[2] <- e.akan_30km.vals[j]
	predictions.m2[i, j] <- 1 / (1 + exp(-x.m2 %*% betas[i,]))
		}
	}
m2b.mean <- apply(X = predictions.m2, MARGIN=2, FUN=mean)
m2b.lower <- apply(predictions.m2, 2, quantile, 0.025)
m2b.upper <- apply(predictions.m2, 2, quantile, 0.975)


### For Mole/Dagbons:
m3 <- glm(vote_ndc_pres ~ akan_30km_l_p  + dev_factor2 + male  +economy_oneyear + poverty + r4, data=dat_mole, family=binomial(link="logit"))
w.dist<-cbind(dat_mole$vote_ndc_pres, dat_mole$akan_30km_l_p, dat_mole$dev_factor2, dat_mole$male, dat_mole$economy_oneyear, dat_mole$poverty, dat_mole$r4, dat_mole$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,8])
mat <- estfun(m3)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m3$rank
df <- (length(unique(w.dist[,8])) / (length(unique(w.dist[,8])) - 1))*((N - 1)/ (N - K))
vcovCL3 <- df*sandwich(m3, meat=crossprod(u)/N)
final<-coeftest(m3, vcovCL3)

m.akan_30km.vals <- seq(quantile(dat_mole$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_mole$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m3<- setx(m3,r4=1,male=1,economy_oneyear=3)
x.m3 <- as.matrix(x.m3)
betas <- mvrnorm(n = 1000, mu=summary(m3)$coefficients[,1], vcovCL3)
predictions.m3 <- matrix(NA, nrow=1000, ncol=length(m.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(m.akan_30km.vals)){
	x.m3[2] <- m.akan_30km.vals[j]
	predictions.m3[i, j] <- 1 / (1 + exp(-x.m3 %*% betas[i,]))
		}
	}
m3b.mean <- apply(X = predictions.m3, MARGIN=2, FUN=mean)
m3b.lower <- apply(predictions.m3, 2, quantile, 0.025)
m3b.upper <- apply(predictions.m3, 2, quantile, 0.975)

#################
#### EVERYONE TOGETHER

m4<-glm(vote_ndc_pres~ akan_30km_l_p + int+ eth_akan + eth_ewe + eth_dagomba + male + economy_oneyear + poverty + dev_factor2 + central + r4, family=binomial(link="logit"), data=data)

### clustered standard errors
w.dist<-cbind(data$vote_ndc_pres, data$akan_30km_l_p ,data$int, data$eth_akan, data$eth_ewe, data$eth_dagomba, data$male, data$economy_oneyear, data$poverty, data$dev_factor2, data$central, data$r4 ,data$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,ncol(w.dist)])
mat <- estfun(m4)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m4$rank
df <- (length(unique(w.dist[,ncol(w.dist)])) / (length(unique(w.dist[,ncol(w.dist)])) - 1))*((N - 1)/ (N - K))
vcovCL4 <- df*sandwich(m4, meat=crossprod(u)/N)
final<-coeftest(m4, vcovCL4)

all.akan_30km.vals <- seq(quantile(data$akan_30km_l_p, probs=c(.1,.9))[1], quantile(data$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
set.seed(02138)
x.m4<- setx(m4,r4=1,male=1,economy_oneyear=3, central=0,int=0)
x.m4 <- as.matrix(x.m4)
betas <- mvrnorm(n = 1000, mu=summary(m4)$coefficients[,1], vcovCL4)
predictions.m4 <- matrix(NA, nrow=1000, ncol=length(all.akan_30km.vals))
for(i in 1:1000){
	for(j in 1:length(all.akan_30km.vals)){
	x.m4[2] <- all.akan_30km.vals[j]
	predictions.m4[i, j] <- 1 / (1 + exp(-x.m4 %*% betas[i,]))
		}
	}
m4b.mean <- apply(X = predictions.m4, MARGIN=2, FUN=mean)
m4b.lower <- apply(predictions.m4, 2, quantile, 0.025)
m4b.upper <- apply(predictions.m4, 2, quantile, 0.975)


#################
## EWE 30 & NPP
#################
### for Akans:
m1 <- glm(vote_npp_pres ~ ewe_30km_l_p + dev_factor2 + male  +economy_oneyear + poverty + central + r4, data=dat_akan, family=binomial(link="logit"))
w.dist<-cbind(dat_akan$vote_npp_pres, dat_akan$ewe_30km_l_p , dat_akan$dev_factor2, dat_akan$male, dat_akan$economy_oneyear, dat_akan$poverty,  dat_akan$central, dat_akan$r4, dat_akan$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,9])
mat <- estfun(m1)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m1$rank
df <- (length(unique(w.dist[,9])) / (length(unique(w.dist[,9])) - 1))*((N - 1)/ (N - K))
vcovCL <- df*sandwich(m1, meat=crossprod(u)/N)

a.ewe_30km.vals <- seq(quantile(dat_akan$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(dat_akan$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m1<- setx(m1,r4=1,male=1,economy_oneyear=3,central=0)
x.m1 <- as.matrix(x.m1)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL)

predictions.m1 <- matrix(NA, nrow=1000, ncol=length(a.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(a.ewe_30km.vals)){
	x.m1[2] <- a.ewe_30km.vals[j]
	predictions.m1[i, j] <- 1 / (1 + exp(-x.m1 %*% betas[i,]))
		}
	}
m1c.mean <- apply(X = predictions.m1, MARGIN=2, FUN=mean)
m1c.lower <- apply(predictions.m1, 2, quantile, 0.025)
m1c.upper <- apply(predictions.m1, 2, quantile, 0.975)


### for Ewes:
m2 <- glm(vote_npp_pres ~ ewe_30km_l_p  + dev_factor2 + male  +economy_oneyear + poverty + central + r4, data=dat_ewe, family=binomial(link="logit"))
w.dist<-cbind(dat_ewe$vote_npp_pres, dat_ewe$ewe_30km_l_p, dat_ewe$dev_factor2, dat_ewe$male, dat_ewe$economy_oneyear, dat_ewe$poverty, dat_ewe$central, dat_ewe$r4, dat_ewe$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,9])
mat <- estfun(m2)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m2$rank
df <- (length(unique(w.dist[,9])) / (length(unique(w.dist[,9])) - 1))*((N - 1)/ (N - K))
vcovCL2 <- df*sandwich(m2, meat=crossprod(u)/N)
final<-coeftest(m2, vcovCL2)

e.ewe_30km.vals <- seq(quantile(dat_ewe$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(dat_ewe$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m2<- setx(m2,r4=1,male=1,economy_oneyear=3,central=0)
x.m2 <- as.matrix(x.m2)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL2)
predictions.m2 <- matrix(NA, nrow=1000, ncol=length(e.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(e.ewe_30km.vals)){
	x.m2[2] <- e.ewe_30km.vals[j]
	predictions.m2[i, j] <- 1 / (1 + exp(-x.m2 %*% betas[i,]))
		}
	}
m2c.mean <- apply(X = predictions.m2, MARGIN=2, FUN=mean)
m2c.lower <- apply(predictions.m2, 2, quantile, 0.025)
m2c.upper <- apply(predictions.m2, 2, quantile, 0.975)

### For Mole/Dagbons:
m3 <- glm(vote_npp_pres ~ ewe_30km_l_p  + dev_factor2 + male  +economy_oneyear + poverty + r4, data=dat_mole, family=binomial(link="logit"))
w.dist<-cbind(dat_mole$vote_npp_pres, dat_mole$ewe_30km_l_p, dat_mole$dev_factor2, dat_mole$male, dat_mole$economy_oneyear, dat_mole$poverty, dat_mole$r4, dat_mole$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,8])
mat <- estfun(m3)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m3$rank
df <- (length(unique(w.dist[,8])) / (length(unique(w.dist[,8])) - 1))*((N - 1)/ (N - K))
vcovCL3 <- df*sandwich(m3, meat=crossprod(u)/N)
final<-coeftest(m3, vcovCL3)

m.ewe_30km.vals <- seq(quantile(dat_mole$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(dat_mole$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m3<- setx(m3,r4=1,male=1,economy_oneyear=3)
x.m3 <- as.matrix(x.m3)
betas <- mvrnorm(n = 1000, mu=summary(m3)$coefficients[,1], vcovCL3)
predictions.m3 <- matrix(NA, nrow=1000, ncol=length(m.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(m.ewe_30km.vals)){
	x.m3[2] <- m.ewe_30km.vals[j]
	predictions.m3[i, j] <- 1 / (1 + exp(-x.m3 %*% betas[i,]))
		}
	}
m3c.mean <- apply(X = predictions.m3, MARGIN=2, FUN=mean)
m3c.lower <- apply(predictions.m3, 2, quantile, 0.025)
m3c.upper <- apply(predictions.m3, 2, quantile, 0.975)


#### EVERYONE TOGETHER

m4<-glm(vote_npp_pres~ ewe_30km_l_p + eth_akan + eth_ewe + eth_dagomba + male + economy_oneyear + poverty + dev_factor2 + central + r4, family=binomial(link="logit"), data=data)

### clustered standard errors
w.dist<-cbind(data$vote_npp_pres, data$ewe_30km_l_p , data$eth_akan, data$eth_ewe, data$eth_dagomba, data$male, data$economy_oneyear, data$poverty, data$dev_factor2, data$central, data$r4 ,data$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,12])
mat <- estfun(m4)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m4$rank
df <- (length(unique(w.dist[,12])) / (length(unique(w.dist[,12])) - 1))*((N - 1)/ (N - K))
vcovCL4 <- df*sandwich(m4, meat=crossprod(u)/N)
final<-coeftest(m4, vcovCL4)

all.ewe_30km.vals <- seq(quantile(data$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(data$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
set.seed(02138)
x.m4<- setx(m4,r4=1,male=1,economy_oneyear=3, central=0)
x.m4 <- as.matrix(x.m4)
betas <- mvrnorm(n = 1000, mu=summary(m4)$coefficients[,1], vcovCL4)
predictions.m4 <- matrix(NA, nrow=1000, ncol=length(all.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(all.ewe_30km.vals)){
	x.m4[2] <- all.ewe_30km.vals[j]
	predictions.m4[i, j] <- 1 / (1 + exp(-x.m4 %*% betas[i,]))
		}
	}
m4c.mean <- apply(X = predictions.m4, MARGIN=2, FUN=mean)
m4c.lower <- apply(predictions.m4, 2, quantile, 0.025)
m4c.upper <- apply(predictions.m4, 2, quantile, 0.975)



#################
## EWE 30 & NDC
#################
### for Akans:
m1 <- glm(vote_ndc_pres ~ ewe_30km_l_p + dev_factor2 + male  +economy_oneyear + poverty + central + r4  , data=dat_akan, family=binomial(link="logit"))
w.dist<-cbind(dat_akan$vote_ndc_pres, dat_akan$ewe_30km_l_p , dat_akan$dev_factor2, dat_akan$male, dat_akan$economy_oneyear, dat_akan$poverty,  dat_akan$central, dat_akan$r4, dat_akan$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,9])
mat <- estfun(m1)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m1$rank
df <- (length(unique(w.dist[,9])) / (length(unique(w.dist[,9])) - 1))*((N - 1)/ (N - K))
vcovCL <- df*sandwich(m1, meat=crossprod(u)/N)

a.ewe_30km.vals <- seq(quantile(dat_akan$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(dat_akan$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m1<- setx(m1,r4=1,male=1,economy_oneyear=3,central=0)
x.m1 <- as.matrix(x.m1)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL)

predictions.m1 <- matrix(NA, nrow=1000, ncol=length(a.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(a.ewe_30km.vals)){
	x.m1[2] <- a.ewe_30km.vals[j]
	predictions.m1[i, j] <- 1 / (1 + exp(-x.m1 %*% betas[i,]))
		}
	}
m1d.mean <- apply(X = predictions.m1, MARGIN=2, FUN=mean)
m1d.lower <- apply(predictions.m1, 2, quantile, 0.025)
m1d.upper <- apply(predictions.m1, 2, quantile, 0.975)


### for Ewes:
m2 <- glm(vote_ndc_pres ~ ewe_30km_l_p  + dev_factor2 + male  +economy_oneyear + poverty + central + r4  , data=dat_ewe, family=binomial(link="logit"))
#summary(m2)
w.dist<-cbind(dat_ewe$vote_ndc_pres, dat_ewe$ewe_30km_l_p, dat_ewe$dev_factor2, dat_ewe$male, dat_ewe$economy_oneyear, dat_ewe$poverty, dat_ewe$central, dat_ewe$r4, dat_ewe$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,9])
mat <- estfun(m2)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m2$rank
df <- (length(unique(w.dist[,9])) / (length(unique(w.dist[,9])) - 1))*((N - 1)/ (N - K))
vcovCL2 <- df*sandwich(m2, meat=crossprod(u)/N)
final<-coeftest(m2, vcovCL2)
#final

e.ewe_30km.vals <- seq(quantile(dat_ewe$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(dat_ewe$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m2<- setx(m2,r4=1,male=1,economy_oneyear=3,central=0)
x.m2 <- as.matrix(x.m2)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL2)
predictions.m2 <- matrix(NA, nrow=1000, ncol=length(e.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(e.ewe_30km.vals)){
	x.m2[2] <- e.ewe_30km.vals[j]
	predictions.m2[i, j] <- 1 / (1 + exp(-x.m2 %*% betas[i,]))
		}
	}
m2d.mean <- apply(X = predictions.m2, MARGIN=2, FUN=mean)
m2d.lower <- apply(predictions.m2, 2, quantile, 0.025)
m2d.upper <- apply(predictions.m2, 2, quantile, 0.975)


### For Mole/Dagbons:
m3 <- glm(vote_ndc_pres ~ ewe_30km_l_p  + dev_factor2 + male  +economy_oneyear + poverty + r4  , data=dat_mole, family=binomial(link="logit"))
w.dist<-cbind(dat_mole$vote_ndc_pres, dat_mole$ewe_30km_l_p, dat_mole$dev_factor2, dat_mole$male, dat_mole$economy_oneyear, dat_mole$poverty, dat_mole$r4, dat_mole$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,8])
mat <- estfun(m3)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m3$rank
df <- (length(unique(w.dist[,8])) / (length(unique(w.dist[,8])) - 1))*((N - 1)/ (N - K))
vcovCL3 <- df*sandwich(m3, meat=crossprod(u)/N)
final<-coeftest(m3, vcovCL3)

m.ewe_30km.vals <- seq(quantile(dat_mole$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(dat_mole$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m3<- setx(m3,r4=1,male=1,economy_oneyear=3)
x.m3 <- as.matrix(x.m3)
betas <- mvrnorm(n = 1000, mu=summary(m3)$coefficients[,1], vcovCL3)
predictions.m3 <- matrix(NA, nrow=1000, ncol=length(m.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(m.ewe_30km.vals)){
	x.m3[2] <- m.ewe_30km.vals[j]
	predictions.m3[i, j] <- 1 / (1 + exp(-x.m3 %*% betas[i,]))
		}
	}
m3d.mean <- apply(X = predictions.m3, MARGIN=2, FUN=mean)
m3d.lower <- apply(predictions.m3, 2, quantile, 0.025)
m3d.upper <- apply(predictions.m3, 2, quantile, 0.975)

#### EVERYONE TOGETHER

m4<-glm(vote_ndc_pres~ ewe_30km_l_p + eth_akan + eth_ewe + eth_dagomba + male + economy_oneyear + poverty + dev_factor2 + central + r4, family=binomial(link="logit"), data=data)

### clustered standard errors
w.dist<-cbind(data$vote_ndc_pres, data$ewe_30km_l_p , data$eth_akan, data$eth_ewe, data$eth_dagomba, data$male, data$economy_oneyear, data$poverty, data$dev_factor2, data$central, data$r4 ,data$fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,12])
mat <- estfun(m4)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m4$rank
df <- (length(unique(w.dist[,12])) / (length(unique(w.dist[,12])) - 1))*((N - 1)/ (N - K))
vcovCL4 <- df*sandwich(m4, meat=crossprod(u)/N)
final<-coeftest(m4, vcovCL4)

all.ewe_30km.vals <- seq(quantile(data$ewe_30km_l_p, probs=c(.1,.9))[1], quantile(data$ewe_30km_l_p, probs=c(.1,.9))[2], 0.001)
set.seed(02138)
x.m4<- setx(m4,r4=1,male=1,economy_oneyear=3, central=0)
x.m4 <- as.matrix(x.m4)
betas <- mvrnorm(n = 1000, mu=summary(m4)$coefficients[,1], vcovCL4)
predictions.m4 <- matrix(NA, nrow=1000, ncol=length(all.ewe_30km.vals))
for(i in 1:1000){
	for(j in 1:length(all.ewe_30km.vals)){
	x.m4[2] <- all.ewe_30km.vals[j]
	predictions.m4[i, j] <- 1 / (1 + exp(-x.m4 %*% betas[i,]))
		}
	}
m4d.mean <- apply(X = predictions.m4, MARGIN=2, FUN=mean)
m4d.lower <- apply(predictions.m4, 2, quantile, 0.025)
m4d.upper <- apply(predictions.m4, 2, quantile, 0.975)


#### FULL PLOT

pdf(file = "Figure4.pdf", height=10, width=10, family="Helvetica")
par(mfrow=c(2,2))

plot(x = 100*a.akan_30km.vals, y= m1.mean, col="dodgerblue4", pch=16, cex=.4, ylab="Predicted Probability of NPP Support", xlab="% Akan in 30km (spatially weighted)", main="(a)", ylim=c(0, 1), xlim=c(0,100), type="l", lwd=2.5, lty="dotdash", cex.lab=1.4)
points(100*a.akan_30km.vals, m1.lower, pch=16, cex=.1, col="gray47")
points(100*a.akan_30km.vals, m1.upper, pch=16, cex=.1, col="gray47")
points(100*e.akan_30km.vals, m2.mean, col="goldenrod", pch=16, cex=.4, type="l", lty="dashed", lwd=2.5)
points(100*e.akan_30km.vals, m2.lower, pch=16, cex=.1, col="gray47")
points(100*e.akan_30km.vals, m2.upper, pch=16, cex=.1, col="gray47")
points(100*m.akan_30km.vals, m3.mean, col="firebrick", pch=16, cex=.4, type="l", lty="longdash", lwd=2.5)
points(100*m.akan_30km.vals, m3.lower, pch=16, cex=.1, col="gray47")
points(100*m.akan_30km.vals, m3.upper, pch=16, cex=.1, col="gray47")
points(100*all.akan_30km.vals, m4.mean, col="darkgreen", pch=16, cex=.4, type="l", lwd=2.5)
points(100*all.akan_30km.vals, m4.lower, pch=16, cex=.1, col="gray47")
points(100*all.akan_30km.vals, m4.upper, pch=16, cex=.1, col="gray47")

plot(x = 100*a.ewe_30km.vals, y= m1c.mean, col="dodgerblue4", pch=16, cex=.4, ylab="Predicted Probability of NPP Support", xlab="% Ewe in 30km (spatially weighted)", main="(b)", ylim=c(0, 1), xlim=c(0,100), type="l", lwd=2.5, lty="dotdash", cex.lab=1.4)
points(100*a.ewe_30km.vals, m1c.lower, pch=16, cex=.1, col="gray47")
points(100*a.ewe_30km.vals, m1c.upper, pch=16, cex=.1, col="gray47")
points(100*e.ewe_30km.vals, m2c.mean, col="goldenrod", pch=16, cex=.4, type="l", lty="dashed", lwd=2.5)
points(100*e.ewe_30km.vals, m2c.lower, pch=16, cex=.1, col="gray47")
points(100*e.ewe_30km.vals, m2c.upper, pch=16, cex=.1, col="gray47")
points(100*all.ewe_30km.vals, m4c.mean, col="darkgreen", pch=16, cex=.4, type="l", lwd=2.5)
points(100*all.ewe_30km.vals, m4c.lower, pch=16, cex=.1, col="gray47")
points(100*all.ewe_30km.vals, m4c.upper, pch=16, cex=.1, col="gray47")

plot(x = 100*a.akan_30km.vals, y= m1b.mean, col="dodgerblue4", pch=16, cex=.4, ylab="Predicted Probability of NDC Support", xlab="% Akan in 30km (spatially weighted)", main="(c)", ylim=c(0, 1), xlim=c(0,100), type="l", lwd=2.5, lty="dotdash", cex.lab=1.4)
points(100*a.akan_30km.vals, m1b.lower, pch=16, cex=.1, col="gray47")
points(100*a.akan_30km.vals, m1b.upper, pch=16, cex=.1, col="gray47")
points(100*e.akan_30km.vals,  m2b.mean, col="goldenrod", pch=16, cex=.4, type="l", lty="dashed", lwd=2.5)
points(100*e.akan_30km.vals, m2b.lower, pch=16, cex=.1, col="gray47")
points(100*e.akan_30km.vals, m2b.upper, pch=16, cex=.1, col="gray47")
points(100*m.akan_30km.vals, m3b.mean, col="firebrick", pch=16, cex=.4, type="l", lwd=2.5, lty="longdash")
points(100*m.akan_30km.vals, m3b.lower, pch=16, cex=.1, col="gray47")
points(100*m.akan_30km.vals, m3b.upper, pch=16, cex=.1, col="gray47")
points(100*all.akan_30km.vals,  m4b.mean, col="darkgreen", pch=16, cex=.4, type="l", lwd=2.5)
points(100*all.akan_30km.vals, m4b.lower, pch=16, cex=.1, col="gray47")
points(100*all.akan_30km.vals, m4b.upper, pch=16, cex=.1, col="gray47")

plot(x = 100*a.ewe_30km.vals, y= m1d.mean, col="dodgerblue4", pch=16, cex=.4, ylab="Predicted Probability of NDC Support", xlab="% Ewe in 30km (spatially weighted)", main="(d)", ylim=c(0, 1), xlim=c(0,100), type="l", lwd=2.5, lty="dotdash", cex.lab=1.4)
points(100*a.ewe_30km.vals, m1d.lower, pch=16, cex=.1, col="gray47")
points(100*a.ewe_30km.vals, m1d.upper, pch=16, cex=.1, col="gray47")
points(100*e.ewe_30km.vals, m2d.mean, col="goldenrod", pch=16, cex=.4, type="l", lty="dashed", lwd=2.5)
points(100*e.ewe_30km.vals, m2d.lower, pch=16, cex=.1, col="gray47")
points(100*e.ewe_30km.vals, m2d.upper, pch=16, cex=.1, col="gray47")
points(100*all.ewe_30km.vals, m4d.mean, col="darkgreen", pch=16, cex=.4, type="l", lwd=2.5)
points(100*all.ewe_30km.vals, m4d.lower, pch=16, cex=.1, col="gray47")
points(100*all.ewe_30km.vals, m4d.upper, pch=16, cex=.1, col="gray47")


dev.off()

##############################################
## FIGURE 5
##############################################
### Data from both rounds


dat1<-read.csv("IchinoNathan_APSR_survey_data.csv")
dat1<-dat1[,c("fid0725","r4","central","urban", "male", "eth_akan", "eth_ewe", "eth_ga", "eth_dagomba", "eth_other", "vote_npp_pres", "vote_ndc_pres", "economy_oneyear", "poverty",  "gov_sentus",   "akan_30km_l_p", "ewe_30km_l_p", "dev_factor2", "popdens5x5")]
dat1 <- na.omit(dat1)
dim(dat1)

dat_akan <- dat1[dat1$eth_akan==1,]
dim(dat_akan)
dat_ewe <- dat1[dat1$eth_ewe==1,]
dim(dat_ewe)
dat_mole <- dat1[dat1$eth_dagomba==1,]
dim(dat_mole)
dat_nonakan <- dat1[dat1$eth_akan==0,]
dim(dat_nonakan)

dat_akan$int<-dat_akan$central*dat_akan$akan_30km_l_p
dat_ewe$int<-dat_ewe$central*dat_ewe$akan_30km_l_p
dat_nonakan$int<-dat_nonakan$central*dat_nonakan$akan_30km_l_p
dat_mole$int<-dat_mole$central*dat_mole$akan_30km_l_p
dat1$int<-dat1$central*dat1$akan_30km_l_p



set.seed(02138)

m1 <- glm(vote_npp_pres ~ akan_30km_l_p + central +gov_sentus + I(gov_sentus*central)+ I(gov_sentus*akan_30km_l_p) + int + I(gov_sentus*int)+male  +economy_oneyear +poverty + dev_factor2, data=dat_nonakan, family=binomial(link="logit"))
w.dist<-cbind(dat_nonakan$vote_npp_pres, dat_nonakan$akan_30km_l_p, dat_nonakan $central, dat_nonakan $gov_sentus, dat_nonakan $int, dat_nonakan $male, dat_nonakan $economy_oneyear, dat_nonakan $poverty, dat_nonakan $dev_factor2,  dat_nonakan $fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,dim(w.dist)[2]])
mat <- estfun(m1)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m1$rank
df <- (length(unique(w.dist[,dim(w.dist)[2]])) / (length(unique(w.dist[,dim(w.dist)[2]])) - 1))*((N - 1)/ (N - K))
vcovCL1 <- df*sandwich(m1, meat=crossprod(u)/N)
final1<-coeftest(m1, vcovCL1)

akan30.vals <- seq(quantile(dat_nonakan$akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_nonakan$akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m1<- setx(m1,male=1,central=0,gov_sentus=0,int=0,economy_oneyear=3)
x.m1 <- as.matrix(x.m1)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL1)

predictions.m1 <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m1[2] <- akan30.vals[j]
	predictions.m1[i, j] <- 1 / (1 + exp(-x.m1 %*% betas[i,]))
		}
	}
m1.mean1 <- apply(X = predictions.m1, MARGIN=2, FUN=mean)
m1.lower1 <- apply(predictions.m1, 2, quantile, 0.025)
m1.upper1 <- apply(predictions.m1, 2, quantile, 0.975)

x.m1g<- setx(m1,male=1,central=0,gov_sentus=1,int=0,economy_oneyear=3)
x.m1g <- as.matrix(x.m1g)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL1)

predictions.m1g <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m1g[2] <- akan30.vals[j]
	predictions.m1g[i, j] <- 1 / (1 + exp(-x.m1g %*% betas[i,]))
		}
	}
m1g.mean1 <- apply(X = predictions.m1g, MARGIN=2, FUN=mean)
m1g.lower1 <- apply(predictions.m1g, 2, quantile, 0.025)
m1g.upper1 <- apply(predictions.m1g, 2, quantile, 0.975)


m2 <- glm(vote_ndc_pres ~ akan_30km_l_p + central +gov_sentus + I(gov_sentus*central)+ I(gov_sentus*akan_30km_l_p) + int + I(gov_sentus*int)+male  +economy_oneyear +poverty + dev_factor2, data=dat_nonakan, family=binomial(link="logit"))
w.dist<-cbind(dat_nonakan $vote_npp_pres, dat_nonakan $akan_30km_l_p, dat_nonakan $central, dat_nonakan $gov_sentus, dat_nonakan $int, dat_nonakan $male, dat_nonakan $economy_oneyear, dat_nonakan $poverty, dat_nonakan $dev_factor2,  dat_nonakan $fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,dim(w.dist)[2]])
mat <- estfun(m2)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m2$rank
df <- (length(unique(w.dist[,dim(w.dist)[2]])) / (length(unique(w.dist[,dim(w.dist)[2]])) - 1))*((N - 1)/ (N - K))
vcovCL2 <- df*sandwich(m2, meat=crossprod(u)/N)
final2<-coeftest(m2, vcovCL2)

akan30.vals <- seq(quantile(dat_nonakan $akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_nonakan $akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m2<- setx(m2,male=1,central=0,gov_sentus=0,int=0,economy_oneyear=3)
x.m2 <- as.matrix(x.m2)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL1)

predictions.m2 <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m2[2] <- akan30.vals[j]
	predictions.m2[i, j] <- 1 / (1 + exp(-x.m2 %*% betas[i,]))
		}
	}
m2.mean1 <- apply(X = predictions.m2, MARGIN=2, FUN=mean)
m2.lower1 <- apply(predictions.m2, 2, quantile, 0.025)
m2.upper1 <- apply(predictions.m2, 2, quantile, 0.975)

x.m2g<- setx(m2,male=1,central=0,gov_sentus=1,int=0,economy_oneyear=3)
x.m2g <- as.matrix(x.m2g)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL1)

predictions.m2g <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m2g[2] <- akan30.vals[j]
	predictions.m2g[i, j] <- 1 / (1 + exp(-x.m2g %*% betas[i,]))
		}
	}
m2g.mean1 <- apply(X = predictions.m2g, MARGIN=2, FUN=mean)
m2g.lower1 <- apply(predictions.m2g, 2, quantile, 0.025)
m2g.upper1 <- apply(predictions.m2g, 2, quantile, 0.975)

akan30.vals1 <- akan30.vals



## Round 3 only

dat1<-read.csv("IchinoNathan_APSR_survey_data.csv")
dat1<-dat1[,c("fid0725","r4","central","urban", "male", "eth_akan", "eth_ewe", "eth_ga", "eth_dagomba", "eth_other", "vote_npp_pres", "vote_ndc_pres", "economy_oneyear", "poverty", "trustother",    "akan_30km_l_p", "ewe_30km_l_p", "dev_factor2", "popdens5x5")]
dat1 <- na.omit(dat1)
dim(dat1)

dat_akan <- dat1[dat1$eth_akan==1,]
dim(dat_akan)
dat_ewe <- dat1[dat1$eth_ewe==1,]
dim(dat_ewe)
dat_mole <- dat1[dat1$eth_dagomba==1,]
dim(dat_mole)
dat_nonakan <- dat1[dat1$eth_akan==0,]
dim(dat_nonakan)
dat_nonewe <- dat1[dat1$eth_ewe==0,]
dim(dat_nonewe)

dat_akan$int<-dat_akan$central*dat_akan$akan_30km_l_p
dat_ewe$int<-dat_ewe$central*dat_ewe$akan_30km_l_p
dat_nonakan$int<-dat_nonakan$central*dat_nonakan$akan_30km_l_p
dat_mole$int<-dat_mole$central*dat_mole$akan_30km_l_p
dat1$int<-dat1$central*dat1$akan_30km_l_p

set.seed(02138)
### Akan 30, interacted with othertrust, for non-Akans 
### support NPP as outcome

m1 <- glm(vote_npp_pres ~ akan_30km_l_p +trustother + I(trustother*akan_30km_l_p) +male  +economy_oneyear +poverty + dev_factor2, data=dat_nonakan, family=binomial(link="logit"))
w.dist<-cbind(dat_nonakan $vote_npp_pres, dat_nonakan $akan_30km_l_p, dat_nonakan$trustother, dat_nonakan$trustother*dat_nonakan$akan_30km_l_p, dat_nonakan $male, dat_nonakan $economy_oneyear, dat_nonakan $poverty, dat_nonakan $dev_factor2,  dat_nonakan $fid0725)
w.dist<-na.omit(w.dist)
dist.fan<-factor(w.dist[,dim(w.dist)[2]])
mat <- estfun(m1)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m1$rank
df <- (length(unique(w.dist[,dim(w.dist)[2]])) / (length(unique(w.dist[,dim(w.dist)[2]])) - 1))*((N - 1)/ (N - K))
vcovCL1 <- df*sandwich(m1, meat=crossprod(u)/N)
final1<-coeftest(m1, vcovCL1)

akan30.vals <- seq(quantile(dat_nonakan $akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_nonakan $akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m1<- setx(m1,male=1,central=0,trustother=0,int=0,economy_oneyear=3)
x.m1 <- as.matrix(x.m1)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL1)

predictions.m1 <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m1[2] <- akan30.vals[j]
	predictions.m1[i, j] <- 1 / (1 + exp(-x.m1 %*% betas[i,]))
		}
	}
m1.mean <- apply(X = predictions.m1, MARGIN=2, FUN=mean)
m1.lower <- apply(predictions.m1, 2, quantile, 0.025)
m1.upper <- apply(predictions.m1, 2, quantile, 0.975)



x.m1g<- setx(m1,male=1,central=0,trustother=1,int=0,economy_oneyear=3)
x.m1g <- as.matrix(x.m1g)
betas <- mvrnorm(n = 1000, mu=summary(m1)$coefficients[,1], vcovCL1)

predictions.m1g <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m1g[2] <- akan30.vals[j]
	predictions.m1g[i, j] <- 1 / (1 + exp(-x.m1g %*% betas[i,]))
		}
	}
m1g.mean <- apply(X = predictions.m1g, MARGIN=2, FUN=mean)
m1g.lower <- apply(predictions.m1g, 2, quantile, 0.025)
m1g.upper <- apply(predictions.m1g, 2, quantile, 0.975)



### Akan 30, interacted with othertrust, for non-Akans 
### support NDC as outcome

m2 <- glm(vote_ndc_pres ~ akan_30km_l_p +trustother + I(trustother*akan_30km_l_p) +male  +economy_oneyear +poverty + dev_factor2, data=dat_nonakan, family=binomial(link="logit"))
w.dist<-cbind(dat_nonakan $vote_ndc_pres, dat_nonakan $akan_30km_l_p, dat_nonakan$trustother, dat_nonakan$trustother*dat_nonakan$akan_30km_l_p, dat_nonakan $male, dat_nonakan $economy_oneyear, dat_nonakan $poverty, dat_nonakan $dev_factor2,  dat_nonakan $fid0725)
dist.fan<-factor(w.dist[,dim(w.dist)[2]])
mat <- estfun(m2)
mat <- na.omit(mat)
N <- nrow(mat)
u <- apply(mat, 2, function(x) tapply(x, dist.fan, sum))
u <- na.omit(u)
K<-m2$rank
df <- (length(unique(w.dist[,dim(w.dist)[2]])) / (length(unique(w.dist[,dim(w.dist)[2]])) - 1))*((N - 1)/ (N - K))
vcovCL1 <- df*sandwich(m2, meat=crossprod(u)/N)
final1<-coeftest(m2, vcovCL1)

akan30.vals <- seq(quantile(dat_nonakan $akan_30km_l_p, probs=c(.1,.9))[1], quantile(dat_nonakan $akan_30km_l_p, probs=c(.1,.9))[2], 0.001)
x.m2<- setx(m2,male=1,central=0,trustother=0,int=0,economy_oneyear=3)
x.m2 <- as.matrix(x.m2)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL1)

predictions.m2 <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m2[2] <- akan30.vals[j]
	predictions.m2[i, j] <- 1 / (1 + exp(-x.m2 %*% betas[i,]))
		}
	}
m2.mean <- apply(X = predictions.m2, MARGIN=2, FUN=mean)
m2.lower <- apply(predictions.m2, 2, quantile, 0.025)
m2.upper <- apply(predictions.m2, 2, quantile, 0.975)



x.m2g<- setx(m2,male=1,central=0,trustother=1,int=0,economy_oneyear=3)
x.m2g <- as.matrix(x.m2g)
betas <- mvrnorm(n = 1000, mu=summary(m2)$coefficients[,1], vcovCL1)

predictions.m2g <- matrix(NA, nrow=1000, ncol=length(akan30.vals))
for(i in 1:1000){
	for(j in 1:length(akan30.vals)){
	x.m2g[2] <- akan30.vals[j]
	predictions.m2g[i, j] <- 1 / (1 + exp(-x.m2g %*% betas[i,]))
		}
	}
m2g.mean <- apply(X = predictions.m2g, MARGIN=2, FUN=mean)
m2g.lower <- apply(predictions.m2g, 2, quantile, 0.025)
m2g.upper <- apply(predictions.m2g, 2, quantile, 0.975)


##### PLOT
pdf(file = "Figure5.pdf", height=10, width=10, family="Helvetica")
#pdf(file = "../Drafts/apsr_rev2/figs/AB_govsentus_trust.pdf", height=10, width=10, family="Helvetica")
par(mfrow=c(2,2))
 
plot(x = 100*akan30.vals1, y= m1.mean1, col="darkgreen", pch=16, cex=.4, cex.lab=1.3, ylab="Pred. Probability of Support for NPP", xlab="(a) % Akan in 30km (spatially weighted)", main="(a) Sent by Government", ylim=c(0, 1), xlim=c(0,80),lty="dashed", type="l", lwd=1.7)
points(100*akan30.vals1, m1.lower1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")
points(100*akan30.vals1, m1.upper1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")

points(100*akan30.vals1, m1g.mean1, col="firebrick", pch=16, type="l", lwd=1.5)
points(100*akan30.vals1, m1g.lower1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")
points(100*akan30.vals1, m1g.upper1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")

plot(x = 100*akan30.vals1, y= m2.mean1, col="darkgreen", pch=16, cex=.4, cex.lab=1.3, ylab="Pred. Probability of Support for NDC", xlab="(b) % Akan in 30km (spatially weighted)", main="(b) Sent by Government", ylim=c(0, 1), xlim=c(0,80),lty="dashed", type="l", lwd=1.7)
points(100*akan30.vals1, m2.lower1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")
points(100*akan30.vals1, m2.upper1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")

points(100*akan30.vals1, m2g.mean1, col="firebrick", pch=16, type="l", lwd=1.5)
points(100*akan30.vals1, m2g.lower1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")
points(100*akan30.vals1, m2g.upper1, pch=16, cex=.3, col="gray47", type="l", lty="dotted")


plot(x = 100*akan30.vals, y= m1.mean, col="darkgreen", pch=16, cex=.4, cex.lab=1.3, ylab="Pred. Probability of Support for NPP", xlab="% Akan in 30km (spatially weighted)", main="(c) Trust Other Groups", ylim=c(0, 1), xlim=c(0,80), type="l", lty="dashed", lwd=1.7)
points(100*akan30.vals, m1.lower, pch=16, cex=.3, col="gray47", type="l", lty="dashed")
points(100*akan30.vals, m1.upper, pch=16, cex=.3, col="gray47", type="l", lty="dashed")
#rug(100*dat_nonakan$akan_30km_l_p)


points(100*akan30.vals, m1g.mean, col="firebrick", pch=16, cex=.4, type="l", lwd=1.7)
points(100*akan30.vals, m1g.lower, pch=16, cex=.3, col="gray47", type="l", lty="dashed")
points(100*akan30.vals, m1g.upper, pch=16, cex=.3, col="gray47", type="l", lty="dashed")

plot(x = 100*akan30.vals, y= m2.mean, col="darkgreen", pch=16, cex=.4, cex.lab=1.3, ylab="Pred. Probability of Support for NDC", xlab="% Akan in 30km (spatially weighted)", main="(d) Trust Other Groups", ylim=c(0, 1), xlim=c(0,80), type="l", lty="dashed", lwd=1.7)
points(100*akan30.vals, m2.lower, pch=16, cex=.3, col="gray47", type="l", lty="dashed")
points(100*akan30.vals, m2.upper, pch=16, cex=.3, col="gray47", type="l", lty="dashed")
#rug(100*dat_nonakan$akan_30km_l_p)

points(100*akan30.vals, m2g.mean, col="firebrick", pch=16, cex=.4, type="l", lty="dashed", lwd=1.7)
points(100*akan30.vals, m2g.lower, pch=16, cex=.3, col="gray47", type="l", lty="dashed")
points(100*akan30.vals, m2g.upper, pch=16, cex=.3, col="gray47", type="l", lty="dashed")


dev.off()
