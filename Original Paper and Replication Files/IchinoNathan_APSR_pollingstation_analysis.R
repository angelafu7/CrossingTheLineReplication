############################################################################################
###  Part 1 of 3: 
###  Replication Code for: 
###  Ichino, Nahomi, and Noah L. Nathan.  2013.  "Crossing the Line: Local Ethnic Geography 
###  and Voting in Ghana," American Political Science Review 107(2): 344-61.
###  http://dx.doi.org/10.1017/S0003055412000664
###
###  IchinoNathan_APSR_pollingstation_analysis.R
###
###  Nahomi Ichino and Noah Nathan
###  Department of Government, Harvard University
###  February 2013
###  
###  Note: 
###  This reproduces Tables 1 and 2 of the published article.
###  This also reproduces Tables 1-3 and Figures 1, 3-6 of the Supplementary Materials 
###
##############
###  Calls: 
###  IchinoNathan_APSR_pollingstation_data.csv
### 
###  Data available at:
###  http://dvn.iq.harvard.edu/dvn/dv/nichino
###  http://hdl.handle.net/1902.1/21461  
############################################################################################


rm(list=ls())
library(foreign)
library(car)
library(arm)
library(xtable)
library(mgcv)



data <- read.csv("IchinoNathan_APSR_pollingstation_data.csv", header=TRUE)
dim(data)

hdat <- data[,c("npp2008ps_pres_p",  "ndc2008ps_pres_p", "akan_p_poly", "ewe_p_poly", "moledagbon_p_poly", "public_semipublic_p", "dev_factor2", "akan_30km_l_p", "akan_20km_l_p", "akan_40km_l_p", "ewe_30km_l_p", "ewe_20km_l_p", "ewe_40km_l_p", "h_30rad_e", "c230_id_h", "otherethn_p_poly", "totalvotes2008ps_pres", "ethfrac_30km_l", "akan_p")]
dim(hdat)
hdat <- na.omit(hdat)
dim(hdat)


#### TABLE 1 of Article ####
data1 <- data[,c("area_sqkm", "sex_total", "akan_p_poly", "gadangbe_p_poly", "ewe_p_poly", "guan_p_poly", "gurma_p_poly", "moledagbon_p_poly", "grusi_p_poly", "mande_p_poly", "other_eth_p_poly", "ethfrac_poly", "english_percent", "public_semipublic_p", "dev_factor2", "akan_30km_l_p", "h_30rad_e")]
table1 <- matrix(NA, nrow=ncol(data1), ncol=4)
for(i in 1:ncol(data1)){
table1[i,1] <- mean(data1[,i])
table1[i,2] <- sd(data1[,i])
table1[i,3] <-min(data1[,i])
table1[i,4] <-max(data1[,i])	
}
colnames(table1)<-c("Mean", "SD", "Min", "Max")
rownames(table1)<-c("area_sqkm", "sex_total", "akan_p_poly", "gadangbe_p_poly", "ewe_p_poly", "guan_p_poly", "gurma_p_poly", "moledagbon_p_poly", "grusi_p_poly", "mande_p_poly", "other_eth_p_poly", "ethfrac_poly", "english_percent", "public_semipublic_p", "dev_factor2", "akan_30km_l_p", "h_30rad_e")
# print result
round(table1, digits=2)



#### TABLE 2 of Article ####
m1 <- lm(npp2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + factor(c230_id_h), data=hdat, weights=totalvotes2008ps_pres)

m2 <- lm(npp2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_30km_l_p + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

m3 <- lm(npp2008ps_pres_p ~ akan_p_poly  + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_30km_l_p + I(akan_30km_l_p*akan_p_poly) + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

# print result
summary(m1)
summary(m2)
summary(m3)



# # # # # # # # # # # # # # # # #
#### SUPPLEMENTAL MATERIALS ####

#### Table 1 ####
m1 <- lm(npp2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + factor(c230_id_h), data=hdat, weights=totalvotes2008ps_pres)

m2 <- lm(npp2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_20km_l_p + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

m3 <- lm(npp2008ps_pres_p ~ akan_p_poly  + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_20km_l_p + I(akan_20km_l_p*akan_p_poly) + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

m4 <- lm(npp2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_40km_l_p + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

m5 <- lm(npp2008ps_pres_p ~ akan_p_poly  + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_40km_l_p + I(akan_40km_l_p*akan_p_poly) + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

# print result
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)


#### Table 2 ####
m1 <- lm(ndc2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + factor(c230_id_h), data=hdat, weights=totalvotes2008ps_pres)

m2 <- lm(ndc2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_30km_l_p + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

m3 <- lm(ndc2008ps_pres_p ~ akan_p_poly  + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_30km_l_p + I(akan_30km_l_p*akan_p_poly) + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

# print result
summary(m1)
summary(m2)
summary(m3)


#### Table 3 ####
m2a <- lm(npp2008ps_pres_p ~ akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + akan_30km_l_p + ethfrac_30km_l + factor(c230_id_h), data=hdat, weights= totalvotes2008ps_pres)

# print result
summary(m2a)

#### Figure 1 ####
pdf(file="SupM_fig1.pdf", family="Helvetica", height=6, width=12)
par(mfrow=c(1,2))
plot(100*hdat$akan_p_poly, 100*hdat$akan_30km_l_p, col="dodgerblue4", pch=1, cex=.8, cex.lab=1.3, xlim=c(0,100), ylim=c(0,100), xlab="% Akan at Polling Station", ylab="% Akan in 30km (Spatially Weighted)")
plot(100*hdat$akan_30km_l_p, 100*hdat$ewe_30km_l_p, col="dodgerblue4", cex=.8,  cex.lab=1.3,   xlim=c(0,100), ylim=c(0,100), xlab="% Akan in 30km (Spatially Weighted)", ylab="% Ewe in 30km (Spatially Weighted)")
dev.off()

#### Figure 3 #### 
gam3<-gam(npp2008ps_pres_p ~ s(akan_30km_l_p) + s(akan_p_poly) + s(moledagbon_p_poly) + s(otherethn_p_poly) + public_semipublic_p + dev_factor2 + factor(c230_id_h), data=hdat, family=gaussian)

pdf(file= "SupM_fig3.pdf", width=6, height=6, family="Helvetica")
plot(gam3,se=TRUE,rug=TRUE, select=1,xlab="% Akan in 30km (spatially weighted)", ylab="", cex.lab=1)
dev.off()

#### Figure 4 ####
pdf(file="SupM_fig4.pdf", height=8, width=3)
par(mfrow=c(3,1))

plot(data$akan_30km_l_p, data$akan_20km_l_p, xlim=c(0,1), ylim=c(0,1), pch=16, col="dodgerblue4", cex=.7, xlab="Akan in 30km (spatially weighted)", ylab="Akan in 20km (spatially weighted)", main="(a)", cex.lab=1.4)
abline(a=0, b=1, col="firebrick", lty="dashed")

plot(data$akan_30km_l_p, data$akan_40km_l_p, xlim=c(0,1), ylim=c(0,1), pch=16, col="dodgerblue4", cex=.7, xlab="Akan in 30km (spatially weighted)", ylab="Akan in 40km (spatially weighted)", main="(b)", cex.lab=1.4)
abline(a=0, b=1, col="firebrick", lty="dashed")

plot(data$akan_30km_l_p, data$akan_30km_lalt_p, xlim=c(0,1), ylim=c(0,1), pch=16, col="dodgerblue4", cex=.7, xlab="Akan in 30km (regular weights)", ylab="Akan in 30km (alternative weights)", main="(c)", cex.lab=1.4)
abline(a=0, b=1, col="firebrick", lty="dashed")

dev.off()


#### Figure 5 ####
hdat$sum <- hdat$npp2008ps_pres_p + hdat$ndc2008ps_pres_p
err <- hdat[hdat$sum > 1,]
dim(err)
head(err)
hdat2 <- hdat[hdat$sum <= 1,]

pdf(file= "SupM_fig5.pdf", width=6, height=6, family="Helvetica")
plot(hdat2$npp2008ps_pres_p, hdat2$ndc2008ps_pres_p, pch=16, col="dodgerblue4", xlim=c(0,1), ylim=c(0,1), xlab="2008 NPP Pres. Vote Share (polling station)", ylab="2008 NDC Pres. Vote Share (polling station)", cex=.6)
abline(a=1, b=-1, col="firebrick", lwd=1.2)
dev.off()


#### Figure 6 ####
model2n <- lmer(npp2008ps_pres_p ~ akan_30km_l_p + akan_p_poly + moledagbon_p_poly + otherethn_p_poly + public_semipublic_p + dev_factor2 + (1 + akan_30km_l_p |c230_id_h), data=hdat)
summary(model2n)
#display() gets masked by loading xtable to make the plot above...

coef(model2n)$c230_id_h[,2]
se.coef(model2n)$c230_id_h[,2]
plotmat <- matrix(c(coef(model2n)$c230_id_h[,2], coef(model2n)$c230_id_h[,2] - 1.96*se.coef(model2n)$c230_id_h[,2], coef(model2n)$c230_id_h[,2] + 1.96*se.coef(model2n)$c230_id_h[,2]), nrow=length(unique(hdat$c230_id_h)), ncol=3, byrow=FALSE)
fixef(model2n)[2]

out <- c()
for(i in 1:length(unique(hdat$c230_id_h))){
	out[i] <- mean(hdat$akan_p[hdat$c230_id_h==unique(hdat$c230_id_h)[i]])	}

pdf(file="SupM_fig6.pdf", family="Helvetica", height=7, width=7)
plot(out, plotmat[,1], col="dodgerblue4", pch=19, cex=.8, ylim=c(0,0.4), main="", ylab="Coefficient on % Akan in 30km (spatially weighted)", xlab="% Akan in Constituency", xlim=c(0,1))
abline(h=fixef(model2n)[2], lty="dashed", lwd=1.5, col="firebrick")
abline(h=0, lty="dotted", lwd=1.5, col="darkgreen")
segments(out, plotmat[,2], out, plotmat[,3], col="gray23", lwd=1.1)
dev.off()

