### En este ejercicio vamos a estudiar las tasas de diversificacion del grupo de gramineas de la subfamilia Pooideae (Poaceae) a lo largo del tiempo. Con 4000 especies distribuidas en todo el mundo este clado es uno de los mas grandes de las angiospermas y esta distribuido por todo el mundo. 
### En este ejercicio vamos a comparar modelos en los que las tasas de diversificacion permanecen constantes, cambian a traves del tiempo, o cambian en función de los cambios de paleotemperaturas durante el Cenozoico. 

###########################################################
##run the analysis
###########################################################

setwd("/Users/ameseguer/Dropbox/Asignatura_Biogegrafia_Master/2020/Ejercicios_Profesor/Ejercicio_Dia 15 Nov-GeoSSE-RPANDA/RPANDA")

#carga las librerías
source("run_diversification_analyses.R")

## corre los modelos. Los resultados apareceran en tu working directory
run_Morlon_models("Pooideae.nexus", sampling_fraction=0.1, number_of_trees=1) #incluye modelos constantes, y time-variable
run_PaleoEnv("Pooideae.nexus", env_data_file="./PaleoEnv/PastTemperature.txt", sampling_fraction=0.1, number_of_trees=1) #incluye paleotemperature models

## organizamos los resultados
# abrimos resultados con time models
read.table("Pooideae_results_Morlon.txt", header=TRUE)->timemodels
#select rows: exclude all linear models 
timemodels[c(1:6),c(1:8)]-> timemodels2
colnames(timemodels2)<-c("Models","NP","logL","AICc","Lambda","Alpha","Mu","Beta")
# abrimos resultados con environmental models
read.table("Pooideae_results_PastTemperature.txt", header=TRUE)->temperaturemodels
#select rows: exclude all linear models 
temperaturemodels[c(1:4),c(1:8)]-> temperaturemodels2
colnames(temperaturemodels2)<-c("Models","NP","logL","AICc","Lambda","Alpha","Mu","Beta")
rbind(timemodels2, temperaturemodels2)->final

## calcule deltaAIC, for model comparisons:
library(qpcR)
final[,9]<-akaike.weights(final[,4])$weights
final[,10]<-akaike.weights(final[,4])$deltaAIC
colnames(final)<-c("Models","NP","logL","AICc","Lambda","Alpha","Mu","Beta","AICw","deltaAIC")
results<-final[order(as.numeric(final[,"deltaAIC"])),]

###########################################################
##plot temperature result
###########################################################

library("picante")
library("pspline")
tree<-read.nexus("Pooideae.nexus")
crown.age<- max(node.age(tree)$ages)

###################
# a) Paleo-temprature
###################

InfTemp<-read.table("./PaleoEnv/PastTemperature.txt", header=T) 
if (crown.age>=100){max.y<-60}
if (crown.age<=100){max.y<-50}
if (crown.age<=66){max.y=22}
if (crown.age<=35){max.y=15}
if (crown.age<=20){max.y=15}
par(mfrow=c(2,1))
plot(-InfTemp[,"Age"],InfTemp[,"Temperature"],xlab="", ylab="Temperature (∞C)", type="p", lwd=0.1, xlim=c(-crown.age,0), ylim=c(0, max.y), col="grey", las=1, cex=0.8, cex.axis=0.8, bty = "n")
temp.spl<-sm.spline(-InfTemp[,"Age"], InfTemp[,"Temperature"],df=100)
lines(temp.spl, col="firebrick1", lwd=2)
legend("topleft", bty="n", c("a) Paleoclimate"), cex=0.7)
abline(v=c(-298.9,-252.2,-201.3,-145,-100.5,-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages


######################################
# d) Diversification as a function of temperature
######################################

load('Pooideae_complete_results_PastTemperature.Rdata')

mean.results<-final_table_tree_file[[1]][order(as.numeric(final_table_tree_file[[1]][,"AICc"])),]
std.err.results<-final_table_tree_file[[2]][order(as.numeric(final_table_tree_file[[2]][,"AICc"])),]

res<-sm.spline(InfTemp[,"Age"],InfTemp[,"Temperature"],df=100) 
Temp_fun<-function(x){predict(res,x)}

# Speciation rate
# Please fill/replace the x in the objects below with the results of your best paleoenvironment-dependent analysis. 
# TempDep_lamb_par1 is lambda, TempDep_lamb_par2 is alpha.
# TempDep_lamb_par1_sd is the standard error of lambda, and TempDep_lamb_par2_sd is the standard error of alpha. 

TempDep_lamb_par1 <- 0.6335106
TempDep_lamb_par2 <- -0.21734576

TempDep_lamb_par1_sd <- 0 # canbe changed if you ran the analyses on a set of trees (then you have the standard errors table)
TempDep_lamb_par2_sd <- 0 # canbe changed if you ran the analyses on a set of trees (then you have the standard errors table)

f.lamb.mean<-function(x){TempDep_lamb_par1*exp(TempDep_lamb_par2*Temp_fun(x))}
f.lamb.low<-function(x){(TempDep_lamb_par1-TempDep_lamb_par1_sd)*exp((TempDep_lamb_par2-TempDep_lamb_par2_sd)*Temp_fun(x))}
f.lamb.high<-function(x){(TempDep_lamb_par1+TempDep_lamb_par1_sd)*exp((TempDep_lamb_par2+TempDep_lamb_par2_sd)*Temp_fun(x))}

max.y<-round(max(f.lamb.high(InfTemp[,"Age"])),1)

plot(-InfTemp[,"Age"], f.lamb.mean(InfTemp[,"Age"]), ty="l",col="chartreuse3",xlim=c(-crown.age,0), ylim=c(0,0.95), lwd=2, yaxt="n", xlab="Time (Myrs ago)",ylab="Speciation (green) and extinction (red) rates", cex.axis=0.8, bty = "n")
axis(2, at = seq(0, 0.95, by = 0.1), las=1, cex.axis=0.6)
legend("topleft", bty="n", c("d) Diversification according to temperature"), cex=0.7)
lines(-InfTemp[,"Age"], f.lamb.low(InfTemp[,"Age"]),ty="l",col="chartreuse3",xlim=c(-crown.age,0),lwd=1,lty="dotted",yaxt="n")
lines(-InfTemp[,"Age"], f.lamb.high(InfTemp[,"Age"]),ty="l",col="chartreuse3",xlim=c(-crown.age,0),lwd=1,lty="dotted",yaxt="n")

# Extinction rate
# Please fill/replace the x in the objects below with the results of your paleoenvironment-dependent analysis. 
# TempDep_mu_par1 is mu, and TempDep_mu_par2 is beta.
# TempDep_mu_par1_sd is the standard error of mu, and TempDep_mu_par2_sd is the standard error of beta. 

TempDep_mu_par1 <- 0
TempDep_mu_par2 <- 0
TempDep_mu_par1_sd <- 0 # canbe changed if you ran the analyses on a set of trees (then you have the standard errors table
TempDep_mu_par2_sd <- 0 # canbe changed if you ran the analyses on a set of trees (then you have the standard errors table

f.mu.mean<-function(x){TempDep_mu_par1*exp(TempDep_mu_par2*Temp_fun(x))}
f.mu.low<-function(x){(TempDep_mu_par1-TempDep_mu_par1_sd)*exp((TempDep_mu_par2-TempDep_mu_par2_sd)*Temp_fun(x))}
f.mu.high<-function(x){(TempDep_mu_par1 + TempDep_mu_par1_sd)*exp((TempDep_mu_par2 + TempDep_mu_par2_sd)*Temp_fun(x))}

lines(-InfTemp[,"Age"], f.mu.mean(InfTemp[,"Age"]), ty="l",col="red",lwd=1, yaxt="n", xlab="Time (Myrs ago)",ylab="Speciation (blue) and extinction (red) rates", cex.axis=0.8, bty = "n")
lines(-InfTemp[,"Age"], f.mu.low(InfTemp[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")
lines(-InfTemp[,"Age"], f.mu.high(InfTemp[,"Age"]),ty="l",col="red",lwd=1,lty="dotted",yaxt="n")

abline(v=c(-298.9,-252.2,-201.3,-145,-100.5,-66,-56,-47.8,-33.9,-28.1,-23.03,-15.97,-11.62,-5.33,-2.58),col="grey",lty="dotted",lwd="1") # add vertical lines to delineate geological periods/epochs/stages



