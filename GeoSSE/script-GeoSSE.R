
# We will be working with the dated phylogeny of Hypericum from Meseguer et al. (2015). We want to test whether diversification rates (extinction and speciation) are linked to the trait "area" using state-dependent speciation extinction (SSE) models (Maddison et al. 2007). Specifically, we will test the Rand Flora hypothesis that extinction rates have been historically higher in Africa than in other regions. We will compare speciation and extinction rates in Africa (AFR) versus speciation and extinction rates in other (non-AFR) regions; the latter comprises all other areas where Hypericum is distributed (Western Palearctic, Nearctic, Eastern Palearctic, South America, and Oceania). We will also test if migration rates have been asymmetric, predominantly from other regions into Africa or the other way around. 

setwd('~/Dropbox/Master-UIMP-Biogeography/2021/Ejercicios_Profesor/Ejercicio-IV-GeoSSE/GEOSSE/')

#  Install the R package diversitree (FitzJohn, 2015) 
install "diversitree"

# Load diversitree:
library(diversitree)

# Upload the dated tree 
#Remember to set first your working directory using the command setwd() and the path to the folder where the input files are located.
tree <- read.tree("Hypericum-newick.tre")

# Read the data file with the distributions for each species: 
# there are three states coded as follows: “1” (AFR), “2” (non-AFR), “0” (widespread: AFR+non-AFR)
data <- read.csv("data.csv", row.names=1, sep=";")

# Transform the data to identify the first column as names and the second as the data.
data.v <- data [,1]
names(data.v) <- row.names(data)
tree$tip.state <-data.v

#We draw the tree with colors assigned to each state for visualization purposes.
statecols <- c("0"="green", "1"="blue", "2"="red")
plot(tree, tip.color=statecols[tree$tip.state+1], cex=0.2)

#  create the GeoSSE likelihood function. 
#The first function we create is the full model in which there is regional dependency of speciation, dispersal or extinction rates. 
lik <- make.geosse (tree,tree$tip.state)

#We start the function by giving some initial values
p <- starting.point.geosse(tree)

# We estimate the Maximum Likelihood Estimator for the function given the parameter values.
fit <- find.mle(lik, p, method="subplex")

#Then, we estimate the log likelihood of the model.
logLik(fit)

#Obtain the MLE estimates of parameters
round(coef(fit), 3)
#Observe that the codes for the states are as following: "A" = 1 (AFR); "B" = 2 (non-AFR); "0" = AB (AFR+non-AFR).

#Can you interpret the results? What does it tell you in terms of speciation, extinction, dispersal rates in AFR vs. non-AFR regions?
 
# We now construct a model with equal dispersal (transition) rates between regions and estimate the likelihood.
lik.d <- constrain(lik, dA ~ dB)
fit.d <- find.mle(lik.d, p[-7])
logLik(fit.d)

#Obtain the MLE estimates of parameters
round(coef(fit.d), 3)

#Can you interpret the results? Are they different from the unconstrained “full” model above (fit)?


# We construct a model in which speciation rates are equal between regions and there is not between-region speciation.
lik.s <- constrain(lik, sA ~ sB, sAB ~ 0)
fit.s <- find.mle(lik.s, p[-c(2,3)])
logLik(fit.s)

#Obtain the MLE estimates of parameters
round(coef(fit.s), 3)

#Can you interpret the results? Are they different from the unconstrained (“full”) model above? And from the “constrained dispersal” model above (lik.d)?

# Finally, we compare the three models using an ANOVA test
anova(fit, equal.d = fit.d, equal.s = fit.s)

######################
### Now we run the same model in a bayesian framework (mcmc)
######################

# como valores iniciales le damos los valores del mejor resultado obtenido por maximum likelihood arriba
p <- coef(fit)

# Prior for all parameters, exponential lamda=0.5 (More probability around less events).
prior <- make.prior.exponential(1/2)

# Random number generator to start markov chain
set.seed(1) 

# preliminary search to identify tunning parameters bellow
tmp <- mcmc(lik, p, nsteps = 100, prior = prior, lower = 0, w = rep(1, 7))

# Sliding window. Tuning parameter for the sampler. w affects how many function evaluations are required between sample updates
w <- diff(sapply(tmp[2:8], quantile, c(0.025, 0.975)))

# run the mcmc chain
mcmc2 <- mcmc(lik, p, nsteps = 3000, prior = prior, lower = 0, w = w)

# a partir de los resultados calculamos las tasas de diversificacion para cada area
mcmc2diff <- with(mcmc2, data.frame(s.diff = sA - sB, x.diff = xA - xB, d.diff = dA - dB, div.A = sA - xA, div.B = sB - xB))
#colMeans(mcmc2diff > 0)

head(mcmc2)
head(mcmc2diff)

### ploteamos los resultados

col1 <- c("red", "orange", "blue", "purple", "black", "gray", "green")
#col2 <- col1[c(1, 3, 5)]

par(mfrow = c(1,4), mar = c(3, 4, 1, 1))
# comparamos las tasas de especiacion
profiles.plot(mcmc2[2:3], col.line = col1, xlab = "", ylab = "", main="speciation")
legend("topright", argnames(lik), col = col1, lty = 1)
# comparamos las tasas de extincion
profiles.plot(mcmc2[5:6], col.line = c("purple","black"), xlab = "", ylab = "", main="extinction")
# comparamos las tasas de dispersion
profiles.plot(mcmc2[7:8], col.line = c("gray", "green"), xlab = "", ylab = "",main="dispersal")
# comparamos las tasas de diversification
profiles.plot(mcmc2diff[4:5], col.line = c("lightblue","lightgreen") , xlab = "", ylab = "", main="diversification")
legend("topleft", c("Africa","Elsewhere"),col = c("lightblue","lightgreen") , lty = 1)
 

# también podemos calcular si las probabilidades posteriores estimadas para cada area son significativamente diferentes
par(mfrow = c(1, 3), mar = c(3, 4, 1, 1))
#son significativas estas diferencias?
profiles.plot(mcmc2diff[1], col.line = "lightblue", xlab = "", ylab = "", main="sA-sB")
abline(v=mean(as.numeric(as.character(mcmc2diff[,1]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
profiles.plot(mcmc2diff[2], col.line = "lightblue", xlab = "", ylab = "", main="xA-xB")
abline(v=mean(as.numeric(as.character(mcmc2diff[,2]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)
profiles.plot(mcmc2diff[3], col.line = "lightblue", xlab = "", ylab = "", main="dAB-dBA")
abline(v=mean(as.numeric(as.character(mcmc2diff[,3]))), col="black", lwd=2)
abline(v=0, col="red", lwd=2)


 