## 
## Topic:  Spatial Regression
## Author: Shaky Sherps
## Date:   November 15 2013
## Course: GTECH 705 - Spatial Data Analysis
## Lab 4
# Here's what we're going to do:
# 1) Define neighbors (k nearest w/point data)
# 2) Create weights matrix
# 3) Moran’s test of observed, Moran scatterplot
# 4) Run OLS regression
# 5) Check residuals for spatial dependence
# 6) Determine which SR model to use w/LM tests
# 7) Run spatial regression model
##


# Clear the workspace
rm(list=ls())

#load required libraries
library(spdep)
library(spgwr)
library(rgdal)

##
## Spatial regression workflow
##

# Load some data
dublin = readOGR("data/dublin_voters.shp","dublin_voters")

# Check what this does
IDs = row.names(dublin)

#Regression model
elec =GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffADD 

# Create neighbors list
dub.gal.nb = poly2nb(dublin)

# Create weights matrix
dub.listw = nb2listw(dub.gal.nb, style="W")

par(mar=c(0,0,0,0))
plot(dublin, dub='lightgrey', border='white')
plot(dub.listw, coords=coordinates(dublin), add=TRUE)

# Test for spatial autocorrelation
moran.test(dublin$GenEl2004, listw=dub.listw, alternative="two.sided")

# Test for local spatial autocorrelation
dev.off()
moran.plot(dublin$GenEl2004, dub.listw, ylab="Spatially lagged Election Turnout", xlab="percentage of population in each ED who voted in 2004 election")
dub.li = localmoran(dublin$GenEl2004, dub.listw)
dublin$localm = dub.li[,4]
spplot(dublin, "localm", main="Local Moran's Ii Z-Scores")

# Run OLS Regression

dub.lm = lm(GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffAdd, data=dublin)
summary(dub.lm)

# Check residuals for spatial dependence
dublin$lmres = residuals(dub.lm) # Grab the residuals (though don't really need to do this)
lm.morantest(dub.lm, dub.listw)

# Perform lagrage multiplier test
# Robust tests used to find a proper alternative
# Only use robust forms when BOTH LMErr and LMLag are significant
lm.LMtests(dub.lm, dub.listw, test="all")

# install.packages("lmtest")
library(lmtest)
bptest(dub.lm)
# Indicates errors are heteroskedastic
# Not surprising since we have spatial dependence!

# Fit spatial regression models

# Spatial lag model(Spatial Simultaneous Autoregressive lag model estimation)
dub.lag = lagsarlm(GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffAdd, data=dublin, listw=dub.listw)
summary(dub.lag)

# Some more diagnostics
bptest.sarlm(dub.lag)#Breusch-Pagan test
# LM test suggests there is no more spatial autocorrelation in the data
# BP test indicates remaining heteroskedasticity in the residuals
#     Most likely due to misspecification

# Spatial error model(Maximum likelihood estimation of spatial simultaneous autoregressive error models of the form:

dub.err = errorsarlm(GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffAdd, data=dublin, listw=dub.listw)
summary(dub.err)
bptest.sarlm(dub.err)


## GWR examples (we'll get to this again later)
##

dub.lm = lm(GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffAdd, data=dublin)
dublin$lmres = residuals(dub.lm) # Extract Model Residuals

#Crossvalidation of bandwidth for GWR
dub.bw = gwr.sel(GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffAdd, data=dublin)
dub.gauss = gwr(GenEl2004 ~ Unempl + LowEduc + Age18_24 + Age25_44 + Age45_64 +SC1 + LARent - DiffAdd, data=dublin, bandwidth=dub.bw, hatmatrix=TRUE)

spplot(dublin, "GenEl2004", main="General Election Turnout in 2004")
spplot(dublin, "Unempl", main="Percentage of the population in each ED who are unemployed")
spplot(dublin, "LowEduc", main="Percentage of the population in each ED who are with little formal education")
spplot(dublin, "Age18_24", main="Percentage of the 18-24 year old population in each ED")
spplot(dublin, "Age25_44", main="Percentage of the 25-44 year old population in each ED")
spplot(dublin, "Age45_64", main="Percentage of the 45-64 year old population in each ED")
spplot(dublin, "SC1", main="percentage of the population in each ED who are social class one (high social class")
#Plot the LM Residuals
spplot(dublin, "lmres", main="Linear Model Residuals")

spplot(dub.gauss$SDF, "Unempl", main="Percentage of the population in each ED who are unemployed")
spplot(dub.gauss$SDF, "LowEduc", main="Percentage of the population in each ED who are with little formal education")

res = dub.gauss$SDF
dublin$gwrres = dublin$GenEl2004 - res$pred
spplot(dublin, "gwrres", main="GWR Model Residuals")

