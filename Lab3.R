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
dub.gal.nb = read.gal(system.file("etc/weights/dublin.gal", package="spdep")[1])

# Create neighbors list
dub.gal.nb2 = poly2nb(dublin)

# Create weights matrix
dub.listw = nb2listw(dub.gal.nb2, style="W")

par(mar=c(0,0,0,0))
plot(dublin, dub='lightgrey', border='white')
plot(dub.listw, coords=coordinates(dublin), add=TRUE)

# Test for spatial autocorrelation
moran.test(dublin$CRIME, listw=dub.listw, alternative="two.sided")

# Test for local spatial autocorrelation
dev.off()
moran.plot(dublin$CRIME, dub.listw, ylab="Spatially lagged CRIME", xlab="CRIME")
dub.li = localmoran(dublin$CRIME, dub.listw)
dublin$localm = dub.li[,4]
spplot(dublin, "localm", main="Local Moran's Ii Z-Scores")

# Run OLS Regression
dub.lm = lm(CRIME ~ INC + HOVAL, data=dublin)
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

# Spatial lag model
dub.lag = lagsarlm(CRIME ~ INC + HOVAL, data=dublin, listw=dub.listw)
summary(dub.lag)

# Some more diagnostics
bptest.sarlm(dub.lag)
# LM test suggests there is no more spatial autocorrelation in the data
# BP test indicates remaining heteroskedasticity in the residuals
#     Most likely due to misspecification

# Spatial error model
dub.err = errorsarlm(CRIME ~ INC + HOVAL, data=dublin, listw=dub.listw)
summary(dub.err)
bptest.sarlm(dub.err)

##
## GWR examples (we'll get to this again later)
##

dub.lm = lm(CRIME ~ INC + HOVAL, data=dublin)
dublin$lmres = residuals(dub.lm)

dub.bw = gwr.sel(CRIME ~ INC + HOVAL, data=dublin)
dub.gauss = gwr(CRIME ~ INC + HOVAL, data=dublin, bandwidth=dub.bw, hatmatrix=TRUE)

spplot(dublin, "CRIME", main="Crimes")
spplot(dublin, "INC", main="Household Income ($1,000s)")
spplot(dublin, "HOVAL", main="Housing Value ($1,000s)")
spplot(dublin, "lmres", main="Linear Model Residuals")

spplot(dub.gauss$SDF, "INC", main="Household Income Params")
spplot(dub.gauss$SDF, "HOVAL", main="Housing Value Params")

spplot(dub.gauss$SDF, "INC_se", main="Household Income Std. Error")
spplot(dub.gauss$SDF, "HOVAL_se", main="Housing Value Std. Error")
res = dub.gauss$SDF
dublin$gwrres = dublin$CRIME - test$pred
spplot(dublin, "gwrres", main="GWR Model Residuals")

