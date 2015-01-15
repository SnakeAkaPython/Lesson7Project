#Snake Aka Pyhton
# 15 January 2015


## libraries
library(raster)


## load data
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/vcfGewata.rda")


vcfGewata


# To avoid values much greater the 100
vcfGewata[vcfGewata > 100] <- NA
plot(vcfGewata, main = "VCF")
summary(vcfGewata)


# Put the 6 bands into a rasterBrick object to summarize together
gewata <- brick(GewataB1,GewataB2, GewataB3, GewataB4, GewataB5, GewataB7)


# Rescale the original reflectance values to their original scale
gewata <- calc(gewata, fun=function(x) x / 10000)


# Make a new raster brick of covariates by adding VCF layer
covs <-  addLayer(gewata, vcfGewata)


# Exploring relationships between raster layers
pairs(covs)
# We conclude that bands 2, 3, 5 and 7 rellate better to the VCF layers


names(covs) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")


# Extract all values into a matrix
valuetable <- getValues(covs)


# Get rid of NA values
valuetable <- na.omit(valuetable)
valuetable <- as.data.frame(valuetable)


LMmodel <- lm(formula = VCF ~ band2 + band3 + band5, data = valuetable)
summary(LMmodel)
# The most important bands in predicting the tree cover is band 2, 3, 5 or 7 

# Predict tree cover 
predictVcf <-predict(covs, model=LMmodel, na.rm=T)

predictVcf[predictVcf < 0] <- 0
names(predictVcf) <- "Pred_VCF"

par(mfrow = c(2, 1))
plot(predictVcf, main = "Predicted VCF")
plot(covs$VCF, main = "Original VCF")

# Calculate the Root maen squared Error between Predicted and the actual tree cover values
Rows <- predictVcf@nrows
Cols <- predictVcf@ncols
RMSE <- sqrt(mean((covs$VCF[Rows, Cols]-predictVcf[Rows, Cols])^2))


# Load the training polygons
load("data/trainingPoly.rda")


# Inspect the data slot of the trainingPoly object
trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
trainingPoly@data


# Assign "code" values to raster cells
classes <- rasterize(trainingPoly, covs$VCF, field='Code')


# Make a new raster brick with the Predicted and the Original VCF
Vcfbrick <- brick(covs$VCF, predictVcf)


# Calculate the mean for every class by the use of zonal()
zonalStat_mean <- zonal(Vcfbrick, classes, fun= 'mean', digits=0, na.rm=T)
zonalStat_mean
zonalstatDF <- as.data.frame(zonalStat_mean)


# Calculate the Root maen squared Error for all of the 3 classes
RMSE_Zonal <- sqrt((zonalstatDF$VCF - zonalstatDF$Pred_VCF)^2)
names(RMSE_Zonal) <- c("Crop", "Forest", "Wetlands")
RMSE_Zonal
# We conclude that there is a small difference between the Predicted and the Original Tree cover