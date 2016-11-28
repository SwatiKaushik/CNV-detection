#implement CBS as discussed in the PSCBS package 
# first 

library(PSCBS)

#load test dataset

data <- read.table("sample_data", header=T)
data <- data[, c("chromosome", "x", "CT")]

# drop outliers 
data <- dropSegmentationOutliers(data)

#removing gaps
gaps <- findLargeGaps(data, minLength = 1e+01)
knownSegments <- gapsToSegments(gaps)

#segmenting 
fit <- segmentByCBS(data, knownSegments = knownSegments, seed = 48879, verbose = -10)
