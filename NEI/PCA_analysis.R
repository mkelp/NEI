

library(RNetCDF)
library((ggplot2))
nc <- open.nc ("spatialResults_mkelp_test_afdust.ncf")
print.nc(nc)

ncread <- read.nc(nc)

#do for each pollutant 
ncprim <- princomp(ncread$`PM25-PRI`, scores = TRUE)

summary(ncprim)
plot(ncprim)

#this ones takes long
biplot(ncprim)




y<- c(0.0002,0.0002,0.0002,0.0003,0.0003,0.0003,0.0004,0.0004,0.0004,0.0005,0.0005,0.0005,0.0006,0.0007,0.0008,0.0009,0.0011,0.0014,0.0023,0.0035,0.0064,0.0168,0.6470)
cy <- cumsum(y)


plot(x=1:length((cy)), y=cy)




