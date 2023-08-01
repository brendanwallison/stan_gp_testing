library(geodata)
source("utils.r")

globe = worldclim_global(var='bio', res=10, path= "./bioclim/globe", version = "2.1")
# 
# climVar <- worldclim_country(country= "US", var="bio", path= "./bioclim/output", version = "2.1")
# plot(climVar, 1)
# 
# climVar <- worldclim_country(country=c("MX", "US","CA"), res=10, var="bio", path= "./Output", version = "2.1")
# plot(climVar, 1)

#test = aggregate(climVar, fun="mean", filename="")
#plot(test, 1)

#globe2 = aggregate(globe, fun="mean", overwrite=FALSE)

test_arr = as.array(globe)
e = ext(-124.848974, -66.885444, 24.396308, 49.384358)
e_expanded = e
e_expanded[c(1,3)] = e_expanded[c(1,3)] - 10.0
e_expanded[c(2,4)] = e_expanded[c(2,4)] + 10.0

cropped_globe = crop(globe, e_expanded)
plot(cropped_globe, 1)
cellSize(cropped_globe)

lats = yFromRow(cropped_globe)
longs = xFromCol(cropped_globe)

lat1 = lats[1:length(lats)-1]
lat2 = lats[2:length(lats)]
lon1 = longs[1:length(longs)-1]
lon2 = longs[2:length(longs)]
lat1 - lat2
lon1 - lon2
# as expected, latitude and longitude are constant differences between neighboring points
# distances and areas will differ

proj_crs = "ESRI:102010"
# for our application, x and y should be the same resolution
# note that it is the distance between cells along each axis
res= c(15000.0, 15000.0) 
method = "cubicspline"
projected_globe = project(cropped_globe, proj_crs, TRUE, method=method, res=res)
plot(projected_globe, 1)
test = crop(projected_globe, e)


lats = yFromRow(projected_globe)
longs = xFromCol(projected_globe)

y1 = lats[1:length(lats)-1]
y2 = lats[2:length(lats)]
x1 = longs[1:length(longs)-1]
x2 = longs[2:length(longs)]
y1 - y2
x1 - x2
# these are now distances between neighbors along a single axis at a time
# As expected from an equidistant projection, steps are all the same distance
# regardless of which axis we are stepping along


xyFromCell(cropped_globe, 1)
xyFromCell(cropped_globe, 2)





allBioVars = as.array(cropped_globe)

BIO1 = allBioVars[,,1] # annual mean temperature
visualize(BIO1)
BIO12 = allBioVars[,,12] # annual precipitation
visualize(BIO12)

BIO2 = allBioVars[,,2] # mean diurnal range
visualize(BIO2)
BIO15 = allBioVars[,,15] # precipitation seasonality
visualize(BIO15)



mean_cellsize = mean(as.array(cellSize(cropped_globe)))
# for numerical stability reasons and/or more interpretable
# distance and dispersal matrices, kernel_units_in_n_kms rescales
# distance matrix from kilometers to ns of kilometers



# 
# 
# 
# #projected_cropped_globe = project(cropped_globe, "EPSG:102003")
# 
# wc <- worldclim_country("Iceland",  "tavg", ".")
# wc[1]
# wc[,,1]
# 
# wc2 <- worldclim_country("Iceland",  "bio", ".")
# wc2
# 
# wc3 <- worldclim_country("Mexico",  "bio", ".")
# wc3
# plot(wc3, 1)
# 
# wc4 <- worldclim_country("USA",  "bio", ".")
# summary(wc4)
# 
# wc4
# plot(wc4, 1)
# 
# 
# subset(wc, 1)
