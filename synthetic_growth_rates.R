library(geodata)

globe = worldclim_global(var='bio', res=10, path= "./bioclim/globe", version = "2.1")
#climVar <- worldclim_country(country= "US", var="bio", path= "./Output", version = "2.1")

e = ext(-124.848974, -66.885444, 24.396308, 49.384358)
e_expanded = e
e_expanded[c(1,3)] = e_expanded[c(1,3)] - 5.0
e_expanded[c(2,4)] = e_expanded[c(2,4)] + 5.0

cropped_globe = crop(globe, e_expanded)

proj_crs = "ESRI:102008" #Albers equal area conic
# for our application, x and y should be the same resolution
# note that it is the distance between cells along each axis
res= c(15000.0, 15000.0) 
method = "cubicspline"
projected_globe = project(cropped_globe, proj_crs, method=method, res=res)
plot(projected_globe, 1)
# group landmasses by 1st layer (all other layers identical)
landmasses = patches(subset(projected_globe,1))
rz <- zonal(cellSize(landmasses), landmasses, sum, as.raster=TRUE)
s <- ifel(rz < 1e12, NA, projected_globe)
plot(s, 1)


# As we flatten the map through projection, NA values come arcing down from the north
# I fix them by including extra values in the north, then cropping 
# However, with With the new projection, our coordinates change from lat/long 
# to meters from a center point. So how do we know how much to crop?
# Right now, it is mostly arbitrary/eyeballed
e_proj = ext(s)
e_proj[c(1,3)] = e_proj[c(1,3)] + abs(e_proj[c(1,3)])*0.25
e_proj[c(2,4)] = e_proj[c(2,4)] - abs(e_proj[c(2,4)])*0.25
bioclim_full = trim(crop(s, e_proj))
plot(bioclim_full, 1)

# standardized each layer
scaled_bioclim_full = scale(bioclim_full)

# calculate growth-rates
# Means
# BIO1 = Annual Mean Temperature
# BIO12 = Annual Precipitation
# Seasonality
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO15 = Precipitation Seasonality (Coefficient of Variation)



allBioVars = as.array(scaled_bioclim_full)

coefs = c(3.5,2.5,1,0.5)
fvi <- function(BIO1, BIO2, BIO3, BIO4, BIO5, BIO6, BIO7, BIO8, BIO9, BIO10, 
                BIO11, BIO12, BIO13, BIO14, BIO15, BIO16, BIO17, BIO18, BIO19){ 
  (3/sum(coefs))*(exp(-coefs[1]*(BIO1-0.5)^2) + exp(-coefs[2]*(BIO12+0.25)^2) + exp(-coefs[3]*(BIO4-0.15)^2) + exp(-coefs[4]*BIO15^2)) - 0.5
  
}
r = lapp(scaled_bioclim_full, fun=fvi)
plot(r)

# Here are a few useful methods
# for coercing spatrasters to data frames and arrays

# Data Frames
bioclim_df = as.data.frame(scaled_bioclim_full)
r_df = as.data.frame(r)
k_df = exp(3*r_df)

# Matrices and arrays
r_mat = as.matrix(r, wide=TRUE) # nxm matrix of growth rates
bioclim_arr = as.array(scaled_bioclim_full) # nxmxb array, where b = 19 bioclim variables

land_mask = r_mat
ocean_indices = is.na(land_mask)
land_mask[ocean_indices] = 0
land_mask[!ocean_indices] = 1

#library("dplyr")
library("ggplot2")
library("reshape2")
# For plotting array as heatmap
visualize <- function(arr) {
  df <- melt(arr)
  ggplot(df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_y_reverse()
}

visualize(land_mask)
visualize(r_mat)

r_df[1:4,]

