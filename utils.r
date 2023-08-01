#devtools::install_github("ropenscilabs/rnaturalearth")
#devtools::install_github("ropenscilabs/rnaturalearthdata")
library("rnaturalearth")
library("rnaturalearthdata")
library("dplyr")
library("sf")
library("reshape2")
library("ggplot2")

### FUNCTIONS ###

shift_centered_array <- function(arr){
  row_shift <- -(dim(arr)[1]%/%2)
  col_shift <- -(dim(arr)[2]%/%2)
  arr <- circshift(arr, c(row_shift, col_shift))
  arr
}

unshift_array <- function(arr){
  row_shift <- (dim(arr)[1]%/%2)
  col_shift <- (dim(arr)[2]%/%2)
  arr <- circshift(arr, c(row_shift, col_shift))
  arr
}

# the value of each cell is r distance from the central cell
# input can be odd or even dimensions
# matrix will be ~4x original size
# scalar is the distance from any given cell to its immediate (queen's) neighbors
synthetic_distance_matrix <- function(basegrid_dim, scalar = 1){
  kernel <- matrix(0, basegrid_dim[1] * 2 - 1, basegrid_dim[2] * 2 - 1)
  center <- basegrid_dim
  for (x in 1:dim(kernel)[1]) {
    for (y in 1:dim(kernel)[2]) {
      kernel[x, y] <- ((x - center[1])^2 + (y - center[2])^2) ^ (1 / 2)
    }
  }
  kernel * scalar
}

# plug base dispersal kernel into this
exp_power_kernel <- function(a, b, distances, scalar = 1) {
  c <- b * exp(-(distances / a) ^ b) / (2 * pi * a^2 * gamma(2 / b))
  return (c * scalar/sum(c))
}

normal_kernel <- function(sigma, distances, scalar=1) {
  c <- exp(-distances / (2*sigma^2)) / (sigma * sqrt(2 * pi))
  return (c * scalar/sum(c))
}

ricker <- function(r, N, K) {
  N * exp(r * (1 - N / K))
}

# spatial feature grid is indexed starting at rowwise
# convert this indexing to row, col indexing
sfgrid_to_rowcol_idx <- function(nrows, ncols, idx){
  r <- nrows - (idx - 1)%/%ncols
  c <- idx%%ncols
  cbind(r, c)
}

next_2 <- function(n){
  round(2^ceiling(log2(n)))
}

# cellsize 10,000 == 100kmx100km
build_sandbox <- function(cellsize, kernel_units_in_n_kms) {
    # get north american map
    na <- ne_countries(continent = 'north america', returnclass = "sf")
    map <- na %>% dplyr::filter(
        sovereignt == "United States of America" | sovereignt == "Canada" )

    # project to North America Equidistant Conic, crop out some islands
    pna <- map %>%
    st_transform(crs = st_crs("ESRI:102010")) %>%
    st_geometry() %>%
    st_crop(xmin = -4000000, ymin = -1800000, xmax = 3800000, ymax = 4936902) %>%
        st_union

    ncols <- ceiling(diff(st_bbox(pna)[c(1, 3)]) / cellsize)
    nrows <- ceiling(diff(st_bbox(pna)[c(2, 4)]) / cellsize)

    contig_grid <- pna %>%
    st_make_grid(cellsize = cellsize, square = TRUE)
    intersections <- st_intersection(contig_grid, pna)
    land_idxes <- attr(intersections, 'idx') [, 1]
    base_filter <- synthetic_distance_matrix(c(nrows, ncols), sqrt(cellsize)/kernel_units_in_n_kms)
    land_idxes_rc <- sfgrid_to_rowcol_idx(nrows, ncols, land_idxes)

    # linear convolution with fft requires zero padding to
    # dim(matrix) + dim(filter) -1
    # fastest fft is when each dim is of the form 2^n
    padded_rownum <- next_2(nrows + dim(base_filter)[1] - 1)
    padded_colnum <- next_2(ncols + dim(base_filter)[2] - 1)

    padded_mask <- array(0, dim = c(padded_rownum, padded_colnum))
    padded_mask[land_idxes_rc] <- 1
    output <- list(pna, base_filter, padded_mask, c(nrows, ncols))
    return(output)
}

visualize <- function(arr) {
  df <- melt(arr)
  ggplot(df, aes(X2, X1, fill = value)) +
    geom_tile() +
    scale_y_reverse()
}
  


### END FUNCTIONS ###
