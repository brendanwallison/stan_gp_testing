library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(dplyr)
library(posterior)
library(bayesplot)
library(ggplot2)

color_scheme_set("brightblue")

file <- file.path("uni_input_multi_output_gp.stan")
mod <- cmdstan_model(file)

# To-Do: Inducing Points
# Currently just takes a random sample of 100 points
N = 100
N_datapoints = dim(r_df)[1]
idx = sample(N_datapoints, N)
X = bioclim_df[idx,1] # Predict from just the first Bioclim variable
Y_all = as.matrix(cbind(r_df, k_df))


Y = Y_all[idx,]

# Visualize our training points
df_vis <- data.frame(X, Y[,1])
df_vis %>%
  ggplot(aes(x=X,y=Y[,1]))+
  geom_point()+
  labs(x="Bioclim 1", y="Growth Rate")

# Visualize our training points
df_vis <- data.frame(X, Y[,2])
df_vis %>%
  ggplot(aes(x=X,y=Y[,2]))+
  geom_point()+
  labs(x="Bioclim 1", y="Carrying Capacity")

# Predictive inference on grid of x values spanning dataset
# In real example, posterior predictions will be over the entire dataset X
xmin = min(bioclim_df[,1])
xmax = max(bioclim_df[,1])
xgrid = seq(xmin, xmax, 0.01)

# names correspond to the data block in the Stan program
data_list <- list(N1 = N, D=2, x1=X, y1=Y, N2 = length(xgrid), x2=xgrid)

fit <- mod$sample(
  data = data_list, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

fit$summary()

# Visualization learned from https://avehtari.github.io/casestudies/Motorcycle/motorcycle_gpcourse.html
odraws_gpcovf <- as_draws_rvars(fit$draws())
subset(odraws_gpcovf, variable=c('sigma','rho','alpha'), regex=TRUE)


Ef=mean(odraws_gpcovf$f)
#Ef=mean(odraws_gpcovf$f2)
#Ef=mean(odraws_gpcovf$f2_with_mo_corr)
sigma=mean(odraws_gpcovf$sigma)
rho=mean(odraws_gpcovf$rho)


ggplot()+
  geom_point(aes(x=X,y=Y[,1]))+
  labs(x="Bioclim 1", y="Growth Rate")+
  geom_line(aes(x=X, y=Ef[,1]))+
  geom_line(aes(x=X, y=Ef[,1]-2*sigma), linetype="dashed")+
  geom_line(aes(x=X, y=Ef[,1]+2*sigma), linetype="dashed")

ggplot()+
  geom_point(aes(x=X,y=Y[,2]))+
  labs(x="Bioclim 1", y="Growth Rate")+
  geom_line(aes(x=X, y=Ef[,2]))+
  geom_line(aes(x=X, y=Ef[,2]-2*sigma), linetype="dashed")+
  geom_line(aes(x=X, y=Ef[,2]+2*sigma), linetype="dashed")

###
Ef=mean(odraws_gpcovf$f2)
ggplot()+
  geom_point(aes(x=X,y=Y[,1]))+
  labs(x="Bioclim 1", y="Growth Rate")+
  geom_line(aes(x=xgrid, y=Ef[,1]))+
  geom_line(aes(x=xgrid, y=Ef[,1]-2*sigma), linetype="dashed")+
  geom_line(aes(x=xgrid, y=Ef[,1]+2*sigma), linetype="dashed")

ggplot()+
  geom_point(aes(x=X,y=Y[,2]))+
  labs(x="Bioclim 1", y="Growth Rate")+
  geom_line(aes(x=xgrid, y=Ef[,2]))+
  geom_line(aes(x=xgrid, y=Ef[,2]-2*sigma), linetype="dashed")+
  geom_line(aes(x=xgrid, y=Ef[,2]+2*sigma), linetype="dashed")

# 
# sq_exponential <- function(x1, x2, sigma, lengthscale) { # create a function with the name my_function
#   distances = outer(x1, x2, "-")
#   distances = sigma^2*exp(-distances^2/(2*lengthscale^2))
# }
# 
# x1x1 = sq_exponential(X, X, rho, 1.0)
# diag(x1x1) = diag(x1x1) + sigma 
# x2x2 = sq_exponential(xgrid, xgrid, rho, 1.0)
# x1x2 = sq_exponential(X, xgrid, rho, 1.0)
# x1x2invx1x1 = t(x1x2)%*%solve(x1x1)
# 
# ustar = x1x2invx1x1%*%Y[,1]
# cstar = x2x2 - x1x2invx1x1%*%x1x2
# 
# library(MASS)
# test= mvrnorm(10, ustar, cstar)
# plot(xgrid, test[2,])
# plot(bioclim_df[,1], Y_all[,1])
