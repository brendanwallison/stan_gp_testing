library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
library(ggplot2)
color_scheme_set("brightblue")

file <- file.path("multi_input_uni_output_gp.stan")
mod <- cmdstan_model(file)

# To-Do: Inducing Points
# Currently just takes a random sample of 100 points
N = 100
N_datapoints = dim(r_df)[1]
idx = sample(N_datapoints, N)
X = as.matrix(bioclim_df[idx,c(1,12)])
Y = r_df[idx,]

# Visualize our training points
ggplot()+
geom_point(aes(x=X[,1],y=Y))+
labs(x="Bioclim 1", y="Growth Rate")

# Visualize our training points
ggplot()+
geom_point(aes(x=X[,2],y=Y))+
labs(x="Bioclim 12", y="Growth Rate")

# Predictive inference on grid of x values spanning dataset
# In real example, posterior predictions will be over the entire dataset X
x1min = min(bioclim_df[,1])
x1max = max(bioclim_df[,1])
x2min = min(bioclim_df[,12])
x2max = max(bioclim_df[,12])
x1vec = seq(x1min, x1max, 0.1)
x2vec = seq(x2min, x2max, 0.1)

xgrid = as.matrix(expand.grid(x1vec, x2vec))

# names correspond to the data block in the Stan program
data_list <- list(N1 = N, input_dim=2, x1=X, y1=Y, N2 = dim(xgrid)[1], x2=xgrid)

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
#subset(odraws_gpcovf, variable=c('sigma','rho','alpha'), regex=TRUE)
Ef=mean(odraws_gpcovf$f2)
sigma=mean(odraws_gpcovf$sigma)

new_df = as.data.frame(xgrid)
new_df$Y = Ef

ggplot(new_df, aes(Var1, Var2, z= Y)) +
  geom_contour_filled()

