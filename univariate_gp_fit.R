library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

file <- file.path("uni_input_uni_output_gp.stan")
mod <- cmdstan_model(file)

# To-Do: Inducing Points
# Currently just takes a random sample of 100 points
N = 100
N_datapoints = dim(r_df)[1]
idx = sample(N_datapoints, N)
X = bioclim_df[idx,1] # Predict from just the first Bioclim variable
Y = r_df[idx,]

# Visualize our training points
df <- data.frame(X, Y)
df %>%
  ggplot(aes(x=X,y=Y))+
  geom_point()+
  labs(x="Bioclim 1", y="Growth Rate")

# Predictive inference on grid of x values spanning dataset
# In real example, posterior predictions will be over the entire dataset X
xmin = min(bioclim_df[,1])
xmax = max(bioclim_df[,1])
xgrid = seq(xmin, xmax, 0.01)

# names correspond to the data block in the Stan program
data_list <- list(N1 = N, x1=X, y1=Y, N2 = length(xgrid), x2=xgrid)

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
Ef=mean(odraws_gpcovf$f2)
sigma=mean(odraws_gpcovf$sigma)

ggplot()+
geom_point(aes(x=X,y=Y))+
labs(x="Bioclim 1", y="Growth Rate")+
geom_line(aes(x=xgrid, y=Ef))+
geom_line(aes(x=xgrid, y=Ef-2*sigma), linetype="dashed")+
geom_line(aes(x=xgrid, y=Ef+2*sigma), linetype="dashed")
