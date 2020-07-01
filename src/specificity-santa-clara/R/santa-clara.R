library("cmdstanr")
library("rstan")
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
spin <- function(x, lower=NULL, upper=NULL, conf=0.95){
  x <- sort(as.vector(x))
  if (!is.null(lower)) {
    if (lower > min(x)) stop("lower bound is not lower than all the data")
    else x <- c(lower, x)
  }
  if (!is.null(upper)) {
    if (upper < max(x)) stop("upper bound is not higher than all the data")
    else x <- c(x, upper)
  }
  n <- length(x)
  gap <- round(conf*n)
  width <- x[(gap+1):n] - x[1:(n-gap)]
  index <- min(which(width==min(width)))
  x[c(index, index + gap)]
}

## Simple model fit using data from Bendavid et al. paper of 11 Apr 2020
sc_model <- cmdstan_model("stan/santa-clara.stan")

fit_1 <- sc_model$sample(data = list(y_sample=50, n_sample=3330, y_spec=369+30, n_spec=371+30, y_sens=25+78, n_sens=37+85), refresh=0, parallel_chains=4, iter_warmup=1e4, iter_sampling=1e4)
print(stanfit(fit_1), digits=3)

# Inference for the population prevalance
draws_1 <- fit_1$draws()
subset <- sample(1e4, 1e3)
x <- as.vector(draws_1[subset,,"spec"])
y <- as.vector(draws_1[subset,,"p"])

pdf("scatter.pdf", height=3.5, width=4.5)
par(mar=c(3,3,0,1), mgp=c(2, .7, 0), tck=-.02)
plot(x, y, xlim=c(min(x), 1), ylim=c(0, max(y)), xaxs="i", yaxs="i", xlab=expression(paste("Specificity, ", gamma)), ylab=expression(paste("Prevalence, ", pi)), bty="l", pch=20, cex=.3)
dev.off()

pdf("hist.pdf", height=3.5, width=5.5)
par(mar=c(3,3,0,1), mgp=c(2, .7, 0), tck=-.02)
hist(y, yaxt="n", yaxs="i", xlab=expression(paste("Prevalence, ", pi)), ylab="", main="")
dev.off()

# Use the shortest posterior interval, which makes more sense than a central interval because of the skewness of the posterior and the hard boundary at 0
print(spin(draws_1[,,"p"], lower=0, upper=1, conf=0.95))

# Simple model fit using pooled specificity and sensitivity data from Bendavid et al. paper of 27 Apr 2020

fit_2 <- sc_model$sample(data = list(y_sample=50, n_sample=3330, y_spec=3308, n_spec=3324, y_sens=130, n_sens=157), refresh=0, parallel_chains=4, iter_warmup=1e4, iter_sampling=1e4)

print(stanfit(fit_2), digits=3)
draws_2 <- fit_2$draws()
print(spin(draws_2[,,"p"], lower=0, upper=1, conf=0.95))

## Hierarchical model allowing sensitivity and specificity to vary across studies, fit using data from Bendavid et al. paper of 27 Apr 2020

sc_model_hierarchical <- cmdstan_model("stan/santa-clara-hierarchical.stan")

santaclara_data = list(y_sample=50, n_sample=3330, J_spec=14, y_spec=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50), n_spec=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52), J_sens=4, y_sens=c(0, 78, 27, 25), n_sens=c(0, 85, 37, 35), logit_spec_prior_scale=1, logit_sens_prior_scale=1)

fit_3a<- sc_model_hierarchical$sample(data=santaclara_data, refresh=0, parallel_chains=4, iter_warmup=1e4, iter_sampling=1e4)

print(stanfit(fit_3a), digits=3)
draws_3a <- fit_3a$draws()
print(spin(draws_3a[,,"p"], lower=0, upper=1, conf=0.95))

# Fit again with stronger priors

santaclara_data$logit_spec_prior_scale <- 0.3
santaclara_data$logit_sens_prior_scale <- 0.3
fit_3b <- sc_model_hierarchical$sample(data=santaclara_data, refresh=0, parallel_chains=4, iter_warmup=1e4, iter_sampling=1e4)

print(stanfit(fit_3b), digits=3)
draws_3b <- fit_3b$draws()
print(spin(draws_3b[,,"p"], lower=0, upper=1, conf=0.95))

# Summarize the fits

print(stanfit(fit_3a), pars=c("p", "spec[1]", "sens[1]", "mu_logit_spec", "mu_logit_sens", "sigma_logit_spec", "sigma_logit_sens"), digits=3)
print(stanfit(fit_3b), pars=c("p", "spec[1]", "sens[1]", "mu_logit_spec", "mu_logit_sens", "sigma_logit_spec", "sigma_logit_sens"), digits=3)


print(spin(draws_3a[,,"p"], lower=0, upper=1, conf=0.95))
print(spin(draws_3a[,,"spec[1]"], lower=0, upper=1, conf=0.95))
print(spin(draws_3a[,,"sens[1]"], lower=0, upper=1, conf=0.95))
print(spin(draws_3a[,,"mu_logit_spec"],  conf=0.95))
print(spin(draws_3a[,,"mu_logit_sens"], conf=0.95))
print(spin(draws_3a[,,"sigma_logit_spec"],  conf=0.95))
print(spin(draws_3a[,,"sigma_logit_sens"], conf=0.95))

print(spin(draws_3b[,,"p"], lower=0, upper=1, conf=0.95))
print(spin(draws_3b[,,"spec[1]"], lower=0, upper=1, conf=0.95))
print(spin(draws_3b[,,"sens[1]"], lower=0, upper=1, conf=0.95))
print(spin(draws_3b[,,"mu_logit_spec"],  conf=0.95))
print(spin(draws_3b[,,"mu_logit_sens"], conf=0.95))
print(spin(draws_3b[,,"sigma_logit_spec"],  conf=0.95))
print(spin(draws_3b[,,"sigma_logit_sens"], conf=0.95))

## MRP model, and allowing prevalence to vary by sex, ethnicity, age category, and zip code.  Model is set up to use the ethnicity, age, and zip categories of Bendavid et al. (2020).

sc_model_hierarchical_mrp <- cmdstan_model("stan/santa-clara-hierarchical-mrp.stan")

# To fit the model, we need individual-level data.  These data are not publcly available, so just to get the program running, we take the existing 50 positive tests and assign them at random to the 3330 people.
N <- 3330
y <- sample(rep(c(0, 1), c(3330 - 50, 50)))
n <- rep(1, 3330)

# Here are the counts of each sex, ethnicity, and age from Bendavid et al. (2020).  We don't have zip code distribution but we looked it up and there are 58 zip codes in Santa Clara County; for simplicity we asssume all zip codes are equally likely.  We then assign these traits to people at random.  This is wrong--actually, these variable are correlated in various ways--but, again,. now we have fake data we can use to fit the model.
male <- sample(rep(c(0,1), c(2101, 1229)))
eth <- sample(rep(1:4, c(2118, 623, 266, 306+17)))
age <- sample(rep(1:4, c(71, 550, 2542, 167)))
N_zip <- 58
zip <- sample(1:N_zip, 3330, replace=TRUE)

# Setting up the zip code level predictor.  In this case we will use a random number with mean 50 and standard deviation 20.  These are arbitrary numbers that we chose just to be able to test the centering and scaling in the model.   In real life we might use %Latino or average income in the zip code
x_zip <- rnorm(N_zip, 50, 20)

# Setting up the poststratification table.  For simplicity we assume there are 1000 people in each cell in the county.  Actually we'd want data from the Census.
J <- 2*4*4*N_zip
N_pop <- rep(NA, J)
count = 1
for (i_zip in 1:N_zip){
  for (i_age in 1:4){
    for (i_eth in 1:4){
      for (i_male in 0:1){
        N_pop[count] <- 1000
        count = count + 1
      }
    }
  }
}

# Put togther the data and fit the model
santaclara_mrp_data <- list(N=N, y=y, male=male, eth=eth, age=age, zip=zip, N_zip=N_zip, x_zip=x_zip,  J_spec=14, y_spec=c(0, 368, 30, 70, 1102, 300, 311, 500, 198, 99, 29, 146, 105, 50), n_spec=c(0, 371, 30, 70, 1102, 300, 311, 500, 200, 99, 31, 150, 108, 52), J_sens=4, y_sens=c(0, 78, 27, 25), n_sens=c(0, 85, 37, 35), logit_spec_prior_scale=0.3, logit_sens_prior_scale=0.3, coef_prior_scale=0.5, J=J, N_pop=N_pop)

fit_4 <- sc_model_hierarchical_mrp$sample(data=santaclara_mrp_data, refresh=0, parallel_chains=4)

# Show inferences for some model parameters. along with p_avg, the population prevalence, also we look at the inferences for the first three poststratification cells just to check that everything makes sense
print(stanfit(fit_4), pars=c("p_avg", "b", "a_age", "a_eth", "sigma_eth", "sigma_age", "sigma_zip", "mu_logit_spec", "sigma_logit_spec",  "mu_logit_sens", "sigma_logit_sens", "p_pop[1]", "p_pop[2]", "p_pop[3]"), digits=3)
draws_4 <- fit_4$draws()
print(spin(draws_4[,,"p_avg"], lower=0, upper=1, conf=0.95))
