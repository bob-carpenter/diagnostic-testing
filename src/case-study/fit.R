library(rstan)
library(ggplot2)

sens_tests <- 3
sens_vals <-
c(78, 85,
  27, 37,
  25, 35)
sens_mat <- matrix(sens_vals, nrow = sens_tests, ncol = 2, byrow = TRUE)
sens_df <- data.frame(pos_tests = sens_mat[ , 1],
                      tests = sens_mat[ , 2],
                      sample_sens = sens_mat[ , 1] / sens_mat[ , 2])

spec_tests <- 13
spec_vals <-
c(368, 371,
   30,    30,
   70,    70,
  1102, 1102,
   300,  300,
   311,  311,
   500,  500,
   198,  200,
    99,   99,
    29,   31,
   146,  150,
   105,  108,
    50,   52)
spec_mat <- matrix(spec_vals, nrow = spec_tests, ncol = 2, byrow = TRUE)
spec_df <- data.frame(neg_tests = spec_mat[ , 1],
                      tests = spec_mat[ , 2],
                      sample_spec = spec_mat[ , 1] / spec_mat[ , 2])

unk_tests <- 1
unk_vals <-
c(50, 3330)
unk_mat <- matrix(unk_vals, nrow = unk_tests, ncol = 2, byrow = TRUE)
unk_df <- data.frame(pos_tests = unk_mat[ , 1],
                     tests = unk_mat[ , 2],
                     sample_prev = unk_mat[ , 1] / unk_mat[ , 2])

print("compiling model")
model <- stan_model("meta.stan")

data <-
list(K_pos = sens_tests,
     N_pos = sens_df$tests,
     n_pos = sens_df$pos_tests,
     K_neg = spec_tests,
     N_neg = spec_df$tests,
     n_neg = spec_df$neg_tests,
     K_unk = unk_tests,
     N_unk = array(unk_df$tests, dim = c(length(unk_df$tests))),
     n_unk = array(unk_df$pos_tests, dim = c(length(unk_df$tests))))

init_fun <- function(chain_id = 1) {
  list(pi = runif(1, 0.01, 0.02),
       mu_logit_sens = rnorm(1, 3, 0.125), mu_logit_spec = rnorm(1, 3, 0.125),
       sigma_logit_sens = runif(1, 0.4, 0.6), sigma_logit_spec = runif(1, 0.4, 0.6),
       logit_sens = array(rnorm(sens_tests, 3, 0.125)),
       logit_spec = array(rnorm(spec_tests, 3, 0.125)),
       logit_sens_unk = array(rnorm(unk_tests, 3, 0.125)),
       logit_spec_unk = array(rnorm(unk_tests, 3, 0.125)))
}

plot_df <- data.frame(sigma_sens = c(), sigma_spec = c(),
                      prevalence = c(), quantile = c())

sigma_senss <- c(0.1, 0.4, 0.7, 1.0)
sigma_specs <- c(0.1, 0.3, 0.5, 0.7, 0.9)
n_sigmas <- length(sigmas)
for (sigma_sens in sigma_senss) {
  for (sigma_spec in sigma_specs) {
    print(c(sigma_sens, sigma_spec))
    data2 <- append(data, list(sigma_sigma_logit_sens = sigma_sens,
                               sigma_sigma_logit_spec = sigma_spec))
    fit <- sampling(model, data = data2, init = init_fun,
                    iter = 10000, seed = 1234,
                    control = list(stepsize = 0.01, adapt_delta = 0.99),
                    open_progress = FALSE, refresh = 0)
    pis <- extract(fit)$pi
    for (q in c(0.05, 0.50, 0.95)) {
      plot_df <- rbind(plot_df, data.frame(sigma_sens = paste("sigma_sens", "=", sigma_sens),
                                           sigma_spec = sigma_spec,
                                           prevalence = quantile(pis, q),
                                           quantile = paste(q, "%", sep = "")))
    }
  }
}

ilogit <- function(v) 1 / (1 + exp(-v))
logit <- function(u) log(u / (1 - u))
brks <- c(1, 5, 9)
labs <- c("0.1", "0.5", "0.9")

plot <- ggplot(plot_df, aes(sigma_sens, sigma_spec, fill = prevalence)) +
  facet_wrap(~ quantile, ncol = 3) +
  coord_equal() +
  geom_tile(color="white", size=0.1) + #  geom_raster() +
  scale_x_continuous(expand = c(0, 0), breaks = brks, labels = labs) +
  scale_y_continuous(expand = c(0, 0), breaks = brks, labels = labs) +
  xlab("sigma_sigma_logit_sens") +
  ylab("sigma_sigma_logit_spec") +
  theme(axis.ticks = element_blank())

plot2 <- ggplot(plot_df, aes(x = sigma_spec, y = prevalence, color = quantile)) +
  facet_wrap(~ sigma_sens, nrow = 1) +
  geom_line(aes(color = quantile)) +
  scale_x_continuous(expand = c(0, 0), breaks = sigma_specs) +
  scale_y_log10(expand = c(0, 0), lim = c(0.001, 0.3), breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3))
