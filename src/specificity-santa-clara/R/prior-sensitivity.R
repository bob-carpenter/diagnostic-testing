library("ggplot2")
library("cmdstanr")

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

print("translating and compiling prior-sensitivity.stan")
model <- cmdstan_model("../stan/prior-sensitivity.stan")

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

ribbon_df <- data.frame(sigma_sens = c(), sigma_spec = c(),
                        prev05 = c(), prev50 = c(), prev95 = c())

sigma_senss <- c(0.01, 0.25, 0.5, 0.75, 1)
sigma_specss <- c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for (sigma_sens in sigma_senss) {
  for (sigma_spec in sigma_specss) {
    print(c(sigma_sens, sigma_spec))
    data2 <- append(data, list(sigma_sigma_logit_sens = sigma_sens,
                               sigma_sigma_logit_spec = sigma_spec))
    fit <-  model$sample(data = data2, 
                    iter_warmup = 5e4, iter_sampling = 5e4, seed = 1234, parallel_chains=4,
                    refresh = 0, adapt_delta=0.9)
    pis <- as.vector(fit$draws()[,,"pi"])
    ribbon_df <- rbind(ribbon_df,
                       data.frame(sigma_sens = paste("sensitivity hyperprior", "=", sigma_sens),
                                  sigma_spec = sigma_spec,
                                  prev05 = quantile(pis, 0.05),
                                  prev50 = quantile(pis, 0.5),
                                  prev95 = quantile(pis, 0.95)))
  }
}

# ribbon version
plot_ribbon <- ggplot(ribbon_df, aes(x = sigma_spec)) +
  facet_wrap(~ sigma_sens, nrow = 1) +
  geom_ribbon(aes(ymin = prev05, ymax = prev95), fill = "gray95") +
  geom_line(aes(y = prev50), size = 0.5) +
  geom_line(aes(y = prev05), color = "darkgray", size = 0.25) +
  geom_line(aes(y = prev95), color = "darkgray", size = 0.25) +
  scale_y_log10(limits = c(0.0009, 1.1), breaks = c(0.001, 0.01, 0.1, 1)) +
  scale_x_continuous(expand = c(0, 0), lim = c(0, 1),
                     breaks = c(0.1, 0.3, 0.5, 0.7, 0.9)) +
  ylab("prevalence") +
  xlab("specificity hyperprior") +
  theme_bw() +
  theme(panel.spacing = unit(0.25, "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

pdf("prior-sensitivity-2.pdf", width=9, height=2.5)
plot_ribbon
dev.off()
