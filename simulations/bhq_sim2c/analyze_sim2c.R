library(tidyverse)
library(latex2exp)
load("simulations/bhq_sim2c/combined_sim2c.rda")

rho_seq <- seq(0, 0.9, by = .1)
res$rho_val <- rep(rho_seq, 200)
summ_res <- res %>% group_by(rho_val) %>%
  summarise(mean_fdr_kp = mean(knock_plus_sdp_fdr),
            mean_fdr_k = mean(knock_sdp_fdr),
            mean_fdr_bh = mean(bhq_base_fdr),
            mean_pwr_kp = mean(knock_plus_sdp_pwr),
            mean_pwr_k = mean(knock_sdp_pwr),
            mean_pwr_bh = mean(bhq_base_pwr))

fdr_res <- summ_res %>% pivot_longer(cols = c("mean_fdr_kp","mean_fdr_k","mean_fdr_bh"), 
                                     names_to = "method") %>%
  rename(fdr = value) %>%
  mutate(method = ifelse(method == "mean_fdr_kp", "Knockoff+", 
                         ifelse(method == "mean_fdr_k", "Knockoff","BHq")))
ggplot(fdr_res, aes(x = rho_val, y = fdr, group = method, color = method, shape = method)) + 
  geom_line() + geom_point() + geom_hline(yintercept=.20, linetype="dashed") +
  xlab(TeX("Feature correlation $\\rho$")) + ylab("FDR") + 
  scale_color_manual(values=c("black", "blue", "red")) + 
  scale_shape_manual(values = c("triangle","square","circle")) +
  scale_y_continuous(breaks = c(0,.05,.10,.15,.20,.25,.30), limits = c(0,.3)) + 
  scale_x_continuous(breaks = c(0,.2,.4,.6,.8), limits = c(0,.9)) + 
  theme_grey(base_size = 20) +
  theme(legend.position = c(0.25, 0.18), legend.title = element_blank())
ggsave("simulations/figs/fig5_fdr.png")

power_res <- summ_res %>% pivot_longer(cols = c("mean_pwr_kp","mean_pwr_k","mean_pwr_bh"), 
                                       names_to = "method") %>%
  rename(power = value) %>%
  mutate(method = ifelse(method == "mean_pwr_kp", "Knockoff+", 
                         ifelse(method == "mean_pwr_k", "Knockoff","BHq")))
ggplot(power_res, aes(x = rho_val, y = power, group = method, color = method, shape = method)) + 
  geom_line() + geom_point() + xlab(TeX("Feature correlation $\\rho$")) + ylab("Power") + 
  scale_color_manual(values=c("black", "blue", "red")) +
  scale_shape_manual(values = c("triangle","square","circle")) +
  scale_y_continuous(breaks = c(0,.2,.4,.6,.8,1), limits = c(0,1)) + 
  scale_x_continuous(breaks = c(0,.2,.4,.6,.8), limits = c(0,.9)) +
  theme_grey(base_size = 20) +
  theme(legend.position = c(0.75, 0.75), legend.title = element_blank())
ggsave("simulations/figs/fig5_power.png")
