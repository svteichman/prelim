library(tidyverse)
load("simulations/bhq_sim1/combined_res.rda")

res1 <- res[2:601,]
n <- nrow(res1)
colMeans(res1)
apply(res1, 2, sd)/sqrt(n)

fdr_names <- names(res1)[1:7]
long_res <- res1 %>% 
  pivot_longer(cols = 1:14, names_to = "method") %>%
  mutate(type = ifelse(method %in% fdr_names, "FDR","power"),
         method = rep(c("Knockoff+ (equi)","Knockoff (equi)", "Knockoff+ (SDP)",
                        "Knockoff (SDP)","BH","BH with log correction",
                        "BH with white noise"), 1200))
long_res$method <- factor(long_res$method, levels = unique(long_res$method))

fdr_res <- long_res %>%
  filter(type == "FDR")
ggplot(fdr_res, aes(x = value)) + 
  geom_histogram() +
  facet_wrap(~method, scales = "free", nrow = 2) + 
  xlab("FDP") + ylab("Count") +
  theme_gray(base_size = 18)
ggsave("simulations/figs/sim1_fdr_hist.png",
       width = 12, height = 4)

power_res <- long_res %>%
  filter(type == "power")
ggplot(power_res, aes(x = value)) + 
  geom_histogram() +
  facet_wrap(~method, scales = "free", nrow = 2) + 
  xlab("Proportion true covariates detected") + ylab("Count") +
  theme_gray(base_size = 18)
ggsave("simulations/figs/sim1_power_hist.png",
       width = 12, height = 4)
