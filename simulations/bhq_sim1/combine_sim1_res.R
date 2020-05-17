base <- "BHq_sim1_res"

res <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
colnames(res) <- c("knock_plus_equi_fdr","knock_equi_fdr","knock_plus_sdp_fdr",
             "knock_sdp_fdr","bhq_base_fdr","bhq_log_fdr",      
             "bhq_white_noise_fdr","knock_plus_equi_pwr","knock_equi_pwr",     
             "knock_plus_sdp_pwr","knock_sdp_pwr","bhq_base_pwr","bhq_log_pwr",
             "bhq_white_noise_pwr")

for (i in 1:600) {
  load(paste0(base, i, ".rda"))
  res <- rbind(res, bhq_sim1_res)
}

save(res, file = paste0("combined_res.rda"))