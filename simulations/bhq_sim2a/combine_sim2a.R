base <- "BHq_sim2a_res"

load(paste0(base, 1, ".rda"))
res <- bhq_sim2_res

for (i in 2:200) {
  load(paste0(base, i, ".rda"))
  res <- rbind(res, bhq_sim2_res)
}

save(res, file = paste0("combined_sim2a.rda"))