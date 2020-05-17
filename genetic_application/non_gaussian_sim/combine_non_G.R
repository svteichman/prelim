base <- "non_G_sim"

load(paste0(base, 1, ".rda"))
res <- sim_res

for (i in 2:1000) {
  load(paste0(base, i, ".rda"))
  res <- rbind(res, sim_res)
}

save(res, file = paste0("combined_non_G.rda"))