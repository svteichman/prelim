base <- "BHq_sim2c"

load(paste0(base, 1, ".rda"))
res <- bhq_sim2c

for (i in 2:200) {
  load(paste0(base, i, ".rda"))
  res <- rbind(res, bhq_sim2c)
}

save(res, file = "combined_sim2c.rda")