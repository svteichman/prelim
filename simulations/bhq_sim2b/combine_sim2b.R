base <- "BHq_sim2b"

load(paste0(base, 1, ".rda"))
res <- bhq_sim2b

for (i in 2:200) {
  load(paste0(base, i, ".rda"))
  res <- rbind(res, bhq_sim2b)
}

save(res, file = "combined_sim2b.rda")