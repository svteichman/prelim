base <- "bhq_res"

load(paste0(base, 1, ".rda"))
res <- bhq_res

for (i in 2:600) {
  load(paste0(base, i, ".rda"))
  res <- rbind(res, bhq_res)
}

save(res, file = paste0("combined_bhq.rda"))