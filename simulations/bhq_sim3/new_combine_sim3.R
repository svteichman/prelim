base <- "BHq_sim3"

load(paste0(base, 1, ".rda"))
res <- bhq_sim3

for (i in 2:1000) {
  if (file.exists(paste0(base,i,".rda"))) {
    load(paste0(base, i, ".rda"))
    res <- rbind(res, bhq_sim3)
  }
}

base1 <- "extra_sim3"

for (i in 1:9) {
  if (file.exists(paste0(base1,i,".rda"))) {
    load(paste0(base1, i, ".rda"))
    res <- rbind(res, bhq_sim3)
  }
}

save(res, file = "combined_sim3.rda")