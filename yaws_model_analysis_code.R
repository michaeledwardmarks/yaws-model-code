library('data.table')

res <- readRDS("/Users/Michael/Documents/Modelling Data/yaws_res.rds")

##Cut the data based on MDA parameters but not estimates of Beta##
  
extinctions <- res[, list(proportion.extinct = round(sum(extinct) / .N * 100)), by = c("tct1", "ttt1", "ttt2", "Rounds of TCT", "Rounds of TTT")]

extinctions[, extinct.groups := cut(proportion.extinct, seq(0, 100, 25), include.lowest = TRUE)]

setnames(extinctions,c("Rounds of TTT","Rounds of TCT"),c("ttt_rounds","tct_rounds"))

##Cut the data based on MDA parameters AND estimates of Beta##
extinctions2 <- res[, list(proportion.extinct = round(sum(extinct) / .N * 100)), by = c("tct1", "ttt1", "ttt2", "Rounds of TCT", "Rounds of TTT", "Beta")]

extinctions2[, extinct.groups := cut(proportion.extinct, seq(0, 100, 25), include.lowest = TRUE)]

setnames(extinctions2,c("Rounds of TTT","Rounds of TCT"),c("ttt_rounds","tct_rounds"))
