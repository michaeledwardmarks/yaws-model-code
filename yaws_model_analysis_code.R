library('data.table')

if (Sys.info()[["user"]] == "seb") {
    data_dir <- path.expand("~/Data/Yaws/12MonthFreq") # Seb
} else {
    data_dir <- path.expand("~/Documents/Modelling Data/") # Michael
}

res <- readRDS(paste(data_dir, "yaws_res.rds", sep = "/"))

##Cut the data based on MDA parameters but not estimates of Beta##
  
extinctions <- res[, list(proportion.extinct = round(sum(extinct) / .N * 100)), by = c("tct1", "ttt1", "ttt2", "Rounds of TCT", "Rounds of TTT")]

extinctions[, extinct.groups := cut(proportion.extinct, seq(0, 100, 25), include.lowest = TRUE)]

setnames(extinctions,c("Rounds of TTT","Rounds of TCT"),c("ttt_rounds","tct_rounds"))

##Cut the data based on MDA parameters AND estimates of Beta##
extinctions2 <- res[, list(proportion.extinct = round(sum(extinct) / .N * 100)), by = c("tct1", "ttt1", "ttt2", "Rounds of TCT", "Rounds of TTT", "Beta")]

extinctions2[, extinct.groups := cut(proportion.extinct, seq(0, 100, 25), include.lowest = TRUE)]

setnames(extinctions2,c("Rounds of TTT","Rounds of TCT"),c("ttt_rounds","tct_rounds"))

## calculate R0
N <-sum(c(S = 10736, I1 = 180, I2 = 180, L = 4996))


dt <- data.table(res)
dt[, R0 := Beta * N * (relapse + latenttreat) /
       ((latent + treat) * latenttreat + relapse * treat)]
dt[, mean(R0)]
## [1] 1.954078

dt[, quantile(R0, c(0.025, 0.5, 0.975))]
##     2.5%      50%    97.5% 
## 1.078005 1.879248 3.324015 

## work out R0 ranges by transmission scenario
Beta_values <- unique(dt[, Beta])

dt[Beta == Beta_values[1], transmission_scenario := "low"]
dt[Beta == Beta_values[2], transmission_scenario := "medium"]
dt[Beta == Beta_values[3], transmission_scenario := "high"]

dt[, list(mean = mean(R0),
          min.95 = quantile(R0, 0.025),
          max.95 = quantile(R0, 0.975)), by = transmission_scenario]

##    transmission_scenario     mean   min.95   max.95
## 1:                   low 1.434656 1.010635 2.137247
## 2:                medium 1.954074 1.376659 2.911574
## 3:                  high 2.473504 1.742236 3.685420
