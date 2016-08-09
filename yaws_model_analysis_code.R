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

## define transmission scenarios in terms of Beta
## Beta_values <- unique(dt[, Beta])

## dt[Beta == Beta_values[1], transmission_scenario := "low"]
## dt[Beta == Beta_values[2], transmission_scenario := "medium"]
## dt[Beta == Beta_values[3], transmission_scenario := "high"]

## dt[, list(mean = mean(R0),
##           min.95 = quantile(R0, 0.025),
##           max.95 = quantile(R0, 0.975)), by = transmission_scenario]

##    transmission_scenario     mean   min.95   max.95
## 1:                   low 1.434656 1.010635 2.137247
## 2:                medium 1.954074 1.376659 2.911574
## 3:                  high 2.473504 1.742236 3.685420

######### define transmission scenarios in terms of R0

## exclude R0 values <= 1 (123226 in total)
dt <- dt[R0 > 1]

## determine low, mid, high third of R0 values
R0_quantiles <- quantile(dt[, R0], seq(0, 1, 1/3))
 ##      0% 33.33333% 66.66667%      100% 
 ## 1.000001  1.635413  2.148715  5.026966 

## split data into low, mid, high according to these quantiles
dt[, R0_quantile := cut(R0, R0_quantiles, include.lowest = TRUE)]

dt[, transmission_scenario :=
       factor(R0_quantile, labels = c("low", "medium", "high"))]
##Cut the data based on MDA parameters and defined transmission scenarios
extinctions3 <- dt[, list(proportion.extinct = round(sum(extinct) / .N * 100)), by = c("tct1", "ttt1", "ttt2", "Rounds of TCT", "Rounds of TTT", "transmission_scenario")]

extinctions3[, extinct.groups := cut(proportion.extinct, seq(0, 100, 25), include.lowest = TRUE)]

setnames(extinctions3,c("Rounds of TTT","Rounds of TCT"),c("ttt_rounds","tct_rounds"))

library('ppcor')

## reorder columns
res <- res[, list(Beta, second, latent, relapse, treat, latenttreat, tct1, ttt1, ttt2, `Rounds of TCT`, `Rounds of TTT`, extinct)]
## calculate PRCC
corr <- pcor(res, method = "spearman")
## create data frame for easier plotting
prcc <- data.table(parameter = colnames(corr$estimate),
                   PRCC = corr$estimate[, "extinct"])
## remove self-correlation
prcc <- prcc[parameter != "extinct"]
prcc[, parameter := factor(parameter, levels = prcc$parameter)]

library('ggplot2')
library('cowplot')

p <- ggplot(prcc[parameter %in% c("Beta", "second", "latent", "relapse", "treat", "latenttreat")], 
            aes(x = parameter, y = PRCC)) +
  geom_bar(stat = "identity") +
  scale_x_discrete("", labels = c('Beta' = expression(beta),
                                  'latent' = expression(eta), 
                                  'latenttreat' = expression(tau[L]),
                                  'relapse' = expression(rho),
                                  'second' = expression(alpha),
                                  'treat' = expression(tau[I]),
                                  'tct1' = "TCT coverage",
                                  'ttt1' = "TTT coverage of active cases",
                                  'ttt2' = "TTT coverage of latent cases")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("prcc_disease_parameters.pdf", p)

p <- ggplot(prcc, aes(x = parameter, y = PRCC)) +
  geom_bar(stat = "identity") +
  scale_x_discrete("", labels = c('Beta' = expression(beta),
                                  'latent' = expression(eta), 
                                  'latenttreat' = expression(tau[L]),
                                  'relapse' = expression(rho),
                                  'second' = expression(alpha),
                                  'treat' = expression(tau[I]),
                                  'tct1' = "TCT coverage",
                                  'ttt1' = "TTT coverage of active cases",
                                  'ttt2' = "TTT coverage of latent cases")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("prcc_all_parameters.pdf", p)
