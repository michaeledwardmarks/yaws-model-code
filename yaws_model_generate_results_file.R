library(data.table)

# get list of files in data directory

files <- list.files(path.expand("/"), pattern = ".*\\.rds", full.names = TRUE)

filenames<- as.data.table(files)

integers <- as.data.table(apply(filenames, 1, function(x) as.integer(substr(x, 47, 48)))) ##character reference will need adjusting based on local file path

nb_events1 <- integers  %/% 6 + 1 ## divide index by 6, keep integer part, add 1
nb_events2 <- integers %% 6 ## divide index by 6, keep modulus 

mdatimes <- cbind(nb_events1,nb_events2)
colnames(mdatimes) <- c("Rounds of TCT","Rounds of TTT")

dt <- NULL

system.time(for (file in files[1:X])
{
  l <- readRDS(file)
  runs <- data.table(t(sapply(l, function(x) { x[["params"]] } )))
  
  runs[, extinct := sapply(l, function(x) {
    last.state <- tail(x[["trajectory"]], 1)
    return(ifelse(last.state[, I1 + I2 + L] == 0, 1, 0))
  })]
  par <- t(sapply(runs, function(x) { tail(x, n = 1) }))
  extinctcount <- t(sum(runs[["extinct"]]))
  colnames(extinctcount) <- c("extinctcount")
  runsum_temp <- merge (extinctcount,par)
  runsum_temp <- runsum_temp[,-(6:11),drop=FALSE] 

    if (is.null(dt))
  {
    dt <- copy(runsum_temp)
  } else
  {
    dt <- rbind(dt, runsum_temp)
  }
})

totaldata <- cbind(dt,mdatimes) 
