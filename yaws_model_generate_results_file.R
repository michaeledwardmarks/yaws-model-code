library(data.table)

if (Sys.info()[["user"]] == "seb") {
    data_dir <- path.expand("~/Data/Yaws/") # Seb
} else {
    data_dir <- path.expand("~/Documents/Modelling Data/") # Michael
}

for (months in c(6, 12))
{
    data_dir_months <- paste0(data_dir, "/", months, "MonthFreq")

    ## get list of files in data directory
    files <- list.files(data_dir_months, pattern = ".*\\.rds", full.names = TRUE)

    pb <- txtProgressBar(min = 0, max = length(files), style = 3)

    system.time(for (file_idx in seq_along(files))
                {
                    l <- readRDS(paste(path, files[file_idx], sep = "/"))
                    runs <- data.table(t(sapply(l, function(x) { x[["params"]] } )))
                    runs[, extinct := sapply(l, function(x) {
                               last.state <- tail(x[["trajectory"]], 1)
                               return(ifelse(last.state[, I1 + I2 + L] == 0, 1, 0))
                           })]

                    integer <- sapply(files[file_idx], function(x) as.integer(substr(x, 10, 11)))
                    nb_events1 <- integer %/% 6 + 1 ## divide index by 6, keep integer part, add 1
                    nb_events2 <- integer %% 6 ## divide index by 6, keep modulus

                    mdatimes <- data.table(cbind(nb_events1, nb_events2))
                    setnames(mdatimes, c("Rounds of TCT","Rounds of TTT"))
                    runs <- cbind(runs, mdatimes)

                    if (is.null(dt))
                    {
                        dt <- copy(runs)
                    } else
                    {
                        dt <- rbind(dt, runs)
                    }
                    setTxtProgressBar(pb, file_idx)
                })

    close(pb)

    saveRDS(dt, paste0("yaws_res_", months, ".rds"))
}
