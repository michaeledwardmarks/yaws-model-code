library(adaptivetau)
library (data.table)

# package for latin hypercube sampling
library (lhs)

# enable compiling of model code -- should speed up the simulations
library(compiler)
enableJIT(1)

# read command line parameter(s)

args <- commandArgs(trailingOnly = TRUE)
filename <- "yaws_runs"

## number of parameter slices
slices <- 5

if (length(args) > 0) {
    sim_nb <- as.integer(args[1])
} else {
    sim_nb <- 1
}

# number of months to use for burn-in
burnin <- 50

# number of months to run the model for
maxtime <- 150

# number of runs per parameter combination
nb_runs <- 1000

#function to run a Stochastic Simulation
stoSI1I2L <- function(y, parms, times, events)
{
    # define transitions
    SI1I2L_transitions <- list(
      c(S=-1,I1=1),# infection
      c(I1=-1,I2=1),# Progression from Primary Disease to Secondary Disease
      c(I1=-1,L=1),# Progression from Primary disease to Latent Disease
      c(I2=-1,L=1),# Progression from Secondary Disease to Latent Disease
      c(L=-1,I2=1),# Relapse from Latent Disease to Secondary Disease
      c(I1=-1,S=1),#Treatment of Primary Disease
      c(I2=-1,S=1),#Treatment of Secondary Disease
      c(L=-1,S=1)#Treatment of latent disease
      
      
      )
    # define rates
    SI1I2L_rateFunc <- function(x, params, t)
    {
      return(c(x["S"] * (params["Beta"]*x["I1"]+params["Beta"]*x["I2"]), 
# infection rate is complex of S, Population with primary and secondary disease and rate of infection - BETA
              +            params["second"]*x["I1"],#Rate of progression from Prim to Second
              +            params["latent"]*x["I1"], # Rate of progression from Prim to latent
              +            params["latent"]*x["I2"],#Rate of progression from Second to latent
              +            params["relapse"]*x["L"], # Rate of relapse from latent to Second
              +            params["treat"]*x["I1"], # Rate of treatment of primary disease
              +            params["treat"]*x["I2"], # Rate of treatment of secondary disease              
              +            params["latenttreat"]*x["L"])) # Rate of treatment of latent disease   
    }

    # we only want events that are before the specified end time
    events$time <- events$time[events$time < max(times)]

    # initialise variable to hold full trajectory
    full.trajectory <- NULL
    # initialise variable to hold current simulation starting time
    current.time <- times[1]
    # convert y to integer
    storage.mode(y) <- "integer"
    # simulate up to event times
    for (event in events[["time"]])
    {
        trajectory <- ssa.adaptivetau(init.values = y,
                                      transitions = SI1I2L_transitions,
                                      rateFunc = SI1I2L_rateFunc,
                                      params = parms,
                                      tf = event - current.time)
        # convert to data frame
        trajectory <- as.data.frame(trajectory)
        # set time to correct value (ssa.adaptivetau automatically starts time at 0)
        trajectory$time <- trajectory$time + current.time
        # apply event to last row of trajectory (except first column,
        # which holds time
        event.y <-
            do.call(events[["func"]],
                    list(t = unlist(trajectory[nrow(trajectory), 1]),
                         y = unlist(trajectory[nrow(trajectory), -1]),
                         parms = parms))
        # assign eventy to y and convert to integer, preserving names
        y[1:length(y)] <- event.y
        storage.mode(y) <- "integer"

        # remove last row (will be first row of next iteration
        trajectory <- trajectory[-nrow(trajectory), ]

        # add to full trajectory, if there is already something in there
        if (is.null(full.trajectory))
        {
            full.trajectory <- trajectory
        } else
        {
            full.trajectory <- rbind(full.trajectory, trajectory)
        }
        current.time <- event
    }

    # simulate from last event to final time
    trajectory <- ssa.adaptivetau(init.values = y,
                                  transitions = SI1I2L_transitions,
                                  rateFunc = SI1I2L_rateFunc,
                                  params = parms,
                                  tf = max(times) - current.time)
    # convert to data frame
    trajectory <- as.data.frame(trajectory)
    # set time to correct value (ssa.adaptivetau automatically starts time at 0)
    trajectory$time <- trajectory$time + current.time
    # add to full trajectory, if there is already something in there
    if (is.null(full.trajectory))
    {
        full.trajectory <- trajectory
    } else
    {
        full.trajectory <- rbind(full.trajectory, trajectory)
    }

    # keep only times that were originally asked for this loops over
    # all columns of "full.trajectory" (via apply), and for each colum
    # interpolates between the times in "full.trajectory" to get the
    # state of the system at the times given to the function as
    # "times" parameter
    result <- cbind(time = times,
                    apply(full.trajectory[, -1], 2,
                          function(col)
                          {
                              approx(x = full.trajectory[, 1],
                                     y = col,
                                     xout = times,
                                     method = "constant")$y}))

    return(result)
}

#define the transmission model

#define MDA Parameters
#Total Community Treatment
eventtct <- function(t, y, parms){
    with(as.list(c(parms,y)),{
        S <- S + tct1*L + tct1*I1 + tct1*I2
        I1 <- I1*(1-tct1)
        I2 <- I2*(1-tct1)
        L <- L*(1-tct1)
        return(c(S,I1,I2,L))
    })
}
#Total Targeted Treatment
#Effect is TTT of I1 and I2 and number of latent cases per clinical cases*efficacy
eventttt <- function(t, y, parms){
    with(as.list(c(parms,y)),{
        S <- S + ttt1*I1 + ttt1*I2 + ttt2*L
        I1 <- I1*(1-ttt1)
        I2 <- I2*(1-ttt1)
        L <-  L*(1-ttt2)
        return(c(S,I1,I2,L))

    })
}

#generate an empty list for the data to be put into

simulationlist_sto <- list()

#MDA Timings

run_nb <- (sim_nb - 1) %% 18
slice_nb <- (sim_nb - 1) %/% 18 + 1

etimes1 <- seq(from = 12, by = 6, length.out = nb_events1)
etimes2 <- seq(from = max(etimes1) + 6, by = 6, length.out = nb_events2)

allevents <- sort(unique(c(etimes1,etimes2)))

#Combine TCT and TTT into a single modelled event
dispatch <- function(t,y,parms){
    ret <- y
    if (t %in% etimes1) ret <- eventtct(t,y,parms)
    if (t %in% etimes2) ret <- eventttt(t,y,parms)
    return (ret)
}

#Set the time frame for the model
dt    <- seq(-burnin, maxtime, 1)

## plist <- cbind(tct1          = runif(nb_runs, min = 0.75, max = 0.95),
##                ttt1          = runif(nb_runs, min = 0.75, max = 0.95),
##                ttt2          = runif(nb_runs, min = 0.09, max = 0.95),
##                Beta          = runif(nb_runs, min = 0.00002762, max = 0.00004603))
# grid search on relevant parameters
param_grid <- expand.grid(tct1          = seq(from = 0.65, to = 0.95, by = 0.05),
                          ttt1          = seq(from = 0.65, to = 0.95, by = 0.05),
                          ttt2          = seq(from = 0.65, to = 0.95, by = 0.05),
                          Beta          = seq(from = 0.00002762, to = 0.00005, by = 0.00001))

## plist <- cbind(second        = runif(nb_runs, min = 0.0278,max = 0.0556),
##                latent        = runif(nb_runs, min = 0.139, max = 0.278),
##                treat         = runif(nb_runs, min = 0.12,  max = 0.30),
##                latenttreat   = runif(nb_runs, min = 0.008, max = 0.042),
##                relapse       = runif(nb_runs, min = 0.012, max = 0.028))
# randomly sample other parameters (latin hypercube sampling)
random_params <- list(second        = c(min = 0.0278, max = 0.0556),
                      latent        = c(min = 0.139, max = 0.278),
                      treat         = c(min = 0.12,  max = 0.30),
                      latenttreat   = c(min = 0.008, max = 0.042),
                      relapse       = c(min = 0.012, max = 0.028))
#Define the spectrum of inital values for population size
## initlist <- cbind(S        = runif(nb_runs, min = 10736, max = 10736),
##                   I1       = runif(nb_runs, min = 180, max=180),
##                   I2       = runif(nb_runs, min = 180, max = 180),
##                   L        = runif(nb_runs, min = 4996, max = 4996))

#fix the inital values for population size
init <- c(S        = 10736,
          I1       = 180,
          I2       = 180,
          L        = 4996)

start <- nrow(param_grid) %/% slices * (slice_nb - 1) + 1
if (slice_nb == slices)
{
    random_plist <- maximinLHS(n = nb_runs, k = 5)
    colnames(random_plist) <- names(random_params)
    plist <- sapply(names(random_params), function(name)
    {
        random_params[[name]][["min"]] +
            random_plist[, name] *
                (random_params[[name]][["max"]] - random_params[[name]][["min"]])
    })
    colnames(plist) <- names(random_params)
    simulationlist_sto[[i]] <- list()
    cat("Exploring parameter set", i, "at", format(Sys.time()), "\n")
    end <- nrow(param_grid)
} else {
    end <- nrow(param_grid) %/% slices * slice_nb
}

cat("Starting to explore parameter sets", start, "to", end, "\n")

for (i in start:end)
{
    {
        #run multiple simulations
        simulationlist_sto[[i]][[j]] <- list()
        simulationlist_sto[[i]][[j]][["params"]] <- c(unlist(param_grid[i, ]), plist[j,])
        simulationlist_sto[[i]][[j]][["trajectory"]] <-
            data.table(stoSI1I2L(y = init,
                                    parms = simulationlist_sto[[i]][[j]][["params"]],
                                    times = dt,
                                    events = list(func = dispatch,
                                    time = allevents)))
        #discard burn-in
        simulationlist_sto[[i]][[j]][["trajectory"]] <-
            simulationlist_sto[[i]][[j]][["trajectory"]][time >= 0]
    }
})

saveRDS(simulationlist_sto, paste0("sims_", run_nb, ".rds"))

