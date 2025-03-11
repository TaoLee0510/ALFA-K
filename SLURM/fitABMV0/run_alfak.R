# Load necessary libraries
print("STARTING JOB")
options(error = function() { print("ERROR DETECTED"); traceback(); q(status=1) })
setwd("~/projects/ALFA-K/")

dataPath <- "data/proc/sweep_inputs.Rds"
outputBasePath <- "data/proc/sweep_results/"

source("utils/ALFA-KQ.R")

args <- commandArgs(trailingOnly=TRUE)
fi_id <- as.numeric(args[1])

y0 <- readRDS(dataPath)
fi <- names(y0)[fi_id]
print(fi_id)
print(fi)
y0 <- y0[[fi]]$data


minobs <- c(20,10,5)
ntp <- c(2,4,8)
pars <- expand.grid(ntp=ntp,minobs=minobs)
tmax <- 1200
pred_iters <- 100
n0 <- 1e5
nb <- 2e6
pm <- 0.00005
num_cores <- 50
nboot <- 50

proc_sweep_input <- function(yi,ntp=8,tmax=1200){
  yi$x <- yi$x[, as.numeric(colnames(yi$x)) < tmax]
  yi$x <- yi$x[, (ncol(yi$x) - ntp + 1):ncol(yi$x)]
  is.absent <- rowSums(yi$x)==0
  yi$clone.fitness <- yi$clone.fitness[!is.absent]
  yi$x <- yi$x[!is.absent,]
  names(yi$clone.fitness) <- row.names(yi$x)
  return(yi)
}

get_prediction_passages <- function(yi,tmax=1200,npassages=10){
  pred_passages <- as.numeric(colnames(yi$x))[as.numeric(colnames(yi$x)) >= tmax]
  pred_passages <- head(pred_passages,npassages)
  return(round(pred_passages*yi$dt))
}

for(i in 1:nrow(pars)){
  tryCatch({
    minobs <- pars$minobs[i]
    ntp <- pars$ntp[i]
    outdir <- paste0(outputBasePath,"minobs_",minobs,"_ntp_",ntp,"/",fi,"/")
    yi <- proc_sweep_input(y0,ntp=ntp,tmax=tmax)
    passage_times <- get_prediction_passages(y0)
    alfak(yi, outdir, passage_times, minobs = minobs,
          nboot = nboot,
          pred_iters = pred_iters,
          n0 = n0,
          nb = nb,
          pred_times = pred_times,
          pm = 0.00005,
          num_cores = num_cores)
  },error=function(e) print(e))
  
}