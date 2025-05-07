# Load necessary libraries
print("STARTING JOB")
options(error = function() { print("ERROR DETECTED"); traceback(); q(status=1) })
setwd("~/projects/ALFA-K/")
inDir <- "data/salehi/alfak_inputs/"
m <- read.csv("data/salehi/metadata.csv")

source("utils/ALFA-KQ.R")

args <- commandArgs(trailingOnly=TRUE)
fi <- args[1]
fi <- tail(unlist(strsplit(fi,split="/")),1)

print(paste("running ",fi))

for(minobs in c(20,10,5)){
  # Parameters
  pred_iters <- 100
  nboot <- 50
  n0 <- 1e5
  nb <- 1e7
  pm <- 0.00005
  num_cores <- NULL
  # Load data
  yi <- readRDS(paste0(inDir, fi))
  ix <- strsplit(fi, split = "_") |> unlist()
  datasetname <- ix[1]
  passage <- ix[2]
  if(grepl("SA535",fi)) {
    datasetname <- paste(ix[1:3],collapse="_")
    passage <- ix[4]
  }
  uid <- m$uid[m$timepoint == passage & m$datasetname == datasetname]
  pred_passages <- c()
  while (sum(uid %in% m$parent) > 0) {
    pred_passages <- c(pred_passages, unique(m$timepoint[m$parent %in% uid & !is.na(m$parent)]))
    uid <- m$uid[m$parent %in% uid & !is.na(m$parent)]
  }
  pred_times <- as.numeric(gsub("X", "", pred_passages))
  pred_times <- yi$dt * pred_times
  passage_times <- min(as.numeric(colnames(yi$x))):max(as.numeric(colnames(yi$x)))
  
  if(length(pred_times)>0) pred_times <- max(passage_times*yi$dt):(max(pred_times))
  
  outdir <- paste0("data/salehi/alfak_outputs_V1a/minobs_",minobs,"/", gsub(".Rds", "", fi))
  print("entering alfak")
  
  tryCatch({
    # Run ALFA-K
    alfak(yi, outdir, passage_times, minobs = minobs,
          nboot = nboot,
          pred_iters = pred_iters,
          n0 = n0,
          nb = nb,
          pred_times = pred_times,
          pm = 0.00005,
          num_cores = num_cores)
  },error=function(e) return(""))

}


