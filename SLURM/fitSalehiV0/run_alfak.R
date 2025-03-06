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

# Parameters
pred_iters <- 100
minobs <- 10
nboot <- 100
n0 <- 1e5
nb <- 1e7
pm <- 0.00005
num_cores <- 50
stop("random error")
# Load data
yi <- readRDS(paste0(inDir, fi))
ix <- strsplit(fi, split = "_") |> unlist()
datasetname <- ix[1]
passage <- ix[2]
uid <- m$uid[m$timepoint == passage & m$datasetname == datasetname]
pred_passages <- c()
while (sum(uid %in% m$parent) > 0) {
  pred_passages <- c(pred_passages, unique(m$timepoint[m$parent %in% uid & !is.na(m$parent)]))
  uid <- m$uid[m$parent %in% uid & !is.na(m$parent)]
}
pred_times <- as.numeric(gsub("X", "", pred_passages))
pred_times <- yi$dt * pred_times

passage_times <- min(as.numeric(colnames(yi$x))):max(as.numeric(colnames(yi$x)))
outdir <- paste0("data/proc/alfak_outputs/minobs_",minobs,"/", gsub(".Rds", "", fi))
print("entering alfak")
# Run ALFA-K
alfak(yi, outdir, passage_times, minobs = minobs,
      nboot = nboot,
      pred_iters = pred_iters,
      n0 = n0,
      nb = nb,
      pred_times = pred_times,
      pm = 0.00005,
      num_cores = num_cores)
