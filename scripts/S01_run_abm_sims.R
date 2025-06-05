setwd("~/projects/ALFA-K/")


setup_and_run_abm <- function(rep_id,wavelength,outDir){
  library(alfakR)
  
  gen_randscape <- function(founder,Nwaves,scalef=NULL,wavelength=0.8){
    if(is.null(scalef)) scalef <- 1/(pi*sqrt(Nwaves))
    f0 <- 0
    ## we want to have simulations where the diploid founder is fit enough to survive.
    ## if scalef <- 1/(pi*sqrt(30))
    ## then this would make the diploid cell in the top 10% fittest
    ## clones:
    while(f0<=0.4){
      pk <- lapply(1:Nwaves, function(i){
        pk <- sample((-10):20,length(founder),replace=T)
      })
      
      d <- sapply(pk,function(ci) {
        sqrt(sum((founder-ci)^2))
      })
      f0=sum(sin(d/wavelength)*scalef)
      
    }
    do.call(rbind,pk)
  }
  
  resample_sim <- function(sim, n_samples) {
    mat <- t(as.matrix(sim[,-1]))
    colnames(mat) <- sim$time
    mat <- apply(mat, 2, function(p) rmultinom(1, n_samples, prob = p / sum(p)))
    keep <- rowSums(mat) > 0
    rownames(mat) <- colnames(sim)[-1]
    mat <- mat[keep, , drop = FALSE]
    mat <- mat[order(rowSums(mat),decreasing = T),]
  }
  
  
  founder <- rep(2,22)
  Nwaves <- 10
  
  l <- gen_randscape(founder,Nwaves,wavelength = wavelength)
  times <- c(0,300)
  x0=c(1)
  names(x0) <- paste(founder,collapse=".")
  pmis <- 0.00005
  sim <- run_abm_simulation_grf(
    centroids = l, lambda = wavelength,
    p = pmis, times = times, x0 = x0,
    abm_pop_size = 5e4,abm_max_pop = 2e6, abm_delta_t = 0.1,
    abm_culling_survival = 0.01,
    abm_record_interval = -1, abm_seed = 42,normalize_freq = F
  )
  
  sim_rs <- resample_sim(sim,1000)
  yi <- list(x=data.frame(sim_rs,check.names = F),dt=1)
  out <- list(abm_output=yi,true_landscape=l)
  
  filePath <- paste0(outDir,
                     "w_",gsub("[.]","p",wavelength),
                     "_m_",gsub("[.]","p",pmis),
                     "_rep_",stringr::str_pad(rep_id,width = 2,pad = 0),".Rds")
  saveRDS(out,filePath)
  return(0)
}

outDir <- "data/raw/ABM/"
dir.create(outDir,showWarnings = F,recursive = TRUE)

reps <- 1:100
w <- c(0.2,0.4,0.8,1.6)

library(parallel)

# Define parameter grid
param_grid <- expand.grid(rep_id = reps, wavelength = w, KEEP.OUT.ATTRS = FALSE)

# Convert to list of argument lists
param_list <- split(param_grid, seq_len(nrow(param_grid)))

# Create cluster
n_cores <- 50
cl <- makeCluster(n_cores)

# Export necessary variables and load packages
clusterExport(cl, varlist = c("setup_and_run_abm", "outDir"))
clusterEvalQ(cl, library(alfakR))

# Run in parallel
parLapplyLB(cl, param_list, function(params) {
  setup_and_run_abm(params$rep_id, params$wavelength, outDir)
})

# Stop cluster
stopCluster(cl)

