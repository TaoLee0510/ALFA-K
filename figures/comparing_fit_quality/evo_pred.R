setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")
source("utils/sim_setup_functions.R")

ids <- c("hTERTa","hTERTb")
Nreps <- 25

for(id in ids){
  x <- readRDS(paste0("example_data/",id,".Rds"))
  x$x <- x$x[,1:(ncol(x$x)-1)]
  fit <- alfak(x)
  
  n0 <- x$x[,ncol(x$x)]
  n0 <- n0[n0>2]
  k <- do.call(rbind,lapply(names(n0),s2v))
  n0 <- round(n0*100000/sum(n0))
  pop <- cbind(k,n0)
  
  sim_dir <- paste0("figures/comparing_fit_quality/evo_pred/",id)
  pars=c(pop_write_freq=100)
  cmd <- setup_abm(sim_dir,pars,pop,fit$fit)
  
  runcmd <- paste("ABM/bin/ABM.exe",cmd)
  for(j in 1:Nreps) system(runcmd)
  
}

x <- readRDS("example_data/hTERTcomb.Rds")
for(i in 1:length(x)){
  x[[i]]$x <- x[[i]]$x[,1:(ncol(x[[i]]$x)-1)]
}
fit <- alfak2(x)

for(i in 1:length(x)){
  n0 <- x[[i]]$x[,ncol(x[[i]]$x)]
  n0 <- n0[n0>2]
  k <- do.call(rbind,lapply(names(n0),s2v))
  n0 <- round(n0*100000/sum(n0))
  pop <- cbind(k,n0)
  
  sim_dir <- paste0("figures/comparing_fit_quality/evo_pred/hTERTcomb_",i)
  pars=c(pop_write_freq=100)
  cmd <- setup_abm(sim_dir,pars,pop,fit$fit)
  
  runcmd <- paste("ABM/bin/ABM.exe",cmd)
  for(j in 1:Nreps) system(runcmd)
}






