
xvalALFAK <- function(alfak_dir = "~/projects/ALFA-K",simDir="data/main/"){
  setwd(alfak_dir)
  conditions <- list.files(simDir)
  ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
  library(parallel)
  cl <- makeCluster(getOption("cl.cores", nCores))
  clusterCall(cl, function() source("utils/ALFA-K.R"))
  clusterExport(cl = cl,c("simDir"),envir=environment())
  x <- do.call(rbind,parLapplyLB(cl=cl,X=conditions, fun=function(fi){
    xi <- tryCatch({
      print(fi)
      di <- paste0(simDir,fi,"/")
      sim_dir <- paste0(di,"train/")
      fit_dir <- paste0(di,"sweep_fits/")
      save_dir <- paste0(di,"loo_res/")
      dir.create(save_dir)
      sim_dir <- paste0(sim_dir,head(list.files(sim_dir),1),"/")
      
      fits <- list.files(fit_dir)
      xi <- do.call(rbind,lapply(fits, function(fit_i){
        ids <- as.numeric(unlist(strsplit(fit_i,split="_"))[c(2,4)])
        fit <- readRDS(paste0(fit_dir,fit_i))
        tt <- seq(100,2000,length.out=ids[2])
        x <- proc_sim(dir=sim_dir,times=tt)
        xfq <- fit$xo[fit$xo$id=="fq",]
        df <- do.call(rbind,lapply(1:nrow(xfq),function(i) optim_loo(i,x,xfq)))
        df$f_est <- xfq$f_est
        df$ntp=ids[2]
        df$minobs <- ids[1]
        df$rep <- fi
        return(df)
      }))
      saveRDS(xi,paste0(di,"loo_summaries.Rds"))
      return(xi)
    },error=function(e) return(NULL))
    return(xi)
  }))
  saveRDS(x,"data/proc/summaries/sweep_xval.Rds")
}




#x <- do.call(rbind,lapply(conditions, function(fi){


