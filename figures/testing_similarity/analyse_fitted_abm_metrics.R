setwd("~/projects/ALFA-K")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

dir <- "data/main/"
conditions <- list.files(dir)
ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
conditions <- conditions[ids=="0.00005"]
dirs <- paste0(dir,conditions,"/")
dirs <- dirs[-c(479,481)] ##IDK why but this hangs! The wasserstein_metric() function call for wl_1p6_rep_78,
##di <- dirs[3], simdir <- dij[3], just never returns. I tryed to catch this with a call to withTImeout
## but the hang evidently occurs in C code so it doesnt work. Cursory inspection of the input data revealed nothing 
## too suspicious. 
dirs <- dirs[length(dirs):1]
wrapper <- function(d0){
  a <- tryCatch({
    print(d0)
    ## first process the training sim
    ff0 <- list.files(paste0(d0,"train/00000/"))
    ff0 <- ff0[!ff0%in%c("summary.txt","log.txt")]
    tt0 <- as.numeric(sapply(ff0,function(fi) unlist(strsplit(fi,split="[.]"))[1]))
    tt_end <- tt0[which.min(abs(tt0-2000))]
    tt0 <- tt0[tt0>=tt_end]
    x0 <- proc_sim(paste0(d0,"train/00000/"),times=c(min(tt0),2500,3000))
    #ntp_x0 <- ncol(x0$x)
    x_0 <- get_mean(x0,t=min(tt0))
    
    
    odir <- paste0(d0,"test_v2/")
    dirs <- paste0(odir,list.files(odir),"/")
    a <- do.call(rbind,lapply(dirs,function(di){
      fij <- list.files(di)
      dij <- paste0(di,fij,"/")
      a <- do.call(rbind,lapply(dij, function(simdir){
        ff <- list.files(simdir)
        ff <- ff[!ff%in%c("summary.txt","log.txt")]
        tt <- as.numeric(sapply(ff,function(fi) unlist(strsplit(fi,split="[.]"))[1]))
        x <- proc_sim(simdir,times=c(0,500,1000))
        colnames(x$x) <- as.numeric(colnames(x$x))+min(as.numeric(colnames(x0$x)))
        eval_times <- c(2500,3000)
        
        a <- sapply(eval_times, function(ti){
          angle_metric(x0,x,t=ti,x0 = x_0)
        })
        
        w <- sapply(eval_times, function(ti){
          wasserstein_metric(x0,x,t=ti)
        })
        
        names(a) <- paste0("a_",eval_times)
        names(w) <- paste0("w_",eval_times)
        return(c(a,w))
      }))
      
      a <- data.frame(a,check.names = F)
      a$id1 <- tail(unlist(strsplit(d0,split="/")),1)
      a$id2 <- tail(unlist(strsplit(di,split="/")),1)
      a$rep <- fij
      return(a)
    }))
    return(a)
  },error=function(e) return(NULL))
  
}

library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))
clusterCall(cl, function() {
  source("utils/ALFA-K.R")
  source("utils/comparison_functions.R")
  #library(R.utils)
  })
#source("utils/ALFA-K.R")
#source("utils/comparison_functions.R")
#a <- do.call(rbind,lapply(dirs,wrapper))
#saveRDS(a,"figures/testing_similarity/fitted_abm_run_metrics.Rds")
a <- do.call(rbind,parLapplyLB(cl=cl,X=dirs, fun=wrapper))
saveRDS(a,"figures/testing_similarity/fitted_abm_run_metrics.Rds")
