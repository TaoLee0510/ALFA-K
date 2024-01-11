setwd("~/projects/ALFA-K")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

dir <- "data/main/"
conditions <- list.files(dir)
ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
conditions <- conditions[ids=="0.00005"]
dirs <- paste0(dir,conditions,"/")

analyse_sims_in_dir <- function(d0){
  ## in future this function shouldn't be here (should be in comparison functions, and be CORRECT as it is here)
  angle_metric <- function(test,ref,t=1200,is.test.multiple.objects=F,is.ref.multiple.objects=F,x0=2){
    xt <- get_mean(test,t,is.multiple.objects = is.test.multiple.objects)-x0
    xr <- get_mean(ref,t,is.multiple.objects = is.ref.multiple.objects)-x0
    getangle(xt,xr)
  }
  a <- tryCatch({
    ff0 <- list.files(paste0(d0,"train/00000/"))
    ff0 <- ff0[!ff0%in%c("summary.txt","log.txt")]
    tt0 <- as.numeric(sapply(ff0,function(fi) unlist(strsplit(fi,split="[.]"))[1]))
    tt_end <- tt0[which.min(abs(tt0-2000))]
    tt0 <- tt0[tt0>=tt_end]
    x0 <- proc_sim(paste0(d0,"train/00000/"),times=tt0)
    ntp_x0 <- ncol(x0$x)
    x_0 <- get_mean(x0,t=min(tt0))
    odir <- paste0(d0,"test/")
    dirs <- paste0(odir,list.files(odir),"/")
    a <- do.call(rbind,lapply(dirs,function(di){
      fij <- list.files(di)
      dij <- paste0(di,fij,"/")
      a <- do.call(rbind,lapply(dij, function(simdir){
        ff <- list.files(simdir)
        ff <- ff[!ff%in%c("summary.txt","log.txt")]
        tt <- as.numeric(sapply(ff,function(fi) unlist(strsplit(fi,split="[.]"))[1]))
        x <- proc_sim(simdir,times=tt)
        
        ## I think instead of using "raw" time, it is better to compare similarity on a passage
        ## by passage basis. That is what we are organizing here:
        ntp_x <- ncol(x$x)
        if(ntp_x>ntp_x0) x$x <- x$x[,1:ncol(x0$x)]
        colnames(x$x) <- colnames(x0$x)[1:ncol(x$x)]
        a <- sapply(as.numeric(colnames(x$x))[-1], function(ti){
          angle_metric(x0,x,t=ti,x0 = x_0)
        })[1:5]
        
        return(a)
      }))
      colnames(a) <- paste0("p",1:ncol(a))
      a <- data.frame(a,check.names = F)
      a <- cbind(a,rep=fij)
      ids <- unlist(strsplit(di,split="/"))
      a$id <- ids[grepl("N_22_w",ids)]
      a$condition <- ids[grepl("minobs_",ids)]
      return(a)
    }))
    return(a)
  },error=function(e) return(NULL))
  return(a)
}

library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))
clusterCall(cl, function() {
  source("utils/ALFA-K.R")
  source("utils/comparison_functions.R")
  })
a <- do.call(rbind,parLapplyLB(cl=cl,X=dirs, fun=analyse_sims_in_dir))
saveRDS(a,"figures/alfak_ABM_tests/fitted_abm_run_angles.Rds")
