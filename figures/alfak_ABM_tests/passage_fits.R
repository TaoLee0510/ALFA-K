mainDir <- "data/main"
simDirs <- list.files(mainDir,full.names = T)

library(parallel)

cl <- makeCluster(getOption("cl.cores", 70))

parLapplyLB(cl=cl,simDirs,function(di){
  subDir <- "train/00000"
  outDir <- "alfak_eq_space_fits"
  
  
  ntp <- c(2,3,4,8)
  min_obs <- c(5,10,20)
  
  get_tt <- function(path){
    tt <- list.files(path)
    tt <- gsub(".csv","",tt)
    tt <- as.numeric(tt[!tt%in%c("log.txt","summary.txt")])
  }
  
  inDir <- file.path(di,subDir)
  times <- get_tt(inDir)
  
  times <- times[times<1200]
  source("utils/ALFA-K.R")
  dir.create(file.path(di,outDir))
  lapply(ntp,function(n){
    lapply(min_obs, function(m){
      tryCatch({
        tt <- tail(times,n)
        x <- proc_sim(inDir,tt)
        fit <- alfak(x,min_obs = m)
        xfq <- fit$xo[fit$xo$id=="fq",]
        df <- do.call(rbind,lapply(1:nrow(xfq),function(i) optim_loo(i,x,xfq)))
        df$f_est <- xfq$f_est
        fobj <- gen_fitness_object(cfig_path = paste0(di,"/config.txt"),lscape_path = paste0(di,"/landscape.txt"))
        
        results <- list(input=x,fit=fit,xval=df,fobj=fobj)
        result_fname <- paste0("minobs_",m,"_ntp_",n,".Rds")
        savePath <- file.path(di,outDir,result_fname)
        saveRDS(results,savePath)
      },error=function(e) print(di))
      
    })
  })
})
  
  






