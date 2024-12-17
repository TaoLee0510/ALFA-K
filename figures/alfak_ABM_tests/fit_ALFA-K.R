setwd("~/projects/ALFA-K")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

dir <- "data/main/"
conditions <- list.files(dir)
ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))
#print(unique(ids))
conditions <- conditions[ids=="0.00005"]
#print(conditions)
conds <- expand.grid(min_obs=c(5,10,20),ntp=c(2,3,4,8))

library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))
clusterCall(cl, function() source("utils/ALFA-K.R"))
clusterExport(cl = cl,c("dir","conds"))


x <- do.call(rbind,parLapplyLB(cl=cl,X=conditions, fun=function(fi){
  print(fi)
  di <- paste0(dir,fi,"/train/")
  fit_dir <- paste0(dir,fi,"/sweep_fits_v2/")
  dir.create(fit_dir)
  rep_id <- list.files(di)
  
  ### I don't think it is necessary to do this for all replicates(!?)
  rep_id <- head(rep_id,1) 
  ###
  
  xi <- do.call(rbind,lapply(rep_id, function(rep_i){
    do.call(rbind,lapply(1:nrow(conds),function(ci){
      ntp <- conds$ntp[ci]
      min_obs <- conds$min_obs[ci]
      dijk <- paste0(di,rep_i)
      x <- proc_sim(dijk,times=seq(0,2000,length.out=ntp))
      fit <- alfak(x,min_obs = min_obs)
      fit_name <- paste0(fit_dir,"minobs_",min_obs,"_ntp_",ntp,"_",rep_i,".Rds")
      saveRDS(fit,fit_name)
      n2 <- gen_all_neighbours(rownames(fit$xo))
      f2 <- predict(fit$fit,n2)
      n2 <- apply(n2,1,paste,collapse=".")
      df2 <- data.frame(f_est=f2,id="n2",row.names = n2)
      df <- rbind(fit$xo[,c("f_est","id")],df2)
      
      fobj <- gen_fitness_object(cfig_path = paste0(dir,fi,"/config.txt"),lscape_path = paste0(dir,fi,"/landscape.txt"))
      f <- sapply(rownames(df), function(ki){
        getf(s2v(ki),fobj)
      })
      df$f <- f
      
      df <- split(df,f=df$id)
      
      res <- do.call(rbind,lapply(df,function(dfi){
        data.frame(id=dfi$id[1],r2=cor(dfi$f_est,dfi$f,use="complete.obs"))
      }))  
      res$cond_id <- fi
      res$rep_id <- rep_i
      res$ntp <- ntp
      res$min_obs <- min_obs
      res
    }))
  }))
  saveRDS(xi,paste0(dir,fi,"/fit_summaries.Rds"))
}))

#saveRDS(x,"figures/alfak_ABM_tests/fit_summaries.Rds")