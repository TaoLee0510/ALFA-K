library("optparse")

option_list = list(
  make_option(c("-c", "--cores"), type="numeric", default=1, 
              help="number of parallel cores [default= %default]"),
  make_option(c("-n", "--ntp"), type="character", default="2,3,4,8", 
              help="number of timepoints to sample[default= %default]"),
  make_option(c("-m", "--minobs"), type="character", default="5,10,20", 
              help="frequent karyotype minimal threshold [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


## https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

root.dir <- gsub("[\\]","/",thisFile())
root.dir <- unlist(strsplit(root.dir,split="/"))
root.dir <- root.dir[1:(length(root.dir)-3)]
root.dir <- paste0(paste(root.dir,collapse="/"),"/")
setwd(root.dir)

dir <- "data/main/"
conditions <- list.files(dir)
ids <- as.character(sapply(conditions, function(i) unlist(strsplit(i,split="_"))[6]))

min_obs <- as.numeric(unlist(strsplit(opt$min_obs,split=",")))
ntp <- as.numeric(unlist(strsplit(opt$ntp,split=",")))

conds <- expand.grid(min_obs=min_obs,ntp=ntp)

library(parallel)
cl <- makeCluster(getOption("cl.cores", opt$cores))
clusterCall(cl, function() source("utils/ALFA-K.R"))
clusterExport(cl = cl,c("dir","conds"))


x <- do.call(rbind,parLapplyLB(cl=cl,X=conditions, fun=function(fi){
  print(fi)
  di <- paste0(dir,fi,"/train/")
  fit_dir <- paste0(dir,fi,"/sweep_fits/")
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
