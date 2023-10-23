target_dir <- "data/main/"
ncores <- 3
library(parallel)
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
script_dir <- paste0(root.dir,"utils/")
source(paste0(script_dir,"comparison_functions.R"))
source(paste0(script_dir,"ALFA-K.R"))
ff <- list.files(target_dir)
df <- expand.grid(test=ff,ref=ff)
df <- df[sample(1:nrow(df),10000),]
df_list <- lapply(1:nrow(df), function(i) df[i,])

cl <- makeCluster(getOption("cl.cores", ncores))
clusterCall(cl, function(script_dir){
  source(paste0(script_dir,"comparison_functions.R"))
  source(paste0(script_dir,"ALFA-K.R"))
},script_dir=script_dir)

eval_metrics <- function(dirs,t_eval=1200){
  test_dir <- paste0(dir1,list.files(dirs[1])[1])
  ref_dirs <- paste0(dir2,list.files(dirs[2])[-1])
  test <- proc_sim(test_dir,times = seq(0,t_eval,t_eval))
  ref <- lapply(ref_dirs,proc_sim,times=seq(0,t_eval,t_eval))
  data.frame(ll=ll_cna(test,ref,t=t_eval),
             angle=angle_metric(test,ref,t=t_eval,is.ref.multiple.objects = T),
             dwass=wasserstein_distance(test,ref,t=t_eval,is.ref.multiple.objects = T))
}

df2<- do.call(rbind,parLapplyLB(cl=cl,x=df_list,eval_metrics))
df <- cbind(df,df2)
saveRDS(df,"figures/testing_similarity/metrics.Rds")

