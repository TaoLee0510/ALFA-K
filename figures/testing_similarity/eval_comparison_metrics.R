target_dir <- "data/main/"
ncores <- 4
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

df_list <- lapply(1:nrow(df), function(i) c(df$test[i],df$ref[i]))

cl <- makeCluster(getOption("cl.cores", ncores))
clusterCall(cl, function(script_dir){
  source(paste0(script_dir,"comparison_functions.R"))
  source(paste0(script_dir,"ALFA-K.R"))
},script_dir=script_dir)

eval_metrics <- function(dirs,target_dir,t_eval=1200){
  out <- tryCatch({
  dirs <- paste0(target_dir,dirs,"/train/")
#  print(dirs)
  test_dir <- paste0(dirs[1],list.files(dirs[1])[1])
  ref_dirs <- paste0(dirs[2],list.files(dirs[2])[-1])
  test <- proc_sim(test_dir,times = seq(0,t_eval,t_eval))
  ref <- lapply(ref_dirs,proc_sim,times=seq(0,t_eval,t_eval))
  data.frame(ll=ll_cna(test,ref,t=t_eval),
             angle=angle_metric(test,ref,t=t_eval,is.ref.multiple.objects = T),
             dwass=wasserstein_distance(test,ref,t=t_eval,is.ref.multiple.objects = T))
},error=function(e) return(data.frame(ll=NaN,angle=NaN,dwass=NaN)))
return(out)
}

df2<- do.call(rbind,parLapplyLB(cl=cl,X=df_list,eval_metrics,target_dir=target_dir))
df <- cbind(df,df2)
saveRDS(df,"figures/testing_similarity/metrics.Rds")

