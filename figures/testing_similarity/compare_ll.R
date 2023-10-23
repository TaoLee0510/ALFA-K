target_dir <- "data/main/"

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

script_dir <- paste0(root.dir,"utils/")
source(paste0(script_dir,"comparison_functions.R"))
source(paste0(script_dir,"ALFA-K.R"))

df <- expand.grid(test=ff,ref=ff)
df$ll <- pbapply::pbsapply(1:nrow(df),function(i){
  dir1 <- paste0(target_dir,df$test[i],"/")
  dir2 <- paste0(target_dir,df$ref[i],"/")
  test_dir <- paste0(dir1,list.files(dir1)[1])
  ref_dirs <- paste0(dir2,list.files(dir2)[-1])
  
  test <- proc_sim(test_dir,times = seq(0,1200,1200))
  ref <- lapply(ref_dirs,proc_sim,times=seq(0,1200,1200))
  
  ll_cna(test,ref,t=1200)
})
saveRDS(df,"figures/testing_similarity/ll_cn.Rds")

