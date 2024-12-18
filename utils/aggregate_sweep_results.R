## various functions for aggregating ABM sweep output.

## this function is used when we have data structure
##  inDir
##    subdir
##      target
## where target is a data frame.
aggregate_fit_summaries <- function(summaryName="fit_summaries.Rds",inDir="data/main/",outPath="data/proc/summaries/"){
  ## make sure directory exists
  dir.create(outPath,recursive = T,showWarnings = F)
  ff <- list.files(inDir)
  ff <- paste0(inDir,ff,"/",summaryName)
  ff <- ff[file.exists(ff)]
  x <- lapply(ff,readRDS)
  saveRDS(x,paste0(outPath,summaryName))
  return(x)
}