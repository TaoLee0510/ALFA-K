##Some of the data we have looked wrong and was excluded. 
## This function removes those elements from the metadata, so they don't show in the
## already cluttered network plots.  
## The premise is that these files are absent from the processed alfak_inputs folder
prune_children <- function(meta,dir="data/processed/salehi/alfak_inputs"){
  filename_patterns <- paste(meta$datasetname,meta$timepoint,sep="_")
  ff <- list.files(dir)
  c1 <- sapply(filename_patterns,function(fi) sum(grepl(fi,ff))>1) 
  ## the other reason these file patterns can be absent from the data is if they're
  ## the initial "root" passages. All these should be kept.
  ## in general just keeping all passages that have a parent is a sound choice:
  c2 <- sapply(meta$uid,function(idi) idi%in%meta$parent)
  meta[c1|c2,]
}