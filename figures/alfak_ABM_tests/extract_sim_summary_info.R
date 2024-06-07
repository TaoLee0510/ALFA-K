setwd("~/projects/ALFA-K")


dir <- "data/main/"
conditions <- list.files(dir)

x <- do.call(rbind,pbapply::pblapply(conditions, function(fi){
  print(fi)
  di <- paste0(dir,fi,"/train/")
  rep_id <- list.files(di)
  do.call(rbind,lapply(rep_id, function(rep_i){
    dijk <- paste0(di,rep_i)
    fijk <- list.files(dijk)
    fijk <- fijk[!fijk%in%c("log.txt","summary.txt")]
    tt <- as.numeric(sapply(fijk, function(ff) unlist(strsplit(ff,split="[.]"))[1]))
    fijk <- paste0(dijk,"/",fijk)
    xijk <- lapply(fijk,read.table,sep=",")
    f <- sapply(xijk, function(xx) sum(xx[,23]*xx[,24])/sum(xx[,23]))
    data.frame(time=tt,fitness=f,cond_id=fi,rep_id=rep_i)
  }))
}))

saveRDS(x,"figures/alfak_ABM_tests/data/sim_summary_info.Rds")