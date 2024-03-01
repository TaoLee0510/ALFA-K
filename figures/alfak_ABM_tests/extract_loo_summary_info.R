setwd("~/projects/ALFA-K")


dir <- "data/main/"
conditions <- list.files(dir)

x <- do.call(rbind,pbapply::pblapply(conditions, function(fi){
  print(fi)
  di <- paste0(dir,fi,"/loo_summaries.Rds")
  if(!file.exists(di)) return(NULL)
  readRDS(di)
}))

saveRDS(x,"figures/alfak_ABM_tests/loo_summary_info.Rds")