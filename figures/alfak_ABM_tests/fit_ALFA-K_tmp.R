setwd("~/projects/ALFA-K")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

indir <- "data/main/"
outdir <- "figures/alfak_ABM_tests/"

ff <- list.files(indir)

x <- sapply(ff, function(fi){
  if(!file.exists(paste0(indir,fi,"/fit_summaries.Rds"))) return(NULL)
  readRDS(paste0(indir,fi,"/fit_summaries.Rds"))
})

names(x) <- ff

saveRDS(x,paste0(outdir,"fit_summaries.Rds"))
