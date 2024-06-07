setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/comparison_functions.R")
source("figures/salehi_data_fitting/scripts/lineage_processing.R")


m <- read.csv("data/salehi/metadata.csv")
lineages <- readRDS("figures/salehi_data_fitting/data/lineages.Rds")
linfo <- do.call(rbind,lapply(lineages,lineage_info))
linfo$filenames <- paste0(rownames(linfo),".Rds")
rownames(linfo) <- linfo$filenames
dir <- "data/salehi/alfak_fits/"
min_obs <- list.files(dir)

x <- pbapply::pblapply(min_obs, function(mo){
  ff <- list.files(paste0(dir,mo))
  ff <- ff[ff%in%linfo$filenames]
  res <- data.frame(do.call(rbind,lapply(ff,get_r2,mo=mo)))
  res <- cbind(res,linfo[ff,])
  
})

x <- data.frame(do.call(rbind,x))
rownames(x) <- NULL
x <- x[,!colnames(x)%in%c("datasetname","timepoint")]
x$r2 <- as.numeric(x$r2)
x$maxf <- as.numeric(x$maxf)

x$r2[x$r2<(-1)] <- -1
x <- x[order(x$r2,decreasing=T),]

saveRDS(x,"figures/salehi_data_fitting/data/fit_summaries.Rds")