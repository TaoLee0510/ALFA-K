setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")
xa <- readRDS("example_data/hTERTa.Rds")
xb <- readRDS("example_data/hTERTb.Rds")
xcomb <- list(xa,xb)
hterta <- alfak(xa,min_obs = 10)
htertb <- alfak(xb,min_obs = 10)
htertcomb <- alfak2(xcomb,min_obs=10)
saveRDS(hterta,file="example_data/hTERTa_fit.Rds")
saveRDS(htertb,file="example_data/hTERTb_fit.Rds")
saveRDS(htertcomb,file="example_data/hTERTcomb_fit.Rds")

wrapLOO <- function(opt,x){
  xfq <- opt$xo[opt$xo$id=="fq",]
  ##leave one out cross validation - use ALFA-K to predict fitness of frequent clones left out of the training data.
  df <- do.call(rbind,lapply(1:nrow(xfq),function(i) optim_loo(i,x,xfq)))
  xfq$loo_pred <- df$pred
  return(xfq)
}

wrapLOO2 <- function(opt,x){
  xfq <- lapply(opt$xo, function(xo) xo[xo$id=="fq",])
  do.call(rbind,lapply(1:length(xfq), function(j){
    do.call(rbind,lapply(1:nrow(xfq[[j]]), function(i){
      optim_loo2(j,i,x,xfq)
    }))
  }))
  
}

x <- wrapLOO2(htertcomb,xcomb)
saveRDS(x,"figures/comparing_fit_quality/loo_xv/hTERTcomb.RDS")

x <- wrapLOO(hterta,xa)
saveRDS(x,"figures/comparing_fit_quality/loo_xv/hTERTa.RDS")

x <- wrapLOO(htertb,xb)
saveRDS(x,"figures/comparing_fit_quality/loo_xv/hTERTb.RDS")

## is there an issue preparing htertcomb for this procedure??

