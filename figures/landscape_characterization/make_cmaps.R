setwd("~/projects/008_birthrateLandscape/ALFA-k")
library(parallel)
source("utils/landscape_functions.R")
source("utils/ALFA-K.R")
dir <- "data/salehi/alfak_fits_minobs_adaptive/"
ff <- list.files(dir)
get_best <- function(x){
  cx <- sapply(x, function(xi){
    if(is.null(xi)) return(-Inf)
    if("loo_pred"%in%colnames(xi$xv_res)){
      return(cor(xi$xv_res$f_est,xi$xv_res$loo_pred)*nrow(xi$xv_res))
    }
    if("f_xv"%in%colnames(xi$xv_res)){
      return(cor(xi$xv_res$f_est,xi$xv_res$f_xv)*nrow(xi$xv_res))
    }
    
  }) 
  return(x[[which.max(cx)]])
}
for(fi in ff){
  x <- get_best(readRDS(paste0(dir,fi)))
  
  xx <- rbind(do.call(rbind,lapply(rownames(x$xo), function(xi) s2v(xi))),
              gen_all_neighbours(rownames(x$xo)))
  
  xs <- apply(xx,1,paste,collapse="")
  
  f <- predict(x$fit,xx)
  cl <- makeCluster(getOption("cl.cores", 3))
  clusterExport(cl, c("f","xx"))
  
  cmap <- parLapplyLB(cl=cl,X = 1:nrow(xx), fun = function(i){
    nn <- which(apply(xx,1,function(xi) mean(xi==(2*xx[i,]))==1 | sum(abs(xi-xx[i,]))==1))
    #if(length(nn)==0) return(NULL)
    list(up=nn[f[i]<f[nn]],down=nn[f[i]>=f[nn]])
  })
  
  pathData <- list(cmap=cmap,xx=xx,f=f)
  saveRDS(pathData,file=paste0("figures/landscape_characterization/cmaps/",fi))
}