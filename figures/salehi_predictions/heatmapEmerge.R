setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")

library(ggplot2)

proc_data <- function(i,x,fit){
  x <- x[!grepl("0",rownames(x)),]
  newk <- rownames(x)[x[,i+1]>0 & x[,i]==0]
  oldk <- rownames(x)[x[,i]>0]
  
  nn <- do.call(rbind,pbapply::pblapply(oldk,function(ki){
    nni <- gen_all_neighbours(ki)
    fi <- predict(fit,nni)
    nni <- apply(nni,1,paste,collapse=".")
    data.frame(k=nni,n=x[ki,i],f=fi,row.names = NULL)
  }))
  
  n2 <- aggregate(list(n=nn$n),by=list(k=nn$k,f=nn$f),sum)
  n2$y <- n2$k%in%newk
  n2$n <- n2$n/sum(x[,i])
  n2$f <- scale(n2$f) 
  return(n2)
}

agg_df <- function(df){
  qtn <- quantile(df$n,probs=seq(0,1,1/3))
  qtf <- quantile(df$f,probs=seq(0,1,1/3))
  
  df$qtn <- sapply(df$n,function(i){
    which.max(diff(i<=qtn))
  })
  
  df$qtf <- sapply(df$f,function(i){
    which.max(diff(i<=qtf))
  })
  
  dfagg <- aggregate(list(y=df$y),by=list(n=df$qtn,
                                          f=df$qtf),mean)
  return(dfagg)
}

wrap_test <- function(ff){
  x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",ff$feval))$x
  fit <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",ff$ftrain))$fit
  
  test <- proc_data(ncol(x)-1,x=x,fit=fit)
  dfagg <- agg_df(test)
  
  
  
  
  
  p <- ggplot(dfagg,aes(x=n,y=f,fill=y))+
    geom_raster()+
    scale_fill_viridis_c(trans="log")+
    ggtitle(ff$feval)
  savename <- head(unlist(strsplit(ff$ftrain,split=".Rds")),1)
  ggsave(paste0("figures/salehi_predictions/heatmaps/",savename,".png"),plot=p)
  return(0)
  
}

xx <- readRDS("figures/salehi_data_fitting/fit_summaries.Rds")
x0 <- xx[xx$r2>0.3&xx$dec1>0,]
lins <- readRDS("figures/salehi_data_fitting/lineages.Rds")

ff <- lapply(x0$filenames,function(fi) {
  id <- head(unlist(strsplit(fi,split=".Rds")),1)
  lii <- lins[[id]]
  dec <- lii$dec1[1]
  xi <- xx[xx$uid==dec,]
  fii <- xi$filenames[which.max(xi$samples)]
  
  list(ftrain = fi,
       feval= fii)
})



df <- lapply(ff,function(ffi)
  tryCatch(wrap_test(ffi),error=function(e) return(NULL))
)
