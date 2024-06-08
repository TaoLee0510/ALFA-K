setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")

library(ggplot2)

proc_data <- function(i,x,fit){
  
  nn <- rownames(fit$xv_res)
  nn <- c(nn,apply(gen_all_neighbours(nn),1,paste,collapse="."))
  nn <- c(nn,apply(gen_all_neighbours(nn),1,paste,collapse="."))
  
  newk <- rownames(x)[x[,i+1]>0 & x[,i]==0]
  oldk <- rownames(x)[x[,i]>0]
  oldv <- do.call(rbind,lapply(oldk,s2v)) 
 
  vn <- do.call(rbind,lapply(nn,s2v))
  
  
  df <- apply(vn,1,function(ni){
    d <- apply(oldv,1,function(oi){
      sum(abs(oi-ni))
    })
    sapply(0:5,function(di){
      sum(x[oldk[d==di],i])
    })
  })
  df <- t(df)
  df <- df/sum(x[,i])
  df <- data.frame(df,row.names=NULL)
  colnames(df) <- paste0("d",0:5)
  df$f <- predict(fit$fit,vn)
  df$y <- nn%in%newk
  df <- df[df[,1]==0,-1]
  
  return(df)

}


wrap_test <- function(ff){
  x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",ff$feval))$x
  fit <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",ff$ftrain))
  
  test <- proc_data(ncol(x)-1,x=x,fit=fit)
  #test$n <- scale(test$n)
  #test$f <- scale(test$f)
  mod <- glm(y~.,data=test,family = binomial)
  res <- summary(mod)$coefficients
  ids <- rownames(res)
  res <- data.frame(res,row.names = NULL)
  res$ids <- ids
  return(res)
}

xx <- readRDS("figures/salehi_data_fitting/data/fit_summaries.Rds")
x0 <- xx[xx$r2>0.3&xx$dec1>0,]
lins <- readRDS("figures/salehi_data_fitting/data/lineages.Rds")

ff <- lapply(x0$filenames,function(fi) {
  id <- head(unlist(strsplit(fi,split=".Rds")),1)
  lii <- lins[[id]]
  dec <- lii$dec1[1]
  xi <- xx[xx$uid==dec,]
  fii <- xi$filenames[which.max(xi$samples)]
  
  list(ftrain = fi,
       feval= fii)
})



df <- pbapply::pblapply(ff,function(ffi)
  tryCatch(wrap_test(ffi),error=function(e) return(NULL))
)

df <- do.call(rbind,df)
df <- data.frame(df)
saveRDS(df,"figures/salehi_predictions/data/multireg_landscape.Rds")
