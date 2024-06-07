setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")

library(ggplot2)

proc_data <- function(i,x,fit,ndist=2,max_size=30000){
  x <- x[!grepl("0",rownames(x)),]
  newk <- rownames(x)[x[,i+1]>0 & x[,i]==0]
  oldk <- rownames(x)[x[,i]>0]
  oldv <- do.call(rbind,lapply(oldk,s2v)) 
  nn <- oldk
  
  #newv <- do.call(rbind,lapply(newk,s2v))
  
  
  for(ni in ndist){
    nn <- sample(nn,min(max_size/44,length(nn)))
    nn <- apply(gen_all_neighbours(nn),1,paste,collapse=".")
  }
  
  vn <- do.call(rbind,lapply(nn,s2v))
  
  df <- apply(vn,1,function(ni){
    d <- apply(oldv,1,function(oi){
      sum(abs(oi-ni))
    })
    data.frame(keep = min(d)==ndist,
               nd =  sum(x[oldk[d==ndist],i]))
  })
  
  df <- do.call(rbind,df)
  df$f <- predict(fit,vn)
  df$y <- nn%in%newk
  df$n <- df$n/sum(x[,i])
  df$f <- df$f-sum(df$f*df$n)/sum(df$n)
  df[keep,]
}


wrap_test <- function(ff,ndist=2,max_size=30000){
  x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",ff$feval))$x
  fit <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",ff$ftrain))$fit
  
  test <- proc_data(ncol(x)-1,x=x,fit=fit,ndist,max_size)
  #test$n <- scale(test$n)
  #test$f <- scale(test$f)
  mod <- glm(y~f+n,data=test)
  res <- summary(mod)$coefficients
  ids <- rownames(res)
  res <- data.frame(res,row.names = NULL)
  res$ids <- ids
  return(res)
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



df <- pbapply::pblapply(ff,function(ffi)
  tryCatch(wrap_test(ffi,ndist=3),error=function(e) return(NULL))
)

df <- do.call(rbind,df)
df <- data.frame(df)
saveRDS(df,"figures/salehi_predictions/multireg_d3.Rds")
