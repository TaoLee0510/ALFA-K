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
  n2$y <- as.numeric(n2$k%in%newk)
  n2$y <- n2$k%in%newk
 # n2$y[n2$y==1] <- x[n2$k[n2$y==1],i+1]/sum(x[,i+1])
  
  n2$n <- n2$n/sum(x[,i])
  n2$f <- n2$f-sum(n2$f*n2$n)/sum(n2$n)
  return(n2)
}


wrap_test <- function(ff){
  x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",ff$feval))$x
  fit <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",ff$ftrain))$fit
  
  test <- proc_data(ncol(x)-1,x=x,fit=fit)
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



df <- lapply(ff,function(ffi)
  tryCatch(wrap_test(ffi),error=function(e) return(NULL))
)

df <- do.call(rbind,df)
df <- data.frame(df)
saveRDS(df,"figures/salehi_predictions/multireg.Rds")
