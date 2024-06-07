setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/ALFA-K.R")

library(caret)

proc_data <- function(i,x,fit){
  for(j in 1:ncol(x)) x[,j] <- x[,j]/sum(x[,j]) 
  newk <- rownames(x)[x[,i+1]>0 & x[,i]==0]
  oldk <- rownames(x)[x[,i]>0]
  nn0 <- gen_all_neighbours(newk)
  nn0 <- nn0[sample(1:nrow(nn0),length(newk)),]
  nn0 <- apply(nn0,1,paste,collapse=".")
  
#  nn <- gen_all_neighbours(newk)
#  nn <- nn[sample(1:nrow(nn),length(newk)),]
 # nn <- apply(nn,1,paste,collapse=".")
  
  k <- unique(c(nn0,newk))
  k <- k[!k%in%oldk]
  ktest <- do.call(rbind,lapply(k,s2v))
  
  y <- k%in%newk
  
  fi <- predict(fit,ktest)
  fi <- fi/median(fi)
  
  ki <- do.call(rbind,lapply(oldk,s2v))
  ni <- x[,i][oldk]
  
  dvec <- do.call(rbind,pbapply::pblapply(1:nrow(ktest),function(i){
    k <- ktest[i,]
    di <- apply(ki,1,function(xj){
      sum(abs(xj-k))
    })
    sapply(1:3, function(j) sum(ni[di==j]))/sum(ni)
  }))
  
  keepers <- rowSums(dvec)>0
  
  y <- cbind(y,fi,dvec)
  colnames(y) <- NULL
  y <- y[keepers,]
  return(y)
}

pfun <- function(pk,dvec){
  apply(dvec,1,function(di){
    fi <- di[1]
    1-prod((1-pk[-1])^(pk[1]*fi+di[-1]))
  })
}

logl <- function(pk,df,use_fitness=F){
  y <- as.logical(df[,1])
  dvec <- df[,-1]
  if(!use_fitness) dvec[,1] <- 1
  pk <- exp(pk)
  pk[1] <- pk[1]*20
  p <- pfun(pk,dvec)
  v <- -(sum(log(p[y]))+sum(log(1-p[!y])))
  if(!is.finite(v)) v <- 1e+10
 # print(v)
  #print(log(pk))
  return(v)
}

wrap_test <- function(ff){
  x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",ff$feval))$x
  l <- readRDS(paste0("data/salehi/alfak_fits/minobs_5/",ff$ftrain))$fit
  
  train <- do.call(rbind,lapply(1:(ncol(x)-2), proc_data,x=x,fit=l ))
  test <- proc_data(ncol(x)-1,x=x,fit=l)
  
  colnames(train) <- c("y","f","n1","n2","n3")
  colnames(test) <- colnames(train) 
  
  tr <- data.frame(train)
  tr$y <- factor(tr$y)
  tst <- data.frame(test)
  tst$y <- factor(tst$y)
  
  fitControl <- trainControl(method = "repeatedcv",number = 10, repeats = 10)
  
  cfit0 <- train(y ~ n1+n2+n3, data = tr, 
                 method = "treebag", 
                 trControl = fitControl)
  pred0 <- predict(cfit0,tst)
  pred0 <- mean(pred0==tst$y)
  cfit1 <- train(y ~ ., data = tr, 
                 method = "treebag", 
                 trControl = fitControl)
  
  pred1 <- predict(cfit1,tst)
  pred1 <- mean(pred1==tst$y)
  
  res <- data.frame(n = nrow(tst),
                    tree_no_fitness=pred0,tree_yes_fitness=pred1)
  print(res)
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

res <- list(df=df,ff=ff)

saveRDS(res,"figures/salehi_predictions/ini_test3.Rds")