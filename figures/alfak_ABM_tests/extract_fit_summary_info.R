setwd("~/projects/ALFA-K")
#setwd("~/projects/008_birthrateLandscape/ALFA-K")

dir <- "data/main/"
conditions <- list.files(dir)


library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))
clusterCall(cl, function() source("utils/ALFA-K.R"))
clusterExport(cl = cl,c("dir"))

x <- do.call(rbind,parLapplyLB(cl=cl,X=conditions, fun=function(fi){
  
  fobj <- gen_fitness_object(cfig_path = paste0(dir,fi,"/config.txt"),lscape_path = paste0(dir,fi,"/landscape.txt"))
  fit_dir <- paste0(dir,fi,"/sweep_fits/")
  fits <- list.files(fit_dir)
  
  do.call(rbind,lapply(fits, function(fij){
    fit <- readRDS(paste0(fit_dir,fij))
    k <- do.call(rbind,lapply(rownames(fit$xo),s2v))
    f <- as.numeric(apply(k,1,getf,fobj=fobj))
    cor_fq <- cor(fit$xo$f_est[fit$xo$id=="fq"],f[fit$xo$id=="fq"],use="complete")
    n_fq <- length(fit$xo$f_est[fit$xo$id=="fq"])
    cor_nn <- cor(fit$xo$f_est[fit$xo$id=="nn"],f[fit$xo$id=="nn"],use="complete")
    n_nn <- length(fit$xo$f_est[fit$xo$id=="nn"])
    n2 <- gen_all_neighbours(rownames(fit$xo))
    f <- as.numeric(apply(n2,1,getf,fobj=fobj))
    fpred <- predict(fit$fit,n2)
    
    cor_n2 <- cor(f,fpred,use="complete")
    n_n2 <- length(f)
    data.frame(cor=c(cor_fq,cor_nn,cor_n2),n=c(n_fq,n_nn,n_n2),id=c("fq","nn","n2"),
               id1=fi,id2=fij)
    
  }))
  
}))

saveRDS(x,"figures/alfak_ABM_tests/data/fit_summary_info.Rds")