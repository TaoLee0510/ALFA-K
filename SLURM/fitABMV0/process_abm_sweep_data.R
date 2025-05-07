setwd("~/projects/ALFA-K")
yin <- readRDS("data/proc/sweep_inputs.Rds")
source("utils/ALFA-KQ.R")
dir <- "data/proc/sweep_results/"
conditions <- list.files(dir)
library(transport)

summarize_results <- function(fi,yi,cj){
  
  repinfo <- strsplit(fi,split="_")|>unlist()
  w <- repinfo[which(repinfo=="w")+1]
  abmrep <- repinfo[which(repinfo=="rep")+1]
  
  fitinfo <- strsplit(cj,split="_")|>unlist()
  minobs <- fitinfo[which(fitinfo=="minobs")+1]
  ntp <- fitinfo[which(fitinfo=="ntp")+1]
  
  info <- data.frame(w,abmrep,minobs,ntp)
  
  fij <- paste0(dir,cj,"/",fi,"/")
 
  if(!file.exists(paste0(fij,"landscape.Rds"))) return(list(info=info))
  lscape <- readRDS(paste0(fij,"landscape.Rds"))
  if(!file.exists(paste0(fij,"xval.Rds"))) return(list(info=info))
  Rxv <- readRDS(paste0(fij,"xval.Rds"))
  k <- do.call(rbind,lapply(lscape$k,s2v))
  
  ## is in the old alfa-k.R
  getf <- function(pk,fobj){
    d <- apply(fobj$peaks,1,function(pki) sqrt(sum((pk-pki)^2)))
    sum(sin(d/fobj$wavelength)*fobj$scale)
  }
  
  
  f <- apply(k,1,function(ki) getf(ki,yi$lscape))
  
  
  summary <- data.frame(r=cor(lscape$mean,f),
                        rfq=cor(lscape$mean[lscape$fq],f[lscape$fq]),
                        rnn=cor(lscape$mean[lscape$nn],f[lscape$nn]),
                        rd2=cor(lscape$mean[!(lscape$fq|lscape$nn)],f[!(lscape$fq|lscape$nn)]),
                        rho=cor(lscape$mean,f,method = "spearman"),
                        rhofq=cor(lscape$mean[lscape$fq],f[lscape$fq],method = "spearman"),
                        rhonn=cor(lscape$mean[lscape$nn],f[lscape$nn],method = "spearman"),
                        rhod2=cor(lscape$mean[!(lscape$fq|lscape$nn)],f[!(lscape$fq|lscape$nn)],method = "spearman"),
                        R=R2R(f,lscape$mean),
                        Rfq=R2R(f[lscape$fq],lscape$mean[lscape$fq]),
                        Rnn=R2R(f[lscape$nn],lscape$mean[lscape$nn]),
                        Rd2=R2R(f[!(lscape$fq|lscape$nn)],lscape$mean[!(lscape$fq|lscape$nn)]),
                        Rxv=Rxv,
                        nfq=sum(lscape$fq))
  
  if(!file.exists(paste0(fij,"predictions.Rds"))) return(list(info=info,summary=summary))
  preds <- readRDS(paste0(fij,"predictions.Rds"))
  for(i in 1:nrow(preds)) {
    preds[i,] <- pmax(0,preds[i,])
    preds[i,] <- preds[i,]/sum(preds[i,])
  }
  
  y <- preds*0
  
  measure_tps <- head(which(as.numeric(colnames(yi$data$x))>=1200),10)
  
  rx <- rownames(yi$data$x)[rownames(yi$data$x)%in%lscape$k]
  
  for(i in 1:length(measure_tps)){
    y[i,rx] <- yi$data$x[rx,measure_tps[i]]/sum(yi$data$x[rx,measure_tps[i]])
  }
  
  y0 <- preds[1,]*0 
  y0[rx] <- yi$data$x[rx,min(measure_tps)-1]/sum(yi$data$x[rx,min(measure_tps)-1])
  
  k0 <- colSums(y0*k)
  
  mag <- function(v) sqrt(sum(v^2))
  getangle <- function(a,b) 180*acos(sum(a*b)/(mag(a)*mag(b)))/pi
  df_a <- do.call(rbind,lapply(1:nrow(y),function(i){
    ky <- colSums(y[i,]*k)
    kp <- colSums(preds[i,]*k)
    
    data.frame(angle=getangle(ky-k0,kp-k0))
    
  }))
  rownames(df_a) <- colnames(yi$data$x)[measure_tps]
  
  df_de <- do.call(rbind,lapply(1:nrow(y),function(i){
    ky <- colSums(y[i,]*k)
    kp <- colSums(preds[i,]*k)
    
    data.frame(de_p=sqrt(sum((kp-ky)^2)),
               de_0=sqrt(sum((k0-ky)^2)))
    
  }))
  rownames(df_de) <- colnames(yi$data$x)[measure_tps]
  
  
  
  s0 <- wpp(k,mass=y0)
  df_dw <- do.call(rbind,lapply(1:nrow(y),function(i){

    tryCatch({
      sy <- wpp(k,mass=y[i,])
      sp <- wpp(k,mass=preds[i,])
      data.frame(dw_p=wasserstein(sy,sp),
                 dw_0=wasserstein(sy,s0))
    },error=function(e){
      data.frame(dw_p=NaN,
                 dw_0=NaN)
    })
    
    
  }))
  rownames(df_dw) <- colnames(yi$data$x)[measure_tps]
  
  df_sc <- do.call(rbind,lapply(1:nrow(y),function(i){
    sc_0 <- sum(y0*y[i,])/prod(sqrt(sum(y0^2)) , sqrt(sum(y[i,]^2)))
    sc_p <- sum(preds[i,]*y[i,])/prod(sqrt(sum(preds[i,]^2)) , sqrt(sum(y[i,]^2)))
    data.frame(sc_p,sc_0)
  }))
  rownames(df_sc) <- colnames(yi$data$x)[measure_tps]
  df_so <- do.call(rbind,lapply(1:nrow(y),function(i){
    so_0 <- sum(pmin(y0, y[i,]))/min(sum(y0), sum(y[i,]))
    so_p <- sum(pmin(preds[i,], y[i,]))/min(sum(preds[i,]), sum(y[i,]))
    data.frame(so_p,so_0)
  }))
  rownames(df_so) <- colnames(yi$data$x)[measure_tps]
  
  
  result <- list(info=info,summary=summary,df_so=df_so,df_sc=df_sc,df_de=df_de,df_dw=df_dw,df_a=df_a)
  saveRDS(result,paste0(fij,"result_summary.Rds"))
  return(result)
}

library(parallel)

cl <- makeCluster(50)

clusterEvalQ(cl, {
  source("utils/ALFA-KQ.R")
  library(transport)
})

clusterExport(cl, varlist = c("dir","conditions","summarize_results"), envir = .GlobalEnv)


ynames <- names(yin)
for(i in 1:length(yin)) yin[[i]]$fi <- ynames[i]




x <- parLapplyLB(cl,yin,fun=function(yi){
  xi <- lapply(conditions,function(cj){
    fi <- yi$fi
    folder_path <- paste0(dir,cj,"/",fi)
    if(!file.exists(folder_path)) return(NULL)
    tryCatch({
      summarize_results(fi,yi,cj)
    },error=function(e) return(NULL))
    
  })
  names(xi) <- conditions
})

stopCluster(cl)

