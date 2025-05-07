#source("SLURM/fitSalehiV1/sweep_metadata.R")
## that will gives us filepaths to load the fits and prediction data
source("utils/ALFA-KQ.R")

proc_data <- function(lscape,x){
  

  ## karyotypes appearing at next timepoint
  newk <- rownames(x)[x[,ncol(x)]>0 & x[,ncol(x)-1]==0]
  ## karyotypes at current timepoint
  oldk <- rownames(x)[x[,ncol(x)]>0]
  oldv <- do.call(rbind,lapply(oldk,s2v)) 
  
  vn <- do.call(rbind,lapply(lscape$k,s2v))
  
  ## feature creating
  df <- apply(vn,1,function(ni){
    d <- apply(oldv,1,function(oi){
      sum(abs(oi-ni))
    })
    sapply(0:5,function(di){
      sum(x[oldk[d==di],ncol(x)-1])
    })
  })
  df <- t(df)
  df <- df/sum(x[,ncol(x)-1])
  df <- data.frame(df,row.names=NULL)
  colnames(df) <- paste0("d",0:5)
  df$f <- lscape$mean
  df$y <- lscape$k%in%newk
  df <- df[df[,1]==0,-1]
  
  return(df)
  
}

df_meta <- readRDS("data/salehi/alfak_outputs_V1a_prediction_metadata.Rds")

xv_thresh <- 0.0
df_meta <- df_meta[df_meta$xv>xv_thresh,]

df <- do.call(rbind,pbapply::pblapply(1:nrow(df_meta),function(i){
  res <- tryCatch({
    lscape <- readRDS(paste0("data/salehi/alfak_outputs_V1a/",df_meta$min_obs[i],"/",df_meta$fi[i],"/landscape.Rds"))
    x <- readRDS(paste0("data/salehi/alfak_inputs/",df_meta$dec_fi[i]))$x
    
    test <- proc_data(lscape,x)
    mod <- glm(y~.,data=test,family = binomial)
    res <- summary(mod)$coefficients
    ids <- rownames(res)
    res <- data.frame(res,row.names = NULL)
    res$ids <- ids
    res <- merge(res,df_meta[i,])
    res
  },error=function(e) return(NULL))
  return(res)
}))


i <- 4
lscape <- readRDS(paste0("data/salehi/alfak_outputs_V1a/",df_meta$min_obs[i],"/",df_meta$fi[i],"/landscape.Rds"))
x <- readRDS(paste0("data/salehi/alfak_inputs/",df_meta$dec_fi[i]))$x
test <- proc_data(lscape,x)
i1 <- unlist(test[which.max(test$d1),1:5])
i2 <- unlist(test[which.min(test$d1),1:5])

dfxmpl <- data.frame(frac=c(i1,i2),d=c(names(i1),names(i2)),
                 id=c(rep("1",length(i1)),rep("2",length(i2))))

res <- list(df=df,dfxmpl=dfxmpl)
saveRDS(res,"figures/misc/data/salehi_novel_kary_plotData.Rds")


