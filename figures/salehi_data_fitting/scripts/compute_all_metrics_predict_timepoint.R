xv5 <- readRDS("data/proc/summaries/salehi_loo_minobs_5.Rds")
xv10 <- readRDS("data/proc/summaries/salehi_loo_minobs_10.Rds")
xv20 <- readRDS("data/proc/summaries/salehi_loo_minobs_20.Rds")
colnames(xv5)[2] <- "R2_5"
colnames(xv10)[2] <- "R2_10"
colnames(xv20)[2] <- "R2_20"
df <- merge(xv5,xv10,all=T)
df <- merge(df,xv20,all=T)

df[is.na(df)]<- (-Inf)

df$best <- c(5,10,20)[apply(df[,2:4],1,which.max)]
df$R2 <- pmax(apply(df[,2:4],1,max),-1)

df <- data.frame(id=gsub(".Rds","",df$file),Rsq=df$R2,best_minobs=df$best)

source("figures/salehi_data_fitting/scripts/process_salehi_data.R")

m <- read.csv("data/salehi/metadata.csv")
lineages <- process_lineages(m)

df <- do.call(rbind,lapply(1:nrow(df),function(i){
  li <- lineages[[df$id[i]]]
  do.call(rbind,lapply(li$dec1,function(dec_i){
    r <- df[i,]
    r$parent <- tail(li$ids,1)
    r$kid <- dec_i
    r$length <- length(li$ids)
    r
  }))
}))

print(df)

## each kid can match to multiple ids, but any will do the job. Each contains at minimum the population of the kid and parent:
kid_lut <- do.call(rbind,lapply(names(lineages),function(nm){
  data.frame(id=nm,kid=tail(lineages[[nm]]$ids,1)) 
}))

kid_lut <- kid_lut[!duplicated(kid_lut$kid),]
tmp <- kid_lut
kid_lut <- tmp$id
names(kid_lut) <- paste0("k",tmp$kid)

source("utils/comparison_functions.R")
preds <- list(p5 = readRDS("data/proc/summaries/salehi_preds_v2_minobs_5.Rds"),
              p10 = readRDS("data/proc/summaries/salehi_preds_v2_minobs_10.Rds"),
              p20 = readRDS("data/proc/summaries/salehi_preds_v2_minobs_20.Rds"))

library(parallel)
cl <- makeCluster(70)
clusterCall(cl, function() {
  source("utils/ALFA-K.R")
  source("utils/comparison_functions.R")
})

# Export variables to the cluster
clusterExport(cl, varlist = c("kid_lut","df","preds"), envir = environment())

r <- parLapplyLB(cl, 1:nrow(df),function(i){
  tryCatch({
    
    xtrain <- readRDS(paste0("data/salehi/alfak_inputs/",df$id[i],".Rds"))
    tt <- tail(as.numeric(colnames(xtrain$x)),2)
    ti <- tt[1]
    tj <- tt[2]
    expect <- data.frame(wasserstein=wasserstein_distance(xtrain,xtrain,ti,tj),
               cosine=compute_cosine_similarity(xtrain,xtrain,ti,tj),
               euclidean=compute_mean_karyotype_distance(xtrain,xtrain,ti,tj),
               overlap=compute_overlap_coefficient(xtrain,xtrain,ti,tj))
    
    input2load <- kid_lut[paste0("k",df$kid[i])]
    xdat <- readRDS(paste0("data/salehi/alfak_inputs/",input2load,".Rds"))
    xpred <- preds[[paste0("p",df$best_minobs[i])]][[df$id[i]]]
    t2 <- as.numeric(tail(colnames(xdat$x),1))
    x0 <- get_mean(xdat,t=as.numeric(head(colnames(xdat$x),1)))
    xt <- get_mean(xdat,t=t2)-x0
    
    r <- do.call(rbind,lapply(xpred,function(xi){
      r <- do.call(rbind,lapply(as.numeric(colnames(xi$x)),function(t1){
        xr <- get_mean(xi,t1)-x0
        data.frame(time=t1,
                   wasserstein=wasserstein_distance(xi,xdat,t1,t2),
                   cosine=compute_cosine_similarity(xi,xdat,t1,t2),
                   euclidean=compute_mean_karyotype_distance(xi,xdat,t1,t2),
                   overlap=compute_overlap_coefficient(xi,xdat,t1,t2),
                   angle=getangle(xt,xr)
        )
      }))
      r <- merge(df[i,],r)
      return(r)
    }))
    

    
    smm <- data.frame(time=r$time,delta=sapply(1:nrow(r),function(j) sum(abs(expect-r[j,colnames(expect)]))))
    smm <- aggregate(list(delta=smm$delta),by=list(time=smm$time),median)
    smm <- smm[!smm$time==0,]
    tpred <- smm$time[smm$delta==min(smm$delta)]
    r[r$time==tpred,]
    
    },error=function(e) return(NULL))
  
})

saveRDS(r,"data/proc/summaries/salehi_all_metrics_predict_timepoint.Rds")