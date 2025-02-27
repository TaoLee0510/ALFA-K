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
preds <- list(p5 = readRDS("data/proc/summaries/salehi_preds_minobs_5.Rds"),
              p10 = readRDS("data/proc/summaries/salehi_preds_minobs_10.Rds"),
              p20 = readRDS("data/proc/summaries/salehi_preds_minobs_20.Rds"))

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
    input2load <- kid_lut[paste0("k",df$kid[i])]
    xdat <- readRDS(paste0("data/salehi/alfak_inputs/",input2load,".Rds"))
    #xdat <- readRDS(paste0("data/proc/alfak_inputs/",input2load,".Rds"))
    xpred <- preds[[paste0("p",df$best_minobs[i])]][[df$id[i]]]
    t2 <- as.numeric(tail(colnames(xdat$x),1))
    x0 <- get_mean(xdat,t=as.numeric(head(colnames(xdat$x),1)))
    xt <- get_mean(xdat,t=t2)
    
    r <- do.call(rbind,lapply(xpred,function(xi){
      r <- do.call(rbind,lapply(as.numeric(colnames(xi$x)),function(t1){
        xr <- get_mean(xi,t1)
        data.frame(time=t1,
                   euc_pred = sqrt(mean((xr-x0)^2)),
                   euc_ref = sqrt(mean((xt-x0)^2)),
                   euc_err = sqrt(mean((xt-xr)^2)))
      }))
      r <- merge(df[i,],r)
      return(r)
    }))},error=function(e) return(NULL))
  
})

saveRDS(r,"data/proc/summaries/salehi_basic_metrics.Rds")