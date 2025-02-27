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

library(parallel)
cl <- makeCluster(70)
clusterCall(cl, function() {
  source("utils/ALFA-K.R")
  source("utils/comparison_functions.R")
})

# Export variables to the cluster
clusterExport(cl, varlist = c("kid_lut","df"), envir = environment())

r <- parSapplyLB(cl, 1:nrow(df),function(i){
  tryCatch({
    input2load <- kid_lut[paste0("k",df$kid[i])]
    xdat <- readRDS(paste0("data/salehi/alfak_inputs/",input2load,".Rds"))
    path2fit <- paste0("data/salehi/alfak_fits/minobs_",df$best_minobs[i],"/",df$id[i],".Rds")
    fit <- readRDS(path2fit)
    fq <- rownames(fit$xo)[fit$xo$id=="fq"]
    nn <- gen_all_neighbours(fq)
    nn <- apply(nn,1,paste,collapse=".")
    n2 <- gen_all_neighbours(nn)
    n2 <- apply(n2,1,paste,collapse=".")
    L <- c(fq,n2,d2)#,n2)
    
    filter <- rownames(xdat$x)%in%L
    
    k <- do.call(rbind,lapply(rownames(xdat$x),s2v))
    f <- predict(fit$fit,k)
    
   # filter <- xdat$x[,1]>5
    
    f1 <- sum(xdat$x[filter,1]*f[filter])/sum(xdat$x[filter,1])
    f2 <- sum(xdat$x[filter,2]*f[filter])/sum(xdat$x[filter,2])
    
    r <- f2-f1
    
    
    
    },error=function(e) return(NaN))
  
})
df$predicted_deltaf <- r

saveRDS(df,"data/proc/summaries/salehi_is_fitness_increasing.Rds")