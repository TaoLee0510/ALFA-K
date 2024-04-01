#setwd("~/projects/008_birthrateLandscape/ALFA-K/")
setwd("~/projects/ALFA-K/")
lineage_info <- function(ln){
  treat_numeric <- c(n=0,y=1)
  has_descendents <- length(ln$dec1)>0
  has_parents <- sum(is.na(m$parent[m$uid%in%ln$ids]))==0
  treatments <- m$on_treatment[m$uid%in%ln$ids]
  nChanges <- sum(abs(diff(treat_numeric[treatments])))
  treatment_status <- "mixed"
  if(mean(treatments=="n")==1) treatment_status <- "off"
  if(mean(treatments=="y")==1) treatment_status <- "on"
  id <- tail(ln$ids,1)
  sid <- paste(gsub(" ","", m$label[id==m$uid],fixed = TRUE),
               m$timepoint[id==m$uid],sep = "_")
  data.frame(id=sid,uid=id,has_descendents,has_parents,nChanges,treatment_status)
  
}

longest_cons <- function(ff){
  ids <- sapply(ff, function(fi){
    paste(head(unlist(strsplit(fi,split="_")),2),collapse="_")
  })
  ff <- split(ff,f=ids)
  
  ff <- sapply(ff, function(fi){
    fi <- fi[order(fi)]
    tail(fi,1)
  })
  as.character(ff)
}

m <- read.csv("data/salehi/metadata.csv")
lineages <- readRDS("figures/salehi_data_fitting/lineages.Rds")
linfo <- do.call(rbind,lapply(lineages,lineage_info))
linfo$filenames <- paste0(rownames(linfo),".Rds")
rownames(linfo) <- linfo$filenames
linfo1 <- linfo[linfo$nChanges==0&!linfo$has_descendents&!linfo$has_parents,]
linfo2 <- linfo[linfo$nChanges==0&!linfo$has_descendents&linfo$treatment_status=="on",]
linfo <- rbind(linfo1,linfo2)
library(parallel)
cl <- makeCluster(getOption("cl.cores", 8))
min_obs <- c(20,10,5)#,10,20)
outdirs <- paste0("data/salehi/alfak_fits/minobs_",min_obs,"/")
sapply(outdirs,dir.create,recursive=T)

ff <- list.files("data/salehi/alfak_inputs_v2/")
#favorites <- readRDS("data/salehi/favorites.Rds")

#ff <- ff[!ff%in%favorites]


#already_fitted <- list.files("data/salehi/alfak_fits/minobs_5/")
#ff <- ff[!ff%in%already_fitted]
#ff <- ff[grepl("SA535",ff)]
#ff <- ff[ff%in%longest_cons(linfo$filenames)]

parLapplyLB(cl=cl,X=ff, function(fi){
  min_obs <- c(20,10,5)#,10,20)
  outdirs <- paste0("data/salehi/alfak_fits/minobs_",min_obs,"/")
  source("utils/ALFA-K.R")
  print(fi)
  tryCatch({
    x <- readRDS(paste0("data/salehi/alfak_inputs_v2/",fi))
    nx <- rowSums(x$x)
    
    fits <- lapply(min_obs,function(ni){
      tryCatch({
        if(sum(nx>ni)<3) stop("insufficient data for XV procedure")
        fit <- alfak(x,min_obs = ni)
        xfq <- fit$xo[fit$xo$id=="fq",]
        fv <- unlist(sapply(1:nrow(xfq),optim_loo,xx=x,xo=xfq))
        xfq$f_xv <- fv
        fit$min_obs <- ni
        fit$xv_res <- xfq
        fit$vx_cor <- tryCatch({cor(xfq$f_est,xfq$f_xv,use="complete")},
                               error=function(e) return(-Inf))
        print(fit$vx_cor)
        fit
      },error=function(e) return(NULL))
      
    })
    for(i in 1:length(fits)){
      if(!is.null(fits[[i]])){
        saveRDS(fits[[i]],paste0(outdirs[i],fi))
      }
    }
  },error=function(e) print(e))
})



