setwd("~/projects/008_birthrateLandscape/ALFA-K/")
outdir <- "data/salehi/alfak_inputs_v2/"
lineages <- readRDS("figures/salehi_data_fitting/data/lineages.Rds")
m <- read.csv("data/salehi/metadata.csv")
cnmat <- readRDS("data/salehi/chrom_level_cn.Rds")
for(i in 1:length(lineages)){
  id <- names(lineages)[i]
  print(i)
  lineage <- lineages[[i]]
  
  
  
  msel <- m[m$uid%in%lineage$ids,]
  pdx_id <- msel$PDX_id[1]
  dt <- 15
  if(pdx_id%in%c("SA039","SA906")) dt <- 5
 
  msel <- do.call(rbind,lapply(1:nrow(msel), function(i){
    lid <- unlist(strsplit(msel$library_ids[i],split=";"))
    data.frame(timepoint=msel$timepoint[i],library_id=lid,row.names = NULL)
  }))
  
  if(mean(msel$library_id%in%names(cnmat))<1) next
  
  msel <- msel[msel$library_id%in%names(cnmat),]
  x <- cnmat[msel$library_id]
  ids <- unlist(lapply(1:nrow(msel), function(i) rep(msel$timepoint[i],nrow(x[[i]]))))
  x <- do.call(rbind,x)
  
  k <- apply(x,1,paste,collapse=".")
  x0 <- unique(k)
  
  x <- split(k,f=ids)
  
  x <- do.call(cbind,lapply(x, function(xi){
    sapply(x0, function(xj) sum(xi==xj))
  }))
  x <- x[order(rowSums(x),decreasing=T),]
  colnames(x) <- gsub("X","",colnames(x))
  x <- x[,order(as.numeric(colnames(x)))]
  x <- list(x=x,pop.fitness=NULL,dt=dt)
  saveRDS(x,paste0(outdir,id,".Rds"))
}


