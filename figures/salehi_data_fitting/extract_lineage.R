setwd("~/projects/008_birthrateLandscape/ALFA-K/data/salehi")

conds <- data.frame(id=c("A96146A","A96149B","A118833A","A110673B","A118862A",
                         "A96219B","A118845B","A98244A","A98256A","A98283A",
                         "A96162B"),
                    dt=c(5,5,15,15,15,15,15,15,15,15,15),
                    mintp=c(1,1,1,4,4,1,1,5,1,4,1))

for(i in 1:nrow(conds)){
  dt <- conds$dt[i]
  mintp <- paste0("X",conds$mintp[i])
  id <- conds$id[i]
  sid <- paste(gsub(" ", "", m$label[grepl(id,m$library_ids)],fixed = TRUE),m$timepoint[grepl(id,m$library_ids)],sep = "_")
  sid <- gsub("/", "", sid,fixed = TRUE)
  m <- read.csv("metadata.csv")
  cnmat <- readRDS("chrom_level_cn.Rds")
  
  uids <- c(m$uid[grepl(id,m$library_ids)],m$parent[grepl(id,m$library_ids)])
  
  while(!is.na(tail(uids,1))){
    uids <- c(uids,m$parent[m$uid==tail(uids,1)])
  }
  
  msel <- m[m$uid%in%uids,]
  
  msel <- do.call(rbind,lapply(1:nrow(msel), function(i){
    lid <- unlist(strsplit(msel$library_ids[i],split=";"))
    data.frame(timepoint=msel$timepoint[i],library_id=lid,row.names = NULL)
  }))
  msel <- msel[msel$timepoint>=mintp,]
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
  saveRDS(x,paste0("alfak_inputs/",sid,".Rds"))
}



