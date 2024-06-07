setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/comparison_functions.R")
source("figures/salehi_data_fitting/scripts/lineage_processing.R")

m <- read.csv("data/salehi/metadata.csv")
lineages <- readRDS("figures/salehi_data_fitting/lineages.Rds")
cnmat <- readRDS("data/salehi/chrom_level_cn.Rds")
x <- readRDS("figures/salehi_data_fitting/fit_summaries.Rds")
xf <- x[x$has_descendents&x$min_obs==5,]
ids <- as.character(sapply(xf$filenames,function(fi){head(unlist(strsplit(fi,split=".Rds")),1)}))
simdir <- "data/salehi/forward_sims_v2/minobs_5/"

tst <- which(ids%in%list.files(simdir))

df <- do.call(rbind,pbapply::pblapply(tst, function(i){
  #print(i)
  #si <- paste0(simdir,ids[i],"/output/00000/")
  odir <- paste0(simdir,ids[i],"/output/")
  s <- paste0(odir,list.files(odir),"/")
  # xs <- proc_sim(si,times=seq(0,500,100))$x
  xs <- lapply(s,function(si) proc_sim(si,times=seq(0,500,100))$x)
  
  xs <- do.call(rbind,xs)
  
  r <- unique(rownames(xs))
  
  xs <- do.call(rbind,lapply(r, function(ri) colSums(xs[rownames(xs)==ri,,drop=F])))
  
  xs <- data.frame(xs,check.names = F,row.names = r)
  
  vs <- do.call(rbind,lapply(rownames(xs),s2v))
  #vs <- lapply(xs, function(xsi) do.call(rbind,lapply(rownames(xsi),s2v)))
  li <- lineages[[ids[i]]]
  
  ## sim output population vectors
  
  #xs <- lapply(1:length(xs),function(z){
  # csx <- colSums(xs[[z]])
  #xs[[z]] <- t(xs[[z]])%*%vs[[z]]
  #for(k in 1:nrow(xs[[z]])) xs[[z]][k,] <- xs[[z]][k,]/csx[k] 
  #return(xs[[z]])
  #})
  
  csx <- colSums(xs)
  xs <- t(xs)%*%vs
  for(k in 1:nrow(xs)) xs[k,] <- xs[k,]/csx[k]
  target_library_ids <- m$library_ids[m$uid==tail(li$ids,1)]
  target_library_ids <- unlist(strsplit(target_library_ids,split=";"))
  x0 <- colMeans(do.call(rbind,lapply(target_library_ids,function(tli) cnmat[[tli]])))
  
  libid_exists <- sapply(li$dec1,function(idi){
    lidi <- m$library_ids[m$uid==idi]
    lidi%in%names(cnmat)
  })
  if(sum(libid_exists)==0) return(NULL)
  li$dec1 <- li$dec1[libid_exists]
  
  dfi <- do.call(rbind,lapply(1:length(li$dec1), function(j){
    target_library_ids <- m$library_ids[m$uid==li$dec1[j]]
    target_library_ids <- unlist(strsplit(target_library_ids,split=";"))
    xd <- colMeans(do.call(rbind,lapply(target_library_ids,function(tli) cnmat[[tli]])))
    
    #df <- do.call(rbind,lapply(1:length(xs),function(z){
    #  angles <- apply(xs[[z]],1,function(xsi){
    #   a <- xsi-x0
    #  b <- xd-x0
    # ai <- getangle(a,b)
    #})
    #names(angles) <- paste0("a",names(angles))
    #angles <- cbind(xf[i,],t(angles))
    #angles
    #}))
    
    angles <- apply(xs,1,function(xsi){
      a <- xsi-x0
      b <- xd-x0
      ai <- getangle(a,b)
    })
    names(angles) <- paste0("a",names(angles))
    df <- cbind(xf[i,],t(angles))
    #df$rep <- 1:nrow(df)
    df$descendent <- li$dec1[j]
    df$declevel <- 1
    return(df)
  }))
  libid_exists <- sapply(li$dec2,function(idi){
    lidi <- m$library_ids[m$uid==idi]
    lidi%in%names(cnmat)
  })
  if(length(libid_exists)>0) li$dec2 <- li$dec2[libid_exists]
  if(length(li$dec2)>0){
    dfii <- do.call(rbind,lapply(1:length(li$dec2), function(j){
      target_library_ids <- m$library_ids[m$uid==li$dec2[j]]
      target_library_ids <- unlist(strsplit(target_library_ids,split=";"))
      xd <- colMeans(do.call(rbind,lapply(target_library_ids,function(tli) cnmat[[tli]])))
      
      angles <- apply(xs,1,function(xsi){
        a <- xsi-x0
        b <- xd-x0
        ai <- getangle(a,b)
      })
      names(angles) <- paste0("a",names(angles))
      df <- cbind(xf[i,],t(angles))
      #df <- do.call(rbind,lapply(1:length(xs),function(z){
      # angles <- apply(xs[[z]],1,function(xsi){
      #  a <- xsi-x0
      # b <- xd-x0
      #ai <- getangle(a,b)
      #})
      #names(angles) <- paste0("a",names(angles))
      #angles <- cbind(xf[i,],t(angles))
      #angles
      #}))
      #df$rep <- 1:nrow(df)
      df$descendent <- li$dec2[j]
      df$declevel <- 2
      return(df)
    }))
    dfi <- rbind(dfi,dfii)
  }
  
  return(dfi)
  
}))

saveRDS(df,"figures/salehi_data_fitting/data/validation_summaries.Rds")