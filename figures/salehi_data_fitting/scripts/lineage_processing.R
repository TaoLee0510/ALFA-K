lineage_wrapper <- function(ids){
  list(ids=ids)
}

lineage_generator <- function(uid){
  ids <- c()
  lineages <- list()
  while(!is.na(uid)){
    ids <- c(uid,ids)
    uid <- m$parent[m$uid==uid]
    if(length(ids)>1) lineages <- c(lineages,lineage_wrapper(ids))
  }
  return(lineages)
}

descendent_checker <- function(ids){
  dec1 <- m$uid[m$parent==tail(ids,1)&!is.na(m$parent)]
  dec2 <- m$uid[m$parent%in%dec1]
  list(ids=ids,dec1=dec1,dec2=dec2)
}

lineage_namer <- function(lineage){
  suffix <- paste0("_l_",length(lineage$ids),"_d1_",length(lineage$dec1),
                   "_d2_",length(lineage$dec2))
  id <- tail(lineage$ids,1)
  sid <- paste(gsub(" ", "", m$datasetname[id==m$uid],fixed = TRUE),m$timepoint[id==m$uid],sep = "_")
  sid <- gsub("/", "", sid,fixed = TRUE)
  paste0(sid,suffix)
}

process_lineages <- function(m){
  lineages <- do.call(c,lapply(m$uid,lineage_generator))
  lineages <- lapply(lineages, descendent_checker)
  
  lnames <- sapply(lineages,lineage_namer)
  names(lineages) <- lnames
}

longest_cons <- function(ff,longest=T){
  ids <- sapply(ff, function(fi){
    fsplit <- unlist(strsplit(fi,split="_"))
    lsplit <- length(fsplit)
    fsplit <- fsplit[1:(lsplit-6)]
    paste(fsplit,collapse="_")
  })
  ff <- split(ff,f=ids)
  
  ff <- sapply(ff, function(fi){
    fi <- fi[order(fi)]
    if(longest) return(tail(fi,1))
    return(head(fi,1))
  })
  as.character(ff)
}
wrap_tree <- function(pdx_id,m){
  x <- m[m$PDX_id==pdx_id,]
  x$x <- as.numeric(gsub("X","",x$x))
  
  xseg <- x[!is.na(x$parent),]
  
  xseg <- cbind(xseg,do.call(rbind,lapply(xseg$parent, function(i){
    tmp <- x[x$uid==i,c("x","y")]
    colnames(tmp) <- paste0(colnames(tmp),"start")
    tmp
  })))
  list(x=x,xseg=xseg)
}

signr <- function(pval){
  if(pval<0.001) return("***")
  if(pval<0.01) return("**")
  if(pval<0.05) return("*")
  return("")
}
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
  cellLine <- m$PDX_id[m$uid==id]
  sid <- paste(gsub(" ","", m$label[id==m$uid],fixed = TRUE),
               m$timepoint[id==m$uid],sep = "_")
  data.frame(id=sid,uid=id,cellLine,has_descendents,has_parents,nChanges,treatment_status)
  
}

longest_cons <- function(ff,longest=T){
  ids <- sapply(ff, function(fi){
    fsplit <- unlist(strsplit(fi,split="_"))
    lsplit <- length(fsplit)
    fsplit <- fsplit[1:(lsplit-6)]
    paste(fsplit,collapse="_")
  })
  ff <- split(ff,f=ids)
  
  ff <- sapply(ff, function(fi){
    fi <- fi[order(fi)]
    if(longest) return(tail(fi,1))
    return(head(fi,1))
  })
  as.character(ff)
}

get_r2 <- function(id,mo){
  x <- readRDS(paste0(dir,mo,"/",id,"."))
  xvr <- x$xv_res
  maxf <- max(xvr$f_est)
  xvr <- xvr[!is.na(xvr$f_xv),]
  r2 <- R2(obs = xvr$f_est,pred=xvr$f_xv)
  id <- unlist(strsplit(id,split=".Rds"))[1]
  id <- unlist(strsplit(id,split="_"))
  if(length(id)>8){
    idt <- tail(id,7)
    idh <- head(id,length(id)-7)
    idh <- paste(idh,collapse="_")
    id <- c(idh,idt)
  }
  mo <- tail(unlist(strsplit(mo,split="_")),1)
  res <- c(id[c(1,2,4,6,8)],mo,maxf,r2)
  names(res) <- c("datasetname","timepoint","samples","dec1","dec2","min_obs","maxf","r2")
  return(res)
}

assign_labels <- function(m){
  root <- m$uid[is.na(m$parent)]
  lineages <- list(c(x=root))
  new_lineages <- list()
  finished_lineages <- list()
  
  
  while(length(lineages)>0){
    for(li in lineages){
      descendents <- m$uid[!is.na(m$parent)&m$parent==tail(li,1)]
      if(length(descendents)==0) finished_lineages <- c(finished_lineages,list(li))
      new_lineages <- c(new_lineages,lapply(descendents, function(di) {
        lij <- c(li,di)
        names(lij)[length(lij)] <- tail(names(li),1)
        if(length(descendents)>1){
          names(lij)[length(lij)] <- paste0(tail(names(li),1),letters[which(descendents==di)])
        }
        return(lij)
      }))
    }
    lineages <- new_lineages
    new_lineages <- list()
    
  }
  
  ids <- unlist(finished_lineages)
  ids <- ids[!duplicated(ids)]
  m <- m[order(m$uid),]
  ids <- ids[order(ids)]
  m$linlab <- names(ids)
  m$linlab <- gsub("x","",m$linlab)
  return(m)
}

adder <- function(xx,dx,dy){
  xx$x$x <- xx$x$x + dx
  xx$xseg$x <- xx$xseg$x + dx
  xx$xseg$xstart <- xx$xseg$xstart + dx
  xx$x$y <- xx$x$y + dy
  xx$xseg$y <- xx$xseg$y + dy
  xx$xseg$ystart <- xx$xseg$ystart + dy
  return(xx)
}

wrap_assign_labels <- function(m){
  m <- split(m,f=m$PDX_id)
  do.call(rbind,lapply(m,assign_labels))
}

