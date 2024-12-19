
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
  return(lineages)
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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

call_cn <- function(x,arms,level="chrom"){
  ##proc cn file
  df <- do.call(rbind,lapply(rownames(x), function(xi) unlist(strsplit(xi,split="_"))))
  df <- data.frame(df[,1:2])
  colnames(df) <- c("chrom","chrompos")
  df$arm <- sapply(1:nrow(df), function(i){
    chrom <- df$chrom[i]
    chrompos <- df$chrompos[i]
    is_q <- arms[arms$chrom==chrom&arms$arm=="q","start"]<chrompos
    is_p <- !is_q
    
    c("q","p")[c(is_q,is_p)]  
    
  })
  
  x <- cbind(df[,c("chrom","arm")],x)
  x <- x[!df$chrom%in%c("X","Y"),]
  xcolnames <- colnames(x)
  if(level=="chrom")  x <- split(x,f=x$chrom)
  if(level=="arm") x <- split(x,f=interaction(x$chrom,x$arm))
  lx <- sapply(x,nrow)
  x <- x[lx>0]
  
  x <- data.frame(do.call(rbind,lapply(x,  function(xi){
    if(level=="chrom")  id <- xi[1,1]
    if(level=="arm") id <- paste(xi[1,1:2],collapse="")
    y <- c(id,as.numeric(apply(xi[,3:ncol(xi)],2,Mode)))
    y <- as.character(y)
  })))
  colnames(x) <- c("chrom",xcolnames[-c(1,2)])
  chrom <- x$chrom
  x <- t(x[,!colnames(x)=="chrom"])
  colnames(x) <- stringr::str_pad(colnames(x),width=4)
  x <- x[,order(colnames(x))]
  return(x)
}

extract_cn_profiles <- function(data_path){
  setwd(data_path)
  arms <- readRDS("../arm_loci.Rds")
  
  ff <- list.files()
  #bins <- lapply(ff, function(fi){
  # if(!file.exists(paste0(fi,"/named_mat.csv"))) return(NULL)
  #x <- data.frame(data.table::fread(paste0(fi,"/named_mat.csv")))
  #x$bin_name
  #})
  
  #bins <- do.call(cbind,bins)
  
  ## for some reason at least one file (SA1035H) has incompletely named bin_names, 
  ## so i took the approach of using the bin names from the first files (which are OK)
  ## this works for all except SA004 which has 6206 bins for some reason.
  ##OK turns out SA004 is very different from the rest (what is it!?) so lets exclude that
  #f0 <- ff[1]
  #bn <- data.frame(data.table::fread(paste0(f0,"/named_mat.csv")))$bin_name
  
  #nbins <- sapply(ff, function(fi) nrow(data.table::fread(paste0(fi,"/named_mat.csv"))))
  ff <- ff[!ff=="SA004"]
  
  ##IF USING POST JUMP DATA UNCOMMENT ABOVE THEN FOLLOW CAPS COMMENTS BELOW
  
  x <- do.call(rbind,lapply(ff, function(fi){
    print(fi)
    if(!file.exists(paste0(fi,"/named_mat.csv"))) return(NULL)
    x <- data.frame(data.table::fread(paste0(fi,"/named_mat.csv")))
    binformat <- unlist(strsplit(as.character(x$bin_name[1]),split="_"))
    if(length(binformat)<3) return(NULL)
    rownames(x) <- x$bin_name ##REPLACE WITH rownames(x) <- bn 
    x <- x[,!colnames(x)=="bin_name"] 
    call_cn(x,arms)
  }))
  
  cx <- colnames(x)
  x <- t(apply(x,1,as.numeric))
  colnames(x) <- cx
  ids <- rownames(x)
  ids <- sapply(ids, function(idi) unlist(strsplit(idi,split="[.]"))[2])
  x <- data.frame(x,check.names=F)
  x <- split(x,f=ids)
  saveRDS(x,file="../chrom_level_cn.Rds")
  
  x <- do.call(rbind,lapply(ff, function(fi){
    print(fi)
    if(!file.exists(paste0(fi,"/named_mat.csv"))) return(NULL)
    x <- data.frame(data.table::fread(paste0(fi,"/named_mat.csv")))
    binformat <- unlist(strsplit(as.character(x$bin_name[1]),split="_"))
    if(length(binformat)<3) return(NULL)
    rownames(x) <- x$bin_name ##REPLACE WITH rownames(x) <- bn 
    x <- x[,!colnames(x)=="bin_name"] 
    call_cn(x,arms,level = "arm")
  }))
  
  cx <- colnames(x)
  x <- t(apply(x,1,as.numeric))
  colnames(x) <- cx
  ids <- rownames(x)
  ids <- sapply(ids, function(idi) unlist(strsplit(idi,split="[.]"))[2])
  x <- data.frame(x,check.names=F)
  x <- split(x,f=ids)
  saveRDS(x,file="../arm_level_cn.Rds")
}

extract_lineages <- function(data_path){
  setwd(data_path)
  m <- read.csv("metadata.csv")
  lineages <- process_lineages(m)
  saveRDS(lineages,"lineages.Rds")
  outdir <- "alfak_inputs/"
  dir.create(outdir,showWarnings = F,recursive = T)
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
}
