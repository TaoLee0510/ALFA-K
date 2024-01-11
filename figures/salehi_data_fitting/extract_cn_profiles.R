setwd("~/projects/008_birthrateLandscape/ALFA-K/data/salehi/raw/")
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

arms <- readRDS("../arm_loci.Rds")

ff <- list.files()

## for some reason at least one file (SA1035H) has incompletely named bin_names, 
## so i took the approach of using the bin names from the first files (which are OK)
## this works for all except SA004 which has 6206 bins for some reason.
##OK turns out SA004 is very different from the rest (what is it!?) so lets exclude that
f0 <- ff[1]
bn <- data.frame(data.table::fread(paste0(f0,"/named_mat.csv")))$bin_name

nbins <- sapply(ff, function(fi) nrow(data.table::fread(paste0(fi,"/named_mat.csv"))))
ff <- ff[!ff=="SA004"]

x <- do.call(rbind,lapply(ff, function(fi){
  print(fi)
  if(!file.exists(paste0(fi,"/named_mat.csv"))) return(NULL)
  x <- data.frame(data.table::fread(paste0(fi,"/named_mat.csv")))
  rownames(x) <- bn
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
  rownames(x) <- bn
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
