library(ggplot2)
base_text_size=16
text_size_theme <- 
  theme(
    text         = element_text(size = base_text_size, family = "sans"),
    axis.title   = element_text(size = base_text_size, family = "sans"),
    axis.text    = element_text(size = base_text_size, family = "sans"),
    legend.title = element_text(size = base_text_size, family = "sans"),
    legend.text  = element_text(size = base_text_size, family = "sans"),
    strip.text   = element_text(size = base_text_size, family = "sans")
  )
setwd("/share/lab_crd/M010_ALFAK_2023/ALFA-K")
base_text_size=16
text_size_theme <- 
  theme(
    text         = element_text(size = base_text_size, family = "sans"),
    axis.title   = element_text(size = base_text_size, family = "sans"),
    axis.text    = element_text(size = base_text_size, family = "sans"),
    legend.title = element_text(size = base_text_size, family = "sans"),
    legend.text  = element_text(size = base_text_size, family = "sans"),
    strip.text   = element_text(size = base_text_size, family = "sans")
  )
source("utils/ALFA-KQ.R")
x <- readRDS("data/salehi/alfak_outputs_proc.Rds")
x <- x[!is.na(x$xv)&x$xv>0&x$test_treat=="NaN",] 
x <- split(x,f=x$fi)
x <- do.call(rbind,lapply(x,function(xi){
  xi[xi$xv==max(xi$xv),]
}))
rownames(x) <- NULL

safe_id <- gsub("SA535_CISPLATIN_Combined","SA535CISPLATINCombined",x$fi)
safe_id <- sapply(safe_id,function(ii){
  paste(unlist(strsplit(ii,split="_"))[1:2],collapse="_")
})
x <- split(x,f=safe_id)

x <- do.call(rbind,lapply(x,function(xi){
  xi[xi$ntrain==max(xi$ntrain),]
}))
rownames(x) <- NULL

m <- read.csv("data/salehi/metadata.csv")

## trying to compare in vivo / PDX:
lut <- m$PDX_id
lut2 <- m$timepoint
names(lut) <- paste(m$datasetname,m$timepoint,sep="_")
names(lut2) <- paste(m$datasetname,m$timepoint,sep="_")
tmp <- sapply(x$fi,function(ffi){
  tmp <- strsplit(ffi,split="_") |> unlist()
  ti <- which(tmp=="l")
  tmp <- paste(tmp[1:(ti-1)],collapse="_")
})

x$pdx <- lut[tmp]
x$tp <- lut2[tmp]
## load each fitted landscape, compute delta-f from fq to neighbours, report as function of ploidy and treatment status

df <- do.call(rbind,pbapply::pblapply(1:nrow(x),function(i){
  print(i)
  inpath <- paste0("data/salehi/alfak_outputs/minobs_",
                   x$min_obs[i],"/",x$fi[i],"/landscape.Rds")
  if(!file.exists(inpath)) return(NULL)
  lscape <- readRDS(inpath)
  
  lut <- lscape$mean
  names(lut) <- lscape$k
  
  #  fq <- lscape$k[lscape$fq|lscape$nn]
  fq <- lscape$k[lscape$fq]#|lscape$nn]
  
  dfi <- do.call(cbind,lapply(fq,function(fqi){
    nn <- apply(gen_all_neighbours(fqi,remove_nullisomes = F),1,paste,collapse=".")
    lut[fqi]-lut[nn]
  }))

  
  dfi <- cor(dfi,use="complete")
  di <- as.matrix(dist(do.call(rbind,lapply(fq,s2v)),method="manh"))
  
  dfi <- data.frame(cor=dfi[upper.tri(dfi)],dist=di[upper.tri(di)])
  
  tmp <-  strsplit(x$train_treat[i],split="") |> unlist()
  dfi$treat_switch <- (tmp |> unique() |> length())>1
  dfi$treated <- grepl("y",x$train_treat[i],fixed=TRUE) 
  dfi$ntrain <- x$ntrain[i]
  dfi$pdx <- x$pdx[i]
  dfi$fi <- x$fi[i]
  return(dfi)
}))

z <- aggregate(list(cor=df$cor),by=list(dist=df$dist),mean)
p0 <- ggplot(df,aes(x=round(dist),y=cor))+
  geom_jitter(height=0,width=0.3)+
  geom_point(data=z,color="red",shape=95,size=10)+
  scale_x_continuous("manhattan distance",limits=c(0,10),breaks=seq(2,8,2))+
  scale_y_continuous("Pearson's r")+
  theme_minimal()+
  text_size_theme

p0
ggsave("figures/misc/figures/cna_fitness_parental.png",p0,width=4,height=4,units="in")

tmp <- do.call(rbind,pbapply::pblapply(1:nrow(x),function(i){
  print(i)
  inpath <- paste0("data/salehi/alfak_outputs/minobs_",
                   x$min_obs[i],"/",x$fi[i],"/landscape.Rds")
  if(!file.exists(inpath)) return(NULL)
  lscape <- readRDS(inpath)
  
  lut <- lscape$mean
  names(lut) <- lscape$k
  
  #  fq <- lscape$k[lscape$fq|lscape$nn]
  fq <- lscape$k[lscape$fq]#|lscape$nn]
  
  dfi <- data.frame(t(do.call(cbind,lapply(fq,function(fqi){
    nn <- apply(gen_all_neighbours(fqi,remove_nullisomes = F),1,paste,collapse=".")
    lut[fqi]-lut[nn]
  }))))
  colnames(dfi) <- NULL
  dfi <- cbind(k=fq,pdx=x$pdx[i],fi=x$fi[i],dfi)
  return(dfi)
}))
#tmp$pdx[tmp$pdx=="SA039"] <- "SA906"
tmat <- tmp[,-c(1,2,3)]
cm <- cor(t(tmat),use="complete")
dm <- as.matrix(dist(do.call(rbind,lapply(tmp$k,s2v)),method="manh"))
sm <- sapply(tmp$pdx,function(a) sapply(tmp$pdx,function(b) a==b))
idm <- sapply(tmp$fi,function(a) sapply(tmp$fi,function(b) paste(a,b)))
km <- sapply(tmp$k,function(a) sapply(tmp$k,function(b) paste(a,b)))

dfx <- data.frame(cm=cm[upper.tri(cm)],dm=dm[upper.tri(dm)],
                  sm=sm[upper.tri(sm)],k=km[upper.tri(km)],
                  id=idm[upper.tri(idm)])
dfx <- dfx[order(dfx$sm,dfx$dm),]

dfy <- aggregate(list(cm=dfx$cm),by=list(dm=dfx$dm,sm=dfx$sm),function(s) c(mean(s)-sd(s),mean(s),mean(s)+sd(s)))
dfy <- cbind(dfy[,1:2],dfy$cm)
colnames(dfy)[3:5] <- c("lo","med","hi")

ninc <- 10
p1 <- ggplot(dfx[!dfx$sm & dfx$dm<=ninc,],aes(x=dm))+
  geom_jitter(aes(y=cm),height=0,width=0.1)+
  geom_point(data=dfy[!dfy$sm&dfy$dm<=ninc,],aes(y=med),color="red",shape=95,size=10)+
  scale_x_continuous("manhattan distance",limits=c(0,10),breaks=seq(2,8,2))+
  scale_y_continuous("Pearson's r")+
  theme_minimal()+
  text_size_theme
p1

dfx <- dfx[order(dfx$sm,dfx$dm),]
