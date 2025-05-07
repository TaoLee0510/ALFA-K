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
setwd("/share/lab_crd/M010_ALFAK_2023/ALFA-K/")

procDir <- "data/salehi/alfak_outputs_V1a/"

minobs <- list.files(procDir)

ff <- unique(list.files(paste0(procDir,minobs)))


x <- do.call(rbind,lapply(ff,function(fi){
  dfi <- do.call(rbind,lapply(minobs,function(mi){
    xval=NaN
    xvpath <- paste0(procDir,mi,"/",fi,"/xval.Rds")
    if(file.exists(xvpath)) xval=readRDS(xvpath)
    data.frame(minobs=gsub("minobs_","",mi),fi=fi,xval=xval)
  }))
  dfi <- dfi[!is.na(dfi$xval),]
  dfi <- dfi[dfi$xval==max(dfi$xval),]
}))

x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- split(x,f=x$clean_id)
x <- do.call(rbind,lapply(x,function(xi) xi[xi$xv==max(xi$xv),]))

m <- read.csv("data/salehi/metadata.csv")
m$clean_id <- paste(m$datasetname,m$timepoint,sep="_")

df <- merge(m,x,all=T)

