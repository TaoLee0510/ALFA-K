library(rootSolve)
library(ggplot2)
library(umap)
library(ggrepel)
library(Rcpp)
source("utils/ALFA-KQ.R")
source("tests/prediction/prediction_functions.R")

base_text_size <- 8
text_size_theme <- theme(
  text         = element_text(size = base_text_size, family = "sans"),
  axis.title   = element_text(size = base_text_size, family = "sans"),
  axis.text    = element_text(size = base_text_size, family = "sans"),
  legend.title = element_text(size = base_text_size, family = "sans"),
  legend.text  = element_text(size = base_text_size, family = "sans"),
  strip.text   = element_text(size = base_text_size, family = "sans")
)
base_theme <- theme_classic() + text_size_theme

sweepDir <- "data/salehi/alfak_outputs_V1_procv0_misseg_sweep/"

paths <- list.files(sweepDir,recursive = T)

info <- data.frame(
  path   = paths,
  run_id = sub("/.*", "", paths),
  minobs = as.numeric(sub(".*minobs_([0-9]+).*", "\\1", paths)),
  raw_p  = sub("_rep.*", "", sub("^.*/p_", "", paths)),
  stringsAsFactors = FALSE
)
info$pval <- as.numeric(gsub("p", ".", info$raw_p, fixed = TRUE)); info$raw_p <- NULL

run_mats <- lapply(split(info, info$run_id), function(df) {
  pvs <- sort(unique(df$pval)); 
  ks  <- colnames(readRDS(file.path(sweepDir, df$path[1])))[-1]
  mat <- pbapply::pbsapply(pvs, function(p) rowMeans(sapply(df$path[df$pval == p], function(f)
    as.numeric(tail(readRDS(file.path(sweepDir, f))[,-1], 1)))))
  rownames(mat) <- ks; colnames(mat) <- as.character(pvs)
  attr(mat, "run_id") <- df$run_id[1]; attr(mat, "minobs") <- unique(df$minobs)
  mat
})

# run_mats is a named list; e.g. run_mats[["SA532_X8_l_6_d1_1_d2_0"]] 
# is your k×pval matrix for that run_id




y <- lapply(run_mats,function(dfi){
  fi <- attr(dfi,"run_id")
  minobs <- attr(dfi,"minobs")
  dfi <- dfi[apply(dfi,1,max)>0.0005,]
  k <- do.call(rbind,lapply(rownames(dfi), s2v))
  print(nrow(k))
  u <- umap::umap(k)$layout
  u <- data.frame(u)
  
  lscape <- readRDS(paste0("data/salehi/alfak_outputs_V1a/minobs_",minobs,"/",fi,"/landscape.Rds"))
  lut <- lscape$mean
  names(lut) <- lscape$k
  
  u$f <- lut[rownames(dfi)]
  
  u <- cbind(u,dfi)
  u <- reshape2::melt(u,id.vars=c("X1","X2","f"))
  
  
  id <- dfi[,2]>dfi[,10]
  ploidy <- apply(k,1,mean)
  
  dfp <- data.frame(ploidy,id)
  dfp$f <- lut[rownames(dfp)]
  
  
  
  tmp <- data.frame(cbind(dfi,id),check.names = F)
  tmp <- reshape2::melt(tmp,id.vars="id")
  sm <- aggregate(list(freq=tmp$value),by=list(id=tmp$id,m=tmp$variable),sum)
  
  sm$lscape <- fi
  dfp$lscape <- fi
  u$lscape <- fi
  
  return(list(sm=sm,dfp=dfp,u=u))
  
})


y[[1]]$u$f <- y[[1]]$u$f-mean(y[[1]]$u$f)
y[[2]]$u$f <- y[[2]]$u$f-mean(y[[2]]$u$f)

u <- rbind(y[[1]]$u,y[[2]]$u)
u <- u[order(u$value,decreasing=F),]

lut <- c("p53 k.o. A\n(scenario A)","SA532\n(scenario B)")
names(lut) <- sel

tmp <- u[u$variable%in%c(0.0005,0.002),]
tmp <- split(tmp,f=interaction(tmp$variable,tmp$lscape))
tmp <- do.call(rbind,lapply(tmp,function(ti) ti[ti$value==max(ti$value),]))
tmp <- tmp[,!colnames(tmp)%in%c("variable","landscape","value")]
tmp$lab <- c("y","x","y","x")

p_freq <- ggplot(u[u$variable%in%c(0.0005,0.002),],aes(x=X1,y=X2))+ ## saturating colors aids visibility
  facet_grid(cols=vars(lut[lscape]),rows=vars(paste0("misrate = ",variable)),scales = "free")+
  geom_point(aes(color=pmin(0.05,pmax(0.000001,value))))+
  scale_color_viridis_c("frequency",breaks=c(0.01,0.001,0.0001,0.00001),trans="log")+
  geom_label_repel(data=tmp,aes(label=lab),box.padding = .5)+
  theme_classic()+
  scale_x_continuous("UMAP 1")+
  scale_y_continuous("UMAP 2")+
  text_size_theme
p_freq


p_fit <- ggplot(u[u$variable%in%c(0.0006,0.005),],aes(x=X1,y=X2))+ ## saturating colors aids visibility
  facet_grid(cols=vars(lut[lscape]),scales = "free")+
  geom_point(aes(color=f))+
  scale_color_viridis_c("fitness",option="magma")+
  geom_label_repel(data=tmp,aes(label=lab),box.padding = .5)+
  theme_classic()+
  scale_x_continuous("UMAP 1")+
  scale_y_continuous("UMAP 2")+
  text_size_theme
p_fit


sm <- rbind(y[[1]]$sm,y[[2]]$sm)

p_fgrp <- ggplot(sm,aes(x=m,y=freq,color=as.character(id),group=as.character(id)))+
  facet_grid(cols=vars(paste0(lut[lscape])))+
  geom_point()+
  geom_line()+
  theme_classic()+
  scale_color_manual("group",labels=c("x","y"),values=c("#CC79A7","#00724D"))+
  scale_x_discrete("misseg. rate")+
  scale_y_continuous("frequency")+
  text_size_theme+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p_fgrp

dfp <- rbind(y[[1]]$dfp,y[[2]]$dfp)

p_ploidy <- ggplot(dfp,aes(x=ploidy,fill=id,group=id))+
  facet_grid(cols=vars(lut[lscape]),scales="free")+
  geom_histogram(bins=12,position="dodge")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  theme_classic()+
  scale_y_sqrt(breaks=c(1,30,100))+
  scale_fill_manual("group",labels=c("x","y"),values=c("#CC79A7","#00724D"))+
  text_size_theme
p_ploidy

p_fit2 <- ggplot(dfp,aes(x=f,fill=id,group=id))+
  facet_grid(cols=vars(lut[lscape]),scales="free")+
  geom_histogram(bins=12,position="dodge")+
  scale_x_continuous("fitness",breaks = scales::pretty_breaks(n = 2))+
  theme_classic()+
  scale_y_sqrt(breaks=c(1,20,60))+
  scale_fill_manual("group",labels=c("x","y"),values=c("#CC79A7","#00724D"))+
  text_size_theme
p_fit2



pleft <- cowplot::plot_grid(p_xmpl,p_freq,labels=LETTERS[c(1,3)],nrow=2)
pright <- cowplot::plot_grid(p_fgrp,p_fit,p_fit2,p_ploidy,labels=LETTERS[c(2,4:6)],ncol=1)

plt <- cowplot::plot_grid(pleft,pright,rel_widths = c(4,3))

ggsave("figures/misc/figures/misseg.png",plt,width=8,height=7,units="in")


