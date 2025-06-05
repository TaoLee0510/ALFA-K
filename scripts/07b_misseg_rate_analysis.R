library(alfakR)
library(ggplot2)
setwd("/share/lab_crd/M010_ALFAK_2023/ALFA-K")

alfak_fit_dir <- "data/processed/salehi/alfak_outputs/"

df <- readRDS("data/processed/salehi/steadyStatePredictions.Rds")

s2v <- function(charvec) as.numeric(unlist(strsplit(charvec,split="[.]")))


## this is a manual screen looking to see if there is a missegregation rate dependent switch in clonal dominance.

pscreen <- lapply(names(df),function(i){
  tryCatch({
    print(i)
    dfi <- df[[i]]
    dfi <- dfi[apply(dfi,1,max)>0.001,]
    k <- do.call(rbind,lapply(rownames(dfi), s2v))
    print(nrow(k))
    u <- umap::umap(k)$layout
    u <- data.frame(u)
    u <- cbind(u,dfi)
    u <- reshape2::melt(u,id.vars=c("X1","X2"))
    p <- ggplot(u,aes(x=X1,y=X2,color=value))+
      facet_wrap(~variable)+
      geom_point()+
      scale_color_viridis_c(trans="log")+
      ggtitle(i)
    p
  },error=function(e) return(NULL))
  
})

run_abm_validation <- function(dfi,Nreps = 3){
  min_obs <- attr(dfi,"min_obs")
  fi <- attr(dfi,"fi")

  lscape <- readRDS(paste0(alfak_fit_dir,fi,"/minobs_",min_obs,"/landscape.Rds"))  
  
  pvec <- as.numeric(colnames(dfi))
  
  pvec
  
  x0 <- rep(1/nrow(lscape),nrow(lscape))
  names(x0) <- lscape$k
  p <- pvec[1]
  out <- alfakR::predict_evo(lscape,p,times=c(0,1000),x0=x0,
                             prediction_type = "ABM",
                             abm_delta_t=0.1,abm_max_pop=1e6,abm_record_interval = 1000)
  
  dfi <- head(dfi,400) ## kind of an arbitrary threshold...
  ploidy <- sapply(rownames(dfi),function(i) mean(as.numeric(s2v(i))))
  lscape <- readRDS(paste0(alfak_fit_dir,fi,"/minobs_",min_obs,"/landscape.Rds"))
  
  ploidy <- sapply(lscape$k,function(i) mean(as.numeric(s2v(i))))
  
  lscape <- lscape[lscape$k%in%rownames(dfi),]
  
  lut <- lscape$mean
  names(lut) <- lscape$k
  fitness <- lut[rownames(dfi)]  
  
  ump <- umap::umap(do.call(rbind,lapply(rownames(dfi),s2v)))
  ump <- data.frame(ump$layout)
  
  cbind(dfi,ump,ploidy,fitness,fi)
}

ids <- c("SA532_X4_l_4_d1_1_d2_1","SA906a_X50_l_5_d1_1_d2_0")
df_list <- df[ids]

x <- do.call(rbind,lapply(df_list,get_plot_data))

xre <- reshape2::melt(x,id.vars=c("X1","X2","ploidy","fitness","fi"))



p <- ggplot(xre,aes(x=X1,y=X2,color=value))+
  facet_grid(cols=vars(variable),rows=vars(fi),scales="free")+
  geom_point()+
  scale_color_viridis_c(trans="log")
p

p <- ggplot(xre,aes(x=X1,y=X2,color=ploidy))+
  facet_grid(cols=vars(variable),rows=vars(fi),scales="free")+
  geom_point()+
  scale_color_viridis_c()
p


