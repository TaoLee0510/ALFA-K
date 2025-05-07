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
## ODE Model equation
chrmod <- function(time,state,parms){
  with(as.list(parms),{
    ds <- state%*%A
    ds <- ds-sum(ds)*state
    return(list(ds))
  })
}
reverse_state_index <- function(index,maxchrom,minchrom=1,ndim){
  ## reset maxchrom to behave as if minchrom was 1,1,...
  mc <- maxchrom-minchrom+1
  if(length(mc)==1 & ndim>1) mc <- rep(mc,ndim)
  ## how many elements does each part of the state represent?
  ## works as prod(numeric(0)) evaluates to 1:
  Nsites <- cumprod(mc)
  cp <- c(1,cumprod(mc)[-length(mc)])
  state <- c()
  for(j in ndim:1){
    ni <- floor((index-1)/cp[j])
    state <- c(ni+1,state)
    index <- index-cp[j]*ni
  }
  state + minchrom-1
}
pij<-function(i, j, beta){
  qij <- 0
  if(abs(i-j)>i){ ## not enough copies for i->j
    return(qij)
  }
  # code fails for j = 0, but we can get this result by noting i->0 always
  ## accompanies i->2i
  if(j==0) j <- 2*i
  s <- seq(abs(i-j), i, by=2 )
  for(z in s){
    qij <- qij + choose(i,z) * beta^z*(1-beta)^(i-z) * 0.5^z*choose(z, (z+i-j)/2)
  }
  ## if there is a mis-segregation: ploidy conservation implies that both daughter cells will emerge simultaneously
  #if(i!=j) qij <- 2*qij
  
  return(qij)
}

get_A <- function(maxchrom,beta){
  
  Ndim <- length(maxchrom)
  Nstates <- prod(maxchrom) 
  A<- matrix(0,nrow=Nstates,ncol=Nstates)
  for(i in 1:nrow(A)){
    state_i <- reverse_state_index(i,maxchrom,minchrom=1,ndim=Ndim)
    for(j in 1:nrow(A)){
      state_j <- reverse_state_index(j,maxchrom,minchrom=1,ndim=Ndim)
      qij <- sapply(1:Ndim, function(k) pij(state_i[k], state_j[k], beta))
      ## joint probability (i1,i2,...)->(j1,j2,...)
      qij <- prod(qij)
      A[i,j] <- 2*qij
      
      # ## case when there is no mis-segregation:
      if(i==j) A[i,j] <- (2*qij-1)
      
    }
  }
  A
}
get_fitness <- function(k,pk1,pk2){
  f1 <- as.numeric(sqrt(sum((k-pk1$centre)^2))<=pk1$rad)*pk1$fitness
  f2 <- as.numeric(sqrt(sum((k-pk2$centre)^2))<=pk2$rad)*pk2$fitness
  max(f1,f2)
}
mean_kary <- function(pm,pk1,pk2,return_df=FALSE){
  A <- get_A(maxchrom=c(8,8),beta=pm)
  state_ids <- do.call(rbind,lapply(1:nrow(A), reverse_state_index,maxchrom=c(8,8),ndim=2))
  f <- apply(state_ids,1,get_fitness,pk1=pk1,pk2=pk2)
  #f <- f/max(f)
  A <- A*f
  parms <- list(A=A)
  out <- runsteady(y=rep(1,length(f))/length(f),func=chrmod,parms=parms)
  df <- data.frame(state_ids)
  colnames(df) <- c("c1","c2")
  df$f <- f
  df$x <- out$y
  df$pm <- pm
  df$deltaf <- abs(pk1$fitness-pk2$fitness)
  if(!return_df) return(sum((df$c1)*df$x))
  if(!return_df) return(sum((df$c1+df$c2)*df$x/2))
  return(df)
  #data.frame(cn=1:8,f=f,pm=pm,deltaf=deltaf,x=out[nrow(out),-1])
}


pk1 <- list(centre=c(6,6),rad=1,fitness=1)
pk2 <- list(centre=c(3,3),rad=1,fitness=0.9)

df1 <- rbind(mean_kary(0.1,pk1,pk2,return_df = T),
             mean_kary(0.0001,pk1,pk2,return_df = T))
df1$id <- "scenario A"

pk3 <- list(centre=c(3,6),rad=0,fitness=1)
pk4 <- list(centre=c(6,3),rad=1.99,fitness=0.9)
df2 <- rbind(mean_kary(0.1,pk3,pk4,return_df = T),
             mean_kary(0.0001,pk3,pk4,return_df = T))
df2$id <-"scenario B"

df <- rbind(df1,df2)

df$pm[df$pm==0.0001] <- "0p0001"

transpchar <- function(x) as.character(gsub("p",".",x))

p_xmpl <- ggplot(df,aes(x=c1,y=c2))+
  facet_grid(rows=vars(paste("misrate:",transpchar(pm))),cols=vars(id))+
  geom_raster(aes(fill=f))+
  geom_point(aes(size=x))+
  scale_color_manual("",values=c("black"))+
  scale_size("karyotype \nfrequency")+
  scale_fill_viridis_c("fitness",begin=0.1)+
  scale_x_discrete("chromosome 1 copy number")+
  scale_y_discrete("chromosome 2 copy number")+
  theme_classic()+
  text_size_theme
p_xmpl
dom_smm <- function(dfi){
  dfi <- dfi[unique(apply(dfi,2,which.max)),,drop=F]
  rownames(dfi) <- LETTERS[1:nrow(dfi)]
  reshape2::melt(dfi)
}

ndom <- function(dfi){
  length(unique(apply(dfi,2,which.max)))
}

# Replace this with your actual list of dataframes
df_list <- readRDS("figures/misc/data/misseg_rate_test.Rds")
tmp <- sapply(df_list,ndom)
table(tmp)
sum(tmp>1)
length(tmp)
mean(tmp>1)

df <- do.call(rbind,lapply(names(df_list),function(i){
  res <- dom_smm(df_list[[i]])
  res$id <- i
  res
}))

p <- ggplot(df,aes(x=Var2,y=value,color=Var1,group=Var1))+
  facet_wrap(~id,scales="free")+
  geom_step()+
  scale_y_log10()
p



if(FALSE){
  pscreen <- lapply(names(df_list),function(i){
    tryCatch({
      print(i)
      dfi <- df_list[[i]]
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
}




procDir <- "data/salehi/alfak_outputs_V1a_procv0/"
files <- list.files(procDir, full.names = TRUE)
# Read all files and keep only those with 13 columns
x_list <- lapply(files, readRDS)
x_list <- x_list[sapply(x_list, ncol) == 13]
x <- do.call(rbind, x_list)
# Keep an original copy for later violin plot
x0 <- x

# Filter data and select observation with maximum xv per file and per clean_id
x <- x[!is.na(x$xv) & x$xv > 0 & x$test_treat != "NaN", ]
x <- do.call(rbind, lapply(split(x, x$fi), function(df) df[df$xv == max(df$xv), ]))
rownames(x) <- NULL
x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- do.call(rbind, lapply(split(x, x$clean_id), function(df) head(df[df$xv == max(df$xv), ],1)))
sel <- names(df_list)[c(27,9)]

y <- lapply(sel,function(i){
  set.seed(42)
  dfi <- df_list[[i]]
  dfi <- dfi[apply(dfi,1,max)>0.001,]
  k <- do.call(rbind,lapply(rownames(dfi), s2v))
  print(nrow(k))
  u <- umap::umap(k)$layout
  u <- data.frame(u)
  
  xi <- x[x$fi==i,]
  lscape <- readRDS(paste0("data/salehi/alfak_outputs_V1a/minobs_",xi$min_obs,"/",xi$fi,"/landscape.Rds"))
  lut <- lscape$mean
  names(lut) <- lscape$k
  
  u$f <- lut[rownames(dfi)]
  
  u <- cbind(u,dfi)
  u <- reshape2::melt(u,id.vars=c("X1","X2","f"))
  
  
  id <- dfi[,1]>dfi[,6]
  ploidy <- apply(k,1,mean)
  
  dfp <- data.frame(ploidy,id)
  dfp$f <- lut[rownames(dfp)]
  

  
  tmp <- data.frame(cbind(dfi,id),check.names = F)
  tmp <- reshape2::melt(tmp,id.vars="id")
  sm <- aggregate(list(freq=tmp$value),by=list(id=tmp$id,m=tmp$variable),sum)
  
  sm$lscape <- i
  dfp$lscape <- i
  u$lscape <- i
 
  return(list(sm=sm,dfp=dfp,u=u))
  
})


check_abm <- function(p_error,si){
  minobs <- x$min_obs[x$fi==si]  
  lscape_path <- paste0("data/salehi/alfak_outputs_V1a/minobs_",minobs,"/",si,"/landscape.Rds")
  lscape <- readRDS(lscape_path)
  
  sim_times=c(0,10000)
  x0 <- rep(100,nrow(lscape))
  x0 <- x0/sum(x0)
  names(x0) <- lscape$k
  abm_res <- predict_evo(lscape = lscape,
                                                    p = p_error,
                                                    times = sim_times,
                                                    x0 = x0,
                                                    prediction_type = "ABM",
                                                    abm_pop_size = 50000,
                                                    abm_delta_t = 1,
                                                    abm_max_pop = 1e6,
                                                    abm_record_interval = 100, # Record every 1 time unit (10 * 0.1 dt)
                                                    rcpp_path = "tests/prediction/abm_core.cpp") 
  
  return(abm_res)
}

if(TRUE){
 
  p_range <- seq(0.001,0.005,0.001)
  abm_res_A <- lapply(p_range,check_abm,si=sel[1])
  p_range <- seq(0.001,0.005,0.001)
  abm_res_B <- lapply(p_range,check_abm,si=sel[2])
  
  abm_res <- list(abm_res_A,abm_res_B)
  names(abm_res) <- sel
  
  saveRDS(abm_res,"figures/misc/data/misseg_rate_abm.Rds")
}else{
  abm_res <- readRDS("figures/misc/data/misseg_rate_abm.Rds")
}

abm_res[[1]] <- do.call(rbind,lapply(abm_res[[1]],function(ai) tail(ai[,-1],1)))
abm_res[[2]] <- do.call(rbind,lapply(abm_res[[2]],function(ai) tail(ai[,-1],1)))

y <- lapply(sel,function(i){
#  set.seed(42)
  dfi <- abm_res[[i]]
  dfi <- t(dfi[,apply(dfi,2,max)>0.0005])
  colnames(dfi) <- seq(0.001,0.005,0.001)
  k <- do.call(rbind,lapply(rownames(dfi), s2v))
  print(nrow(k))
  u <- umap::umap(k)$layout
  u <- data.frame(u)
  
  xi <- x[x$fi==i,]
  lscape <- readRDS(paste0("data/salehi/alfak_outputs_V1a/minobs_",xi$min_obs,"/",xi$fi,"/landscape.Rds"))
  lut <- lscape$mean
  names(lut) <- lscape$k
  
  u$f <- lut[rownames(dfi)]
  
  u <- cbind(u,dfi)
  u <- reshape2::melt(u,id.vars=c("X1","X2","f"))
  
  
  id <- dfi[,1]>dfi[,5]
  ploidy <- apply(k,1,mean)
  
  dfp <- data.frame(ploidy,id)
  dfp$f <- lut[rownames(dfp)]
  
  
  
  tmp <- data.frame(cbind(dfi,id),check.names = F)
  tmp <- reshape2::melt(tmp,id.vars="id")
  sm <- aggregate(list(freq=tmp$value),by=list(id=tmp$id,m=tmp$variable),sum)
  
  sm$lscape <- i
  dfp$lscape <- i
  u$lscape <- i
  
  return(list(sm=sm,dfp=dfp,u=u))
  
})


y[[1]]$u$f <- y[[1]]$u$f-mean(y[[1]]$u$f)
y[[2]]$u$f <- y[[2]]$u$f-mean(y[[2]]$u$f)

u <- rbind(y[[1]]$u,y[[2]]$u)
u <- u[order(u$value,decreasing=F),]

lut <- c("p53 k.o. A\n(scenario A)","SA532\n(scenario B)")
names(lut) <- sel

tmp <- u[u$variable%in%c(0.001,0.005),]
tmp <- split(tmp,f=interaction(tmp$variable,tmp$lscape))
tmp <- do.call(rbind,lapply(tmp,function(ti) ti[ti$value==max(ti$value),]))
tmp <- tmp[,!colnames(tmp)%in%c("variable","landscape","value")]
tmp$lab <- c("y","x","y","x")

p_freq <- ggplot(u[u$variable%in%c(0.001,0.005),],aes(x=X1,y=X2))+ ## saturating colors aids visibility
  facet_grid(cols=vars(lut[lscape]),rows=vars(paste0("misrate = ",variable)),scales = "free")+
  geom_point(aes(color=pmin(0.05,pmax(0.000001,value))))+
  scale_color_viridis_c("frequency",breaks=c(0.01,0.001,0.0001,0.00001),trans="log")+
  geom_label_repel(data=tmp,aes(label=lab),box.padding = .5)+
  theme_classic()+
  scale_x_continuous("UMAP 1")+
  scale_y_continuous("UMAP 2")+
  text_size_theme
p_freq


p_fit <- ggplot(u[u$variable%in%c(0.001,0.005),],aes(x=X1,y=X2))+ ## saturating colors aids visibility
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


