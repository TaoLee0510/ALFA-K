setwd("~/projects/008_birthrateLandscape/ALFA-K/figures/alfak_ABM_tests/")
x <- readRDS("fit_summaries_tres.Rds")
x <- x[!sapply(x,is.null)]
x <- x[sapply(x,nrow)==27]
x <- do.call(rbind,lapply(x,function(xi){
  xi$tres <- c(rep(1/3,9),rep(2/3,9),rep(1,9))
  xi
}))
x$w <- sapply(x$cond_id,function(i){
  unlist(strsplit(i,split="_"))[4]
})

p <- ggplot(x[x$min_obs==10,],aes(x=tres,y=r2,))+
  facet_grid(rows=vars(id),cols=vars(w))+
  #geom_jitter(height=0,width=0.05)+
  geom_violin(aes(group=tres))+
  scale_x_continuous("temporal range")+
  scale_y_continuous("Pearson coef.")
p
