setwd("~/projects/008_birthrateLandscape/ALFA-K/")
library(ggplot2)
source("utils/landscape_functions.R")

set.seed(42)
nchrom <- 2
Nwaves <- 10
wavelength <- 1
pk <- gen_rf_landscape(founder = sample(1:10,nchrom,replace = T),Nwaves = Nwaves,wavelength = wavelength)

df <- do.call(rbind,lapply(1:3, function(i){
  pki <- pk[i]
  dfi <- expand.grid(x=seq(0,10,0.25),y=seq(0,10,0.25))
  dfi$f <- apply(dfi,1,function(xi) get_rf_fitness(k=xi,pk=pki,wavelength=wavelength))
  dfi$waveid <- i
  dfi
}))

dfi <- expand.grid(x=seq(0,10,0.25),y=seq(0,10,0.25))
dfi$f <- apply(dfi,1,function(xi) get_rf_fitness(k=xi,pk=pk,wavelength=wavelength))
dfi$waveid <- 10

df <- rbind(df,dfi)

wvnm <- c(" wave 1 ", " wave 2 ", " wave 3 ",rep(), "sum 10\nwaves")
names(wvnm) <- c("1","2","3","10")
df$id <- wvnm[as.character(df$waveid)]

x <- df[df$waveid==10,]

y <- do.call(rbind,pbapply::pblapply(1:nrow(x),function(i){
  dx <- x$x[i]-x$x
  dy <- x$y[i]-x$y
  
  d <- sqrt(dx^2+dy^2)
  pr <- x$f[i]*x$f
  
  data.frame(d,pr)
}))

y <- aggregate(list(pr=y$pr),by=list(d=y$d),mean)

p_grf_demo <- ggplot(df,aes(x=x,y=y,fill=f))+
  facet_grid(cols=vars(id))+
  geom_raster()+
  scale_fill_viridis_c("fitness")+
  scale_x_continuous("chromosome 1",breaks=0:10)+
  scale_y_continuous("chromosome 2",breaks=0:10)+
  theme_classic(base_size = 8)+
  labs(tag="B")+
  guides(fill = guide_colorbar(title = "fitness",
                               frame.colour = "black",
                               barwidth = 0.5,
                               barheight = 3)) 
p_grf_demo