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

procDir <- "data/salehi/alfak_outputs_V1a_proc/"
fIn <- paste0(procDir,list.files(procDir))

x <- lapply(fIn,readRDS)
x <- do.call(rbind,x[sapply(x,ncol)==13])
x0 <- x

mean(x$base>x$pred)
x <- x[!is.na(x$xv)&x$xv>0&x$test_treat!="NaN",] 
mean(x$base>x$pred)
x <- split(x,f=x$fi)
x <- do.call(rbind,lapply(x,function(xi){
  xi[xi$xv==max(xi$xv),]
}))
rownames(x) <- NULL
mean(x$base>x$pred)

x$clean_id <- sub("_l_\\d+_d1_\\d+_d2_\\d+$", "", x$fi)
x <- split(x,f=x$clean_id)
x <- do.call(rbind,lapply(x,function(xi) xi[xi$xv==max(xi$xv),]))
mean(x$base>x$pred)

m <- read.csv("data/salehi/metadata.csv")

lut <- m$datasetname
names(lut) <- paste(m$datasetname,m$timepoint,sep="_")

x$pdx <- lut[sapply(x$fi,function(i) strsplit(i,split="_")|>unlist() |>head(-6) |>paste(collapse="_"))]

sample_lineage_map <- c(
  SA1035T   = "SA1035",
  SA1035U   = "SA1035",
  SA532     = "SA532",
  SA609     = "SA609",
  SA609UnBU = "SA609",
  SA906a    = "p53 k.o",
  SA906b    = "p53 k.o",
  SA039U="p53 w.t",
  SA535_CISPLATIN_CombinedT="SA535",
  SA535_CISPLATIN_CombinedU="SA535"
)



x$type <- sample_lineage_map[x$pdx]

x <- split(x,f=x$fi)

x <- lapply(x,function(xi){
  maxp <- 15
  if(xi$type[1]%in%c("p53 k.o","p53 w.t.")) maxp <- 25
  xi <- xi[xi$tt%in%c(1:maxp),]
  xi$eval_tp=FALSE
  xi$eval_tp[xi$tt==maxp] <- TRUE
  xi
})

x <- do.call(rbind,x)
rownames(x) <- NULL

x$win <- x$base>x$pred
x$win[x$metric%in%c("overlap","cosine")] <- (x$base<x$pred)[x$metric%in%c("overlap","cosine")] 
#x[as.numeric(gsub("X","",x$t_train))<7,]

z <- aggregate(list(fwin=x$win),by=list(lineage=x$type,metric=x$metric,day=x$tt),mean)
z2 <- aggregate(list(fwin=x$win),by=list(metric=x$metric,day=x$tt),mean)

tmp <- x[grepl("SA039",x$fi),]
xx <- readRDS("/share/lab_crd/M010_ALFAK_2023/ALFA-K/data/salehi/alfak_inputs/SA039U_X51_l_3_d1_1_d2_0.Rds")$x
yy <- readRDS("/share/lab_crd/M010_ALFAK_2023/ALFA-K/data/salehi/alfak_outputs_V1a/minobs_5/SA039U_X30_l_2_d1_1_d2_1/predictions.Rds")
yy <- yy[,rownames(head(xx))[rownames(head(xx))%in%colnames(yy)]]
xx <- xx[rownames(xx)%in%colnames(yy),]
yy <- yy[-1,]
yy <- yy[1:(21*5),]
yy <- data.frame(yy,check.names = F)
yy$tt <- 30+(1:nrow(yy))/5
yy <- reshape2::melt(yy,id.vars="tt")

for(i in 1:ncol(xx)) xx[,i] <- xx[,i]/sum(xx[,i])
xx <- data.frame(xx,check.names = F)
xx$k <- rownames(xx)
xx <- reshape2::melt(xx,id.vars="k")
colnames(xx)[1:2] <- c("variable","tt")
xx$tt <- as.numeric(as.character(xx$tt))



pp <- ggplot(z,aes(x=day,y=fwin))+
  facet_grid(rows=vars(metric))+
  geom_line(aes(color=lineage,group=lineage))+
  scale_y_continuous("fraction beating baseline")+
  scale_x_continuous("time (days)")+
  theme_minimal()+
  text_size_theme
pp
ggsave("figures/misc/figures/salehi_validation_pp.png",pp,width=5,height=6,units="in")

x <- x[x$eval_tp,]
z <- aggregate(list(fwin=x$win),by=list(lineage=x$type,metric=x$metric),mean)
z$n <- aggregate(list(n=x$win),by=list(lineage=x$type,metric=x$metric),length)$n

# Add fill aesthetic
p0 <- ggplot(z, aes(x = lineage, y = fwin, fill = n)) +
  facet_grid(rows=vars(metric))+
  geom_col() +
  coord_flip() +
  scale_y_continuous("fraction beating baseline",breaks=c(0,0.5,1)) +
  scale_x_discrete("") +
  scale_fill_viridis_c("num.\nsub-lins.",option="magma",breaks=c(2,4,6,8))+
  theme_minimal()+
  text_size_theme

p0
ggsave("figures/misc/figures/salehi_validation_p0.png",p0,width=5,height=8,units="in")


x0 <- x0[!is.na(x0$xv),]
p1 <- ggplot(x0[x0$metric=="overlap"&x0$tt==1,],aes(x=pmax(-1,xv),y=ntrain,group=ntrain))+
  geom_violin()+
  geom_jitter(width=0,height=0.2)+
  scale_y_continuous("training samples",breaks=2:9)+
  scale_x_continuous("CV score")+
  theme_minimal()+
  text_size_theme
p1

ggsave("figures/misc/figures/salehi_validation_p1.png",p1,width=5,height=5,units="in")

p2 <- ggplot(x[x$metric=="wasserstein",],aes(x=base,y=pred0))+
  geom_point()+
  scale_x_continuous("actual distance")+
  scale_y_continuous("predicted distance")+
  theme_minimal()+
  text_size_theme
p2
ggsave("figures/misc/figures/salehi_validation_p2.png",p2,width=4,height=4,units="in")

z <- x[x$metric=="euclidean",]
dSphereAngle <- function(theta, N) {
  coef <- integrate(function(t) sin(t)^(N - 2), lower = 0, upper = pi)$value
  sin(theta)^(N - 2) / coef
}
cSphereAngle <- function(theta, N) {
  if (!is.finite(theta)) return(1)
  integrate(function(t) dSphereAngle(t, N), lower = 0, upper = theta)$value
}
nulldf <- data.frame(angle = 0:180)
nulldf$rads <- pi * nulldf$angle / 180
nulldf$CDF <- sapply(nulldf$rads, cSphereAngle, N = 22)
ggplot(z, aes(x = angle)) +
  stat_ecdf(geom = "step", size = 1, color = "blue") +        # Empirical ECDF
  geom_line(data = nulldf, aes(x = angle, y = CDF),            # Null CDF
            color = "red", size = 1, linetype = "dashed") +
  labs(x = "Angle (degrees)",
       y = "Cumulative Probability") +
  theme_minimal()


# Convert empirical angles from degrees to radians.
z$rads <- z$angle * pi / 180

# Define the null CDF function that takes a vector of radians.
F_null <- function(x) sapply(x, cSphereAngle, N = 22)

# Run the KS test.
# Here, alternative="greater" tests the hypothesis that
# F_empirical(x) > F_null(x) for some x, indicating the empirical values are lower.
ks_test_result <- ks.test(z$rads, F_null, alternative = "greater")
print(ks_test_result)
