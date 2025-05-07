setwd("~/projects/ALFA-K/")
source("utils/ALFA-K.R")

lscape <- readRDS("data/proc/sweep_results/minobs_20_ntp_8/N_22_w_1p6_m_0.00005_rep_99/landscape.Rds")

predictions <- readRDS("data/proc/sweep_results/minobs_20_ntp_8/N_22_w_1p6_m_0.00005_rep_99/predictions.Rds")

xin<- readRDS("data/proc/sweep_inputs.Rds")[["N_22_w_1p6_m_0.00005_rep_99"]]
xin <- xin$data$x
tmp <- head(which(as.numeric(colnames(xin))>1200),10)
tmp <- c(min(tmp)-1:8,tmp)
tmp <- tmp[order(tmp)]
xin <- xin[,tmp]

kpred <- do.call(rbind,lapply(colnames(predictions),s2v))
kin <- do.call(rbind,lapply(rownames(xin),s2v))

for(i in 1:ncol(xin)) xin[,i] <- xin[,i]/sum(xin[,i])
vin <- t(apply(xin,2,function(xi) colSums(xi*kin)))

vpred <- t(apply(predictions,1,function(xi) colSums(xi*kpred)))

vpred <- rbind(vin[8,],vpred)


combined <- rbind(vin, vpred)

pca      <- prcomp(combined, center=TRUE)
coords   <- pca$x[,1:2]
colnames(coords) <- c("PC1","PC2")
n        <- nrow(vin)

# Extract all centroids
time      <- 0:(n-1)
df_actual <- data.frame(coords[1:n, ], time, group="Actual")
df_pred   <- data.frame(coords[(n+1):nrow(coords), ], time=((n+1):nrow(coords))-1, group="Predicted")


# Plot
ggplot() +
  # full trajectories
  geom_path(data=df_actual, aes(PC1, PC2), size=1, color="black") +
  geom_point(data=df_actual, aes(PC1, PC2,color=time), size=2) +
  geom_path(data=df_pred, aes(PC1, PC2), size=1, linetype="dashed", color="red") +
  geom_point(data=df_pred, aes(PC1, PC2,color=time), size=2) +
  theme_void() +
  theme(
    plot.background=element_rect(fill="white", color=NA)
  )


