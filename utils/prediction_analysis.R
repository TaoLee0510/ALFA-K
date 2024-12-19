
compute_probabilities <- function(counts, alpha) {
  total_counts <- sum(counts) # Total population size
  unique_karyotypes <- length(counts) # Number of unique karyotypes
  probabilities <- (counts + alpha) / (total_counts + alpha * unique_karyotypes)
  return(log(probabilities))
}
# Function to compute weighted total probability
compute_total_probability <- function(x_col, y_col, x_initial, alpha, log_transform = FALSE) {
  # Compute smoothed probabilities for y[,i] and x[,0]
  
  tmp <- names(x_col)[!names(x_col)%in%names(y_col)]
  y_extra <- rep(0,length(tmp))
  names(y_extra) <- tmp
  y_col <- c(y_col,y_extra)
  
  tmp <- names(x_col)[!names(x_col)%in%names(x_initial)]
  x_extra <- rep(0,length(tmp))
  names(x_extra) <- tmp
  x_initial <- c(x_initial,x_extra)
  
  p_y <- compute_probabilities(y_col, alpha)
  p_x0 <- compute_probabilities(x_initial, alpha)
  
  
  # Extract probabilities for observed karyotypes
  prob_y <- p_y[names(x_col)]
  prob_x0 <- p_x0[names(x_col)]

  
  # Compute weighted total probabilities
  total_prob_y <- sum(x_col * prob_y)
  total_prob_x0 <- sum(x_col * prob_x0)
  
  # Return results
  return(data.frame(
    alpha=alpha,
    Total_Prob_Y = total_prob_y,
    Total_Prob_X0 = total_prob_x0
  ))
}

setwd("~/projects/008_birthrateLandscape/ALFA-K/")
source("utils/comparison_functions.R")
source("utils/ALFA-K.R")

x <- proc_sim("data/main/N_22_w_1p6_m_0.00005_rep_08/train/00000/",times=seq(2000,3000,length.out=11))$x
y <- proc_sim("data/main/N_22_w_1p6_m_0.00005_rep_08/test_v2/minobs_20_ntp_8_00000/00000/",times=seq(0,1000,length.out=11))$x

x <- x[order(x[,1],decreasing = T),]
y <- y[order(y[,1],decreasing = T),]

# Example usage
# Assuming `x` and `y` are matrices or data frames with karyotype counts
alpha <- 1 # Smoothing parameter
i <- 3 # Choose a specific timepoint
log_transform <- TRUE # Set to TRUE if you want to log-transform probabilities


# Compute total probabilities
total_probs <- do.call(rbind,lapply(2:10,function(i){
  do.call(rbind,lapply(exp(-5:5),function(ai){
    df <- compute_total_probability(x[x[,i]>0, i], y[y[,i]>0, i], x[x[,1]>0, 1], ai, log_transform)
    df$timepoint <- i
    return(df)
  }))
}))

z <- reshape2::melt(total_probs,id.vars=c("alpha","timepoint"))

p <- ggplot(z,aes(x=alpha,y=value,color=variable,group=variable))+
  geom_line()+
  facet_wrap(~timepoint,scales="free")+
  scale_x_log10()
p

