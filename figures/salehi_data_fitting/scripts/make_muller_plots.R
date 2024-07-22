# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(viridis)

# Define the data manually
data <- data.frame(
  Time = rep(c(1, 2, 3, 4), times = 5),
  Clone = rep(c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5"), each = 4),
  Frequency = c(10, 5, 3, 1,
                8, 4, 2, 1,
                4, 2, 0, 0,
                0, 2, 4, 8,
                0, 0, 1, 3)
)

# Normalize frequencies at each timepoint to sum to 1
data <- data %>%
  group_by(Time) %>%
  mutate(Frequency = Frequency / sum(Frequency))

# Add group information
data <- data %>%
  mutate(Group = case_when(
    Clone %in% c("Clone1", "Clone2") ~ "Group A",
    Clone == "Clone3" ~ "Group B",
    Clone %in% c("Clone4", "Clone5") ~ "Group C"
  ))

clone_names <- c("Clone1", "Clone2", "Clone3", "Clone4", "Clone5")
zeta <- clone_names[1:3]
psi <- clone_names[c(1:2,4:5)]
psiNzeta <- clone_names[c(4:5)]

groups <- list(zeta=zeta,
               psi=psi,
               psiNzeta=psiNzeta)
tmp <- do.call(rbind,lapply(unique(data$Group),function(di) {
  tmp <- data
  tmp$Group2 <- di
  tmp
  }))

tmp$Clone2 <- FALSE
tmp$Clone2[tmp$Group2=="Group A" & tmp$Clone%in%groups[[1]]] <- TRUE
tmp$Clone2[tmp$Group2=="Group B" & tmp$Clone%in%groups[[2]]] <- TRUE
tmp$Clone2[tmp$Group2=="Group C" & tmp$Clone%in%groups[[3]]] <- TRUE

lablr <- c('Group B' = " \u03A8 ",
           'Group A' = "  \u03B6  ",
           'Group C' = "\u0060\u03B6\u2229\u03A8")

tmp$Group2 <- lablr[tmp$Group2]

# Create a Muller plot using ggplot2
ggplot(tmp, aes(x = Time, y = Frequency, group=Clone,fill = Clone,alpha=Clone2)) +
  geom_area(color="black") +
  scale_fill_viridis_d("karyotype",labels=letters[1:5])+
  guides(alpha="none")+
  scale_alpha_manual(values=c(0,1))+
  facet_wrap(~ Group2, ncol = 3) +
  theme_classic()+
  scale_x_continuous("",breaks=c(1,4),labels=c(expression(S[0],S[t])))+
  scale_y_continuous("clone frequency")

