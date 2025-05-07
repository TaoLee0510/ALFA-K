library(tidyverse); library(igraph); library(tidygraph); library(ggraph); library(dplyr)
library(ggrepel)  # ensure ggrepel is loaded for geom_text_repel

base_text_size = 16
text_size_theme <- 
  theme(
    text         = element_text(size = base_text_size, family = "sans"),
    axis.title   = element_text(size = base_text_size, family = "sans"),
    axis.text    = element_text(size = base_text_size, family = "sans"),
    legend.title = element_text(size = base_text_size, family = "sans"),
    legend.text  = element_text(size = base_text_size, family = "sans"),
    strip.text   = element_text(size = base_text_size, family = "sans")
  )

# Load libraries
setwd("/share/lab_crd/M010_ALFAK_2023/ALFA-K")
#df <- read.csv("data/salehi/metadata.csv", stringsAsFactors = FALSE)  # Load metadata

prune_children <- function(df) {
  while(sum(is.na(df$xval) & !df$uid %in% df$parent) > 0){
    df <- df[!(is.na(df$xval) & !df$uid %in% df$parent),]
  }
  df
}

sample_lineage_map <- c(
  SA1035   = "SA1035",
  SA532    = "SA532",
  SA609    = "SA609",
  SA906    = "p53 k.o",
  SA039    = "p53 w.t",
  SA535    = "SA535"
)

df <- readRDS("figures/misc/data/annotated_metadata.Rds")
df <- prune_children(df)
df$linlab <- sample_lineage_map[df$PDX_id]

# Create dummy root and combine with data; fix missing parent values
dummy_uid <- "ROOT"; 
dummy_row <- df[1,]
dummy_row$uid <- dummy_uid; dummy_row$datasetname = "ROOT"
df_fixed <- rbind(df, dummy_row); 

df_fixed$parent[is.na(df_fixed$parent) | df_fixed$parent == ""] <- dummy_uid

# Assign lineage: specific mapping; others get "TNBC PDX"; mark dummies
lm <- c("SA039U" = "p53 wt", "SA906a" = "p53 ko", "SA906b" = "p53 ko")
df_fixed$lineage <- ifelse(df_fixed$datasetname %in% names(lm), lm[df_fixed$datasetname], "TNBC PDX")
df_fixed$lineage[grepl("x", df_fixed$uid)] <- "dummy"
df_fixed$lineage[df_fixed$uid == "ROOT"] <- "dummy"
df_fixed <- cbind(df_fixed[, c("uid", "parent")], df_fixed[, !colnames(df_fixed) %in% c("uid", "parent")])

# Build edge list (flag edges from dummy) and create graph
edges_df <- df_fixed %>% 
  filter(uid != dummy_uid) %>% 
  mutate(from_dummy = (parent == dummy_uid)) %>% 
  select(from = parent, to = uid, treatment = on_treatment, from_dummy)
g <- graph_from_data_frame(d = edges_df, vertices = df_fixed, directed = TRUE)

# Convert to tidygraph; mark dummy nodes and classify node types
tg <- as_tbl_graph(g) %>% 
  mutate(dummy = (lineage == "dummy"),
         node_type = ifelse(dummy, "dummy", ifelse(parent == dummy_uid, "root", "non-root")))

# Define color palette for lineage (white for dummy)
lineage_colors <- c("dummy" = "white", "p53 wt" = "#E69F00", "p53 ko" = "#56B4E9", "TNBC PDX" = "#009E73")
set.seed(42) 
layout <- create_layout(tg, layout = "dendrogram", circular = TRUE, 
                        height = -node_distance_to(1, mode = "all"))

tmp <- layout$x
layout$x <- layout$y
layout$y <- tmp
layout$y <- -layout$y
# <<--- Minimal change: new lineage label logic (filter points with parent equal to dummy_uid) --->
# Compute label positions from the layout by filtering for nodes with parent == "ROOT"
# Define a named vector with multipliers for each label (linlab)
radius_multipliers <- c(
  "p53 w.t" = 5,
  "SA1035"  = 6,
  "SA532"   = 10,
  "SA535"   = 8,
  "SA609"   = 5,
  "p53 k.o" = 8
)

# Compute label positions from the layout (only the nodes with parent equal to dummy_uid)
lineage_labels <- layout %>% 
  filter(parent == dummy_uid) %>% 
  mutate(
    r = sqrt(x^2 + y^2),            # current radial distance from center
    theta = atan2(y, x),            # angle in radians for each node
    r_multiplier = radius_multipliers[linlab],  # look up the multiplier based on linlab text
    r_label = r * r_multiplier,     # new radial distance for labels
    x_label = r_label * cos(theta), # new x coordinate for labels
    y_label = r_label * sin(theta)  # new y coordinate for labels
  )




# Plot with added connector lines and adjustable label positions
# Adding connector lines and the repelled text labels:
pp <- ggraph(layout) + 
  
  # Connector lines from node to label
  geom_segment(data = lineage_labels,
               aes(x = x, y = y, xend = x_label, yend = y_label),
               color = "gray80", linetype = "dashed") +
  
  geom_edge_link(aes(color = ifelse(from_dummy, "dummy", treatment)), edge_width = 1) +
  
  # Plotting dummy nodes
  geom_node_point(data = filter(layout, dummy),
                  aes(x = x, y = y), color = "white") +
  
  # Plotting non-dummy nodes (using xval color scale)
  geom_node_point(data = filter(layout, !dummy),
                  aes(x = x, y = y, color = xval), size = 4) +
  scale_edge_color_manual(name = "",
                          values = c("y" = "red", "n" = "black", "dummy" = "white"),
                          labels = c("y" = "Treatment ON", "n" = "Treatment OFF", "dummy" = "")) +
  scale_color_viridis_c("CV\nscore") +
  
  
  # Repelled text labels at the computed positions
  geom_text(data = lineage_labels,
            aes(x = x_label, y = y_label, label = linlab),
            fontface = "bold") +
  
  theme(legend.position = c(0.01, 0.01),
        legend.justification = c("left", "bottom"),
        legend.box = "vertical",
        legend.spacing.y = unit(0, "cm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text         = element_text(size = base_text_size, family = "sans"),
        legend.title = element_text(size = base_text_size, family = "sans"),
        legend.text  = element_text(size = base_text_size, family = "sans"),
        strip.text   = element_text(size = base_text_size, family = "sans"),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA))

ggsave("figures/misc/figures/ABM_lineages.png", plot = pp, width = 6, height = 6, units = "in")
