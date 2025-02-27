melt_for_plotting <- function(data_obj, nclones=5, karyotypes=NULL, fit_obj=NULL){
  if(!is.null(fit_obj) && is.null(karyotypes)){
    nclones <- min(nclones, sum(fit_obj$xo$id == "fq"))
  }
  
  dt <- data_obj$dt
  x <- data_obj$x
  
  # If karyotypes are provided, subset by them
  if(!is.null(karyotypes)){
    x <- x[rownames(x) %in% karyotypes, , drop=FALSE]
  } else { 
    # Otherwise, default to selecting top nclones
    x <- x[order(rowSums(x), decreasing=TRUE), , drop=FALSE]
    x <- head(x, nclones)
  }
  
  # Normalize frequencies
  for(i in 1:ncol(x)) x[,i] <- x[,i] / sum(x[,i])
  
  # Reshape data for plotting
  x <- reshape2::melt(x)
  colnames(x) <- c("karyotype", "time", "frequency")
  x$time <- x$time * dt
  
  out <- list(data=x, fit=NULL)
  
  # Handle fit_obj if provided
  if(!is.null(fit_obj)){
    y <- fit_obj$xo
    y <- y[unique(as.character(x$karyotype)), , drop=FALSE]
    
    t <- seq(min(x$time), max(x$time), length.out=100)
    
    z <- y$u0 + y$f_est %*% t(t)
    z <- apply(z, 2, function(zi) exp(zi) / sum(exp(zi)))
    colnames(z) <- t
    rownames(z) <- rownames(y)
    
    z <- reshape2::melt(z)
    colnames(z) <- c("karyotype", "time", "frequency")    
    out$fit <- z
  }
  return(out)
}

treemaker <- function(z,colorby="Rsq",edgeColorby=NULL,margin = 2){
  library(ggplot2)
  library(dplyr)
  library(viridis)
  
  z$visCol <- z[,colorby]
  
  # Ensure uid and parent are character strings
  z <- z %>% 
    mutate(uid = as.character(uid),
           parent = as.character(parent))
  
  # We'll compute a layout for each phylogeny (each unique value in pdx)
  phylo_ids <- unique(z$pdx)
  
  # Prepare lists to collect node layouts and edge layouts
  layout_list <- list()
  edge_list_all <- list()
  
  # Global horizontal offset for packing trees
  offset_x <- 0
  # horizontal margin between trees
  
  # Process each phylogeny individually
  for (p in phylo_ids) {
    
    # Subset the data for this phylogeny.
    # (z may not list a tree’s root because it appears only as a parent.
    # Therefore, we need the union of the parent and uid values.)
    df <- z %>% filter(pdx == p)
    
    # All nodes in this tree are those that appear as a child or as a parent.
    nodes_in_tree <- unique(c(df$uid, df$parent))
    
    # Build a children map: for every node, list its children.
    # (Nodes with no children will be assigned an empty vector.)
    children_map <- split(df$uid, df$parent)
    for (n in nodes_in_tree) {
      if (!(n %in% names(children_map))) {
        children_map[[n]] <- character(0)
      }
    }
    
    # Identify the root: assume it’s a node that appears as a parent but never as a child.
    root_candidates <- setdiff(unique(df$parent), df$uid)
    if (length(root_candidates) == 0) {
      # Fall back in the (unlikely) event no candidate is found:
      root_candidates <- setdiff(nodes_in_tree, df$uid)
      if (length(root_candidates) == 0)
        root_candidates <- setdiff(nodes_in_tree, df$parent)
    }
    root <- root_candidates[1]
    
    # --- Custom layout using a recursive algorithm ---
    # We use a (local) counter to assign x–positions for leaves.
    x_counter <- 0
    # We’ll collect coordinates in a named list: one entry per node.
    result_coords <- list()
    
    # The recursive function: for a given node and depth,
    # if a leaf then assign the next x value;
    # if not, layout its children and assign x as the mean of its children’s x.
    layout_recursive <- function(node, depth) {
      if (length(children_map[[node]]) == 0) {
        x_counter <<- x_counter + 1
        x_val <- x_counter
      } else {
        # Recurse over all children; note that sapply returns a numeric vector.
        child_xs <- sapply(children_map[[node]], function(child) layout_recursive(child, depth + 1))
        x_val <- mean(child_xs)
      }
      result_coords[[node]] <<- c(x = x_val, y = -depth)
      return(x_val)
    }
    
    # Start recursion from the identified root
    layout_recursive(root, 0)
    
    # In case some nodes were not reached (should not happen in a connected tree),
    # run the recursion for them at depth 0.
    reached <- names(result_coords)
    not_reached <- setdiff(nodes_in_tree, reached)
    if (length(not_reached) > 0) {
      for (n in not_reached) {
        layout_recursive(n, 0)
      }
    }
    
    # Convert the results into a data frame.
    tree_layout <- do.call(rbind, result_coords) %>% as.data.frame()
    tree_layout$uid <- rownames(tree_layout)
    
    # Add the phylogeny’s label (pdx) for later annotation.
    tree_layout$pdx <- p
    
    # Merge in the visCol value from z (nodes that appear only as a parent will have NA)
    tree_layout <- left_join(tree_layout, z %>% select(uid, visCol), by = "uid")
    
    # Mark the root node for labelling (we already computed the root)
    tree_layout$is_root <- tree_layout$uid == root
    
    # Normalize the tree’s x–coordinates so that the leftmost node is at 0.
    tree_layout <- tree_layout %>% mutate(x = x - min(x))
    
    # Record the width of this tree.
    tree_width <- max(tree_layout$x)
    
    # Offset the tree’s x–coordinates by the current global offset.
    tree_layout <- tree_layout %>% mutate(x = x + offset_x)
    
    # Save the layout for this phylogeny.
    layout_list[[p]] <- tree_layout
    
    # --- Build edge data ---
    # For every row in df (an edge from parent to child) we want the coordinates.
    edges_df <- df %>%
      select(parent, uid) %>%
      mutate(parent = as.character(parent),
             uid = as.character(uid)) %>%
      left_join(tree_layout %>% select(uid, x, y), by = c("parent" = "uid")) %>%
      rename(x_parent = x, y_parent = y) %>%
      left_join(tree_layout %>% select(uid, x, y), by = c("uid" = "uid")) %>%
      rename(x_child = x, y_child = y)
    
    edge_list_all[[p]] <- edges_df
    
    # Update the global offset: add the width of this tree and some margin.
    offset_x <- offset_x + tree_width + margin
  }
  
  # Combine all trees’ nodes and edges into single data frames.
  nodes_all <- do.call(rbind, layout_list)
  edges_all <- do.call(rbind, edge_list_all)
  
  if(!is.null(edgeColorby)){
    lut <- z[,edgeColorby]
    names(lut) <- paste0("x",z$uid)
    edges_all$edgeColor <- lut[paste0("x",edges_all$uid)]
  }
  
  # --- Plot with ggplot2 ---
  p <- ggplot() 
    if(is.null(edgeColorby)){
      p <- p+geom_segment(data = edges_all,
                          aes(x = x_parent, y = y_parent, xend = x_child, yend = y_child),
                          color = "grey50")
    }else{
      p <- p+geom_segment(data = edges_all,size=1,
                          aes(x = x_parent, y = y_parent, xend = x_child, 
                              yend = y_child,color=edgeColor))+
        scale_color_discrete(edgeColorby)
    }
    # Draw the edges (branches)
    p <- p +
    # Draw the nodes as circles, filled by visCol using the viridis palette.
    geom_point(data = nodes_all,
               aes(x = x, y = y, fill = visCol),
               shape = 21, size = 3) 
    # Add a legend for visCol (nodes with missing visCol appear in grey50)
    if(is.numeric(nodes_all$visCol)){
      p <- p+scale_fill_viridis_c(colorby,option = "viridis", na.value = "grey50") 
    }else{
      p <- p+scale_fill_viridis_d(colorby,option = "viridis", na.value = "grey50") 
    }
    # Label the root node of each tree with its pdx.
    p <- p+ geom_text(data = nodes_all %>% filter(is_root),
              aes(x = x, y = y, label = pdx),
              nudge_y = 0.5, fontface = "bold", size = 4) +
    theme_minimal() +
      theme(
        legend.box = "horizontal",
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0, 0),  # Bottom-left corner
        legend.justification = c(0, 0),  # Anchor legend to the bottom-left
        legend.direction = "horizontal"  # Arrange legend items horizontally
      )
  return(p)
}

