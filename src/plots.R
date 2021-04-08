library(ggplot2)
library(wesanderson)
library(colorspace)

f1 <- function(x) {
  return(x^2.5)
}

custom_color_scale <- function(low, mid, high, f, n = 10) {
  lows = rev(half_scale(mid, low, f, n))
  highs = half_scale(mid, high, f, n)
  return(c(lows, hex(mid), highs))
}

half_scale <- function(c1, c2, f, n) {
  colors <- character(n)
  for (i in 1:n) {
    colors[i] <- hex(mixcolor(f(i / n), c1, c2))
  }
  return(colors)
}



# load_colors <- custom_color_scale("#0000FF", "#808080", "#FF0000", f1)
load_colors <- custom_color_scale(RGB(0, 0, 1), RGB(0.5, 0.5, 0.5), RGB(1, 0, 0), f1)



# colors <- c("#0000FF", "#1515EA", "#2B2BD5", "#4040C0", "#5555AA", "#6B6B95", "#808080", "#956B6B", )

## Loadings plot

plot_loadings <- function(csv_path, plot_name, labels_path = NA, height_scale=1.0){
  df  <- read.csv(csv_path, header=TRUE, encoding="UTF-8")
  
  if (!is.na(labels_path)) {
    labels_df <- read.csv(labels_path, header=TRUE, encoding="UTF-8")
    trait_levs <- labels_df$pretty
    cat_levs <- unique(labels_df$cat)
    
    cat_dict <- labels_df$cat
    names(cat_dict) <- labels_df$trait
    df$cat <- cat_dict[df$trait]

    names(trait_levs) <- labels_df$trait

    df$trait <- trait_levs[df$trait]
  } else {
    trait_levs <- unique(df$trait)
    df$cat <- rep("NA", nrow(df))
    cat_levs <- c("NA")
  }

  wes_pal <- wes_palette("GrandBudapest2")
  pal <- c(wes_pal[1], "grey", wes_pal[4])
  df$trait <- factor(df$trait, levels=trait_levs)
  df$cat <- factor(df$cat, levels=cat_levs)
  df$L <- sapply(df$L, as.numeric)
  df$sign_perc <- df$perc
  for (i in 1:length(df$perc)) {
    if (df$L[i] < 0.0) {
      df$sign_perc[i] <- 1.0 - df$sign_perc[i]
    }
  }
  
  df$sign <- sapply(df$L, sign)
  n <- dim(df)[[1]]
  for (i in 1:n) {
    if (df$hpdl[i] <= 0 && df$hpdu[i] >= 0) {
      df$sign[i] = 0
    }
  }
  df$sign <- factor(df$sign, levels=c(1,0, -1))
  ymin = min(df$hpdl)
  ymax = max(df$hpdu)
  absmax = max(abs(ymin), abs(ymax))
  
  k <- max(df$factor)
  facet_labels <- character(k)
  facet_names <- integer(k)
  for (i in 1:k) {
    facet_labels[i] = paste("factor", i)
    facet_names[i] = i
  }
  names(facet_labels) <- facet_names
  
  n_groups <- length(cat_levs)
  group_labels <- character(n_groups)
  if (length(n_groups) == 1) {
    names(group_labels) <- cat_levs
  }
  
  # ps <- c()
  # for (i in 1:max(df$row)) {
  #   ps <- c(ps, plot_single_row(df, 1, ymin, ymax))
  # }

  
  p <- ggplot(df) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_point(aes(y=trait, x=L, color=sign_perc), size=1.5) +
    geom_errorbarh(aes(y=trait, xmin=hpdl, xmax=hpdu, color=sign_perc), height=0.0, size=1) +
    scale_color_gradientn(colors = load_colors, limits=c(0, 1), name="probability > 0") +
    scale_y_discrete(limits=rev) +
    # scale_color_gradient2(low="orange", mid="white", high="purple", limits=c(-1, 1), name="L") +
    #facet_grid(~ cat, scales="free_x", space="free_x") +
    #geom_tile() +
    #scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
    #scale_x_discrete(position = "top") +
    # scale_y_discrete() +
    labs(y="", x="loadings value") +
    theme_minimal() +
    theme(#axis.text.x = element_text(angle=0, hjust=1),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
      # legend.position = "none"
    ) +
    xlim(-absmax, absmax) +
    facet_grid(rows=vars(cat),
               cols=vars(factor),
               labeller = labeller(factor=facet_labels),
               scales="free_y",
               space="free_y",
               switch = "y")
  
  if (length(cat_levs) == 1) {
    p <- p + theme(strip.text.y = element_blank())
  }
  
  # axis.title.y = element_text())
  n_traits <- length(trait_levs)
  ggsave(plot_name, width=k * 2 + 2, height= height_scale * (n_traits * 0.15 + 1), units="in", limitsize=FALSE)
  gc()
}

# plot_loadings("nwmUpScaled_loadingsStatistics.csv", "test.pdf")

## Factor plot


library(ggtree)
library(phytools)
library(tidytree)
library(ggplot2)
library(aplot)
library(RColorBrewer)
library(ggnewscale)
library(phyclust)
require(treeio)

prep_trait <- function(tree, trait){
  fit <- phytools::fastAnc(tree, trait, vars=TRUE, CI=TRUE)
  td <- data.frame(node = nodeid(tree, names(trait)),
                   trait = trait)
  nd <- data.frame(node = names(fit$ace), trait = fit$ace)
  d <- rbind(td, nd)
  d$node <- as.numeric(d$node)
  tree <- full_join(tree, d, by = 'node')
  return(tree)
}

plot_tree <- function(tree, colors, 
                      border=FALSE,
                      line_width=0.1,
                      tip_labels=FALSE,
                      layout="rectangular",
                      color_tree=TRUE,
                      fan.angle=15,
                      new_labels=NA) {
  p <- ggtree(tree, layout=layout, open.angle=fan.angle, size = 0)
  if (border){
    # p <- p + geom_tree(size=line_width)
  } else {
    #   p <- p + ggtree(tree)
  }
  if (color_tree) {
    p <- p + geom_tree(aes(color=trait),
                       continuous = TRUE, size=line_width) +
      scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]]) +
      labs(color="factor value")
  } else {
    p <- p + geom_tree(size=line_width)
  }
  
  if (tip_labels) {
    if (!is.na(new_labels)) {
      p <- p %<+% new_labels + geom_tiplab(aes(label=genus), offset=0.02)
    } else {
      p <- p + geom_tiplab(aligh=TRUE)
    }
  } else {
    # p <- p + geom_tiplab(aes(label=character(1)), align=TRUE)
  }
  
  # p <- ggtree(tree, size=2.8) + 
  #   geom_tree(aes(color=trait),
  #             continuous = TRUE, size=2) +
  #   scale_color_gradient2(midpoint = 0, low='purple', mid='white', high='orange')
  return(p)
}

plot_factor_tree <- function(name, tree_path, factors_path, factors = NA,
                             classes = NA,
                             height=20,
                             width=20,
                             border=FALSE,
                             line_width=1.0,
                             tip_labels=TRUE,
                             layout="rectangular",
                             class_palette=NA,
                             combined=TRUE,
                             fan.angle=30.0,
                             new_labels=NA
) {

  x <- as.matrix(read.csv(factors_path, header=TRUE))
  k <- ncol(x) - 1
  if (is.na(factors)) {
    factors <- 1:k
  }
  
  taxa <- x[, 1]
  x <- x[, factors + 1]
  x <- apply(as.matrix(x), 2, as.numeric)
  rownames(x) <- taxa
  n <- length(taxa)
  x <- scale(x)
  
  combined_layout = layout
  if (layout=="circular") {
    combined_layout = "fan"
  }
  
  include_class = !all(is.na(classes))
  class_colors = discrete_scale("fill", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set3")), name=colnames(classes)[1]) #scale_fill_brewer(palette="Set2")
  if (!all(is.na(class_palette))) {
    class_colors = scale_fill_manual(values=class_palette, name=colnames(classes)[1])
  }
  
  n_factors = length(factors)
  fac_names <- character(n_factors)
  for (i in 1:n_factors) {
    fac_names[i] = paste("factor", i)
  }
  colnames(x) <- fac_names
  base_tree <- read.tree(tree_path)
  max_height <- max(node.depth.edgelength(base_tree))
  base_tree$edge.length <- base_tree$edge.length / max_height
  
  heat_width=0.05
  
  
  if (combined) {
    p1 <- plot_tree(base_tree, "", border=border, line_width=line_width, tip_labels = tip_labels, layout=combined_layout, color_tree=FALSE, new_labels=new_labels)
    pname <- paste(name, "_factors.svg", sep="")
    p2 <- p1
    if (include_class) {
      p2 <- gheatmap(p1, classes, offset=0, width=heat_width, colnames_angle=90, colnames=TRUE, legend_title=colnames(classes)[1]) + class_colors
    }
    
    p3 <- p2 + new_scale_fill()
    p4 <- gheatmap(p3, x, offset=heat_width * 2, width=heat_width * n_factors, colnames_angle=90, colnames_offset_y = 0, legend_title="factor value") +
      scale_fill_gradient2(midpoint = 0.0, low="blue", mid='grey', high="red") + labs(fill="factor value")
    svg(pname, height=height, width=width)
    print(p4)
    dev.off()
    
  } else {
    
    for (i in 1:length(factors)) {
      k <- factors[i]
      trait <- x[,k]
      tree <- prep_trait(base_tree, trait)
      
      p <- plot_tree(tree, c("blue", "red"), border=border, line_width=line_width, tip_labels=tip_labels, layout=layout, new_labels=new_labels)
      # p <- gheatmap(p, odf2, offset=0, width=0.05, legend_title="origin") + scale_fill_brewer(palette = "Set2") + scale_x_ggtree() + scale_y_continuous(expand=c(0, 0.3))
      # p <- p + scale_x_continuous(expand = c(.1, .1))
      pname <- paste(name, k, ".svg", sep="")
      
      if (include_class) {
        p <- gheatmap(p, classes, offset=0, width=heat_width, colnames=FALSE, legend_title=colnames(classes)[1]) + class_colors
        # scale_fill_manual(values=class_palette, name=colnames(classes)[1])
        # scale_x_ggtree() +
        # scale_y_continuous(expand=c(0, 0.3))
      }
      
      print(paste("saving", pname))
      # p <- p_origins %>% insert_left(p, width=8)
      
      
      svg(pname, height=height, width=width)
      print(p)
      dev.off()
    }
  }
  gc()
}



