library(ggplot2)
library(wesanderson)




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
  # df$sign_perc <- df$perc
  # for (i in 1:length(df$perc)) {
  #   if (df$L[i] < 0.0) {
  #     df$sign_perc[i] <- -df$sign_perc[i]
  #   }
  # }
  
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
    geom_point(aes(y=trait, x=L, color=sign), size=1.5) +
    geom_errorbarh(aes(y=trait, xmin=hpdl, xmax=hpdu, color=sign), height=0.0, size=1) +
    scale_color_manual(values=pal) + 
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
      panel.grid.minor.x = element_blank(),
      legend.position = "none"
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
  ggsave(plot_name, width=k * 2 + 2, height= height_scale * (n_traits * 0.15 + 1), units="in")
  gc()
}

