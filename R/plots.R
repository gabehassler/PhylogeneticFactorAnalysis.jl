library(ggplot2)
library(wesanderson)
library(colorspace)
library(tidyr)

f1 <- function(x) {
  return(x^2.5)
}

f2 <- function(x) {
  return(x^0.25)
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

x <- 255
c11 <- 133
c12 <- 212
c13 <- 227
c21 <- 244
c22 <- 181
c23 <- 189

load_colors <- custom_color_scale(sRGB(c11/x, c12/x, c13/x), sRGB(0.75, 0.75, 0.75), sRGB(c21/x, c22/x, c23/x), f1)
fac_colors <- custom_color_scale(sRGB(1, 0, 0), sRGB(1, 1, 1), sRGB(0, 0, 1), f2)


# colors <- c("#0000FF", "#1515EA", "#2B2BD5", "#4040C0", "#5555AA", "#6B6B95", "#808080", "#956B6B", )

## Loadings plot

plot_loadings <- function(csv_path, plot_name, labels_path = NA, factors = NA,
                          height_scale=1.0, width_scale=1.0,
                          verbose = FALSE,
                          lims = NA){

  if (verbose) {print("Starting loadings plot")}
  if (verbose) {print("Reading data")}

  df  <- read.csv(csv_path, header=TRUE, encoding="UTF-8")
  if (!all(is.na(factors))) {
    if (verbose) {print("Subsetting factors")}
    df <- df[df$factor %in% factors,]
  }

  trait_levs <- unique(df$trait)
  trait_labels <- c()

  if (verbose) {print("X")}

  if (!is.na(labels_path)) {
    labels_df <- read.csv(labels_path, header=TRUE, fileEncoding="UTF-8-BOM")
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[labels_df$trait[i]] = labels_df$pretty[i]
    }

    # trait_levs <- labels_df$pretty
    cat_levs <- unique(labels_df$cat)

    cat_dict <- labels_df$cat
    names(cat_dict) <- labels_df$trait
    df$cat <- cat_dict[df$trait]

    # names(trait_levs) <- labels_df$trait

    # df$trait <- trait_levs[df$trait]
  } else {
    for (i in 1:length(trait_levs)[[1]]) {
      trait_labels[trait_levs[i]] = trait_levs[i]
    }
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
    if (df$hpdu[i] == 0 && df$hpdl[i] == 0) {
      df$sign_perc[i] = NA
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
  if (is.na(lims)) {
      absmax = max(abs(ymin), abs(ymax))
      lims <- c(-absmax, absmax)
  }

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
    scale_color_gradientn(colors = load_colors, limits=c(0, 1), name="probability > 0", na.value="grey75") +
    scale_y_discrete(limits=rev, labels=trait_labels) +
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
    xlim(lims) +
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
  ggsave(plot_name, width=width_scale * (k * 2 + 2), height= height_scale * (n_traits * 0.15 + 1), units="in", limitsize=FALSE)
  gc()
  return(p)
}

# mammals_dir = "C:\\Users\\gabeh\\OneDrive\\SD_storage\\many_traits\\public\\integrated_factors\\data\\mammalian_life_history"
# labels_path = file.path(mammals_dir, "mammals_labels.csv")
# stats_path = file.path(mammals_dir, "mammals", "mammals_loadingsStatistics.csv")
# plot_loadings(stats_path, "test.pdf", labels_path = labels_path)

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
                      new_labels=NA,
                      limits=NA,
                      labels_offset=0) {
  p <- ggtree(tree, layout=layout, open.angle=fan.angle, size = 0)
  if (border){
    # p <- p + geom_tree(size=line_width)
  } else {
    #   p <- p + ggtree(tree)
  }
  if (color_tree) {
    p <- p + geom_tree(aes(color=trait),
                       continuous = TRUE, size=line_width)
    if (is.na(limits)) {
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]])
    } else {
      p <- p + scale_color_gradient2(midpoint = 0.0, low=colors[[1]], mid='grey', high=colors[[2]], limits=limits)

    }

    p <- p + labs(color="factor value")

  } else {
    p <- p + geom_tree(size=line_width)
  }

  if (tip_labels) {
    if (!is.na(new_labels)) {
      p <- p %<+% new_labels + geom_tiplab(aes(label=genus), offset=labels_offset)
    } else {
      p <- p + geom_tiplab(align=FALSE, offset=labels_offset)
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
                             class_path = NA,
                             height=NA,
                             width=NA,
                             border=FALSE,
                             line_width=1.0,
                             tip_labels=TRUE,
                             layout="rectangular",
                             class_palette=NA,
                             combined=TRUE,
                             fan.angle=30.0,
                             new_labels=NA,
                             common_scale=FALSE,
                             scale=TRUE,
                             extra_offset=0,
                             labels_offset=0.02,
                             fac_names=NA,
                             factor_fill = scale_fill_gradient2(midpoint = 0.0, low="blue", mid='white', high="red"),
                             factor_color = scale_color_gradient2(midpoint = 0.0, low="blue", mid='white', high="red"),
                             relabel = NA,
                             include_only = NA
) {

  x <- as.matrix(read.csv(factors_path, header=TRUE))
  k <- ncol(x) - 1
  if (all(is.na(factors))) {
    factors <- 1:k
  }

  base_tree <- read.tree(tree_path)
  taxa <- x[, 1]

  include_class = !is.na(class_path)
  if (include_class) {
    class_df = read.csv(class_path)
  }

  if (!all(is.na(include_only))) {
    drop_taxa = setdiff(taxa, include_only)
    keep_rows = match(include_only, taxa)
    x <- x[keep_rows, ]
    for (taxon in drop_taxa) {
      base_tree <- drop.tip(base_tree, taxon)
    }

    if (include_class) {
      keep_rows <- match(include_only, class_df$taxon)
      class_df <- class_df[keep_rows, ]
    }
  }

  taxa <- x[, 1]

  if (!all(is.na(relabel))) {
    matched_rows <- match(relabel$original, taxa)
    taxa[matched_rows] <- relabel$new

    matched_rows <- match(relabel$original, base_tree$tip.label)
    base_tree$tip.label[matched_rows] <- relabel$new

    if (include_class) {
      matched_rows <- match(relabel$original, class_df$taxon)
      class_df$taxon[matched_rows] <- relabel$new
    }

  }

  x <- x[, factors + 1]
  x <- apply(as.matrix(x), 2, as.numeric)
  rownames(x) <- taxa
  n <- length(taxa)
  if (scale) {
    x <- scale(x)
  }

  limits <- NA
  if (common_scale) {
    x_min = min(x)
    x_max = max(x)
    limits = c(x_min, x_max)
  }






  combined_layout = layout
  if (is.na(height)) {
    height <- 20
    if (layout == "rectangular") {
      height <- 0.2 * n
    }
  }

  if (is.na(width)) {
    width <- 20
  }
  if (layout=="circular") {
    combined_layout = "fan"
  }

  # include_class = !is.na(class_path)
  if (include_class) {

    # stop()
    # class_df = read.csv(class_path)
    class_df[,2] <- as.factor(class_df[,2])
    classes <- data.frame(x = class_df[,2])
    colnames(classes) <- c(colnames(class_df)[2])
    row.names(classes) <- class_df$taxon

    class_fills = discrete_scale("fill", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set3")), name=colnames(classes)[1]) #scale_fill_brewer(palette="Set2")
    class_colors = discrete_scale("color", "manual", colorRampPalette(brewer.pal(min(nrow(unique(classes)), 12), "Set3")), name=colnames(classes)[1])

    if (!all(is.na(class_palette))) {
      class_fills = scale_fill_manual(values=class_palette, name=colnames(classes)[1])
      class_colors = scale_color_manual(values=class_palette, name=colnames(classes)[1])
    }
  }

  n_factors = length(factors)

  if (all(is.na(fac_names))) {
    fac_names <- character(n_factors)
    for (i in 1:n_factors) {
      fac_names[i] = paste("factor", i)
    }
  }


  colnames(x) <- fac_names
  max_height <- max(node.depth.edgelength(base_tree))
  base_tree$edge.length <- base_tree$edge.length / max_height

  heat_width=0.05


  if (combined) {
    p1 <- plot_tree(base_tree, "", border=border, line_width=line_width, tip_labels = tip_labels, layout=combined_layout, color_tree=FALSE, new_labels=new_labels, labels_offset=labels_offset)
    pname <- paste(name, "_factors.svg", sep="")
    p2 <- p1
    if (include_class) {
      p2 <- my_gheatmap(p1, classes, offset=extra_offset, width=heat_width, colnames_angle=90, colnames=TRUE, legend_title=colnames(classes)[1]) + class_fills + class_colors
    }

    p3 <- p2 + new_scale_fill() + new_scale_color()
    p4 <- my_gheatmap(p3, x, offset=extra_offset + heat_width * 2, width=heat_width * n_factors, colnames_angle=90, colnames_offset_y = 0, legend_title="factor value") +
      factor_fill + factor_color +
      labs(fill="factor value") +
      guides(color=FALSE)

    svg(pname, height=height, width=width)
    print(p4)
    dev.off()
    gc()
    return(p4)

  } else {

    for (i in 1:length(factors)) {
      k <- factors[i]
      trait <- x[,k]
      tree <- prep_trait(base_tree, trait)

      p1 <- plot_tree(tree, c("blue", "red"), border=border, line_width=line_width, tip_labels=tip_labels, layout=layout, new_labels=new_labels, limits=limits, labels_offset=labels_offset)
      # p <- gheatmap(p, odf2, offset=0, width=0.05, legend_title="origin") + scale_fill_brewer(palette = "Set2") + scale_x_ggtree() + scale_y_continuous(expand=c(0, 0.3))
      # p <- p + scale_x_continuous(expand = c(.1, .1))
      p2 <- p1 # + new_scale_fill() + new_scale_color()

      if (include_class) {
        # p3 <- my_gheatmap(p2, classes, offset=0, width=heat_width, colnames_angle=90, colnames=TRUE, legend_title=colnames(classes)[1]) + class_fills # + class_colors
        p3 <- gheatmap(p2, classes, offset=extra_offset, width=heat_width, colnames=FALSE, legend_title=colnames(classes)[1]) + class_fills # + class_colors
        # scale_fill_manual(values=class_palette, name=colnames(classes)[1])
        # scale_x_ggtree() +
        # scale_y_continuous(expand=c(0, 0.3))
      } else {
        p3 <- p2
      }

      pname <- paste(name, k, ".svg", sep="")
      print(paste("saving", pname))
      # p <- p_origins %>% insert_left(p, width=8)


      svg(pname, height=height, width=width)
      print(p3)
      dev.off()
      gc()
      if (length(factors) == 1) {
        return(p3)
      }
    }
  }
  gc()
}


my_gheatmap <- function(p, data, offset=0, width=1, low="green", high="red", color="white",
                     colnames=TRUE, colnames_position="bottom", colnames_angle=0, colnames_level=NULL,
                     colnames_offset_x = 0, colnames_offset_y = 0, font.size=4, family="", hjust=0.5, legend_title = "value") {

  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL

  ## if (is.null(width)) {
  ##     width <- (p$data$x %>% range %>% diff)/30
  ## }

  ## convert width to width of each cell
  width <- width * (p$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)

  isTip <- x <- y <- variable <- value <- from <- to <- NULL

  ## handle the display of heatmap on collapsed nodes
  ## https://github.com/GuangchuangYu/ggtree/issues/242
  ## extract data on leaves (& on collapsed internal nodes)
  ## (the latter is extracted only when the input data has data on collapsed
  ## internal nodes)
  df <- p$data
  nodeCo <- intersect(df %>% filter(is.na(x)) %>%
                        select(.data$parent, .data$node) %>% unlist(),
                      df %>% filter(!is.na(x)) %>%
                        select(.data$parent, .data$node) %>% unlist())
  labCo <- df %>% filter(.data$node %in% nodeCo) %>%
    select(.data$label) %>% unlist()
  selCo <- intersect(labCo, rownames(data))
  isSel <- df$label %in% selCo

  df <- df[df$isTip | isSel, ]
  start <- max(df$x, na.rm=TRUE) + offset

  dd <- as.data.frame(data)
  ## dd$lab <- rownames(dd)
  i <- order(df$y)

  ## handle collapsed tree
  ## https://github.com/GuangchuangYu/ggtree/issues/137
  i <- i[!is.na(df$y[i])]

  lab <- df$label[i]
  ## dd <- dd[lab, , drop=FALSE]
  ## https://github.com/GuangchuangYu/ggtree/issues/182
  dd <- dd[match(lab, rownames(dd)), , drop = FALSE]


  dd$y <- sort(df$y)
  dd$lab <- lab
  ## dd <- melt(dd, id=c("lab", "y"))
  dd <- gather(dd, variable, value, -c(lab, y))

  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels=colnames(data))
  } else {
    dd$variable <- factor(dd$variable, levels=colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from=dd$variable, to=V2)
  mapping <- unique(mapping)

  dd$x <- V2
  dd$width <- width
  dd[[".panel"]] <- factor("Tree")
  height <- 1.00
  size <- 0.1
  if (is.null(color)) {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), width=width, height=height, size=size, inherit.aes=FALSE)
  } else {
    p2 <- p + geom_tile(data=dd, aes(x, y, fill=value, color=value), width=width, height=height, size=size, inherit.aes=FALSE)
  }
  if (is(dd$value,"numeric")) {
    p2 <- p2 + scale_fill_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
    p2 <- p2 + scale_color_gradient(low=low, high=high, na.value=NA, name = legend_title) # "white")
  } else {
    p2 <- p2 + scale_fill_discrete(na.value=NA, name = legend_title) #"white")
    p2 <- p2 + scale_color_discrete(na.value=NA, name = legend_title) #"white")
  }

  if (colnames) {
    if (colnames_position == "bottom") {
      y <- 0
    } else {
      y <- max(p$data$y) + 1
    }
    mapping$y <- y
    mapping[[".panel"]] <- factor("Tree")
    p2 <- p2 + geom_text(data=mapping, aes(x=to, y = y, label=from), size=font.size, family=family, inherit.aes = FALSE,
                         angle=colnames_angle, nudge_x=colnames_offset_x, nudge_y = colnames_offset_y, hjust=hjust)
  }

  p2 <- p2 + theme(legend.position="right")
  ## p2 <- p2 + guides(fill = guide_legend(override.aes = list(colour = NULL)))

  if (!colnames) {
    ## https://github.com/GuangchuangYu/ggtree/issues/204
    p2 <- p2 + scale_y_continuous(expand = c(0,0))
  }

  attr(p2, "mapping") <- mapping
  return(p2)
}

# aqui_dir = "C:\\Users\\gabeh\\OneDrive\\SD_storage\\many_traits\\public\\integrated_factors\\data\\mammalian_life_history"
# newick = file.path(aqui_dir, "mammals_newick_processed.txt")
# facs = file.path(aqui_dir, "mammals", "mammals_factorMeans.csv")
# class_path <- file.path(aqui_dir, "mammals_classification.csv")
# # aqui_dir = "C:\\Users\\gabeh\\OneDrive\\SD_storage\\many_traits\\public\\integrated_factors\\data\\aquilegia"
# # newick = file.path(aqui_dir, "aquilegia_newick.txt")
# # facs = file.path(aqui_dir, "aquilegiaBinary", "aquilegiaBinary_factorMeans.csv")
# # class_path = file.path(aqui_dir, "aquilegia_pollinators - Copy.csv")
#
# plot_factor_tree("test", newick, facs, tip_labels = FALSE, layout="circular", height=40, width=40, class_path = class_path)
#
# class_df = read.csv(class_path)
# class_df[,2] <- as.factor(class_df[,2])
# classes <- data.frame(x = class_df[,2])
# colnames(classes) <- c(colnames(class_df)[2])
# # row.names(classes) <- class_df$taxon