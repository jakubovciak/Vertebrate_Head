#### Markos et al 2023
#### functions used in transitions workflow


#### Modified URD::createURD() function to use integrated data

createURD_nn<-function (count.data,norm.data, meta = NULL, min.cells = 1, min.genes = 1, 
          min.counts = 0, gene.max.cut = Inf, max.genes.in.ram = 5000, 
          verbose = T) 
{
  if (verbose) 
    message(paste0(Sys.time(), ": Filtering cells by number of genes."))
  num.genes <- apply(count.data, 2, function(x) sum(x > 0))
  names(num.genes) <- colnames(count.data)
  cells.enough.genes <- names(num.genes[which(num.genes > 
                                                min.genes)])
  shhhh <- gc()
  if (verbose) 
    message(paste0(Sys.time(), ": Filtering genes by number of cells."))
  num.cells <- apply(count.data[, cells.enough.genes], 1, 
                     function(x) sum(x > 0))
  genes.enough.cells <- names(num.cells[which(num.cells > 
                                                min.cells)])
  shhhh <- gc()
  if (verbose) 
    message(paste0(Sys.time(), ": Filtering genes by number of counts across entire data."))
  num.counts <- apply(count.data[, cells.enough.genes], 1, 
                      sum)
  genes.enough.counts <- names(num.counts[which(num.counts > 
                                                  min.counts)])
  shhhh <- gc()
  if (verbose) 
    message(paste0(Sys.time(), ": Filtering genes by maximum observed expression."))
  maxgene <- apply(count.data[, cells.enough.genes], 1, max)
  genes.not.above.max <- rownames(count.data)[maxgene <= gene.max.cut]
  shhhh <- gc()
  genes.use <- intersect(intersect(genes.enough.cells, genes.enough.counts), 
                         genes.not.above.max)
  if (verbose) 
    message(paste0(Sys.time(), ": Creating URD object."))
  object <- methods::new("URD", count.data = as(count.data[genes.use, 
                                                           cells.enough.genes], "dgCMatrix"))
  shhhh <- gc()
  # if (verbose) 
  #   message(paste0(Sys.time(), ": Determining normalization factors."))
  # cs <- apply(object@count.data, 2, sum)
  # norm_factors <- (10^ceiling(log10(median(cs))))/cs
  # shhhh <- gc()
  # if (verbose) 
  #   message(paste0(Sys.time(), ": Normalizing and log-transforming the data."))
  # n.chunks <- ceiling(ncol(object@count.data)/max.genes.in.ram)
  # n.cells <- ncol(object@count.data)
  # lognorm.chunks <- lapply(1:n.chunks, function(chunk) {
  #   shhhh <- gc()
  #   i <- (chunk - 1) * max.genes.in.ram + 1
  #   j <- min((chunk * max.genes.in.ram), n.cells)
  #   as(round(log2(sweep(object@count.data[, i:j], 2, norm_factors[i:j], 
  #                       "*") + 1), digits = 2), "dgCMatrix")
  # })
  # shhhh <- gc()
  #object@logupx.data <- do.call(what = "cbind", lognorm.chunks)
  object@logupx.data <- as(norm.data[genes.use,cells.enough.genes], "dgCMatrix")
  if (verbose) 
    message(paste0(Sys.time(), ": Finishing setup of the URD object."))
  initial.group <- factor(unlist(lapply(colnames(object@logupx.data), 
                                        function(x) strsplit(x, "_|-")[[1]][1])))
  names(initial.group) <- colnames(object@logupx.data)
  object@group.ids <- data.frame(initial.group)
  names(object@group.ids) <- "init"
  object@meta <- data.frame(n.Genes = num.genes[colnames(object@count.data)])
  object@meta[, "n.Trans"] <- apply(object@count.data, 2, 
                                    sum)
  if (!is.null(meta)) {
    object@meta <- cbind(object@meta, meta[rownames(object@meta), 
    ])
  }
  if (verbose) 
    message(paste0(Sys.time(), ": All done."))
  return(object)
}

#### Modified URD::plotDim() function to allow UMAP as reduced dimensionality coordinates


plotDim_mod<-function (object, label, label.type = "search", reduction.use = c("umap",
                                                                               "pca", "dm"), dim.x = 1, dim.y = 2, colors = NULL, discrete.colors = NULL,
                       point.size = 1, alpha = 1, point.shapes = F, plot.title = label,
                       legend = T, legend.title = "", legend.point.size = 3 * point.size,
                       label.clusters = F, cells = NULL, x.lim = NULL, y.lim = NULL,
                       color.lim = NULL, na.rm = F, transitions.plot = 0, transitions.alpha = 0.5,
                       transitions.df = NULL)
{
  if (length(reduction.use) > 1)
    reduction.use <- reduction.use[1]
  if (tolower(reduction.use) == "umap") {
    data.plot <- object@tsne.y
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2])
      stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("UMAP", dim.x)
    dim.y <- paste0("UMAP", dim.y)
  }
  else if (tolower(reduction.use) == "pca") {
    data.plot <- object@pca.scores
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2])
      stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("PC", dim.x)
    dim.y <- paste0("PC", dim.y)
    data.plot <- data.plot[, c(dim.x, dim.y)]
  }
  else if (tolower(reduction.use) == "dm") {
    data.plot <- object@dm@eigenvectors
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2])
      stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("DC", dim.x)
    dim.y <- paste0("DC", dim.y)
    data.plot <- as.data.frame(data.plot[, c(dim.x, dim.y)])
  }
  else {
    stop("The reduction provided is invalid.")
  }
  sig.score <- data.for.plot(object, label = label, label.type = label.type,
                             as.color = F, as.discrete.list = T)
  data.plot$SIG <- sig.score[[2]][rownames(data.plot)]
  if (na.rm) {
    data.plot <- data.plot[complete.cases(data.plot), ]
  }
  if (!is.null(cells)) {
    cells <- intersect(cells, rownames(data.plot))
    data.plot <- data.plot[cells, ]
  }
  if (is.null(transitions.plot) || transitions.plot > 0 ||
      !is.null(transitions.df)) {
    if (is.null(transitions.df))
      transitions.df <- edgesFromDM(object, cells = rownames(data.plot),
                                    edges.return = transitions.plot)
    transitions.df$x1 <- data.plot[transitions.df$from,
                                   dim.x]
    transitions.df$x2 <- data.plot[transitions.df$to, dim.x]
    transitions.df$y1 <- data.plot[transitions.df$from,
                                   dim.y]
    transitions.df$y2 <- data.plot[transitions.df$to, dim.y]
    transitions.df$alpha <- transitions.df$weight/max(transitions.df$weight) *
      transitions.alpha
  }
  this.plot <- ggplot(data = data.plot, aes_string(x = dim.x,
                                                   y = dim.y))
  if (!is.null(transitions.df))
    this.plot <- this.plot + geom_segment(inherit.aes = F,
                                          data = transitions.df, aes(x = x1, y = y1, xend = x2,
                                                                     yend = y2, alpha = alpha))
  if (sig.score[[1]]) {
    if (point.shapes) {
      shape.rep <- ceiling(length(unique(data.plot$SIG))/4) +
        1
      this.plot <- this.plot + geom_point(aes(color = SIG,
                                              shape = SIG), size = point.size, alpha = alpha) +
        scale_shape_manual(values = rep(c(0, 2, 8, 9),
                                        shape.rep))
    }
    else {
      this.plot <- this.plot + geom_point(aes(color = SIG),
                                          size = point.size, alpha = alpha, stroke = 0)
    }
    if (!is.null(discrete.colors)) {
      this.plot <- this.plot + scale_color_manual(values = discrete.colors)
    }
  }
  else {
    if (is.null(colors))
      colors <- defaultURDContinuousColors()
    this.plot <- this.plot + geom_point(aes(color = SIG),
                                        size = point.size) + scale_color_gradientn(colors = colors,
                                                                                   limits = color.lim)
  }
  this.plot <- this.plot + labs(title = plot.title, color = legend.title,
                                shape = legend.title)
  this.plot <- this.plot + theme_bw() + theme(panel.grid.minor = element_blank(),
                                              panel.grid.major = element_blank(), plot.title = element_text(face = "bold"))
  if (label.clusters && sig.score[[1]]) {
    data.plot$CLUSTER <- data.plot$SIG
    k.centers <- aggregate(data.plot[, 1:2], by = list(data.plot$CLUSTER),
                           FUN = "mean")
    this.plot <- this.plot + geom_label(data = k.centers,
                                        aes_string(x = dim.x, y = dim.y, label = "Group.1"),
                                        color = "black", alpha = 0.6, show.legend = F)
  }
  if (!legend) {
    this.plot <- this.plot + guides(color = "none", shape = "none")
  }
  else if (sig.score[[1]]) {
    this.plot <- this.plot + guides(color = guide_legend(override.aes = list(size = legend.point.size)))
  }
  this.plot <- this.plot + guides(alpha = "none")
  if (!is.null(x.lim))
    this.plot <- this.plot + xlim(x.lim[1], x.lim[2])
  if (!is.null(y.lim))
    this.plot <- this.plot + ylim(y.lim[1], y.lim[2])
  return(this.plot)
}

#### Modified URD::calcDM() function to allow setting dm_pcs=50 parameters for destiny::DiffusionMap() function

calcDM_mod<-function (object, genes.use = object@var.genes, cells.use = NULL,
          knn = NULL,
          sigma.use = NULL,
          n_local = 5:7,
          distance = c("euclidean","cosine", "rankcor"),
          density.norm = T,
          dcs.store = 100,
          verbose = T,
          dm_pcs=50)
{
  require(destiny,quietly = TRUE)
  if (is.null(genes.use) || length(genes.use) == 0)
    genes.use <- rownames(object@logupx.data)
  if (is.null(cells.use))
    cells.use <- colnames(object@logupx.data)
  data.use <- t(object@logupx.data[genes.use, cells.use])
  rownames(data.use) <- NULL
  colnames(data.use) <- NULL
  data.use <- as.matrix(data.use)
  if (is.null(sigma.use)) {
    sigma.use <- find_sigmas(data.use, steps = 25, verbose = F)@optimal_sigma
    if (verbose)
      print(paste("destiny determined an optimal global sigma of",
                  round(sigma.use, digits = 3)))
  }
  else if (is.numeric(sigma.use)) {
    if (verbose)
      print(paste("Using provided global sigma", sigma.use))
  }
  else if (sigma.use == "local") {
    if (verbose)
      print(paste("Using local sigma."))
  }
  else {
    stop("sigma must either be NULL, 'local', or numeric.")
  }
  if (is.null(knn)) {
    knn <- find_dm_k(length(cells.use))
    if (verbose)
      print(paste("destiny will use", knn, "nearest neighbors."))
  }
  dm <- DiffusionMap(data.use, sigma = sigma.use, k = knn,
                     n_eigs = dcs.store, density_norm = density.norm,
                     distance = distance[1],
                     n_pcs=dm_pcs)
  rownames(dm@eigenvectors) <- cells.use
  rownames(dm@transitions) <- cells.use
  colnames(dm@transitions) <- cells.use
  object <- importDM(object, dm)
  return(object)
}

#### filtering and visualisation of transition matrices in igraph format using ggplot2
#### for arrow parameters see https://github.com/mdhall272/ggarchery#position_attractsegment-allows-you-to-automatically-shave-the-ends-of-arrow-segments

igraph2ggplot <- function(layout_df = NULL, # graph layout object constructed using igraph::layout_(graph_o) or similar
                          graph_o = NULL, # igraph object
                          width_limit_min = 0, # remove edges with lower weight than width_limit_min
                          scale_by_node = FALSE, # scale edge width by each node
                          keep_max_to = FALSE, # prevent disconnecting the graph
                          filter_skip = FALSE, # remove edges skipping a timepoint
                          return_data = FALSE, # return data frame describing the graph
                          scale_coords = TRUE, # scale the node coordinates to c(-1,1)
                          arrow_length = 0.07, # length of edge arrow
                          arrow_offset = 0.375, # space between node coordinate and the arrow tip
                          max_width = 4, # upper limit of edge width scale
                          node_size = 10, # node point size
                          label_col = 'gold4', # color of the node label
                          lab_frame = TRUE, # use geom_text_repel() or geom_label_repel() for node label
                          lab_size = 5 # size of the node label text
                          ) {
  require(igraph)
  require(ggplot2)
  require(ggrepel)
  require(ggarchery)
  
  # scale node coordinates to -1;1 scale (optional)
  
  if (scale_coords) {
    layout_df <-
      apply(layout_df, 2, function(x) {
        scales::rescale(x, to = c(-1, 1))
      })
  }
  
  # reformat node coordinates data frame
  g_layout_df <- as.data.frame(layout_df) %>%
    mutate(node = V(graph_o)$name) %>%
    mutate(time = factor(
      gsub(".*_", '', node),
      levels = c('G4', 'N0', 'N2', 'N5'),
      ordered = TRUE
    ))
  
  if (!'node_name' %in% colnames(g_layout_df)) {
    g_layout_df <- mutate(g_layout_df, node_name = node)
  }
  
  # extract values of the edges (transitions) from the igraph object and reformat to receive node-centered data frame
  
  graph_df <- get.data.frame(graph_o, 'edges')
  
  graph_df$from.x <-
    g_layout_df$V1[match(graph_df$from, g_layout_df$node)]
  graph_df$from.y <-
    g_layout_df$V2[match(graph_df$from, g_layout_df$node)]
  graph_df$to.x <-
    g_layout_df$V1[match(graph_df$to, g_layout_df$node)]
  graph_df$to.y <-
    g_layout_df$V2[match(graph_df$to, g_layout_df$node)]
  
  graph_df <- mutate(
    graph_df,
    from_time = factor(
      gsub(".*_", '', from),
      levels = c('G4', 'N0', 'N2', 'N5'),
      ordered = TRUE
    ),
    to_time = factor(
      gsub(".*_", '', to),
      levels = c('G4', 'N0', 'N2', 'N5'),
      ordered = TRUE
    )
  )
  
  rownames(graph_df) <- paste0(graph_df$from, '_to_', graph_df$to)
  
  # keep edges pointing to the next timepoint
  graph_df <- filter(graph_df, to_time > from_time)
  
  # remove edges skipping a timepoint (optional)
  if (filter_skip) {
    graph_df <-
      filter(graph_df, (as.numeric(to_time) - as.numeric(from_time)) == 1)
  }
  
  # scale transition probabilities for each node separately (optional)
  
  if (scale_by_node) {
    df_f <- data.frame()
    for (from_ct in unique(graph_df$from)) {
      df_sub <- graph_df[graph_df$from == from_ct, , drop = FALSE] %>%
        mutate(weight = scales::rescale(weight))
      df_f <- rbind(df_f, df_sub)
    }
    graph_df <- df_f
  }
  
  # identify edges with highest probability per node to prevent disconnecting the graph
  
  top_to_df <- group_by(graph_df, to) %>%
    filter(weight == max(weight)) %>%
    as.data.frame()
  
  rownames(top_to_df) <- paste0(top_to_df$from, '_to_', top_to_df$to)
  
  # remove edges with values lower than selected threshold
  
  graph_df <- filter(graph_df, weight >= width_limit_min)
  
  # keep edges with highest probability per node to prevent disconnecting the graph (optional)
  
  if (keep_max_to) {
    graph_df <-
      rbind(graph_df, top_to_df[!rownames(top_to_df) %in% rownames(graph_df), ])
  }
  
  # create plot object
  
  graph_plot <- ggplot() +
    geom_segment(
      data = graph_df,
      aes(
        x = from.x,
        xend = to.x,
        y = from.y,
        yend = to.y,
        size = weight
      ),
      arrow = arrow(type = "closed",
                    length = unit(arrow_length, "inches")),
      position = position_attractsegment(
        start_shave = 0,
        end_shave = arrow_offset,
        type_shave = 'distance'
      )
    ) +
    scale_size_continuous(limits = c(0, NA), range = c(0.25, max_width)) +
    geom_point(data = g_layout_df,
               aes(x = V1, y = V2),
               size = node_size + 1,
               colour = "black") +
    geom_point(data = g_layout_df,
               aes(x = V1, y = V2, color = time),
               size = node_size) +
    scale_x_continuous(expand = expansion(add = 0.2)) +
    scale_y_continuous(expand = expansion(add = 0.2)) +
    coord_fixed() +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )
  
  # add labels in textbox or plain text format
  
  if (lab_frame) {
    graph_plot <- graph_plot +
      geom_label_repel(
        data = g_layout_df,
        aes(x = V1, y = V2, label = node_name),
        alpha = 0.75,
        color = label_col,
        size = lab_size
      )
  } else {
    graph_plot <- graph_plot +
      geom_text_repel(
        data = g_layout_df,
        aes(x = V1, y = V2, label = node_name),
        alpha = 0.75,
        color = label_col,
        size = lab_size
      )
  }
  
  # return format
  
  if (return_data==TRUE) {
    return(graph_df)
  } else {
    return(graph_plot)
  }
}
