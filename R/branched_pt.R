#' Branched pt function
#'
#' This function reconstructs transcriptional trajectories between a developmental root and the tips.
#' @param env An enviroment produced by the oposSOM pipeline
#' @param root A single group label defined as the developmental source
#' @param tips A vector of group labels defined as the developmental sinks
#' @keywords urd
#' @import URD oposSOM
branched_pt <- function(env, root, tips) {
  filename <- file.path(
    paste(env$files.name, "- Results"),
    "Pseudotemporal Analysis",
    "Cell State Trajectories",
    "Report.pdf"
  )
  util.info("Writing:", filename)
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)

  discrete.colors <- env$groupwise.group.colors[sort(names(env$groupwise.group.colors))]
  obj <- new("URD")
  obj@count.data <- as(as.matrix(env$metadata), "dgCMatrix")
  obj@logupx.data <- as(as.matrix(env$metadata), "dgCMatrix")
  rownames(obj@logupx.data) <- as.character(1:nrow(obj@logupx.data))
  obj@group.ids <- data.frame(
    group.labels = env$group.labels,
    pat.labels = env$pat.labels
  )

  
  var.genes <- sapply(env$spot.list.group.overexpression$spots, function(x){x$metagenes})
  var.genes <- as.character(sort(unique(unlist(var.genes))))
  obj@var.genes <- var.genes
  obj <- calcPCA(obj, verbose = FALSE, do.print = FALSE)
  plot(plotDimArray(
    obj,
    label = "group.labels",
    reduction.use = "pca",
    dims.to.plot = 1:18,
    discrete.colors = discrete.colors,
    legend = FALSE,
    outer.title = "PCA (group labels)",
    plot.title = ""
  ))

  obj <- calcKNN(obj)
  knn <- round(sqrt(ncol(obj@logupx.data)))
  sigma.use <- "local"
  obj <- calcDM(obj, sigma.use = sigma.use, knn = knn, verbose = FALSE)
  plot(plotDimArray(
    obj,
    label = "group.labels",
    reduction.use = "dm",
    dims.to.plot = 1:18,
    discrete.colors = discrete.colors,
    legend = FALSE,
    outer.title = "Diffusion map (group labels)",
    plot.title = ""
  ))
  
  obj <- calcTsne(
    object = obj, 
    perplexity = max(30,ncol(obj@logupx.data)/10),
    which.dims = c(1:18), 
    theta = 0.5,
    max_iter = 3000
  )

    plot(plotDim(
    obj,
    label = "group.labels",
    reduction.use = "tsne",
    discrete.colors = discrete.colors,
    legend = FALSE,
    plot.title = "tSNE (group labels)"
  ))

  root.cells <- rownames(obj@group.ids)[obj@group.ids$group.labels == root]
  flood.result <- floodPseudotime(obj, root.cells = root.cells)
  obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")
  
  
  plot(plotDim(
    obj,
    label = "pseudotime",
    reduction.use = "tsne",
    plot.title = "tSNE (pseudotime)"
  ))

  tip.clusters <- obj@group.ids$pat.labels[obj@group.ids$group.labels %in% tips]
  names(tip.clusters) <- rownames(obj@group.ids)[obj@group.ids$group.labels %in% tips]
  tip.cluster.number <- seq_along(unique(tip.clusters))
  names(tip.cluster.number) <- unique(tip.clusters)
  tip.cluster.number <- tip.cluster.number[match(tip.clusters, names(tip.cluster.number))]
  obj@group.ids[, "tip.clusters"] <- NA
  obj@group.ids[names(tip.clusters), "tip.clusters"] <- as.character(tip.cluster.number)

  obj.ptlogistic <- pseudotimeDetermineLogistic(
    obj,
    "pseudotime",
    optimal.cells.forward = 40,
    max.cells.back = 80,
    do.plot = FALSE,
    print.values = FALSE
  )

  obj.biased.tm <- as.matrix(
    pseudotimeWeightTransitionMatrix(
      obj,
      "pseudotime",
      logistic.params = obj.ptlogistic
    )
  )

  obj@group.ids <- obj@group.ids[colnames(obj.biased.tm), ]

  obj.walks <- simulateRandomWalksFromTips(
    obj,
    verbose = F,
    tip.group.id = "tip.clusters",
    root.cells = root.cells,
    transition.matrix = obj.biased.tm
  )

  obj <- processRandomWalksFromTips(obj, obj.walks, verbose = F)

  obj.tree <- loadTipCells(obj, "tip.clusters")

  obj.tree <- buildTree(obj.tree,
    pseudotime = "pseudotime",
    tips.use = unique(na.omit(obj@group.ids$"tip.clusters")),
    divergence.method = "preference",
    verbose = FALSE
  )

  plot(plotTree(obj.tree, "group.labels",
    title = "Branching tree (group labels)", 
    cell.size = 1, label.segments = T,
    discrete.colors = discrete.colors,
    cell.alpha = 1,
    legend = FALSE
  ))

  plot(plotTree(obj.tree, "pseudotime",
                title = "Branching tree (pseudotime)", 
                cell.size = 1, label.segments = T,
                cell.alpha = 1,
                legend = FALSE
  ))

  plot(plotTree(obj.tree, "segments",
                title = "Branching tree (tree segments)", 
                cell.size = 1, label.segments = T,
                cell.alpha = 1,
                legend = FALSE
  ))
  

  plot(plotDim(
    obj.tree,
    label = "segments",
    reduction.use = "tsne",
    plot.title = "tSNE (tree segments)"
  ))
  
    
  visitation <- data.frame(
    cell = rownames(obj.tree@diff.data),
    seg = obj.tree@diff.data$segment,
    stringsAsFactors = F,
    row.names = rownames(obj.tree@diff.data)
  )

  visitation$visit <- log10(
    apply(
      visitation, 1, function(cr) {
        obj.tree@diff.data[as.character(cr["cell"]), paste0("visitfreq.raw.", as.character(cr["seg"]))]
      }
    )
    + 1
  )

  robustly.visited.cells <- visitation[visitation$visit >= 0.5, "cell"]
  final.tips <- segTerminal(obj.tree)
  obj.tree <- treeForceDirectedLayout(
    obj.tree,
    method = "fr",
    cells.to.do = robustly.visited.cells,
    tips = final.tips,
    verbose = FALSE
  )

  plot(plotTreeForce2D(
    obj.tree,
    "group.labels",
    colors = discrete.colors,
    title = "Force directed layout (group labels)"))
  dev.off()
}

