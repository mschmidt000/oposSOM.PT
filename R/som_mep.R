#' SOM Maximum Expression Path (MEP) Function
#'
#' This function connects root and tips along the SOM clusters providing the maximum cumulative weights
#' and then defines the gene state trajectory ensuring maximum cumulative expression between source and
#' sink
#' @param env An enviroment produced by the oposSOM pipeline.
#' @param root A single group label defined as the developmental source.
#' @param tips A vector of group labels defined as the developmental sinks.
#' @param lineages A list of vectors of group labels sorted by pseudotime.
#' @import raster landscapemetrics igraph oposSOM
som_mep <- function(env, root = c(), tips = c(), lineages = c()) {
  filename <- file.path(
    paste(env$files.name, "- Results"),
    "Pseudotemporal Analysis",
    "Gene State Trajectories",
    "SOM Path Portraits.pdf"
  )
  util.info("Writing:", filename)
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
  par(mfrow = c(2, 3))

  spot.list <- env$spot.list.kmeans
  group.metadata <- do.call(
    cbind,
    by(t(env$metadata), env$group.labels, colMeans)
  )[, unique(env$group.labels)]
  group.spotdata <- do.call(
    cbind,
    by(t(spot.list$spotdata), env$group.labels, colMeans)
  )[, unique(env$group.labels)]

  loglog.group.metadata <- apply(group.metadata, 2, function(x) {
    meta.sign <- sign(x)
    meta <- log10(abs(x))
    meta <- meta - min(meta, na.rm = TRUE)
    return(meta * meta.sign)
  })

  r <- raster(
    matrix(
      spot.list$overview.mask,
      env$preferences$dim.1stLvlSom,
      env$preferences$dim.1stLvlSom
    )
  )
  crs(r) <- CRS("+proj=robin +datum=WGS84")

  adj <- get_adjacencies(r, neighbourhood = 4, what = "unlike", upper = FALSE)
  adj <- adj[[1]]
  adj[adj > 0] <- 1
  rownames(adj) <- colnames(adj) <- names(spot.list$spots)

  adj.graph <- graph.adjacency(adj, mode = "undirected")

  if (!length(tips) == 0) {
    start.spot <- names(which.max(group.spotdata[, root]))
    graph.list <- list()

    for (i in seq_along(tips)) {
      
      tip <- tips[i]
      graph.list[[tip]]$stop.spot <- names(which.max(group.spotdata[, tip]))
      graph.list[[tip]]$labels <- tip

      add.data <- rep(0, nrow(group.metadata))
      for (e in c(root,graph.list[[tip]]$labels)) {
        act.metagenes <- loglog.group.metadata[, e]
        add.data <- add.data + act.metagenes
      }

      add.mean <- sapply(spot.list$spots, function(x) {
        mean(add.data[x$metagenes])
      })
      stop.spot <- graph.list[[tip]]$stop.spot
      capacities <- sapply(seq_along(E(adj.graph)), function(x) {
        sum(
          add.mean[ends(adj.graph, x)[, 1]],
          add.mean[ends(adj.graph, x)[, 2]]
        )
      })
      names(capacities) <- E(adj.graph)
      for (e in 1:length(capacities)) {
        capacities[i] <- (add.mean[ends(adj.graph, i)[, 1]] + add.mean[ends(adj.graph, i)[, 2]])
      }

      capacities <- capacities + abs(min(capacities)) + 1

      max.flow <- max_flow(adj.graph, start.spot, stop.spot, capacities)
      flow <- which(!max.flow$flow == 0)
      flow.graph <- graph_from_data_frame(as_edgelist(adj.graph)[flow, ], vertices = names(spot.list$spots))

      V(flow.graph)$label <- names(spot.list$spots)
      graph.list[[tip]]$graph <- flow.graph
      graph.list[[tip]]$max.flow <- max.flow
      graph.list[[tip]]$cutoff.flow <- flow
      col <- colorRampPalette(c("blue2", "white", "red2"))(1000)
      image(
        matrix(add.data, env$preferences$dim.1stLvlSom),
        axes = TRUE, cex.main = 1.5, col = col
      )
      par(new = TRUE)
      plot(flow.graph,
        layout = t(sapply(spot.list$spots, function(x) {
          x$position
        })),
        vertex.label.color = "black",
        vertex.label.cex = 1,
        vertex.shape = "none",
        xlim = c(-0.9, 0.9),
        ylim = c(-0.9, 0.9),
        edge.width = abs(max.flow$flow[flow]) * 0.3,
        edge.color = "gray20",
        edge.arrow.mode = 0,
        main = tip
      )
    }
    dev.off()
  } else if (!length(lineages) == 0) {
    for (i in seq_along(lineages)) {
      root <- head(lineages[[i]], 1)
      tip <- tail(lineages[[i]], 1)

      start.spot <- names(which.max(group.spotdata[, root]))
      graph.list <- list()
      graph.list[[tip]]$stop.spot <- names(which.max(group.spotdata[, tip]))
      graph.list[[tip]]$labels <- lineages[[i]]

      add.data <- rep(0, nrow(group.metadata))
      for (e in lineages[[i]]) {
        act.metagenes <- loglog.group.metadata[, e]
        add.data <- add.data + act.metagenes
      }

      add.mean <- sapply(spot.list$spots, function(x) {
        mean(add.data[x$metagenes])
      })
      stop.spot <- graph.list[[tip]]$stop.spot
      capacities <- sapply(seq_along(E(adj.graph)), function(x) {
        sum(
          add.mean[ends(adj.graph, x)[, 1]],
          add.mean[ends(adj.graph, x)[, 2]]
        )
      })
      names(capacities) <- E(adj.graph)
      capacities <- capacities + abs(min(capacities)) + 1

      max.flow <- max_flow(adj.graph, start.spot, stop.spot, capacities)
      flow <- which(!max.flow$flow == 0)
      flow.graph <- graph_from_data_frame(as_edgelist(adj.graph)[flow, ], vertices = names(spot.list$spots))

      V(flow.graph)$label <- names(spot.list$spots)
      graph.list[[tip]]$graph <- flow.graph
      graph.list[[tip]]$max.flow <- max.flow
      graph.list[[tip]]$cutoff.flow <- flow
      col <- colorRampPalette(c("blue2", "white", "red2"))(1000)
      image(
        matrix(add.data, env$preferences$dim.1stLvlSom),
        axes = TRUE, cex.main = 1.5, col = col
      )
      par(new = TRUE)
      plot(
        flow.graph,
        layout = t(sapply(spot.list$spots, function(x) {
          x$position
        })),
        vertex.label.color = "black",
        vertex.label.cex = 1,
        vertex.shape = "none",
        xlim = c(-0.9, 0.9),
        ylim = c(-0.9, 0.9),
        edge.width = abs(max.flow$flow[flow]) * 0.3,
        edge.color = "gray20",
        edge.arrow.mode = 0,
        main = tip
      )
    }
    dev.off()
  }
}
# cutoff = sd(abs(max.flow$flow))
# flow = which(abs(max.flow$flow)>cutoff)
# arr <- max.flow$flow[flow]
# arr[which(max.flow$flow[flow] < 0)] <- 2
# arr[which(max.flow$flow[flow] > 0)] <- 1
