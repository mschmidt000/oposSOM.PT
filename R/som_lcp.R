#' SOM Least Cost Path (LCP) function
#'
#' This function calculate the difference portraits between the root and the tips
#' to generate gradient trajectories along topocrafic least cost paths
#' @param env An enviroment produced by the oposSOM pipeline.
#' @param root A single group label defined as the developmental source.
#' @param tips A vector of group labels defined as the developmental sinks.
#' @param lineages A list of vectors of group labels sorted by pseudotime.
#' @keywords topografic distance LCP
#' @import raster topoDistance oposSOM
som_lcp <- function(env, root = c(), tips = c(), lineages = c()) {
  
  filename <- file.path(
    paste(env$files.name, "- Results"),
    "Pseudotemporal Analysis",
    "Gene State Trajectories",
    "Topological Distance Portraits.pdf"
  )
  util.info("Writing:", filename)
  pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
  
  spot.list <- env$spot.list.kmeans
  group.metadata <- do.call(
    cbind, 
    by(t(env$metadata), env$group.labels, colMeans)
  )[, unique(env$group.labels)]
  group.spotdata <- do.call(
    cbind, 
    by(t(spot.list$spotdata), env$group.labels, colMeans)
  )[, unique(env$group.labels)]
  
  if (!length(tips) == 0) {
    
    start.pos <- spot.list$spots[[names(which.max(group.spotdata[, root]))]]$position
    
    for (i in seq_along(tips)) {
      
      tip <- tips[i]
      stop.pos <- spot.list$spots[[names(which.max(group.spotdata[, tip]))]]$position
      diff.data <- group.metadata[, root] - group.metadata[, tip]
      
      m <- matrix(
        diff.data,
        env$preferences$dim.1stLvlSom,
        env$preferences$dim.1stLvlSom
      )
      m <- m - min(m)
      
      r <- raster(
        t(m[, ncol(m):1]),
        xmn = 0.5,
        xmx = nrow(m) + .5,
        ymn = 0.5,
        ymx = ncol(m) + 0.5
      )
      projection(r) <- CRS("+init=epsg:27700")
      
      xy <- matrix(
        ncol = 2, byrow = TRUE,
        c(
          start.pos,
          stop.pos
        )
      )
      colnames(xy) <- c("longitude", "latitude")

      tLCP <- topoLCP(r, costSurface = r, pts = xy, paths = TRUE)
      topoPathMap(
        r,
        xy,
        topoPaths = tLCP, type = "hillshade", pathWidth = 4, cex = 2, bg = "blue",
        main = paste0("Least Cost Path", " (", tip, ")")
      )
      tryCatch(
        suppressWarnings({
        topoProfile(r, topoPaths = tLCP)
      }), error = function(w) {
        print(paste(tip, "is not a tip"))
      })
    }
    
    dev.off()
    
  } else if (!length(lineages) == 0) {
    
    for (i in seq_along(lineages)) {
      
      root <- head(lineages[[i]], 1)
      tip <- tail(lineages[[i]], 1)
      start.pos <- spot.list$spots[[names(which.max(group.spotdata[, root]))]]$position
      stop.pos <- spot.list$spots[[names(which.max(group.spotdata[, tip]))]]$position
      diff.data <- group.metadata[, root] - group.metadata[, tip]
      
      m <- matrix(
        diff.data,
        env$preferences$dim.1stLvlSom,
        env$preferences$dim.1stLvlSom
      )
      m <- m - min(m)
      
      r <- raster(
        t(m[, ncol(m):1]),
        xmn = 0.5,
        xmx = nrow(m) + .5,
        ymn = 0.5,
        ymx = ncol(m) + 0.5
      )
      projection(r) <- CRS("+init=epsg:27700")

      xy <- matrix(
        ncol = 2, byrow = TRUE,
        c(
          start.pos,
          stop.pos
        )
      )
      colnames(xy) <- c("longitude", "latitude")

      tLCP <- topoLCP(r, costSurface = r, pts = xy, paths = TRUE)
      topoPathMap(r, xy,
        topoPaths = tLCP, type = "hillshade", pathWidth = 4, cex = 2, bg = "blue",
        main = paste0("Least Cost Path", " (", names(lineages)[i], ")")
      )
      tryCatch(suppressWarnings({
        topoProfile(r, topoPaths = tLCP)
      }), error = function(w) {
        print(paste(tip, "is not a tip"))
      })
    }
    dev.off()
  }
}  
# nrz <- nrow(m)
# zfacet <- m[-1, -1] + m[-1, -nrz] + m[-nrz, -1] + m[-nrz, -nrz]
# facetcol <- cut(zfacet, 1000)
# coords <- tLCP[[2]]@lines[[1]]@Lines[[1]]@coords
# # for(a in 1:nrow(coords)){
# #  m[coords[a,1], coords[a,2]] <- 1
# # }
# mycol <- env$color.palette.portraits(1000)[facetcol]
# mycol[apply(coords, 1, function(x) {
#   which(env$som.result$node$x == x[1] & env$som.result$node$y == x[2])
# }) + 61] <- "black"
# 
# persp(r,
#   axes = FALSE, border = NA, expand = 0.5, col = mycol,
#   phi = 45, theta = -5, xlab = "", ylab = "", zlab = "", box = TRUE
# )
