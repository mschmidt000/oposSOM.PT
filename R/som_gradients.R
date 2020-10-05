#' SOM gradients function
#'
#' This function provides the vector field (gradient field) and the corresponding stream lines in expression portraits.
#' @param env An enviroment produced by the oposSOM pipeline.
#' @param root A single group label defined as the developmental source.
#' @param tips A vector of group labels defined as the developmental sinks.
#' @param lineages A list of vectors of group labels sorted by pseudotime.
#' @import rasterVis raster RColorBrewer oposSOM
som_gradients <- function(env, root = c(), tips = c(), lineages = c()) {
  if (!length(tips) == 0) {

    filename <- file.path(
      paste(env$files.name, "- Results"),
      "Pseudotemporal Analysis",
      "Gene State Trajectories",
      "Vector and Stream Portraits.pdf"
    )
    util.info("Writing:", filename)
    pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)
    

    par(mfrow = c(5, 6))
    par(mar = c(0.5, 2.5, 4.5, 1.5))

    group.metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans))[, unique(env$group.labels)]
    stack.list <- list()
    stack.names <- c()
    for (i in seq_along(tips)) {
      diff.data <- group.metadata[, root] - group.metadata[, tips[i]]
      m <- matrix(diff.data, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
      r <- raster(t(m[, ncol(m):1]), xmn = 0.5, xmx = nrow(m) + .5, ymn = 0.5, ymx = ncol(m) + 0.5)
      projection(r) <- CRS("+init=epsg:27700")
      stack.list[[i]] <- r
      stack.names[i] <- paste0(root, " - ", tips[i])
    }
    names(stack.list) <- stack.names
    stacks <- stack(stack.list)
    for (a in seq_along(tips)) {
      plot(vectorplot(stack.list[[a]],
        scaleSlope = TRUE, aspX = 0.5,
        # narrows = 500,
        margin = list(FUN = "median"),
        # lwd.arrows = 0.6,
        contour = TRUE,
        par.settings = rasterTheme(
          region = rev(brewer.pal(9, "RdBu"))
        ),
        main = paste0("Vector field (", stack.names[a], ")")
      ))

      ### Stream Maps ###
      plot(streamplot(stack.list[[a]],
        contour = TRUE,
        region = TRUE,
        par.settings = streamTheme(
          region = rev(brewer.pal(9, "RdBu")),
          symbol = rev(brewer.pal(5, "Greys")),
        ),
        parallel = TRUE,
        mc.cores = 40,
        droplet = list(pc = 3), streamlet = list(L = 20),
        main = paste0("Streams (", stack.names[a], ")")
      ))
    }

    dev.off()
  } else if (!length(lineages) == 0) {
    for (i in seq_along(lineages)) {
      filename <- file.path(paste(env$files.name, "- Results"), "Pseudotemporal Analysis", "Gene State Trajectories", paste0("Vector and Stream Portraits ", names(lineages)[i], ".pdf"))
      util.info("Writing:", filename)
      pdf(filename, 29.7 / 2.54, 21 / 2.54, useDingbats = FALSE)

      par(mfrow = c(5, 6))
      par(mar = c(0.5, 2.5, 4.5, 1.5))

      group.metadata <- do.call(cbind, by(t(env$metadata), env$group.labels, colMeans))[, unique(env$group.labels)]
      stack.list <- list()
      stack.names <- c()
      for (a in 1:(length(lineages[[i]]) - 1)) {
        diff.data <- group.metadata[, lineages[[i]][a]] - group.metadata[, lineages[[i]][a + 1]]
        m <- matrix(diff.data, env$preferences$dim.1stLvlSom, env$preferences$dim.1stLvlSom)
        r <- raster(t(m[, ncol(m):1]), xmn = 0.5, xmx = nrow(m) + .5, ymn = 0.5, ymx = ncol(m) + 0.5)
        projection(r) <- CRS("+init=epsg:27700")
        stack.list[[a]] <- r
        stack.names[a] <- paste0(lineages[[i]][a], " - ", lineages[[i]][a + 1])
      }
      names(stack.list) <- stack.names
      stacks <- stack(stack.list)
      for (a in 1:(length(lineages[[i]]) - 1)) {
        plot(vectorplot(stack.list[[a]],
          scaleSlope = TRUE, aspX = 0.5,
          # narrows = 500,
          margin = list(FUN = "median"),
          # lwd.arrows = 0.6,
          contour = TRUE,
          par.settings = rasterTheme(
            region = rev(brewer.pal(9, "RdBu"))
          ),
          main = paste0("Vector field (", stack.names[a], ")")
        ))

        ### Stream Maps ###
        plot(streamplot(stack.list[[a]],
          contour = TRUE,
          region = TRUE,
          par.settings = streamTheme(
            region = brewer.pal(9, "RdBu"),
            symbol = rev(brewer.pal(5, "Greys")),
          ),
          parallel = TRUE,
          mc.cores = 40,
          droplet = list(pc = 3), streamlet = list(L = 20),
          main = paste0("Streams (", stack.names[a], ")")
        ))
      }
      dev.off()
    }
  }
}
# ### Vector Maps ###
# plot(vectorplot(stacks, scaleSlope = TRUE,aspX=0.5,
#                 #narrows = 500,
#                 #margin = TRUE,
#                 #lwd.arrows = 0.6,
#                 contour=TRUE,
#                 par.settings = RdBuTheme(),
#                 main = "Vector Portraits"))
#
#
#
# ### Stream Maps ###
# plot(streamplot(stacks, contour=TRUE,
#                 region = TRUE,
#                 par.settings = streamTheme(
#                   region =  brewer.pal(9, 'RdBu'),
#                   symbol =  rev(brewer.pal(5, 'Greys')),
#                 ),
#                 parallel = TRUE,
#                 mc.cores = 40,
#                 droplet = list(pc = 3), streamlet = list(L = 10),
#                 main = "Stream Portraits"))
