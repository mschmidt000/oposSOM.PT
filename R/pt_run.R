#' oposSOM PT run Function
#'
#' This function allows you to run the oposSOM pseudotemporal analysis.
#' @param env An enviroment produced by the oposSOM pipeline.
#' @param root A single group label defined as the developmental source.
#' @param tips A vector of group labels defined as the developmental sinks.
#' @param lineages A list of vectors of group labels sorted by pseudotime.
#' @export
pt_run <- function(env, root = c(), tips = c(), lineages = list(), loom_path = c()) {
  if (length(lineages) == 0 & length(tips) == 0) {
    "Nether tips nor lineages specified."
  } else if (!length(lineages) == 0 & !length(tips) == 0) {
    "Please specify either tips or lineages."
  } else if (!length(lineages) == 0 & length(tips) == 0) {
    
    dir.create(paste(env$files.name, "- Results"), showWarnings = FALSE)
    dir.create(
      file.path(
        paste(env$files.name, "- Results"), "Pseudotemporal Analysis"
      ),
      showWarnings = FALSE
    )
    dir.create(
      file.path(
        paste(env$files.name, "- Results"), "Pseudotemporal Analysis", "Cell State Trajectories"
      ), 
      showWarnings = FALSE
    )
    dir.create(
      file.path(
        paste(env$files.name, "- Results"), "Pseudotemporal Analysis", "Gene State Trajectories"
      ), 
    showWarnings = FALSE
    )

    som_gradients(env, lineages = lineages)
    som_lcp(env, lineages = lineages)
    som_mep(env, lineages = lineages)
    # som_arrows(env, lineages = lineages)
    
  } else if (length(lineages) == 0 & !length(tips) == 0) {
    
    dir.create(paste(env$files.name, "- Results"), showWarnings = FALSE)
    dir.create(
      file.path(
        paste(env$files.name, "- Results"), "Pseudotemporal Analysis"
      ),
      showWarnings = FALSE
    )
    dir.create(
      file.path(
        paste(env$files.name, "- Results"), "Pseudotemporal Analysis", "Cell State Trajectories"
      ), 
      showWarnings = FALSE
    )
    dir.create(
      file.path(
        paste(env$files.name, "- Results"), "Pseudotemporal Analysis", "Gene State Trajectories"
      ), 
      showWarnings = FALSE
    )
    
    som_gradients(env, root, tips)
    som_lcp(env, root, tips)
    som_mep(env, root, tips)
    # som_arrows(env, root, tips)
    branched_pt(env, root, tips)
  }
  
  # if (!length(loom_path) == 0){
  #   som_velocity(env, loom_path = loom_path)
  # }
}
