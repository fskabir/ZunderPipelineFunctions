#' Pull out markers for transformation
#'
#' Finds markers in panel file marked for either clustering or plotting
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.markers <- function(panel){
  transform.markers <- panel$Metal[panel$Clustering == 1 | panel$Plotting == 1]
  return(transform.markers)
}

#' Get marker names in more legibile format for Concat_Transformed file
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.annotate <- function(panel){
  transform.markers.annotate <- paste0(panel$Antigen[panel$Clustering == 1 | panel$Plotting == 1],"_",panel$Metal[panel$Clustering == 1 | panel$Plotting == 1])
  return(transform.markers.annotate)
}

#' Get clustering variables in legible format
#' #' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.clustering.annotate <- function(panel){
  clustering.markers.annotate <- paste0(panel$Antigen[panel$Clustering == 1],"_",panel$Metal[panel$Clustering == 1])
  return(clustering.markers.annotate)
}
