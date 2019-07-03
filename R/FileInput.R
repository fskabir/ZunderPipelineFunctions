#' Read metadata
#'
#' Reads metadata file and produces format friendly to scripts in Zunder lab pipeline. Metadata
#' strategy inspired by Nowicka et al., 2017, F1000 Research
#' @param input.folder directory containing metadata file
#' @param md.filename metadata filename, defaults to metadata.csv as outputted by generation script
#' @export
read.metadata <- function(input.folder,md.filename = "metadata.csv"){
  md <- read.csv(paste0(input.folder,md.filename))
  #make md into df
  md <- as.data.frame(md)
  #make filenames character vectors
  md$file_name <- as.character(md$file_name)
  #testing: print file_name data type
  print(typeof(md$file_name))
  return(md)
}

#' Read panel
#'
#' Reads panel file and produces format friendly to scripts in Zunder lab pipeline
#' @param input.folder directory containing metadata file
#' @param panel.filename panel filename, defaults to metadata.csv as outputted by generation script
#' @export
read.panel <- function(input.folder,panel.filename = "panel.csv"){
  panel <- read.csv(paste0(input.folder,panel.filename))
  panel$Antigen <- gsub("-", "_", panel$Antigen)
  return(panel)
}

#' Pull out markers for transformation
#'
#' Finds markers in panel file marked for either clustering or plotting
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.markers <- function(panel){
  transform.markers <- panel$Metal[panel$Clustering == 1 || panel$Plotting == 1]
  return(transform.markers)
}

#' Get marker names in more legibile format for Concat_Transformed file
#' @param panel panel formatted following generation by getMetadata.R and input by read.panel
#' function
#' @export
get.transform.annotate <- function(panel){
  transform.markers.annotate <- paste0(panel$Antigen[panel$Clustering == 1 || panel$Plotting == 1],"_",panel$Metal[panel$Clustering == 1 || panel$Plotting == 1])
  return(transform.markers.annotate)
}
