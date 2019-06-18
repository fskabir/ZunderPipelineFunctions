#' Read metadata
#'
#' Reads metadata file and produces format friendly to scripts in Zunder lab pipeline. Metadata
#' strategy inspired by Nowicka et al., 2017, F1000 Research
#' @param input.folder directory containing metadata file
#' @param md.filename metadata filename, defaults to metadata.csv as outputted by generation script
#' @export
read.metadata <- function(input.folder,md.filename = "metadata.csv"){
  md <- read.csv(paste0(input.folder,md.filename))
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
