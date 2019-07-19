#' Metadata file input
#'
#' Reads metadata file and produces format friendly to scripts in Zunder lab pipeline. Metadata
#' strategy inspired by Nowicka et al., 2017, F1000 Research
#' @param input.folder directory containing metadata file
#' @param md.filename metadata filename, defaults to metadata.csv as outputted by generation script
#' @export
read.metadata <- function(input.folder,md.filename = "metadata.csv"){
  md <- read.csv(paste0(input.folder,md.filename))
  #make filenames character vectors
  md$file_name <- as.character(md$file_name)
  return(md)
}

#' Panel file input
#'
#' Reads panel file and produces format friendly to scripts in Zunder lab pipeline
#' @param input.folder directory containing metadata file
#' @param panel.filename panel filename, defaults to panel.csv as outputted by generation script
#' @export
read.panel <- function(input.folder,panel.filename = "panel.csv"){
  panel <- read.csv(paste0(input.folder,panel.filename))
  panel$Antigen <- gsub("-", "_", panel$Antigen)
  return(panel)
}

#' Concat_Transformed file input
#'
#' Reads Concat_Transformed file for use in Zunder lab pipeline
#' @param input.folder directory containing concat_transformed file
#' @param concat.transformed.filename concat_transformed filename, defaults to
#' Concat_Transformed.csv as outputted by generation script
#' @export
read.concat.transformed <- function(input.folder,
                                    concat.transformed.filename = "Concat_Tranformed.csv"){
  concat.transformed <- read.csv(paste0(input.folder,concat.transformed.filename))
  return(concat.transformed)
}

#' Layout file input
#'
#' Reads uMAP or other 2D layout for plotting
#' @param input.folder directory containing layout file
#' @param layout.filename layout filename, defaults to UMAP_layout.csv as outputted by generation
#'  script
#' @export
read.layout <- function(input.folder,layout.filename = "UMAP_layout.csv"){
  layout <- read.csv(paste0(input.folder,layout.filename))
  return(layout)
}

#' Cluster file input
#'
#' Reads clusters from pipeline
#' @param input.folder directory containing layout file
#' @param cluster.filename clustering replicate filename, defaults to clusters.csv as outputted by
#'  generation script
#' @importFrom data.table fread
#' @export
read.clusters <- function(input.folder,clusters.filename = "clusters.csv"){
  #Add 1 since clustering done in Python begins indexing at zero
  clusters <- t(fread(paste0(input.folder,clusters.filename))) + 1
  return(clusters)
}

#' Clustering replicate file input
#'
#' Reads clustering replicates for consensus clustering and stability analysis
#' @param input.folder directory containing layout file
#' @param clustering.rep.filename clustering replicate filename, defaults to
#' ClusterStabilityCheck.csv as outputted by generation script
#' @export
read.clustering.rep <- function(input.folder,cluster.rep.filename = "ClusterStabilityCheck.csv"){
  #Add 1 since clustering done in Python begins indexing at zero
  clustering.rep <- read.csv(paste0(input.folder,cluster.rep.filename),header = FALSE) + 1
  return(clustering.rep)
}
