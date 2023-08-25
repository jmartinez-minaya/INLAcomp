#' Arabidopsis thaliana in the Iberian Peninsula
#'
#' Collection of 301 accessions of the annual plant \textit{Arabidopsis thaliana} on the Iberian Peninsula. For each
#' accession, the probability of belonging to each of the 4 genetic clusters (GC) inferred in \cite{martinez-minaya2019}, namely, GC1, GC2, GC3 and GC4 represents the Compositional Data.
#'
#' @format ## `arabidopsis`
#' A data frame with 301 rows and 15 columns:
#' \describe{
#'   \item{acronym}{Country name}
#'   \item{x, y}{Coordinates}
#'   \item{gc1, gc2, gc3, gc4}{Probability for belonging to each of the four genetic clusters: GC1, GC2, GC3 and GC4.}
#'   \item{bio1, bio2, bio3, bio4, bio8, bio12, bio15, bio18}{Bioclimatic variables corresponding with the definition in Worldclim}
#'   ...
#' }
#' @source <https://zenodo.org/record/2552025>
"arabidopsis"


#' Spatial Polygons of the Iberian Peninsula
#'
#' Spatial polygons for the Iberian Peninsula
#'
#' @format ## `polygon_IP`
#' SpatialPolygonsDataFrame
#' \describe{
#'   \item{GID_0}{Country name}
#'   \item{NAME_0}{Country name}
#'   ...
#' }
#' @source GADM
"polygon_IP"

#' Spatial Polygons of the Iberian Peninsula low resolution
#'
#' Spatial polygons for the Iberian Peninsula low resolution
#'
#' @format ## `polygon_IP_low`
#' SpatialPolygonsDataFrame
#' \describe{
#'   \item{GID_0}{Country name}
#'   \item{NAME_0}{Country name}
#'   ...
#' }
#' @source GADM
"polygon_IP_low"


#' Raster for predicting spatial term and future predictions
#'
#' Raster for predicting spatial term and future predictions
#'
#' @format ## `raster_predict`
#' raster_predict
#' \describe{
#'   ...
#' }
#' @source GADM
"raster_predict"
