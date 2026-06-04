#' Create a character vector by reading the $BeadsType keyword from a list of FCS files in a cytoSet object
#' @param data A cytoSet object
#' @param x index of the file within the cytoSet from which to retrieve the $DATE keyword
getBeads <- function(x, data){
    beadsType <- flowCore::description(data[[x]])$BeadsType
    return(beadsType)
}