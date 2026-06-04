#' Create a vector of dates by reading the $DATE keyword from a list of FCS files in a cytoSet object
#' @param data A cytoSet object
#' @param x index of the file within the cytoSet from which to retrieve the $DATE keyword
getDate <- function(x, data){
    originalDate <- flowCore::description(data[[x]])$`$DATE`
    newDate <- as.character(originalDate)
    return(newDate)
}