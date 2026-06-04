#' This function will return the MFI of a specified parameter in a sample from a specified gate.
#' 
#' @param gs a gating set object from which to retrieve 
#' @param x the index number of a file within gs from which to retrieve the MFI
#' @param subset the name of a subpopulation that exists in gs for which to retrieve the MFI, default is "root"
#' @param parameter the measurement parameter from which the median should be calculated
#' @param inverse.transform TRUE or FALSE whether to revert the transformed data back to the raw values, default is TRUE
#' 
retrieveMFI <- function(x, gs, subset="root", parameter, inverse.transform=TRUE){
    MFIData <- flowWorkspace::gs_pop_get_data(gs[[x]], y=subset, inverse.transform=inverse.transform)
    Exprs <- MFIData[[1]] |> flowCore::exprs() |> as.data.frame()
    column <- Exprs[[parameter]] 
    MFI <- median(column)
    return(MFI)
    }