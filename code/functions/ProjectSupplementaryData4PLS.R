#' Title Project supplementary rows onto PLSC dimension
#'
#' @param pls.res results from tepPLS
#' @param new.data1 supplementary rows for Data1
#' @param new.data2 supplementary rows for Data2
#'
#' @return lx and ly
#' @export
#'
#' @examples
ProjectSupplementaryData4PLS <- function(pls.res, new.data1, new.data2){
  
  ## Check data dimensions
  if (ncol(new.data1) != nrow(pls.res$TExPosition.Data$fi)) {
    stop("The number of columns of your new.data1 do not match the number of variables in the original data1!")
  }
  
  if (ncol(new.data2) != nrow(pls.res$TExPosition.Data$fj)) {
    stop("The number of columns of your new.data2 do not match the number of variables in the original data2!")
  }
  
  ## Check pls.res
  if (inherits(pls.res, "tepPLS")) 
    stop("This function only takes pls.res from tepPLS output of the `TExPosition` package!")
  
  ## get SVD results
  res <- pls.res$TExPosition.Data
  
  ## preprocessed new data tables
  data1.processed <- scale(new.data1, 
                           center = res$data1.norm$center,
                           scale = res$data1.norm$scale)
  data2.processed <- scale(new.data2, 
                           center = res$data2.norm$center,
                           scale = res$data2.norm$scale)
  
  ## projection
  lx <- data1.processed %*% res$pdq$p
  ly <- data2.processed %*% res$pdq$q
  
  ## return
  return(list(lx = lx, ly = ly))
}