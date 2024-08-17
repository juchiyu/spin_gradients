#' Title k-fold cross validation for PLSC
#'
#' @param data1 Data matrix 1 that entered tepPLS
#' @param data2 Data matrix 2 that entered tepPLS
#' @param pls.res results from tepPLS
#' @param center1 (Default = TRUE) center or not for Data1 in tepPLS
#' @param center2 (Default = TRUE) center or not for Data2 in tepPLS
#' @param scale1 (Default = "SS1") how to scale Data1 in tepPLS
#' @param scale2 (Default = "SS1") how to scale Data2 in tepPLS
#' @param k (Default = 10) number of folds
#' @param DESIGN (Default = NULL) the vector that describes the design 
#' that should be kept when creating the folds
#'
#' @return
#' @export
#'
#' @examples
PLS.kFoldCV <- function(data1, data2, pls.res, 
                        center1 = TRUE, center2 = TRUE, 
                        scale1 = "SS1", scale2 = "SS1", 
                        k = 10, DESIGN = NULL){
  
  ## Check data dimensions
  if (nrow(data1) != nrow(data2)) {
    stop("The number of rows of your two data tables do not match!")
  }
  ## Check pls.res
  if (inherits(pls.res, "tepPLS")) 
    stop("This function only takes pls.res from tepPLS output of the `TExPosition` package!")
  
  ## get index for participants
  ID <- seq(nrow(data1))
  n.ID <- length(ID)
  
  ## Create k folds
  if (is.null(DESIGN)){ ## if resample with no structure
    ID.order <- sample(ID, size = length(ID), replace = FALSE)
    fold.idx <- setNames(
      cut(seq(n.ID), unique(quantile(seq(n.ID), probs = seq(0, 1, length = k+1))), 
          include.lowest = TRUE, 
          labels = paste0("Fold", seq(k))), 
      ID.order)
    kfold.design <- fold.idx[order(ID.order)]
  }else{ ## resample with design
    grp.ID <- split(ID, DESIGN)
    grp.n.ID <- lapply(grp.ID, length)
    grp.ID.order <- lapply(grp.ID, function(x) sample(x, size = length(x), replace = FALSE))
    grp.fold.idx <- lapply(grp.n.ID, function(x){
      cut(seq(x), unique(quantile(seq(x), probs = seq(0, 1, length = k+1))), 
          include.lowest = TRUE, 
          labels = paste0("Fold", seq(k)))
    })
    grp.kfold.design <- purrr::map2(grp.fold.idx, grp.ID.order, setNames)
    fold.idx <- setNames(unlist(grp.kfold.design, use.names = FALSE), unlist(grp.ID.order, use.names = FALSE))
    kfold.design <- fold.idx[order(unlist(grp.ID.order, use.names = FALSE))]
  }
  
  ## Extracts folds
  fold.lev <- levels(kfold.design)
  
  ## create empty matrices for 10 folds
  pls.kfcv <- list(lx.hat = matrix(NA, 
                                   nrow = nrow(pls.res$TExPosition.Data$lx),
                                   ncol = ncol(pls.res$TExPosition.Data$lx),
                                   dimnames = dimnames(pls.res$TExPosition.Data$lx)),
                   ly.hat = matrix(NA, 
                                   nrow = nrow(pls.res$TExPosition.Data$ly),
                                   ncol = ncol(pls.res$TExPosition.Data$ly),
                                   dimnames = dimnames(pls.res$TExPosition.Data$ly)),
                   lxly.cor = list(),
                   p = list(),
                   q = list(),
                   ci = list(),
                   cj = list(),
                   Dv = list(),
                   eig = list())
  
  ## run PLS
  for (i in 1:k){
    ## select target fold
    tg.lev <- fold.lev[i]
    ## get train and test sets
    train.data1 <- data1[kfold.design != tg.lev, ]
    test.data1 <- data1[kfold.design == tg.lev, ]
    train.data2 <- data2[kfold.design != tg.lev, ]
    test.data2 <- data2[kfold.design == tg.lev, ]
    ## run PLS with train set
    train.pls <- tepPLS(train.data1, train.data2,
           center1 = center1, center2 = center2,
           scale1 = scale1, scale2 = scale2, graphs = FALSE)
    
    ## Flip the sign if the correlation to the original loadings are negative
    flip <- diag(cor(pls.res$TExPosition.Data$pdq$q, train.pls$TExPosition.Data$pdq$q)) < 0
    train.pls$TExPosition.Data$pdq$p[,flip] <- train.pls$TExPosition.Data$pdq$p[,flip]*-1
    train.pls$TExPosition.Data$pdq$q[,flip] <- train.pls$TExPosition.Data$pdq$q[,flip]*-1
    
    ## project test sets
    test.pls <- ProjectSupplementaryData4PLS(train.pls, test.data1, test.data2)
    
    ## save results from the train set
    pls.kfcv$p[[tg.lev]] <- train.pls$TExPosition.Data$pdq$p
    pls.kfcv$q[[tg.lev]] <- train.pls$TExPosition.Data$pdq$q
    pls.kfcv$ci[[tg.lev]] <- train.pls$TExPosition.Data$pdq$p^2
    pls.kfcv$cj[[tg.lev]] <- train.pls$TExPosition.Data$pdq$q^2
    pls.kfcv$Dv[[tg.lev]] <- train.pls$TExPosition.Data$pdq$Dv
    pls.kfcv$eig[[tg.lev]] <- train.pls$TExPosition.Data$pdq$eigs
    
    ## save results from the test set
    pls.kfcv$lx.hat[rownames(test.data1),] <- test.pls$lx
    pls.kfcv$ly.hat[rownames(test.data2),] <- test.pls$ly
  }
  
  pls.kfcv$lxly.cor[["cor.lxhatlyhat"]] <- diag(cor(pls.kfcv$lx.hat, pls.kfcv$ly.hat))
  pls.kfcv$lxly.cor[["cor.lxlxhat"]] <- diag(cor(pls.res$TExPosition.Data$lx, pls.kfcv$lx.hat))
  pls.kfcv$lxly.cor[["cor.lylyhat"]] <- diag(cor(pls.res$TExPosition.Data$ly, pls.kfcv$ly.hat))
  
  return(list(resample.grpvec = kfold.design,
              cross.validation.res = pls.kfcv))
  
}
