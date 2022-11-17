# run permutation for second eigenvalue

lx1 <- as.matrix(pls.res$TExPosition.Data$lx[,1])
ly1 <- pls.res$TExPosition.Data$ly[,1]

X.proc <- expo.scale(behav.reg, center = TRUE, scale = "SS1")
Y.proc <- expo.scale(grad.reg, center = TRUE, scale = "SS1")

# # subtract from X and Y
# Xhat <- (diag(length(lx1)) - lx1 %*% solve(crossprod(lx1)) %*% t(lx1)) %*% X.proc
# Yhat <- (diag(length(ly1)) - ly1 %*% solve(crossprod(ly1)) %*% t(ly1)) %*% Y.proc

Xhat <- lm(X.proc~lx1)$residuals
Yhat <- lm(Y.proc~ly1)$residuals

pls.inf2 <- data4PCCAR::perm4PLSC(Xhat, Yhat, 
                                 center1 = FALSE, scale1 = FALSE, 
                                 center2 = FALSE, scale2 = FALSE, 
                                 nIter = 1000)

pls.inf2$pOmnibus

PlotScree(pls.inf2$fixedEigenvalues, pls.inf2$pEigenvalues)

# # I think this is incorrect
# pls.inf2 <- data4PCCAR::perm4PLSC(Xhat, Yhat, 
#                                   center1 = TRUE, scale1 = "SS1", 
#                                   center2 = TRUE, scale2 = "SS1", 
#                                   nIter = 1000)
