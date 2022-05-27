## Plot the rectangular matrix into a symmetrical square matrix
## Note: Map factor scores to a square matrix that has the same rows and columns as the correlation matrix that is analyzed 
## Created by Ju-Chi Yu, 02/27/2019
##====================================================
## Function name: vec2sqmat
## Input:
##      (1) data2map: n x 1 vector; 
##                a vector with matching lengths of n, and n is the number of intersections in the upper triangle of the 
##                square matrix (diagonal included)
## Output:
##      An m x m square matrix with the matching length of dimensions
## Note:
##      m is found by solving the quadratic equation: m^2 - m - 2*n = 0, because m = sqrt(2*n - m)
##      a = 1, b = -1, c = -2*n, m = (-b + sqrt(b^2 - 4ac))/2a, only the positive solution is considered
#=====================================================
vec2sqmat <- function(data2map){
  # find the length of the input vector
  nVec <- length(data2map)
  # find the matching dimension
  nDim <- (1 + sqrt(1+8*nVec))/2
  # build empty square matrix
  res_mat <- matrix(0, nDim, nDim)
  # map the vector back to the upper triangle
  res_mat[upper.tri(res_mat)] <- data2map
  # copy and paste the upper triangle to the lower one
  res_mat[lower.tri(res_mat)] <- t(res_mat)[lower.tri(res_mat)]
  return(sqmat = res_mat)
}