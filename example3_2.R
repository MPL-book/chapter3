# R program for profile MLE

mpleltrc <- function(Y, X, bin_num, beta0, 
                     maxiter=2000, tol=1e-5) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  if(missing(beta0)) {beta0 <- matrix(1,p,1)}
  else {beta0 <- matrix(beta0,p,1)}
  
  if(missing(bin_num)) {bin_num <- floor(sqrt(sum(Y[,3])))}
  
  num_boundary <- bin_num+1
  boundary <- unname(quantile(unique(Y[Y[,3]==1,2]), 
                              probs = seq(0, 100, length.out=num_boundary)/100))
  boundary[1] <- min(Y[,2]) 
  boundary[num_boundary] <- max(Y[,2])
  upper_b <- boundary[2:num_boundary]
  lower_b <- boundary[1:(num_boundary-1)]
  binwid <- upper_b-lower_b
  
  h <- hist(Y[Y[,3]==1,2], breaks = boundary, right = F, 
            plot = F)
  binCount <- matrix(h$counts,bin_num,1)
  bin_ind <- t(sapply(1:n, function(i){Y[i,1]<=upper_b & Y[i,2]>lower_b}))
  
  
  Xbeta <- X%*%beta0
  A_mat <- bin_ind*matrix(rep(exp(Xbeta),bin_num),n,bin_num)
  sumA_mat <- apply(A_mat,2,sum)
  sumA_mat[sumA_mat==0] <- 1e-10
  weights <- A_mat/matrix(rep(sumA_mat,n),n,bin_num, 
                          byrow = T)
  allpl <- NULL
  
  betaold <- beta0
  for(k in 1:maxiter) {
    proflik <- Y[,3]%*%Xbeta-log(sumA_mat)%*%binCount
    allpl <- c(allpl, proflik)
    gradient <- t(X)%*%(Y[,3]-weights%*%binCount)
    Xbar_mat <- t(X)%*%weights
    xtwx <- matrix(0,p,p)
    for(u in 1:bin_num) {
      xtmp <- t(X)%*%(matrix(rep(weights[,u],p),n,p)*X)
      xtwx <- xtwx+binCount[u,]*xtmp
    }
    profhess <- xtwx-Xbar_mat%*% diag(as.numeric(binCount))%*%t(Xbar_mat)
    
    betanew <- betaold + solve(profhess)%*%gradient
    
    Xbeta <- X%*%betanew
    A_mat <- bin_ind*matrix(rep(exp(Xbeta),bin_num),n,bin_num)
    sumA_mat <- apply(A_mat,2,sum)
    sumA_mat[sumA_mat==0] <- 1e-10
    weights <- A_mat/matrix(rep(sumA_mat,n),n,bin_num, byrow = T)
    if(all(abs(betaold-betanew)<tol)) {break}
    else{betaold <- betanew}
    
  }
  theta <- as.numeric(binCount)/binwid*sumA_mat
  retheta <- theta/binwid
  
  return(list(betanew, retheta, gradient, allpl, 
              sqrt(diag(solve(profhess)))))
  
}
