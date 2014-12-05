my_kernel <- rbfdot(sigma=0.1)
my_matrix <- as.matrix(mydat)
my_kernel(my_matrix[1, ], my_matrix[2, ]) # Distance using my_kernel function
sigma=0.1 #Set sigma value
##Check it through Gaussian function
exp(sigma * (2 * crossprod(my_matrix[1, ], my_matrix[2, ]) -
               crossprod(my_matrix[1, ]) - crossprod(my_matrix[2, ])))
class(my_kernel)
getClass('rbfkernel')
K <- kernelMatrix(kernel=my_kernel, x=as.matrix(mydat))
dim(K)
class(K)
C <- 1
my_c <- (-1)*diag(K)
my_H <- (2)*K
my_A <- rep(1, nrow(my_matrix))
my_b <- 1
my_l <- rep(0, nrow(my_matrix))
my_u <- rep(C, nrow(my_matrix))
my_r <- 0
my_solution <- ipop(c=my_c, H=my_H, A=my_A, b=my_b, l=my_l, u=my_u, r=my_r, maxiter=300, margin=1E-6)
my_alphas <- my_solution@primal # use @ symbol to access s4 slot
my_alphas <- round(my_alphas, abs(log10(1E-6)))
sum(my_alphas)
any(my_alphas < 0)
svdd <- function(x, k, cost, tol=1E-6){
  K <- kernelMatrix(kernel=k, x=x)
  n <- nrow(x)
  solution <- ipop(c=(-1)*diag(K),
                   H=2*K,
                   A=rep(1, n),
                   b=1,
                   l=rep(0, n),
                   u=rep(cost, n),
                   r=0,
                   maxiter=300,
                   margin=tol)
  alphas <- solution@primal
  alphas <- round(alphas, abs(log10(tol))) # round based on tol
  svs <- x[alphas > 0, ] # these are the support vectors
  return(new('SVDD',
             data=x,
             cost=as.numeric(cost),
             alphas=alphas,
             support_vectors=svs,
             kernel=k,
             tolerance=tol,
             fun_call=match.call()
  ))
}
# supplying no value to tol implies it should use the default value tol = 1E-6
euclidean_solution <- svdd(x=my_matrix, k=my_kernel, cost=5)
summary(euclidean_solution)
setMethod(
  f='plot',
  signature='SVDD',
  definition=function(x, y){
    if(class(x@kernel)=='rbfkernel'){
      rotated_data <- prcomp(x@data, center=F, scale=F)$x
    }else{
      rotated_data <- prcomp(x@data, center=T, scale=T)$x
    }
    SV <- ifelse(x@alphas > 0, 19, 1)
    legend_x <- quantile(rotated_data[, 'PC1'], 0.02)
    legend_y <- quantile(rotated_data[, 'PC2'], 1)
    plot(rotated_data[, 'PC1'], rotated_data[, 'PC2'], type='p', pch=SV, xlab='PC 1', ylab='PC 2')
    legend(x=legend_x, y=legend_y, pch=c(19, 1), legend=c('Support Vectors', 'Non-support Vectors'))
    title(paste('Cost =', x@cost))
  })
plot(euclidean_solution)
setMethod(
  f='predict',
  signature='SVDD',
  definition=function(object, newdata, type=c('radii', 'boolean')){
    # isolate non-zero alpha values (corresponding to support vectors)
    sv_alphas <- object@alphas[object@alphas > 0]
    # calculate k(x_s, x_s) for each support vector
    term1 <- apply(object@support_vectors, 1, function(x) object@kernel(x, x))
    # kernelMult computes (-2) sum(alpha_i k(x_s, x_i)) for each support vector x_s
    term2 <- (-2)*kernelMult(kernel=object@kernel, x=object@support_vectors, z=sv_alphas)
    # kernelPol computes z_i z_j k(x_i, x_j)
    term3 <- sum(kernelPol(kernel=object@kernel, x=object@support_vectors, z=sv_alphas))
    radii <- term1+term2+term3 # term1 and term2 vectors, term3 scalar
    r2 <- mean(radii)
    z1 <- apply(newdata, 1, function(x) object@kernel(x, x))
    z2 <- (-2)*kernelMult(kernel=object@kernel, x=newdata, y=object@support_vectors, z=sv_alphas)
    z_rad <- z1+z2+term3
    if(type=='radii'){
      return(list(mod_r2=r2, newdata_r2=c(z_rad)))
    }else{
      return(as.numeric(z_rad>r2))
    }
  })
predict(euclidean_solution, newdata=as.matrix(mydat2), type='boolean')
predict(euclidean_solution, newdata=as.matrix(mydat2), type='radii')
detect_outliers <- function(SVDD_obj, newdata, plot=T){
  x_mat <- as.matrix(newdata) # in case newdata is a data.frame
  pred <- predict(SVDD_obj, newdata=x_mat, type='radii')
  rad <- pred$newdata_r2
  r2 <- pred$mod_r2
  outlier_ind <- rad>r2
  num_outliers <- sum(outlier_ind)
  if(plot==T){
    plot(rad, type='b')
    abline(h=r2, col='red')
  }
  msg <- paste(num_outliers, 'outliers out of', nrow(x_mat), 'observations detected. \n')
  cat(msg)
  return('outliers'=which(outlier_ind))
}
detect_outliers(SVDD_obj=euclidean_solution, newdata=mydat2)
malhalanobis <- function(covmat){
  rval <- function(x, y = NULL) {
    if (!is(x, 'vector'))
      stop('x must be a vector')
    if (!is(y, 'vector') && !is.null(y))
      stop('y must be a vector')
    if (is(x, 'vector') && is.null(y)) {
      t(x)%*%solve(covmat)%*%x
    }
    if (is(x, 'vector') && is(y, 'vector')) {
      if (!length(x) == length(y))
        stop('number of dimension must be the same on both data points')
      t(x)%*%solve(covmat)%*%y
    }
  }
  return(new('malhalanobis', .Data=rval, kpar=list(covmat=covmat)))
}
my_mal_kernel <- malhalanobis(covmat=cov(mydat))
class(my_mal_kernel)
my_mal_kernel(c(1, 1, 1, 1), c(1, 1, 1, 1)) # test case
mal_solution <- svdd(as.matrix(mydat), k=my_mal_kernel, cost=1)
summary(mal_solution)
plot(mal_solution)
predict(mal_solution, newdata=as.matrix(mydat2), type='boolean')
detect_outliers(SVDD_obj=mal_solution, newdata=mydat2)