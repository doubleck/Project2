

# Clear working environment
rm(list=ls())
library(ggplot2) # for plots
library(kernlab)

# Options for document compilation
knitr::opts_chunk$set(warning=FALSE, message=FALSE, comment=NA, fig.width=6, fig.height=5)



##Robert's Path
mydat <- read.table('/home/robert/cloud/Classes/STA6106 Stat Computing/Project2/Project3/training dataset.txt')
mydat2 <- read.table('/home/robert/cloud/Classes/STA6106 Stat Computing/Project2/Project3/data_project3.txt')

##Jung-Han's Path
# mydat <- read.table("E:/Cloud Storage/Dropbox/Life long study/Ph.D/Lecture/2014 Fall/Statistical Computing/Project 2/Project2/Project3/training dataset.txt")
# mydat2 <- read.table("E:/Cloud Storage/Dropbox/Life long study/Ph.D/Lecture/2014 Fall/Statistical Computing/Project 2/Project2/Project3/data_project3.txt")



my_kernel <- vanilladot()



my_matrix <- as.matrix(mydat)
my_kernel(my_matrix[1, ], my_matrix[2, ]) # dot prod using kernel function
crossprod(my_matrix[1, ], my_matrix[2, ]) # dot prod using base R function
my_matrix[1, ] %*% my_matrix[2, ] # old school matrix multiplication operator



class(my_kernel)



getClass('vanillakernel')



K <- kernelMatrix(kernel=my_kernel, x=as.matrix(mydat))



dim(K)
class(K)



my_c <- (-1)*diag(K)
my_H <- (2)*K
my_A <- rep(1, nrow(my_matrix))
my_b <- 1
my_l <- rep(0, nrow(my_matrix))
my_u <- rep(1, nrow(my_matrix))
my_r <- 0
my_solution <- ipop(c=my_c, H=my_H, A=my_A, b=my_b, l=my_l, u=my_u, r=my_r, maxiter=300, margin=1E-6)
my_alphas <- my_solution@primal # use @ symbol to access s4 slot



my_alphas <- round(my_alphas, abs(log10(1E-6)))



sum(my_alphas)



any(my_alphas < 0)



# define an S4 class "SVDD"
SVDD <- setClass(
  # set the name of the formal S4 class being defined
  'SVDD',
  # Define the slots and their types
  slots=c(
    data='matrix',
    alphas='numeric',
    support_vectors='matrix',
    kernel='kernel',
    tolerance='numeric',
    fun_call='call'),
  # A function to see if the object is a valid SVDD object
  validity=function(object){
    if(class(object@data)!='matrix'){
      return('data should be of class matrix.')
      }else if(class(object@alphas)!='numeric'){
        return('alphas should be of class numeric.')
        }else if(any(object@alphas < 0)){
          return('alphas should all be non-negative.')
          }else if(class(object@support_vectors)!='matrix'){
            return('support_vectors should be of class matrix.')
            }else if(!inherits(object@kernel, 'kernel')){
              return('kernel should be a subclass of the kernel class.')
              }else if(class(object@tolerance)!='numeric' || object@tolerance > 1){
                return('tolerance should be a numeric value less than one.')
                }else if(class(object@fun_call)!='call'){
                  return('call should be of class function')
                  }else{
                    return(TRUE)
                    }
    }
  )

# create a summary method for class SVDD
setGeneric('summary', # reserve the method name "summary"
           def=function(object){
             standardGeneric('summary')
             }
           )

setMethod( # define the method's behavior
  f='summary',
  signature='SVDD',
  definition=function(object){
    k <- class(object@kernel)
    num_obs <- nrow(object@data)
    num_svs <- nrow(object@support_vectors)
    msg <- paste0('SVDD object using kernel function ', k, ':\n', num_svs, ' support vectors identified out of ', num_obs, ' observations.\n\n')
    cat(msg)
    return(object@support_vectors)
    })



svdd <- function(x, k, tol=1E-6){
  K <- kernelMatrix(kernel=k, x=x)
  n <- nrow(x)
  solution <- ipop(c=(-1)*diag(K), 
                   H=2*K, 
                   A=rep(1, n), 
                   b=1, 
                   l=rep(0, n), 
                   u=rep(1, n), 
                   r=0, 
                   maxiter=300, 
                   margin=tol)
  
  alphas <- solution@primal
  alphas <- round(alphas, abs(log10(tol))) # round based on tol
  
  svs <- x[alphas > 0, ] # these are the support vectors
  
  return(new('SVDD',
             data=x,
             alphas=alphas, 
             support_vectors=svs, 
             kernel=k, 
             tolerance=tol, 
             fun_call=match.call()
             ))
  }



# supplying no value to tol implies it should use the default value tol = 1E-6
euclidean_solution <- svdd(x=my_matrix, k=my_kernel)
summary(euclidean_solution)



setMethod(
  f='plot',
  signature='SVDD',
  definition=function(x, y){
    if(class(x@kernel)=='vanillakernel'){
      rotated_data <- prcomp(x@data, center=F, scale=F)$x
    }else{
     rotated_data <- prcomp(x@data, center=T, scale=T)$x 
    }
    SV <- ifelse(x@alphas > 0, 19, 1)
    legend_x <- quantile(rotated_data[, 'PC1'], 0.02)
    legend_y <- quantile(rotated_data[, 'PC2'], 1)
    plot(rotated_data[, 'PC1'], rotated_data[, 'PC2'], type='p', pch=SV, xlab='PC 1', ylab='PC 2')
    legend(x=legend_x, y=legend_y, pch=c(19, 1), legend=c('Support Vectors', 'Non-support Vectors'))
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



setClass(
  # name the class
  'malhalanobis',
  # define the slots
  slots=c(.Data='function',
          kpar='list'),
  validity=function(object){
    if(class(object@.Data)!='function'){
      return('.Data should be of class function.')
      }else if(class(object@kpar)!='list'){
        return('kpar should be of class list.')
        }else if(length(kpar)>1){
          return('kpar should only have one element; covmat.')
          }else{
            return(TRUE)
            }
    },
  contains='kernel'
  )



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



mal_solution <-  svdd(as.matrix(mydat), k=my_mal_kernel)
summary(mal_solution)
plot(mal_solution)



predict(mal_solution, newdata=as.matrix(mydat2), type='boolean')



detect_outliers(SVDD_obj=mal_solution, newdata=mydat2)





