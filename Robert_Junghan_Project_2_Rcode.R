# Clear working environment
rm(list=ls())
library(glmnet) # for fitting elastic net models
library(RCurl) # for getting data from the web
library(ggplot2) # for plots
# Options for document compilation
knitr::opts_chunk$set(warning=FALSE, message=FALSE, comment=NA, fig.width=4, fig.height=3)
# #### Begin Problem 1 ####
deriv_g <- deriv3(g~4*x*y+(x+y^2)^2, c('x', 'y'), c('x', 'y'))
my_deriv <- deriv_g(x=1, y=1) # test the function deriv_g()
attr(my_deriv, 'gradient') # checks out
attr(my_deriv, 'hessian') # checks out
newton_update <- function(x, y){
  my_deriv <- deriv_g(x=x, y=y) # evaluate gradient and hessian
  gradient <- matrix(attr(my_deriv, 'gradient'), nrow=2) # extract gradient
  hessian <- matrix(attr(my_deriv, 'hessian'), nrow=2) # extract hessian
  update <- c(x, y) - c(solve(hessian)%*%gradient) # calculate update
  names(update) <- c('x', 'y')
  return(update)
}
dist_type <- 'euclidean' # see ?dist for other options
epsilon <- 1E-9 # convergence criterion
x_current <- c('x'=1, 'y'=1) # inital values
x_prev <- c('x'=Inf, 'y'=Inf) # no previous values yet!
counter <- 0 # for counting iterations
verbose <- TRUE # do you want to get updates?
while(dist(rbind(x_current, x_prev), method=dist_type) > epsilon){
  counter <- counter+1
  x_prev <- x_current
  x_current <- newton_update(x=x_current['x'], y=x_current['y'])
  if(verbose==T){
    msg <- paste0('Iteration ', counter, ': x = (', x_current[1], ', ', x_current[2], ')')
    cat(msg)
    cat('\n')
  }
}
gsection = function(ftn, x.l, x.r, x.m, tol = 1e-9) {
  # applies the golden-section algorithm to minimize ftn
  # we assume that ftn is a function of a single variable
  # and that x.l < x.m < x.r and ftn(x.l), ftn(x.r) >= ftn(x.m)
  # the algorithm iteratively refines x.l, x.r, and x.m and terminates
  # when x.r - x.l <= tol, then returns x.m
  # golden ratio plus one
  gr1 = 1 + (1 + sqrt(5))/2
  # successively refine x.l, x.r, and x.m
  f.l = ftn(x.l)
  f.r = ftn(x.r)
  f.m = ftn(x.m)
  while ((x.r - x.l) > tol) {
    if ((x.r - x.m) > (x.m - x.l)) {
      y = x.m + (x.r - x.m)/gr1
      f.y = ftn(y)
      if (f.y <= f.m) {
        x.l = x.m
        f.l = f.m
        x.m = y
        f.m = f.y
      } else {
        x.r = y
        f.r = f.y
      }
    } else {
      y = x.m - (x.m - x.l)/gr1
      f.y = ftn(y)
      if (f.y <= f.m) {
        x.r = x.m
        f.r = f.m
        x.m = y
        f.m = f.y
      } else {
        x.l = y
        f.l = f.y
      }
    }
  }
  return(x.m)
}
f <- function(x) {
  y <- 4*x[1]*x[2]+(x[1]+x[2]^2)^2
  return(y)
}
gradf<- function (x)
{##Calculate First Derivative of Function to x[1] (f1)
  f1<-4*x[2]+2*(x[1] + x[2]^2)
  ##Calculate First Derivative of Function to x[2] (f2)
  f2<-4 * x[1] + 2 * (2 * x[2] * (x[1] + x[2]^2))
  return(c(f1, f2))
}
line.search <- function(f, x, gradf, tol = 1e-9, a.max = 100) {
  # x and gradf are vectors of length d
  # g(a) =f(x +a*gradf) hasa local minumum at a,
  # within a tolerance
  # if no local minimum is found then we use 0 or a.max for a
  # the value returned is x + a*y
  if (sum(abs(gradf)) == 0) return(x) # g(a) constant
  g <- function(a) return(f(x - a*gradf))
  # find a.l < a.m < a.r such that
  # g(a.m) >=g(a.l) and g(a.m) >= g(a.r)
  # a.l
  a.l <- 0
  g.l <- g(a.l)
  # a.m
  a.m <- 1
  g.m <- g(a.m)
  while ((g.m > g.l) & (a.m > tol)) {
    a.m <- a.m/2
    g.m <- g(a.m)
  }
  # if a suitable a.m was not found then use 0 for a
  if ((a.m <= tol) & (g.m >= g.l)) return(x)
  # a.r
  a.r <- 2*a.m
  g.r <- g(a.r)
  while ((g.m >= g.r) & (a.r < a.max)) {
    a.m <- a.r
    g.m <- g.r
    a.r <- 2*a.m
    g.r <- g(a.r)
  }
  # if a suitable a.r was not found then use a.max for a
  if ((a.r >= a.max) & (g.m > g.r)) return(x - a.max*gradf)
  # apply golden-section algorithm to g to find a
  a <- gsection(g, a.l, a.r, a.m)
  return(x - a*gradf)
}
descent <- function(f,gradf, x0, tol = 1e-9, n.max = 100) {
  # steepest descent algorithm
  # find a local minimum of f starting at x0
  # function gradf is the gradient of f
  x <- x0
  x.old <- x
  x <- line.search(f, x, gradf(x))
  n <- 1
  while (f(x.old)-(f(x)> tol) & (n < n.max)) {
    x.old <- x
    x <- line.search(f, x, gradf(x))
    n <- n + 1
  }
  return(x)
}
descent(f,gradf,c(1,0) )
# #### Begin Problem 2 ####
# load challenger data
library(mcsm)
data(challenger)
my_deriv_func <- deriv3(l~y*log(1/(1+exp(-b0-b1*x)))+(1-y)*log(1/(1+exp(b0+b1*x))), namevec=c('b0', 'b1'), function.arg=c('x', 'y', 'b0', 'b1'))
deriv_all <- function(b0, b1, deriv_func){
  # create gradient and hessian object for each (x, y) pair
  deriv_list <- mapply(deriv_func,
                       x=challenger$temp, y=challenger$oring,
                       MoreArgs=list(b0=b0, b1=b1),
                       SIMPLIFY=FALSE
  )
  # extract all gradient components
  gradient_list <- lapply(deriv_list, function(g) matrix(attr(g, 'gradient'), nrow=2))
  # and sum them to get overall gradient at (b0, b1)
  gradient <- Reduce('+', gradient_list)
  # extract all hessian components
  hessian_list <- lapply(deriv_list, function(g) matrix(attr(g, 'hessian'), nrow=2))
  # and sum them to get overall hessian at (b0, b1)
  hessian <- Reduce('+', hessian_list)
  return(list('gradient'=gradient, 'hessian'=hessian))
}
deriv_all(b0=0.5, b1=0, deriv_func=my_deriv_func)
newton_update <- function(b0, b1){
  # get gradient and hessian
  my_deriv <- deriv_all(b0=b0, b1=b1, deriv_func=my_deriv_func)
  # calculate update
  update <- c(b0, b1) - c(solve(my_deriv$hessian)%*%my_deriv$gradient)
  # retun update
  names(update) <- c('b0', 'b1')
  return(update)
}
dist_type <- 'euclidean' # see ?dist for other options
epsilon <- 1E-9 # convergence criterion
##x_current <- c('b0'=0, 'b1'=0) # inital values
x_current <- c('b0'=0, 'b1'=0) # inital values
x_prev <- c('b0'=Inf, 'b1'=Inf) # no previous values yet!
counter <- 0 # for counting iterations
verbose <- TRUE # do you want to get updates?
while(dist(rbind(x_current, x_prev), method=dist_type) > epsilon){
  counter <- counter+1
  x_prev <- x_current
  x_current <- newton_update(b0=x_current['b0'], b1=x_current['b1'])
  if(verbose==T){
    msg <- paste0('Iteration ', counter, ': b = (',
                  paste(signif(x_current, 5), collapse=', '),
                  ')')
    cat(msg)
    cat('\n')
  }
}
glm_mod <- glm(oring~temp, data=challenger, family=binomial)
glm_mod$coefficients
y <- challenger$oring
x <- challenger$temp
Z <- cbind(rep(1, nrow(challenger)), challenger$temp)
find_pis <- function(beta, x){
  # beta a vector of parameter estimates, x a vector of independent variable observations
  pis <- exp(beta%*%t(Z))/(1+exp(beta%*%t(Z)))
  return(pis)
}
find_W <- function(pis){
  W <- diag(c(pis*(1-pis)))
  return(W)
}
find_update <- function(beta, Z, W, y, pis){
  numerator <- t(Z)%*%(y-t(pis))
  denominator <- t(Z)%*%W%*%Z
  update <- solve(denominator)%*%numerator
  return(beta+c(update))
}
beta <- c(0, 0)
old_beta <- c(NA, NA)
counter <- 0
decimal_places_of_precision <- 7
verbose <- TRUE
while(!identical(beta, old_beta)){
  pis <- find_pis(beta, x)
  W <- find_W(pis)
  # set old beta equal to current beta before updating
  old_beta <- beta
  # replace beta with updated estimates
  beta <- find_update(beta, Z, W, y, pis)
  # round to the number of decimal places desired
  beta <- round(beta, decimal_places_of_precision)
  counter <- counter+1
  # show us the current value of beta if "verbose" is TRUE
  if(verbose==T){
    msg <- paste0('Iteration ', counter, ': Beta = (', beta[1], ', ', beta[2], ')')
    cat(msg)
    cat('\n')
  }
}
# see our estimates
print(beta)
# see the glm estimates
glm(oring~temp, data=challenger, family=binomial)$coefficients
# #### Begin Problem 3 ####
myfile <- getURL("http://www-bcf.usc.edu/~gareth/ISL/Credit.csv") # grab the content of this page
credit <- read.csv(textConnection(myfile), header=T) # read the page content as a .csv
credit$X <- NULL # get rid of row ID column
# find only the column with a numeric class (or integer)
num_cols <- names(credit)[sapply(credit, class)%in%c('numeric', 'integer')]
# now subset the data set to only take the numeric columns
credit <- credit[, num_cols]
k <- 10 # for 10-fold cross-validation
n <- nrow(credit) # sample size
seg_size <- floor(n/k) # approx size of each segment
leftovers <- n%%k # n mod k
seg_assignments <- rep(1:k, each=seg_size) # create segment assignments
if(leftovers > 0){
  seg_assignments <- c(seg_assignments, 1:leftovers) # add "leftover" observations if necessary
}
set.seed(8675309) # set random seed (is it bad that I always use 867-5309?)
seg_assignments <- sample(seg_assignments) # randomly permute segment assignments
get_val_errs <- function(alpha, lambdas, train_data, train_response, test_data, test_response){
  # fit the elastic net for all lambda values
  fit <- glmnet(train_data, train_response, alpha=alpha, lambda=lambdas)
  # use the models fit to predict test data
  preds <- predict(fit, newx=test_data)
  # calculate the mean squared error in the test data for each value of lambda
  MSEs <- apply(preds, 2, function(x) mean((test_response-x)^2))
  # return a data.frame with alpha, lambdas, and corresponding MSEs
  return(data.frame(alpha=rep(alpha, length(MSEs)), lambda=lambdas, MSE=MSEs))
}
lambda_vec <- seq(40, 0, by=-1) # lambdas to try for each alpha
alphas <- seq(0, 1, by=0.01) # alphas to try
# preallocate an empty list with a slot for each of the k validation folds
results <- vector(mode='list', length=k)
for(i in 1:k){ # k validation folds
  # define training data matrix
  x <- model.matrix(Balance~., subset(credit, seg_assignments!=i))[,-1]
  # and the training response vector
  y <- credit$Balance[seg_assignments!=i]
  # define the test data matrix
  x_test <- model.matrix(Balance~., subset(credit, seg_assignments==i))[,-1]
  # and the test response vector
  y_test <- credit$Balance[seg_assignments==i]
  # succinctly loop over alphas, finding MSE for each lambda as we go
  result_list <- sapply(alphas, get_val_errs, lambdas=lambda_vec, train_data=x, train_response=y, test_data=x_test, test_response=y_test, simplify=F)
  # convert the list into a data.frame using rbind()
  result_df <- Reduce(rbind, result_list)
  # assign a field for k
  result_df$k <- i
  # deposit results into preallocated list
  results[[i]] <- result_df
}
# average the MSEs from all k folds
results_all <- Reduce(rbind, results)
results_agg <- aggregate(MSE~alpha+lambda, data=results_all, mean)
bestMSE <- results_agg[which.min(results_agg$MSE),]
bestMSE
# plot a heat map of the MSEs by alpha and lambda
ggplot(results_agg, aes(x=alpha, y=lambda, z=MSE))+
  geom_tile(aes(fill=MSE))+
  stat_contour(bins=12)+
  geom_point(data=bestMSE, color='red', size=3)+
  geom_text(data=bestMSE, aes(label=paste0('(', alpha, ', ', lambda, ')')), color='white', vjust=-1)+
  theme_bw(base_size=10)
# calculate std dev of MSEs
results_SEs <- aggregate(MSE~alpha+lambda, data=results_all, sd)
# divide by sqrt(k)
results_SEs$MSE <- results_SEs$MSE/sqrt(k)
# rename column for less confusion
names(results_SEs) <- gsub('MSE', 'StdErr', names(results_SEs))
SE_best_model <- results_SEs[results_SEs$alpha==bestMSE$alpha & results_SEs$lambda==bestMSE$lambda, ]
SE_best_model
MSE_max_tolerable <- bestMSE$MSE+SE_best_model$StdErr
candidate_lambdas <- results_agg[results_agg$MSE<MSE_max_tolerable, c('alpha', 'lambda')]
# plot a heat map of the MSEs by alpha and lambda
ggplot(results_agg, aes(x=alpha, y=lambda, z=MSE))+
  geom_tile(aes(fill=MSE))+
  stat_contour(bins=12)+
  geom_point(data=bestMSE, color='red', size=3)+
  geom_text(data=bestMSE, aes(label=paste0('(', alpha, ', ', lambda, ')')), color='white', vjust=-1)+
  geom_tile(data=candidate_lambdas, aes(fill=NULL, z=NULL), fill='white', alpha=0.5)+
  theme_bw(base_size=10)
# define data matrix
x <- model.matrix(Balance~., credit)[,-1]
# and the response vector
y <- credit$Balance
# define a function
find_nonzero_parameters <- function(alpha, lambda, data, response){
  fit <- glmnet(data, response, alpha=alpha, lambda=c(lambda))
  betas <- fit$beta@x # @ accesses a slot in an s4 object
  num_nonzero <- sum(abs(betas)>0.001)
  return(num_nonzero)
}
# use mapply() to iterate over alpha and lambda together
candidate_lambdas$num_params <- mapply(find_nonzero_parameters,
                                       alpha=candidate_lambdas$alpha,
                                       lambda=candidate_lambdas$lambda,
                                       MoreArgs=list(data=x, response=y))
# how many params does the smallest model have?
min(candidate_lambdas$num_params)
final_candidate_lambdas <- candidate_lambdas[candidate_lambdas$num_params==4,]
final_candidate_lambdas
# merge with MSEs by alpha and lambda
final_candidate_lambdas <- merge(final_candidate_lambdas, results_agg, by=c('alpha', 'lambda'), all.x=T, all.y=F)
# best, simplest model using 1-SE rule
final_candidate_lambdas[which.min(final_candidate_lambdas$MSE),]