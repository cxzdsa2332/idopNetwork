#' @title calculate log-sum-exp values
#' @param v numeric vector
#' @return log-sum-exp values
#' @examples logsumexp(c(100, 1000, 10000))
#' @export
logsumexp <- function(v) {
  vm = max(v)
  log(sum(exp(v-vm))) + vm
}

#' @title min-max normalization
#' @param x numeric vector
#' @param z scalar add minimum value to avoid 0
#' @return normalized vector
#' @examples normalization(runif(100,min = -100, max = 100))
#' @export
normalization <- function(x, z = 0.2){(x-min(x))/(max(x)-min(x))+z}

#' @title curve fit with modified logistic function
#' @param mu_par vector with five number
#' @param times vector of time point
#' @return numeric vector with the same length to times
#' @examples get_mu(mu_par = 1:5, times = 1:14)
#' @export
get_mu <- function(mu_par,times){
  mu <- mu_par[1]/(1 + mu_par[2] * exp(-mu_par[3] * times)) - (mu_par[4] * exp(-mu_par[5] * times))
  return(mu)
}

#' @title generate mean vectors with ck and stress condition
#' @param par vector with ten number, first five for ck and the rest for stress
#' @param times vector of time point
#' @return numeric vector with the double length to times
#' @examples get_mu2(par = 1:10, times = 1:14)
#' @export
get_mu2 <- function(par,times){
  mu2 <- c(get_mu(par[1:5],times),get_mu(par[6:10],times))
  return(mu2)
}

#' @title generate standard SAD1 covariance matrix
#' @param par vector with two number for SAD1 covariance matrix
#' @param n scalar indicate length of time d
#' @return SAD1 covariance matrix
#' @examples get_SAD1_covmatrix(par = c(2,0.5), n = 14)
#' @export
get_SAD1_covmatrix <- function(par,n){
  phi <- par[1]; gamma <- par[2];
  sigma <- array(dim=c(n,n))
  #formula 1, diag element
  diag(sigma) <- sapply(1:n, function(c)(1-phi^(2*c))/(1-phi^2) )
  #formula 2, non-diag element
  sigma[lower.tri(sigma)] <- do.call(c,lapply(1:(n-1),function(c)phi^seq(1:(n-c))*diag(sigma)[c]))
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(t(sigma))]
  return(gamma^2*sigma)
}

#' @title generate biSAD1 covariance matrix
#' @param par vector with four number, first two for ck and the rest for stress
#' @param n1 scalar indicate length of time1
#' @param n2 scalar indicate length of time2
#' @return biSAD1 covariance matrix
#' @examples get_biSAD1(par=c(2,0.5,2,0.1),n1=4, n2 = 5)
#' @export
get_biSAD1 <- function(par, n1, n2){
  sig1 <- get_SAD1_covmatrix(par[1:2],n1)
  sig2 <- get_SAD1_covmatrix(par[3:4],n2)
  n = n1+n2
  sig = array(0, dim=c(n,n))
  sig[1:n1,1:n1] = sig1
  sig[(n1+1):(n1+n2),(n1+1):(n1+n2)] = sig2

  return(sig)
}

#' @title generate legendre matrix
#' @importFrom orthopolynom legendre.polynomials polynomial.values scaleX
#' @param legendre_order the order of legendre polynomials
#' @param x vector equal to the x value for legendre polynomials(in this case times)
#' @return the polynomials value of each order
#' @examples get_legendre_matrix(1:14,4)
#' @export
get_legendre_matrix <- function(x,legendre_order){
  legendre_coef <- legendre.polynomials(n = legendre_order, normalized=F)
  legendre_matrix <- as.matrix(as.data.frame(polynomial.values(
    polynomials = legendre_coef, x = scaleX(x, u = -1, v = 1))))
  colnames(legendre_matrix) <- paste0("legendre_",0:legendre_order)
  return(legendre_matrix[,2:(legendre_order+1)])
}

#' @title use legendre polynomials to fit a given data
#' @importFrom stats lm coef
#' @param legendre_order scalar of legendre polynomials
#' @param x vector equal to the x value for legendre polynomials(in this case times)
#' @param y vector equal to the y observed data(in this case generic effect)
#' @return the polynomials coefficients
#' @examples get_legendre_par(14:1,4,1:14)
#' @export
get_legendre_par <- function(y,legendre_order,x){
  #lm_method
  legendre_par <- as.numeric(coef(lm(y~get_legendre_matrix(x,legendre_order))))
  return(legendre_par)
}

#' @title generate curve based on legendre polynomials
#' @importFrom orthopolynom legendre.polynomials polynomial.values scaleX
#' @param par vector of legendre polynomials coefficients
#' @param x vector equal to the x value for legendre polynomials(in this case times)
#' @return the polynomials value
#' @examples legendre_fit(rep(1,5),1:14)
#' @export
legendre_fit <- function(par,x){
  legendre_order = length(par)
  fit <- sapply(1:length(par), function(c)
    par[c]*legendre.polynomials(n=legendre_order, normalized=F)[[c]])
  legendre_fit <- as.matrix(as.data.frame(polynomial.values(
    polynomials = fit, x = scaleX(x, u = -1, v = 1))))
  x_interpolation <- rowSums(legendre_fit)
  return(x_interpolation)
}

#' @title make color more dark
#' @importFrom grDevices rgb col2rgb
#' @param color hex color code
#' @param factor scalar for darken level
#' @return darkened hex color code
#' @examples darken("#FF0000")
#' @export
darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
