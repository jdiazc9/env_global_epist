# Functions for polynomial approxaimation of a density first for recovering the density and sampling
# in the approximated density by the accept-reject algorithm then, fun specific to the simulated
# example of the paper

############################################################################
## 1. General fun for approximated distribution and sampling in it

# now by simulation, accept reject algo for sampling in the distribution accept-reject algo that uses
# polynomial approximation simulate r.v. in the approximation
# accept_reject_copy=function(poly_approx,a,b,N_sim=1000){

accept_reject = function(poly, weight, a, b, N_sim) { #  = 1000
  X = rbeta(N_sim, a, b)
  U = runif(N_sim)
#   most_support = (weight>10^(-4)*sum(weight))
M = max(poly)*gamma(a)*gamma(b)/gamma(a+b)
#   M = max(poly*most_support)*gamma(a)*gamma(b)/gamma(a+b)
#   M = max(poly*most_support)#*gamma(a)*gamma(b)/gamma(a+b)
  Ngrid = length(poly)
  # since the poly function is discretized on a grid of size Ngrid, I recover the right evaluation of
  # the polynomial by using the index floor(X*Ngrid)+1
  Y = X[(M * U) < (poly[floor(X * Ngrid) + 1])]
  Y = Y[!is.na(Y)]
  Y = sample(x = Y, size = N_sim, replace = TRUE)
  Y
}

# computes the two parameters of the weight function of G-Jacobi polynomials
first_weight_fun = function(mu1, mu2) {
  (mu1 - mu2)/(mu2 - mu1^2) - 1
}

second_weight_fun = function(mu1, mu2) {
  mu1 * (mu1 - mu2)/(mu2 - mu1^2)
}

# gives the polynomial P in approx f = w * P, plus a and b, coef of w seen as beta(a,b)
compute_Jacobi_unweighted = function(moments, N_moments = length(moments), xgrid) {
    # compute moments for the weight function, warning starts at 2
  first_weight = first_weight_fun(moments[1], moments[2])  # aka p
  second_weight = second_weight_fun(moments[1], moments[2])  # aka q
  
  # compute the weight function
  weight_fun = function(x) {
    jacobi.g.weight(x, first_weight, second_weight)
  }
  
  # define the first n G-Jacobi polynomials that are normalized, ie in the weighted scalar product, norm
  # is 1 with unnormalized poly, the vector of the norms obtained by jacobi.g.inner.products has 0's...
  Jacobi_list <- jacobi.g.polynomials(N_moments, first_weight, second_weight, normalized = TRUE)
  # warnings: contains moment of order zero as first element
  # Jacobi_squared_norm <- jacobi.g.inner.products(n, first_weight, second_weight)
  
  # first Jacobi polynomial, WARNING : unweighted
  series_unweighted_Jacobi = coef(Jacobi_list[[1]])^2
  series_Jacobi = coef(Jacobi_list[[1]])^2 * weight_fun(xgrid)
  
  # then, n following increments
  for (i in 1:N_moments) {
    G_i = Jacobi_list[[i+1]]
    lambda_i = sum(coef(G_i) * c(1,moments[1:i]))
    G_i = as.function(G_i)
    P_i = function(x) {
      lambda_i * G_i(x)
    }
    series_unweighted_Jacobi = P_i(xgrid) + series_unweighted_Jacobi
    series_Jacobi = weight_fun(xgrid) * P_i(xgrid) + series_Jacobi
  }
  
  # sum them with cumsum to get the right polynomials
  a = second_weight
  b = first_weight - second_weight + 1
  # series_unweighted_Jacobi=apply()
  poly_a_b = list(poly = series_unweighted_Jacobi, 
                  poly_weighted = series_Jacobi, 
                  a = a, 
                  b = b, 
                  weight = weight_fun(xgrid),
                  xgrid = xgrid)
  poly_a_b
  # accept_reject(series_unweighted_Jacobi,a,b,N_sim=N_sim)
}

# create object of class "momentified.res"
momentify = function(moments, N_moments = length(moments), N_sim = 1000, xgrid = seq(0,1,length = 200)) {
  poly_a_b = compute_Jacobi_unweighted(moments, N_moments, xgrid = xgrid)
  psample = accept_reject(poly_a_b$poly, poly_a_b$weight, 
                     poly_a_b$a, poly_a_b$b, N_sim)
  res = list(psample = psample,
             xgrid = xgrid,
             approx_density = poly_a_b$poly_weighted
  )
  class(res) = "momentify"
  res
}

# A few graphical functions

plot.momentify = function(res, ...){
  plot(x = res$xgrid, y = res$approx_density, type = 'l', ...)
}

hist.momentify = function(res, ...){
  hist(x = res$psample, ...)
}