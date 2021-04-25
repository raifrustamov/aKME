library(spatstat)
library(harmonicmeanp)
library(BayesFactor)
library(polynom)

myintegral <- function(p, ddim){
  cc = coef(p)
  n = length(cc)
  ss = (0:(n-1))/2
  sum(cc*2^(ss)*gamma(ss+ddim/2)/gamma(ddim/2))
}

QuadratureRoots <- function(ddim, qorder) {
  #compute quadrature roots
  #ddim dimension of space
  #qorder quadrature qorderer
  a = rep(NA, qorder)
  b = rep(NA, qorder)
  xpoly = polynomial(coef = c(0,1)) #f(x)=x
  plist = polylist(polynomial(coef = c(1)))
  
  a[2] = myintegral(xpoly*plist[[1]]*plist[[1]], ddim)
  pnext = xpoly-a[2]
  plist[[2]] = pnext
  
  for (j in 3:(qorder+1)){
    a[j] = myintegral(xpoly*plist[[j-1]]*plist[[j-1]], ddim)/ myintegral(plist[[j-1]]*plist[[j-1]], ddim)
    b[j] = myintegral(plist[[j-1]]*plist[[j-1]], ddim)/myintegral(plist[[j-2]]*plist[[j-2]], ddim)
    pnext = (xpoly-a[j])*plist[[j-1]] - b[j]*plist[[j-2]]
    plist[[j]] = pnext
  }
  
  #roots/knots; note that we do not need the weights, they get cancelled out in t-test
  Roots = solve(plist[[qorder+1]])
  Roots
}

#from https://github.com/maremun/quffka
rnsimp <- function(m){
    S = matrix(0.0, nrow=m, ncol=m + 1)
    mp = m + 1
    for (i in 1:m){
        rv = sqrt(mp / ((m - i + 1.) * m * (m - i + 2.)))
        S[i, i] = (m - i + 1) * rv
        S[i, (i+1):(m+1)] = -rv
    }
    S
}

#generate rays emanating from the origin on which we will project data
proj_mat <- function(ddim){
  #simplex in R^ddim
  simp = rnsimp(ddim)
  
  
  #random rotation https://nhigham.com/2020/04/22/what-is-a-random-orthogonal-matrix/
  G = matrix(rnorm(ddim*ddim), ncol=ddim)
  tmp = qr(G)
  Q = qr.Q(tmp)%*%diag(sign(diag(qr.R(tmp))))
  
  #rotated simplex -- has n_radial=ddim+1 rays emanating from origin
  t(Q%*%simp)
}

#maps every point in the point pattern using aKME process
prep_aKME <- function(ddim, sigmas, n_linear=4, seed=1234) {
  #quadrature roots and weights
  Roots = QuadratureRoots(n_linear, n_linear)
  
  #radial projections -- generate one random frame per sigma
  set.seed(seed)
  W = do.call('rbind',
    lapply(seq_along(sigmas), function(i){
       do.call('rbind', lapply(seq_along(Roots), 
                               function(j) Roots[j]/sigmas[i] * 
                                 proj_mat(ddim)))
    })
  )
  
  #return function that computes aKME embedding
  function(p){
    #project the points onto the lines, c.f. Figure 1
    proj = W%*%t(as.matrix(p))
    
    #embedding, one row for each point of the pattern
    #dim = (#points,  2 x #sigmas x (ddim+1) x n_linear)
    emb = t(rbind(cos(proj), sin(proj)))
    emb
  }
}

#data dimensionality
ddim = 3 
#we use several Gaussian widths at once
#one can repeat the same sigma to get multiple random radial frames for the same sigma
sigma_list = c(1 / 16, 1 / 8, 1 / 4)
get_aKME = prep_aKME(ddim, sigma_list)


#some data -- null case
d1 = rpoispp3(100)
d2 = rpoispp3(100)

#embedding
aKME1 = get_aKME(d1)
aKME2 = get_aKME(d2)

#p-values for t-tests on each dimension of the embeddings
p_ttests = sapply(seq_len(ncol(aKME1)), function(i)
  t.test(aKME1[, i], aKME2[, i])$p.value)

#overall p-value using harmonic mean combo
p_harmonic = p.hmp(p_ttests, L = length(p_ttests), multilevel=FALSE)
p_harmonic

#overall p-value using Cauchy combo
p_cauchy= 0.5 - atan(mean(tan((0.5-p_ttests)*pi)))/pi
p_cauchy

#mean Bayes factor
bf10_mean = mean(sapply(seq_len(ncol(aKME1)), function(i)
  extractBF(ttestBF(aKME1[, i], aKME2[, i]))$bf))
bf10_mean
