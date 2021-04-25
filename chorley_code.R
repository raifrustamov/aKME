library(spatstat)
library(harmonicmeanp)
library(BayesFactor)

radialGaussHermiteData <- function(k) {
	if (k == 4) {
		Roots4 = c(
			0.3961205684809970482046989062,
			1.1767346714185921750900701015,
			2.2010676629189581946543092946,
			3.4836109980286615716953100856
		)
		Weights4 = c(
			0.2279981086730596270303131898,
			0.538464759426994159094177725,
			0.2215779133653168055615931299,
			0.0119592185346294083139159552
		)
	    list(x = Roots4, w = Weights4)
	} else {
		stop(
			"Please refer to the values in
			http://www.jaeckel.org/RootsAndWeightsForRadialGaussHermiteQuadrature.cpp"
		)
	} 
}

#maps every point in the point pattern using aKME process
get_aKME <- function(p, sigmas, n_radial = 4, n_linear = 4) {
	#radial unit vectors
	theta = seq(0, pi, length.out = n_radial + 1)[-(n_radial + 1)]
	W0 = cbind(cos(theta), sin(theta))
	
	#quadrature roots and weights
	temp = radialGaussHermiteData(n_linear)
	Roots = temp$x
	Weights = temp$w
	
	W = do.call('rbind', lapply(seq_along(Roots), function(i) Roots[i] * W0))
	coeff = rep(sqrt(Weights / n_radial), each=nrow(W0))
	
	#project the points onto the lines, c.f. Figure 1
	proj = W%*%rbind(p$x, p$y)
	
	#embedding: for each sigma wrap on the circle of corresponding radius
	emb = lapply(seq_along(sigmas), function(i)
		t(rbind(cos(proj / sigmas[i]) * coeff, sin(proj / sigmas[i]) * coeff)))
	
	#concatenated embedding, one row for each point of the pattern
	#dim = (#points,  #sigmas x 2 x # n_radial x n_linear)
	do.call('cbind', emb)
}

#load data
d = chorley
d1 = split(d)$lung
d2 = split(d)$larynx

#we use several Gaussian widths at once
#these values are for [0,1]x[0,1] window
sigma_list0 = c(1 / 16, 1 / 8, 1 / 4) 
sigma_list = diameter(Window(d)) * sigma_list0 / sqrt(2)

aKME1 = get_aKME(d1, sigma_list)
aKME2 = get_aKME(d2, sigma_list)

#p-values for t-tests on each dimension of the embeddings
p_ttests = sapply(seq_len(ncol(aKME1)), function(i)
	t.test(aKME1[, i], aKME2[, i])$p.value)

#overall p-value using harmonic mean combo
p_harmonic = p.hmp(p_ttests, L = length(p_ttests), multilevel=FALSE)
#0.9654329

#overall p-value using Cauchy combo
p_cauchy = 0.5 - atan(mean(tan((0.5-p_ttests)*pi)))/pi
#0.6727789

#mean Bayes factor
bf10_mean = mean(sapply(seq_len(ncol(aKME1)), function(i)
	extractBF(ttestBF(aKME1[, i], aKME2[, i]))$bf))
#0.2443478