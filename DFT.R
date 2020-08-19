# Version 1.0.1
# Copyright 2011 Neil Stewart
# The program is distributed under the terms of the GNU General Public License
#
# DFT.R is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# DFT.R is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# <http://www.gnu.org/licenses/>
#
# Version history
# 1.0.1 Fixed bug in make.S to include abs()
#       Fixed bug in DFT.choices to allow parallel processing
#       Set seed in DFT.choices to be a random integer instead of based on clock time      
#

###############################################################################
#
# Code for multiattribute decision field theory
#
###############################################################################

dyn.load("DFT.so")
dyn.load("DFT_parallel.so")

DFT <- function(C, M, W, S, T, P, sigma, iterations, seed, parallel=FALSE) {
#
# Calls the C library (either the normal single-threadversion or the multi-threaded parallel version)
# C is a square contrast matrix 
# M is the option attribute matrix, with rows for each option and colums for each dimension
# W is a vector of weights for each dimension. NOTE: W is passed to the library in cumulative form
# S is the inhibition matrix
# T is a vector of accumulator thresholds
# P is a vector of accumulator starting points
# iterations in an integer giving the number of times to repeat the simulation
# seed is the seed for the random number generator. NOTE: The DFT function sets this using the system time. NOTE: Because of the dynamic nature of the threading, using the same seed will give different results in the parallel version of the model
# Set parallel=TRUE to use the parallel version
#
	counter <- rep(0, nrow(M))
	if(parallel)
		model.version <- "simulate_DFT_parallel"
	else
		model.version <- "simulate_DFT"	

	.C(model.version,
		as.double(t(C)), # t() to pass vector ordered by rows then columns, as C expects
		as.double(t(M)), 
		as.double(cumsum(W)), # W is passed in cumulative form
		as.double(t(S)),
		as.double(T),
		as.double(P),
		as.double(sigma),
		as.integer(nrow(M)),
		as.integer(ncol(M)),
		as.integer(iterations),
		as.integer(seed),
		as.integer(counter)
	) [[12]]
}


###############################################################################
#
# Extra code for representing gambles and choices, and converting these into
# a form for the DFT code
#
###############################################################################

gamble <- function(x, p) {
#
# gamble(x=c(x1, x2, ...), p=c(p1, p2, ...)) returns a matrix with rows for
# each branch (pi chance of xi) and columns labeled "x" and "p" for each
# attribute
#
    g <- matrix(c(x,p), ncol=2)
    colnames(g) <- c("x", "p")
    g
}


attributes.from.gamble <- function(gamble, attribute.type) {
#
# attributes.from.gamble(gamble=g, attribute.type="x"|"p") returns either
# the amounts or the probabilities from a gamble in a one column matrix
#
    gamble[ ,attribute.type, drop=F]
}


matrix.rep <- function(Z, n) {
#
# Repeats matrix Z to the right of itself n times
# Used in make.states()
#
	matrix(Z, byrow=FALSE, ncol=ncol(Z)*n, nrow=nrow(Z))
}


make.states <- function(choice, independent=TRUE) {
#
# Convert two gambles into a states of the world representation (an M matrix
# and W vector). Non-independent gambles must have the same branch
# probabilities---this is not checked
#
	n <- length(choice)
    if(independent) {
		W <- attributes.from.gamble(choice[[1]], "p")
		M <- matrix(attributes.from.gamble(choice[[1]], "x"), ncol=nrow(choice[[1]]))
		for (i in 2:n) {
			W <- rep(W, nrow(choice[[i]])) * 
				rep(attributes.from.gamble(choice[[i]], "p"), each=length(W))
			M <- rbind(
				matrix.rep(M, nrow(choice[[i]])), 
				rep(attributes.from.gamble(choice[[i]], "x"), each=ncol(M))
			)
		}
    } else {
		W <- c(attributes.from.gamble(choice[[1]], "p"))
		M <- c(attributes.from.gamble(choice[[1]], "x"))
		for (i in 2:n)
			M <- rbind(M, c(attributes.from.gamble(choice[[i]], "x")))
	}
	rownames(M) <- NULL
    list(M=M, W=W)
}


make.C <- function(n) {
#
# Make default contrast matrix C where each option is compared with the average of the others
#
	C <- matrix(-1/(n-1), nrow=n, ncol=n)
	diag(C) <- 1
	C
}


make.S <- function(M, W, alpha, beta, q=2, r=2) {
#
# Make the inhibition matrix based on the distances between the outcome
# vectors, weighted by the state probabilities. beta gives level of
# inhibition between identical options. 1-beta is the level of self
# excitation. alpha is in units of 1/outcome, with larger values giving less
# inhibition between options. q gives the form of the generalisation
# gradient. r gives the Minkowski metric. This function is my addition, not
# based on the Roe et al. (2001) paper, in which elements of S are free
# parameters. See Hotaling, Busemeyer, and Li (2010) for an alternative
# method for generating S in multiattribute choice
#
	n <- nrow(M)
	S <- matrix(nrow=n, ncol=n)
	for (i in 1:n) {
		for(j in 1:n) {
			S[i,j] <- (W %*% (abs(M[i,] - M[j,])^r))^(1/r)
		}
	}
	diag(n) - beta * exp(-alpha * S^q)
}

DFT.choices <- function(choices, independent=TRUE, alpha=NULL, beta=NULL, q=NULL, r=NULL, sigma,
	threshold, S=NULL, iterations=100, parallel=TRUE) {
#
# For each choice in choices = list(c1, c2, ...) call the simulate function iterations times
# Returns the proportion of times each option won the accumulations for each choice
#
# Either: 
# (a) a pre-defined inhibition matrix S should be passed, or 
# (b) values for alpha, beta, q, and r should be passed. These and are used
# by make.S() to make an inhibition matrix for each choice based on the
# distance between options
#
# sigma is the standard deviation of the normally distributed noise added
#   each step in the accumulation
# threshold is the value the accumulators race to
# iterations in the number of times to repeat the accumulation process for each choice
# Set parallel=TRUE to use the multi-thread code or FALSE for the single-thread code
#
	predictions <- c()
	for (i in 1:length(choices)) {
		n <- length(choices[[i]])
		# Make the contrast matrix C based on the number of gambles in the choice
		C <- make.C(n)
		# Make M and W based on the choice
		states <- make.states(choices[[i]], independent=independent)
		# Inhibition matrix S: Either use predefined matrix or make one based on distances between options
		if(is.null(S))
			S <- make.S(states$M, states$W, alpha, beta, q, r)
		# Set equal thresholds at threshold for each accumulator
		T <- rep(threshold, n)
		# Unbiased starting position
		P <- rep(0, n)
		# Run the simulation. Store the proportion of wins for each option in predicitons[[i]]
		predictions[[i]] <- DFT(C, states$M, states$W, S, T, P, sigma, iterations=iterations, 
			seed=as.integer(runif(1) * .Machine$integer.max), parallel=parallel
		) / iterations
		# Copy the names of each gamble for output
		names(predictions[[i]]) <- names(choices[[i]])

	}
	# Copy the names of each choice for output
	names(predictions) <- names(choices)
	predictions
}
