# This is the simplest code for applying the DFT model using the DFT.R
# library and the DFT.so and DFT_parallel.so C libraries to the similarity,
# attraction, and compromise effects described in Roe, Busemeyer, and
# Townsend (2001). Results should match the matlab code at
# http://mypage.iu.edu/~jbusemey/Dec_lect/Dec_prg/sim_mdf.txt (in
# particular, the first, fourth, and sixth problems). The single-thread
# version of the code is about 300 times faster than the matlab version. The
# multi-threaded version is, obviously, faster still.

# R CMD BATCH big_three_example.R should produce example output

rm(list=ls())

# Load the DFT library
# Note, you will need the DFT.so and DFT_parallel.so libraries in the same directory
source("DFT.R")


iterations <- 100000
(C <- matrix(c(1.0, -0.5, -0.5, 
              -0.5,  1.0, -0.5, 
              -0.5, -0.5,  1.0), nrow=3, byrow=TRUE))
(W <- c(.5, .5))
(T <- rep(17.5, 3))
(P <- rep(0.0, 3))
(sigma <- 1.0)

# Similarity effect
(M <- matrix(c(1.0, 3.0, 
               3.0, 1.0, 
               0.9, 3.1), nrow=3, byrow=T))
(S <- matrix(c(0.9500000, -0.0122316, -0.04999996, 
              -0.0122316,  0.9500000, -0.0090303, 
              -0.0499996, -0.0090303, 0.9500000),
    nrow=3, byrow=TRUE))

result <- DFT(C, M, W, S, T, P, sigma, iterations=iterations, seed=1, parallel=TRUE)
result/iterations



# Compromise effect
(M <- matrix(c(1.0, 3.0, 
               3.0, 1.0, 
               2.0, 2.0), nrow=3, byrow=T))
(S <- matrix(c(0.950000, -0.012232, -0.045788, 
              -0.012232,  0.950000, -0.045788, 
              -0.045788, -0.045788,  0.950000),
    nrow=3, byrow=TRUE))

result <- DFT(C, M, W, S, T, P, sigma, iterations=iterations, seed=1, parallel=TRUE)
result/iterations



# Attraction effect
(M <- matrix(c(1.0, 3.0, 
               3.0, 1.0, 
               0.5, 2.5), nrow=3, byrow=T))
(S <- matrix(c(0.95, -1.2232e-02, -2.2647e-02, 
              -1.2232e-02, 0.95, -6.7034e-04,
              -2.2647e-02, -6.7034e-04, 0.95),
    nrow=3, byrow=TRUE))

result <- DFT(C, M, W, S, T, P, sigma, iterations=iterations, seed=1, parallel=TRUE)
result/iterations

