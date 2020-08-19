// Version 1.0.0
// Copyright 2011 Neil Stewart
// The program is distributed under the terms of the GNU General Public License
//
// DFT.c is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// DFT.c is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// <http://www.gnu.org/licenses/>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <R.h>

/* Compile time optimization */
/* Generator algorithm for Gaussian random variates */
//#define GSL_RAN_GAUSSIAN gsl_ran_gaussian
#define GSL_RAN_GAUSSIAN gsl_ran_gaussian_ziggurat	// 64% faster than gsl_ran_gaussian
/* Use uniformly distributed noise for speed */
//#define UNIFORM_APPROXIMATION // 25% faster than ziggurat
/* Algorithm for random number generator */
#define RANDOM_NUMBER_GENERATOR gsl_rng_taus2

gsl_rng *g;

inline int sample_dimension(double *w)
{
/*
  sample_dimension(w) returns an integer in the range 0 - n-1 with
  probabilities p0, p1, ... p_n-1. w is a 1D array of cumulative
  probabilities (p0, p0+p1, ..., 1.0). The last element in w must be 1.0,
  otherwise the function may never return.
*/
	int i = 0;
	double r;

	r = gsl_rng_uniform(g);

	while (1) {
		if (r < w[i])
			return i;
		i++;
	}
}

inline void gaussian_random_vector(gsl_vector * X, int n, double *sigma)
{
/*
  Fill X with Gaussian random variates mean 0 standard deviation *sigma. If
  UNIFORM_APPROXIMATION is defined, use uniform deviates with same standard
  deviation
*/
	int i;

	for (i = 0; i < n; i++) {
#ifdef UNIFORM_APPROXIMATION
		gsl_vector_set(X, i,
		  (gsl_rng_uniform(g) * 2 * sqrt(3) - sqrt(3)) * (*sigma));
#else
		gsl_vector_set(X, i, GSL_RAN_GAUSSIAN(g, *sigma));
#endif
	}
}

void
simulate_DFT(double *Cin, double *Min, double *Win, double *Sin,
  double *Tin, double *Pstartin, double *sigma, int *nin,
  int *din, int *iter, int *seed, int *counter)
/*
    Cin, the contrast matrix, returned unchanged
    Min, the attribute value matrix, returned unchanged
    Win, the dimension weights, must be cumulative, returned unchanged
    Sin, the inhibition matrix, returned unchanged
    Tin, the threshold vector, normally all equal, returned unchanged
    Pstartin, the initial preference vector, normally all zero for no bias, returned unchanged
    sigma, the standard deviation of the normally distributed noise added each step, returned unchanged
    nin, the number of options, returned unchanged
    din, the number of dimensions, returned unchanged
    iter, the number of interations, returned unchanged
    seed, the seed for the random number generator, returned unchanged
    counter, should be all zero to start, returns a vector giving the number of times each option is selected, should sum to inter
*/
{
	int n = *nin;
	int d = *din;
	int i;
	const gsl_rng_type *G;
	gsl_vector **V;
	gsl_matrix_view C = gsl_matrix_view_array(Cin, n, n);
	gsl_matrix_view M = gsl_matrix_view_array(Min, n, d);
	gsl_matrix_view S = gsl_matrix_view_array(Sin, n, n);
	gsl_vector_view T = gsl_vector_view_array(Tin, n);
	gsl_vector_view Pstart = gsl_vector_view_array(Pstartin, n);
	gsl_vector *P = gsl_vector_alloc(n);
	gsl_vector *SP = gsl_vector_alloc(n);
	gsl_vector *diff = gsl_vector_alloc(n);
	gsl_matrix *CM = gsl_matrix_alloc(n, d);
	gsl_vector *epsilon = gsl_vector_alloc(n);

	// Set up random number generator
	gsl_rng_env_setup();
	G = RANDOM_NUMBER_GENERATOR;	/* Override any command line option */
	gsl_rng_default_seed = *seed;	/* (e.g., GSL_RNG_TYPE="rand" GSL_RNG_SEED=0) */
	g = gsl_rng_alloc(RANDOM_NUMBER_GENERATOR);

	// Precompute V for attention to each possible dimension
	// V <- CM * W 
	// So V[0] is V if first dimension is attended, V[1] for second dimension, and so on
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &C.matrix, &M.matrix, 0.0, CM);	// CM <- C * M
	V = (gsl_vector **) malloc(sizeof(gsl_vector *) * d);
	for (i = 0; i < d; i++) {
		gsl_vector *Wbasis = gsl_vector_alloc(d);
		gsl_vector_set_basis(Wbasis, i);
		V[i] = gsl_vector_alloc(n);
		gsl_blas_dgemv(CblasNoTrans, 1.0, CM, Wbasis, 0.0, V[i]);
		gsl_vector_free(Wbasis);
	}
	gsl_matrix_free(CM);

	{
		for (i = 0; i < *iter; i++) {
			// Repeat the accumulation iter times
			gsl_vector_memcpy(P, &Pstart.vector);
			do {
				// Race preference state P to threshold T

				// P <- SP + V + C * noise
				gsl_blas_dgemv(CblasNoTrans, 1.0, &S.matrix, P, 0.0, SP);	// SP <- S * P
				gsl_vector_memcpy(P, SP);	// Put SP in P
				gsl_vector_add(P, V[sample_dimension(Win)]);	// Add a _randomly_selected_ V to P
				gaussian_random_vector(epsilon, n, sigma);	// Make random vector epsilon
				gsl_blas_dgemv(CblasNoTrans, 1.0, &C.matrix, epsilon, 1.0, P);	// P <- P + C * epsilon

				// Test to see if threshold is reached
				gsl_vector_memcpy(diff, &T.vector);
				gsl_vector_sub(diff, P);
			}
			while (gsl_vector_isnonneg(diff));
			counter[gsl_vector_min_index(diff)]++;	// Get the index of the element in P which exceeds threshold
		}
	}

	// Free up memory
	for (i = 0; i < d; i++)
		gsl_vector_free(V[i]);
	free(V);
	gsl_vector_free(SP);
	gsl_vector_free(P);
	gsl_vector_free(diff);
	gsl_vector_free(epsilon);
}

int main(void)
/* Useful for debugging, not used by the R call */
{
	double C[] = { 1.0, -0.5, -0.5, -0.5, 1.0, -0.5, -0.5, -0.5, 1.0 };
	double W[] = { 0.5, 1.0 };	// Use cumulative form
	double T[] = { 17.5, 17.5, 17.5 };
	double Pstart[] = { 0.0, 0.0, 0.0 };
	double sigma = 1.0;
	int n = 3;
	int d = 2;
	int iter = 1000000;
	int seed = 1;
	int counter[] = { 0, 0, 0 };
	int i;
	// Attraction effect
	double M[] = { 1.0, 3.0, 3.0, 1.0, 0.5, 2.5 };
	double S[] = {
		0.95, -1.2232e-02, -2.2647e-02, -1.2232e-02, 0.95, -6.7034e-04,
		-2.2647e-02, -6.7034e-04, 0.95
	};

	simulate_DFT(C, M, W, S, T, Pstart, &sigma, &n, &d, &iter, &seed, counter);

	printf("counter[] = { ");
	for (i = 0; i < n; i++)
		printf("%d ", counter[i]);
	printf("}\n");

	return 0;
}
