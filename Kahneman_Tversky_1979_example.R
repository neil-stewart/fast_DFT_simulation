# This is the simplest code for applying the DFT model using the DFT.R
# library and the DFT.so and DFT_parallel.so C libraries to the Kahneman and
# Tversky (1979) data.
#
# R CMD BATCH Kahneman_Tversky_1979_example.R should produce example output

rm(list=ls())

# Load the DFT library
# Note, you will need the DFT.so and DFT_parallel.so libraries in the same directory
source("DFT.R")

# Load the Kahneman and Tversky (1979) choices
source("Kahneman_Tversky_1979_data.R")
# See the choices
choices
# See the choice proportion data
proportion.A.chosen.data

predictions <- DFT.choices(choices=choices, independent=T, alpha=.1,
	beta=.1, q=2, r=2, sigma=.1, threshold=10, iterations=10000, parallel=TRUE)
predictions

# Run the DFT model on the choice, using a distance-dependent inhibition matrix for each choice
predictions.parallel <- DFT.choices(choices=choices, independent=T, alpha=.1, beta=.1, q=2, r=2, sigma=.1, threshold=10, iterations=100000, parallel=TRUE)
predictions.serial<- DFT.choices(choices=choices, independent=T, alpha=.1, beta=.1, q=2, r=2, sigma=.1, threshold=10, iterations=100000, parallel=FALSE)
predictions

# Extract the probability of choosing the first (A) option for each choice
probability.A.chosen.predictions.parallel <- sapply(predictions.parallel,function(z) {z[[1]]})
probability.A.chosen.predictions.serial <- sapply(predictions.serial,function(z) {z[[1]]})
probability.B.chosen.predictions.parallel <- sapply(predictions.parallel,function(z) {z[[2]]})
probability.B.chosen.predictions.serial <- sapply(predictions.serial,function(z) {z[[2]]})


( probability.total.chosen.predictions.parallel <- sapply(predictions.parallel,function(z) {z[[1]]+z[[2]]}) )
( probability.total.chosen.predictions.serial <- sapply(predictions.serial,function(z) {z[[1]]+z[[2]]}) )


# Plot model predictions against choice proportions
plot(proportion.A.chosen.data, probability.A.chosen.predictions, 
	xlab="Proportion of A Choices (Data)", ylab="Probability of an A Choice (Model)",
	xlim=c(0,1), ylim=c(0,1)
)
# The fit is not very good, but then the parameters were chosen quite arbitrarily.

# Or you can run the DFT model using a single predefined S matrix for all choices
S <- matrix(c(0.95, -0.05, -0.05, -0.05, 0.95, -0.05, -0.05, -0.05, 0.95), nrow=3, byrow=T)
predictions <- DFT.choices(choices=choices, independent=T, S=S,
	sigma=.1, threshold=10, iterations=10000, parallel=TRUE)
predictions
