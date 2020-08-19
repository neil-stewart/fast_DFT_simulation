# Fast DFT Simulation

This source code implements decision field theory (DFT) as reported in Roe, Busemeyer, and Townsend (2001).

DFT.R is the library of functions used in the model. The code is written in version 2.13.0 of the R programming language.

The library calls optimised C code for speed, and you will need to compile this code before you can use the library. DFT.c is a single-thread version of the model and DFT_parallel.c is a multi-threaded version that will use all of the CPU power available.

## Compiling the C libraries

The makefile will compile DFT.c and DFT_parallel.c to give DFT.so and DFT_parallel.so. You will need the GNU Scientific Library installed. I used version 1.14 of the GNU Scientific Library. I compiled the code using gcc version 4.4.5 from the GNU Compiler Collection.

If you are using Debian or a derivative, installing the packages r-base, libgsl0-dev, gcc, and make should give everything you need. Otherwise, install the corresponding packages for your distribution. If you are using Windows, check out Cygwin, which includes the GNU Compiler Collection and GNU Scientific Library.

Put DFT.c, DFT_parallel.c, and makefile in the same direct and run the command "make". This will create DFT.so and DFT_parallel.so, which are the shared libraries used by the R code and should be in the same directory at DFT.R.

## Big Three

big_three_example.R shows how to use the DFT() function from the R library to simulate the similarity, attraction, and compromise effects. See big_three_example.Rout for output.

## Kahneman and Tverksky (1979) Choices

Kahneman_Tversky_1979_example.R fits DFT to the choices between prospects (gambles) listed in Kahneman and Tversky (1979). Kahneman_Tversky_1979_data.R contains the choices. The function DFT.choices() simulates a series of risky choices. DFT is applied to risky prospects using a states-of-the-world representation. For example, a choice between (A) a 25% chance of 6,000 otherwise nothing and (B) a 25% chance of 4,000 otherwise a 25% chance of 2,000 otherwise nothing is represented by

```{r}
$M
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,] 6000    0 6000    0 6000    0
[2,] 4000 4000 2000 2000    0    0

$W
[1] 0.0625 0.1875 0.0625 0.1875 0.1250 0.3750
```

Each column of M represents the outcomes in different states of the world (e.g., win A and win B, lose A and win B, ...), and each element of W gives the probability of that state. The function make.states() creates an option attribute matrix M and a attribute weighting vector W from a risky choice, and is used by DFT.choices().

A second issue is how to make the inhibition matrix S from the choice attributes. The DFT.choices() function gives two choices. Either you supply a pre-defined S matrix to apply to all choices, or you supply the parameters for an algorithm that makes a matrix for each choice based on the distance between the gambles in outcome state space, with each dimension weighted by the probability of that state. The Kahneman_Tversky_1979_example.R code contains examples of each approach.

## Running the Scripts

You can run the R scripts by starting the R program and copying and pasting in the commands one-by-one or in chunks from the scripts big_three_example.R or Kahneman_Tversky_1979_example.R.

Whole .R scripts can be run from the command line in batch mode with "R CMD BATCH filename.R". The .Rout files show the result of running the .R files.

Please send bugs or comments to Neil.Stewart@warwick.ac.uk.


