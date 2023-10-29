# Numerical_integration

Code related to the thesis "On galaxy rotation curves and alternative gravity theories" for Bachelor Degree in Physics, October 2023, at Università degli Studi di Milano,
by Alessandro Musmeci.

This C++ code calculates the rotation curves for a separable disk with exponential surface density and different thickness profiles in newtonian physics. Only adimensional formulas are calculated. 
The code produces a machine-readable file with 3 columns: the reduced distance $\tilde{r}$, the rotation velocity $\tilde{v}$, and the error on $\tilde{v}^2$.

## Code structure

The random number generator is from the Numerical Simulation Laboratory at Dipartimento di Fisica, Unimi, by Prof. D. Galli. Correlated files are random.cpp, random.h, Primes, seed.in.

File functions.h contains two functions: one initializes the random number generator, the other calculates the error on the integrals.

The file classes.h contains 4 classes: The class f and distro contain the integrands, factorized in the product of a Bessel function and a probability distribution. In the class Metropolis a standard Metropolis algorithm is implemented for the sampling of the probability distributions. The parameters have been fixed in order to have an acceptance probability of about 50%. In the Blockaverage data blocking for the calculation of the integrals is implemented.

The main file is integral.cpp. It depends on the input.in file, where the following parameters have to be specified: the thickness profile (0 for slab, 1 for exponential, 2 for sech2 profile), the scale length of the disk (examples from 3 galaxies are already in the file) and the distance step h. The code calculates an integral for each value of $\tilde{r}$, at steps of h, until $\tilde{r} = 15$.

Lastly, for the makefile, in the CFLAGS variable, one has to specify the location of the boost library for the use of Bessel functions on the machine in the following way:
-I <path_to_library>/boost/<version>/include

