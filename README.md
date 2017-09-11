# Numerical analysis algorithms in C++

This repository contains some basic numerical analysis algorithms implemented in C++. More specifically, the following techniques were implemented:

-   derivative of a single-variable function;
-   partial derivatives of two-variable function;
-   Newton-Raphson method for finding roots of single-variable functions;
-   Gradient descent method for finding (local) minima of functions;
-   Numerical integration using the Newton-Cotes formulae (rectangle, trapezoidal and Simpson's functions);
-   Adaptive quadrature, implemented according to Numerical Recipes 3rd edition;
-   Monte Carlo integration for single variable functions and for the approximation of the volume and center of mass of a tridimensional region.

Numerical integration methods have built-in support for OpenMP. Compile the package with the `-fopenmp` flag in order to use it.