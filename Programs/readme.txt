List of programs and files:

besselzero.m: generates a desired number of zeros of a Bessel function of desired order. The function is taken from here:
Jason Nicholson (2022). Bessel Zero Solver (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver), MATLAB Central File Exchange. Retrieved February 3, 2022. 

importdiskfunctions.m: For input N and reso, checks if the necessary files containing the base functions and the coefficients for the initial conditions are already present. If they are not, they are generated in the Data folder by calling the functions CreateBesselArray.m and BumpCoefficients.m

DiskWavenumbers2000.mat: Contains the lowest 2000 dirichlet wavenumbers for the unit disk.

simps.m: An implementation of Simpson's rule for numerical integration. The code is taken from here:
Damien Garcia (2022). Simpson's rule for numerical integration (https://www.mathworks.com/matlabcentral/fileexchange/25754-simpson-s-rule-for-numerical-integration), MATLAB Central File Exchange. Retrieved February 3, 2022. 

polygonize.m: Accepts an array and a contour (an array with two rows for x and y) as inputs and returns the same array with all entries set to NaN which are outside the contour. Used on arrays which need to be plotted. In practice, it is faster to call this function only once for a given shape and create a "NaNarray", and henceforth use that array for comparison instead of calling polygonize.m repeatedly.
