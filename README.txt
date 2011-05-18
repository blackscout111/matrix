################################################################################
README: (MATRIX)
################################################################################
This is a set of c++ classes and header files that allow for a MATLAB-like
object called "matrix" to be manipulated.  This project was origially designed
as a means of performing image manipulation.

To use this project simply include the "matrix.h" file into your code (be sure
to add it to the working directory of your code as well.  To use the
"matrixmath" programs copy "matrixmath.h" and all of the "matrixmathXX.h"
files into your code's working directory and include "matrixmath." into your
code (This file will include the rest of the matrixmath files.).



################################################################################
TO DO:
################################################################################
[0] Fix variance and standard deviation functions in the matrixmath header
	files.  They currently will overflow if the matrix is too large.
[1] Implement Fourier transform function
[2]	Implement Gaussian Reduction function into matrixmath header files

