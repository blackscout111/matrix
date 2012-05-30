################################################################################
README: (MATRIX)
################################################################################
This is a set of c++ classes and header files that allow for a MATLAB-like
object called "matrix" to be manipulated.  This project was originally designed
as a means of performing image manipulation.  This project is not intended for
high performance computation.  It was made so that I could perform some linear
algebra techniques and matrix manipulations without having to install a more
advanced linear algebra library such as boost or LAPACK.

To use this project simply include the "matrix.h" file into your code (be sure
to add it to the working directory of your code as well.  To use the
"matrixmath" programs copy "matrixmath.h" and all of the "matrixmathXX.h"
files into your code's working directory and include "matrixmath.h" into your
code. (This file will include the rest of the matrixmath files.) The
"matrixmath" functions are inside of the "matrixmath" namespace. The "matrix"
object is not inside of any namespace.

The "test.cpp" file in this directory is only for testing purposes.  It is used
to test the project code and its functions as it is being developed.

Everyone is free to use this code!  Please just reference me in your code as
"Bryan A. Clifford". Feel free to email at "blackscout111@gmail.com" with
constructive feedback/comments.


