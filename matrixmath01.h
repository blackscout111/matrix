#ifndef matrixmath01_h
#define matrixmath01_h

#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "matrix.h"
#include "matrixmath00.h"   // Needed for checking the type of matrix

using namespace std;

namespace matrixmath
{

/*******************************************************************************
* This file normally included from matrixmat.h
*
* Section 01:	Operator overloads (+,=+,-,-=,*)
*******************************************************************************/
//______________________________________________________________________________
// Overloads the += opperator
template <class type>
matrix<type>& operator += (matrix<type>& a, const matrix<type>& b)
{
	// Checks to see if the matrices can do math together
	// (Ends the program if either of them can't)
	continueIfMathType(a);
	continueIfMathType(b);


	// The two matricies must have the same dimensions
	if ((a.width() == b.width()) &&
		(a.height() == b.height()))
	{

		// Fills the a(i,j) with a(i,j) + b(i,j)
		for (_DIM i = 0; i < a.height(); i++)
		{
			for (_DIM j = 0; j < a.width(); j++)
			{
				a(i,j) += b(i,j);
			}
		}

		return a;
	}
	else
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'matrix<type>"
				<< "matrix<type>::operator"
				<< " += (matrix<type>& a, matrix<type>& b)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "Dimension mis-match!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}



//______________________________________________________________________________
// Overloads the + opperator
template <class type>
const matrix<type> operator + (const matrix<type>& a, const matrix<type>& b)
{
	matrix<type> result(a);
	
	result += b;

	return result;
}



//______________________________________________________________________________
// Overloads the - opperator
template <class type>
matrix<type>& operator -= (matrix<type>& a, const matrix<type>& b)
{
	// Checks to see if the matrices can do math together
	// (Ends the program if either of them can't)
	continueIfMathType(a);
	continueIfMathType(b);

	// The two matricies must have the same dimensions
	if ((a.width() == b.width()) &&
		(a.height() == b.height()))
	{
		// Fills the a(i,j) with a(i,j) - b(i,j)
		for (_DIM i = 0; i < a.height(); i++)
		{
			for (_DIM j = 0; j < a.width(); j++)
			{
				a(i,j) -= b(i,j);
			}
		}

		return a;
	}
	else
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'matrix<type>"
				<< "matrix<type>::operator"
				<< " -= (matrix<type>& a, matrix<type>& b)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "Dimension mis-match!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}



//______________________________________________________________________________
// Overloads the - opperator
template <class type>
const matrix<type> operator - (const matrix<type>& a, const matrix<type>& b)
{
	matrix<type> result(a);

	result -= b;

	return result;
}



//______________________________________________________________________________
// Overloads the *= opperator when a matrix is multiplied by a number
template <class type>
matrix<type>& operator *= (matrix<type>& a, const type& b)
{
	for (_DIM i = 0; i < a.height(); i++)
	{
		for (_DIM j = 0; j < a.width(); j++)
		{
			a(i,j) *= b;
		}
	}

	return a;
}



//______________________________________________________________________________
// Overloads the *= opperator when a number is multiplied by a matrix
template <class type>
matrix<type>& operator *= (const type& a, matrix<type>& b)
{
	for (_DIM i = 0; i < b.height(); i++)
	{
		for (_DIM j = 0; j < b.width(); j++)
		{
			b(i,j) *= a;
		}
	}

	return b;
}



//______________________________________________________________________________
// Overloads the * opperator when a matrix is multiplied by a number
template <class type>
const matrix<type> operator * (const matrix<type>& a, const type& b)
{
	matrix<type> result(a);

	result *= b;

	return result;
}



//______________________________________________________________________________
// Overloads the * opperator when a matrix is multiplied by a number
template <class type>
const matrix<type> operator * (const type& a, const matrix<type>& b)
{
	matrix<type> result(b);

	result *= a;

	return result;
}



//______________________________________________________________________________
// Overloads the * opperator when a matrix is multiplied by another matrix
template <class type>
const matrix<type> operator * (const matrix<type>& a, const matrix<type>& b)
{
	// Checks to see if the matrices can do math together
	// (Ends the program if either of them can't)
	continueIfMathType(a);
	continueIfMathType(b);

	// The width of a must equal the height of b
	if (a.width() == b.height())
	{
		// Creates a matrix that is the same height as a and the same width as b
		matrix<type> innerProduct(a.height(), b.width());

		// Fills the innerProduct with a(i,j) - b(i,j)
		for (_DIM i = 0; i < innerProduct.height(); i++)
		{
			for (_DIM j = 0; j < innerProduct.width(); j++)
			{
				// Temporary summing variable
				type sum = 0;

				// Dot row i of a with column j of b
				for (_DIM k = 0; k < a.width(); k++)
				{
					sum += a(i,k) * b(k,j);
				}

				innerProduct(i,j) = sum;
			}
		}

		return innerProduct;
	}
	else
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'matrix<type>"
				<< "matrix<type>::operator"
				<< " * (matrix<type>& a, matrix<type>& b)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "Dimension mis-match!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}



//______________________________________________________________________________
// Overloads the *= opperator when a matrix is multiplied by another matrix
template <class type>
matrix<type>& operator *= (matrix<type>& a, const matrix<type>& b)
{
	matrix<type> result(a * b);

	a = result;

	return a;
}

};// namespace matrixmath

#endif
