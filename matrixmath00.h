#ifndef matrixmath00_h
#define matrixmath00_h

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "matrix.h"

using namespace std;

/*******************************************************************************
* This file normally included from matrixmath.h
*
* Section 00:	Check to see if a matrix has appropriate type of elements
*******************************************************************************/
//______________________________________________________________________________
// Checks to see if a matrix's entries are of an appropriate type
// (e.g. int, double, float, char...)
template <class type>
bool isMathType(matrix<type> mat)
{
	// Its a type that can have math done with it
	if (
			(typeid(mat) == typeid(matrix<long>)) ||
			(typeid(mat) == typeid(matrix<long long>)) ||

			(typeid(mat) == typeid(matrix<short>)) ||

			(typeid(mat) == typeid(matrix<int>)) ||
			(typeid(mat) == typeid(matrix<short int>)) ||
			(typeid(mat) == typeid(matrix<unsigned short int>)) ||
			(typeid(mat) == typeid(matrix<long int>)) ||
			(typeid(mat) == typeid(matrix<long int>)) ||
			(typeid(mat) == typeid(matrix<long long int>)) ||
			(typeid(mat) == typeid(matrix<unsigned int>)) ||
			(typeid(mat) == typeid(matrix<unsigned long int>)) ||
			(typeid(mat) == typeid(matrix<unsigned long long int>)) ||

			(typeid(mat) == typeid(matrix<char>)) ||
			(typeid(mat) == typeid(matrix<unsigned char>)) ||

			(typeid(mat) == typeid(matrix<double>)) ||
			(typeid(mat) == typeid(matrix<long double>)) ||

			(typeid(mat) == typeid(matrix<float>))
	   )
	{
		return true;
	}
	else // Not a type that can have math done with it.
	{
		return false;
	}
}



//______________________________________________________________________________
// Terminates the program and displays error message if a matrix's entries are
// not of an appropriate type.
// (Determined by isMathType(...))
template <class type>
void continueIfMathType(matrix<type> mat)
{
	// If the matrix can have math done with it
	if (isMathType(mat))
	{
		// Do nothing
	}
	else
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'template <class type> void continueIfMathType(matrix<type>"
				<< " mat)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "This type of matrix can't have math done with it!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}



#endif
