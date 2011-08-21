#ifndef matrixmath02_h
#define matrixmath02_h

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <typeinfo>

#include "matrix.h"
#include "matrixmath00.h" // Needed for checking the type of matrix

using namespace std;

namespace matrixmath
{

/*******************************************************************************
* This file normally included from matrixmath.h
*
* Section 02:	Special matrix functions (e.g. transpose, variance,...)
*******************************************************************************/
//______________________________________________________________________________
// Returns the transpose of a matrix
template <class type>
const matrix<type> trans(const matrix<type>& mat)
{
	matrix<type> trans(mat.width(), mat.height());

	for (size_t i = 0; i < trans.height(); i++)
	{
		for (size_t j = 0; j < trans.width(); j++)
		{
			trans(i,j) = mat(j,i);
		}
	}

	return trans;
}


//______________________________________________________________________________
// Returns the trace of a matrix
template <class type>
const type trace(const matrix<type>& mat)
{
	// Check to see if the matrix contains elements of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// The eventual trace of the matrix
	type trace = 0;

	// Check to see if the matrix is square
	if (mat.height() != mat.width())
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'template <class type> const type trace(const matrix<type>&"
				<< " mat)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "Can only take the trace of a square matrix!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}

	// Adds up the values along the main diagonal
	for (size_t i = 0; i < mat.height(); i++)
	{
		trace += mat(i,i);
	}

	return trace;
}


//______________________________________________________________________________
// Returns a column matrix with the same height as 'mat' and with elements
// corresponding to the maximum values of the elements in the rows of 'mat'
template <class type>
const matrix<type> max(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	matrix<type> max(mat.height(), 1);

	// The current maximum value
	type maxVal;

	for (size_t i = 0; i < mat.height(); i++)
	{
		maxVal = mat(i,0);

		for (size_t j = 1; j < mat.width(); j++)
		{
			// A new maximum has been found
			if (maxVal < mat(i,j))
			{
				maxVal = mat(i,j);
			}
		}

		max(i,0) = maxVal;
	}

	return max;
}


//______________________________________________________________________________
// Returns the maximum value of all the elements in 'mat'
template <class type>
const type max2d(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// The current maximum value
	type maxVal;

	for (size_t i = 0; i < mat.height(); i++)
	{
		for (size_t j = 0; j < mat.width(); j++)
		{
			// The first value checked
			if ((i == 0) && (j == 0))
			{
				maxVal = mat(i,j);
			}

			// A new maximum has been found
			else if (maxVal < mat(i,j))
			{
				maxVal = mat(i,j);
			}
		}
	}

	return maxVal;
}


//______________________________________________________________________________
// Returns a column matrix with the same height as 'mat' and with elements
// corresponding to the minimum values of the elements in the rows of 'mat'
template <class type>
const matrix<type> min(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	matrix<type> min(mat.height(), 1);

	// The current minimum value
	type minVal;

	for (size_t i = 0; i < mat.height(); i++)
	{
		minVal = mat(i,0);

		for (size_t j = 0; j < mat.width(); j++)
		{
			// A new minimum has been found
			if (minVal > mat(i,j))
			{
				minVal = mat(i,j);
			}
		}

		min(i,0) = minVal;
	}

	return min;
}


//______________________________________________________________________________
// Returns the minimum value of all the elements in 'mat'
template <class type>
const type min2d(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// The current minimum value
	type minVal;

	for (size_t i = 0; i < mat.height(); i++)
	{
		for (size_t j = 0; j < mat.width(); j++)
		{
			// The first value checked
			if ((i == 0) && (j == 0))
			{
				minVal = mat(i,j);
			}

			// A new minimum has been found
			else if (minVal > mat(i,j))
			{
				minVal = mat(i,j);
			}
		}
	}

	return minVal;
}


//______________________________________________________________________________
// Returns a column matrix with the same height as 'mat' and with elements
// corresponding to the sum of the values of the elements in the rows of 'mat'
template <class type>
const matrix<type> sum(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	matrix<type> sum(mat.height(), 1);

	// The sum of all of the elements in a row
	type rowSum = 0;
	for (size_t i = 0; i < mat.height(); i++)
	{

		for (size_t j = 0; j < mat.width(); j++)
		{
			// Add up the sum
			rowSum += mat(i,j);
		}

		sum(i,0) = rowSum;
		rowSum = 0;
	}

	return sum;
}


//______________________________________________________________________________
// Returns the average value of all the elements in 'mat'
template <class type>
const type sum2d(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// The sum of all of the elements in the matrix
	type sum = 0;
	for (size_t i = 0; i < mat.height(); i++)
	{
		for (size_t j = 0; j < mat.width(); j++)
		{
			// Add up the sum
			sum += mat(i,j);
		}
	}

	return sum;
}


//______________________________________________________________________________
// Returns a column matrix with the same height as 'mat' and with elements
// corresponding to the average values of the elements in the rows of 'mat'
template <class type>
const matrix<double> mean(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	matrix<double> mean(mat.height(), 1);

	double inv_n = 1/(double)mat.width();

	mean.fill(0);
	for (size_t i = 0; i < mean.height(); ++i)
	{
		for (size_t j = 0; j < mat.width(); ++j)
		{
			mean(i,0) += (double)mat(i,j);
		}
		mean(i,0) *= inv_n;
	}

	return mean;
}


//______________________________________________________________________________
// Returns the average value of all the elements in 'mat'
template <class type>
const double mean2d(const matrix<type>& mat)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	return ((double)sum2d(mat))/((double)mat.height()*(double)mat.width());
}


//______________________________________________________________________________
// Returns a column matrix with the same height as 'mat' and with elements
// corresponding to the variance of the values of the elements in the
// rows of 'mat'.
// *	If 'isPop' = true then the population variance is calculated
// *	If 'isPop' = false (default) then the sample variance is calculated
template <class type>
const matrix<double> vari(const matrix<type>& mat, bool isPop= false)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	matrix<double>	variance(mat.height(), 1);

	double	n = mat.width(),
			N = n - 1,
			invn = 1/n,
			invN = 0,
			sum = 0,	// sum of row elements
			sum2 = 0;	// sum of square of row elements

	if (isPop) N += 1;
	if (N == 0)
	{
		// Display error message and terminates the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'template <class type>const double vari(const matrix<type>&"
				<< " mat, bool isPop= false)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "cannot compute variance for a sample size of 1!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
	invN = 1/N;

	for (size_t row = 0; row <  mat.height(); ++row)
	{
		for (size_t col = 0; col < mat.width(); ++col)
		{
			sum += mat(row,col);
			sum2 += mat(row,col)*mat(row,col);

#if _MTXSAFE
			if (sum2 > 0.999*DBL_MAX)
				// Display overflow warning
				cerr	<< endl
						<< "matrixmath.h: In function "
						<< "'template <class type>const double vari(const "
						<< "matrix<type>& mat, bool isPop= false)'"
						<< endl
						<< "matrixmath.h: warning: "
						<< "possible overflow"
						<< endl << endl;
#endif
		}
		variance(row,0) = invN*sum2 - invN*invn*sum*sum;
		sum = 0;
		sum2 = 0;
	}

	return variance;
}


//______________________________________________________________________________
// Returns the population variance of the values of all the elements in 'mat'
// *	If 'isPop' = true then the population variance is calculated
// *	If 'isPop' = false (default) then the sample variance is calculated
template <class type>
const double vari2d(const matrix<type>& mat, bool isPop= false)
{
	// Checks to see that the elements are of a type that can do math
		// (Exits the program if they can't)
		continueIfMathType(mat);

		double	n = mat.width()*mat.height(),
				N = n - 1,
				invn = 1/n,
				invN = 0,
				sum = 0,	// sum of row elements
				sum2 = 0;	// sum of square of row elements

		if (isPop) N += 1;
		if (N == 0)
		{
			// Display error message and terminates the program
			cerr	<< endl
					<< "matrixmath.h: In function "
					<< "'template <class type>const double vari2d(const "
					<< "matrix<type>& mat, bool isPop= false)'"
					<< endl
					<< "matrixmath.h: error: "
					<< "cannot compute variance for a sample size of 1!"
					<< endl << endl;

			// Exits program
			std::exit(EXIT_FAILURE);
		}
		invN = 1/N;

		for (size_t row = 0; row <  mat.height(); ++row)
		{
			for (size_t col = 0; col < mat.width(); ++col)
			{
				sum += mat(row,col);
				sum2 += mat(row,col)*mat(row,col);

#if _MTXSAFE
				if (sum2 > 0.999*DBL_MAX)
					// Display overflow warning
					cerr	<< endl
							<< "matrixmath.h: In function "
							<< "'template <class type>const double "
							<< "vari2d(const matrix<type>& mat, "
							<< "bool isPop= false)'"
							<< endl
							<< "matrixmath.h: warning: "
							<< "possible overflow"
							<< endl << endl;
#endif
			}
		}

		return (invN*sum2 - invN*invn*sum*sum);
}


//______________________________________________________________________________
// Returns a column matrix with the same height as 'mat' and with elements
// corresponding to the standard deviation of the values of the elements in the
// rows of 'mat'.
// *	If 'isPop' = true then the population variance is calculated
// *	If 'isPop' = false (default) then the sample variance is calculated
template <class type>
const matrix<double> stdev(const matrix<type>& mat, bool isPop= false)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// Holds the standard deviation values
	matrix<double>	stdeviation;

	// Calculate standard deviations (the square root of the variance)
	stdeviation = vari(mat, isPop);
	for (size_t i = 0; i < mat.height(); i++)
	{
		stdeviation(i,0) = std::sqrt(stdeviation(i,0));
	}

	return stdeviation;
}


//______________________________________________________________________________
// Returns the population standard deviation of the values of all the elements
// in 'mat'
// *	If 'bias' = 0 (default) then the population standard deviation is
//		calculated
// *	If 'bias' = 1 then the sample standard deviation is calculated
template <class type>
const double stdev2d(const matrix<type>& mat, bool isPop= false)
{
	// Checks to see that the elements are of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// Calculate standard deviation (the square root of the variance)
	return std::sqrt(vari2d(mat, isPop));
}


//______________________________________________________________________________
// Returns the determinant of the matrix 'mat'
// * Uses method of cofactors to calulate determinent
// * If dim = 1 then will calculate cofactors along row 'idx'
// * If dim = 2 then will calculate cofactors along column 'idx'
template <class type>
const double det(const matrix<type>& mat, unsigned char dim= 1, size_t idx= 0)
{
	// Check to see if the matrix contains elements of a type that can do math
	// (Exits the program if they can't)
	continueIfMathType(mat);

	// Check to see if the matrix is square
	if (mat.height() != mat.width())
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'template <class type>const double det(const matrix<type>&"
				<< " mat, unsigned char dim= 1, size_t idx= 0)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "Can only take the determinant of a square matrix!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}

	// Value of the determinant 'sum' and the sign of of the minors
	double detSum = 0;

	// Cannot take the determinant of a 1x1 matrix
	if (mat.height() == 1)
	{
		// Display error message and terminate the program
		cerr	<< endl
				<< "matrixmath.h: In function "
				<< "'template <class type>const double det(const matrix<type>&"
				<< " mat, unsigned char dim= 1, size_t idx= 0)'"
				<< endl
				<< "matrixmath.h: error: "
				<< "Cannot take the determinant of a 1x1 matrix!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
	// If 2x2 matrix then take the determinant
	else if (mat.height() == 2)
	{
		detSum = mat(0,0)*mat(1,1) - mat(0,1)*mat(1,0);
	}

	// Use method of cofactors
	else
	{
		switch(dim)
		{
			////////////////////////////////////////////////////////////////////
			// Calculate cofactors along the row 'idx'
			case 1:
				  for(size_t curCol = 0; curCol < mat.width(); curCol++)
				  {
					  if (mat(idx,curCol) != 0)
					  {
						  // Cofactor matrix
						  matrix <double> cofac(mat.height() - 1, mat.width() -1);
						  size_t idxCnt = 0;

						  // Fill Cofactor matrix
						  for(size_t row = 0; row < mat.height(); row++)
						  {
							  if (row != idx)
							  {
								  for(size_t col = 0; col < mat.width(); col++)
								  {
									if (col != curCol)
									{
										cofac(idxCnt/cofac.width(),
											  idxCnt%cofac.width())
											= mat(row,col);
										idxCnt++;
									}
								  }
							  }
						  }

						  // The determinant of 'mat' is the sum of its signed
						  // minors.
						  detSum += (1 - 2*(int)(curCol%2))*(1 - 2*(int)(idx%2))
							  		*mat(idx,curCol)*det(cofac);
					  }
				  }
				break;
			////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////
			// Calculate cofactors along the column 'idx'
			case 2:
				  for(size_t curRow = 0; curRow < mat.height(); curRow++)
				  {
					  if( mat(curRow,idx) != 0)
					  {
						  // Cofactor matrix
						  matrix <double> cofac(mat.height() - 1, mat.width() -1);
						  size_t idxCnt = 0;

						  // Fill Cofactor matrix
						  for(size_t row = 0; row < mat.height(); row++)
						  {
							  if (row != curRow)
							  {
								  for(size_t col = 0; col < mat.width(); col++)
								  {
									if (col != idx)
									{
										cofac(idxCnt/cofac.width(),
											  idxCnt%cofac.width())
											= mat(row,col);
										idxCnt++;
									}
								  }
							  }
						  }

						  // The determinant of 'mat' is the sum of its signed
						  // minors.
						  detSum += (1 - 2*(int)(curRow%2))*(1 - 2*(int)(idx%2))
							  		*mat(curRow,idx)*det(cofac);
					  }
				  }
				  break;
			////////////////////////////////////////////////////////////////////


			////////////////////////////////////////////////////////////////////
			default: 
				// Display error message and terminates the program
				cerr	<< endl
						<< "matrixmath.h: In function "
						<< "'template <class type>const double det(const "
						<< "matrix<type> mat, "
						<< "unsigned char dim= 1, size_t idx= 0)'"
						<< endl
						<< "matrixmath.h: error: "
						<< "dim= " << (unsigned int)dim << " is out of bounds!"
						<< endl
						<< "dim can be 1 or 2"
						<< endl << endl;

				// Exits program
				std::exit(EXIT_FAILURE);

				break;
			////////////////////////////////////////////////////////////////////
		}
	}

	return detSum;
}


//______________________________________________________________________________
// Reduces the matrix 'mat' to its row-echelon form


//______________________________________________________________________________
// Reduces the matrix 'mat' to its reduced row-echelon form


} // namespace matrixmath

#endif
