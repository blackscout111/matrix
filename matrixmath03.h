#ifndef matrixmath03_h
#define matrixmath03_h

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

#include "matrix.h"
#include "matrixmath00.h" // Needed for checking the type of matrix

using namespace std;

/*******************************************************************************
* This file normally included from matrixmath.h
*
* Section 03:	Functions that generate special matricies (e.g. identity matrix)
*******************************************************************************/
//______________________________________________________________________________
// Returns an identity matrix with dimensions of 'height' and 'width'
template <class type>
const matrix<type> idnty(_DIM height, _DIM width)
{
	matrix<type> mat(height,width);

	for (_DIM i = 0; i < height; i++)
	{
		for (_DIM j = 0; j < width; j++)
		{
                    if (i == j)
                    {
                        mat(i,j) = (type)1;
                    }
                    else
                    {
                        mat(i,j) = (type)0;
                    }
		}
	}

	return mat;
}



//______________________________________________________________________________
// Returns a square identity matrix with dimension of 'dim'
template <class type>
const matrix<type> idnty(_DIM dim=1)
{
    return idnty<type>(dim,dim);
}



//______________________________________________________________________________
// Returns a matrix filled with ones with dimensions of 'height' and 'width'
template <class type>
const matrix<type> ones(_DIM height, _DIM width)
{
    matrix<type> mat(height,width);

    // Fill matrix with ones
    mat.fill(1);

    return mat;
}



//______________________________________________________________________________
// Returns a square matrix filled with ones width dimension of 'dim'
template <class type>
const matrix<type> ones(_DIM dim=1){return ones<type>(dim,dim);}


//______________________________________________________________________________
// Returns a matrix filled with zeros width dimensions of 'height' and 'width'
template <class type>
const matrix<type> zeros(_DIM height, _DIM width)
{
    matrix<type> mat(height,width);

    // Fill matrix with zeros
    mat.fill(0);

    return mat;
}



//______________________________________________________________________________
// Returns a square matrix filled with zeros width dimension of 'dim'
template <class type>
const matrix<type> zeros(_DIM dim=1){return zeros<type>(dim,dim);}



#endif
