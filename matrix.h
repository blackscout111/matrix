#ifndef matrix_h
#define matrix_h

#include <iostream>
#include <cstdio>
#include <cstdlib>


using namespace std;

namespace matrix
{

// matrix is a 2-D array with elements of type 'type'.
template <class type>
class matrix
{
	public:
		// Default constructor
		// (Allocates memory for matrix given numRows and numCols)
		matrix(size_t height = 1, size_t width = 1);

		// Returns numRows
		size_t width() const;

		// Returns numCols
		size_t height() const;

		// Copy constructor
		matrix(const matrix<type>& other);

		// Destructor
		~matrix();

		// Overloads the () operator
		const type& operator ()(size_t row, size_t col) const;
		type& operator ()(size_t row, size_t col);
		matrix<type> operator ()(size_t firstRow,
								 size_t lastRow,
								 size_t firstCol,
								 size_t lastCol);

		// Sets all entries to "val"
		void fill(const type& val);

		// Swaps row 'i' with row 'j'
		void rowSwap(size_t i, size_t j);

		// Swaps column 'i' with row 'j'
		void colSwap(size_t i, size_t j);
		
		// Crops the matrix
		void crop(size_t minrow,
			  	  size_t maxrow,
			  	  size_t mincol,
			  	  size_t maxcol);

		// Increases the size of the matrix
		void grow(size_t above,
			  	  size_t below,
			  	  size_t left,
			  	  size_t right);

		// Increases the size of the matrix
		// (New spaces initialized to "val")
		void grow(size_t above,
			  	  size_t below,
			  	  size_t left,
			  	  size_t right,
			  	  type val);

		// Overloads the = operator
		matrix<type>& operator =(const matrix<type>& other);

	protected:
		size_t	numRows,	// # of rows of elem
				numCols;	// # of columns of elem

		type	*elem;		// The elements matrix
};

/**************************************************************************
* Function Definitions													  *
**************************************************************************/
//_________________________________________________________________________
template <class type>
matrix<type>::matrix(size_t height, size_t width)
{
	// Sets dimensions of matrix
	numRows = height;
	numCols = width;

	// Allocates memory
	elem = new type [numRows * numCols];
}


//_________________________________________________________________________
template <class type>
size_t matrix<type>::width()	const
{
	return numCols;
}


//_________________________________________________________________________
template <class type>
size_t matrix<type>::height()	const
{
	return numRows;
}
	

//_________________________________________________________________________
template <class type>
matrix<type>::matrix(const matrix<type>& other)
{
	// Copies height and width values
	numRows = other.height();
	numCols = other.width();

	// Allocates memory
	elem = new type [numRows * numCols];

	// Copies data from other into element array
	for (size_t i = 0; i < numRows; ++i)
	{
		for (size_t j = 0; j < numCols; ++j)
		{
			(*this)(i,j) = other(i,j);
		}
	}
}


//_________________________________________________________________________
template <class type>
matrix<type>::~matrix()
{
	delete[] elem;
	elem = NULL;
}


//_________________________________________________________________________
template <class type>
const type& matrix<type>::operator ()(size_t row,
									  size_t col) const
{
	if (row >= 0 && row < numRows && col >= 0 && col < numCols)
	{
		return elem[numCols * row + col];
	}
	else
	{
		// Display error message and exit the program
		cerr	<< endl
				<< "matrix.h: In member function "
				<< "'template <class type>matrix<type>&"
				<< "matrix<type>::operator"
				<< " ()(size_t row, size_t col)'"
				<< endl
				<< "matrix.h: error: "
				<< "Index ("
				<< row << "," << col
				<< ") is out of bounds!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}


//_________________________________________________________________________
template <class type>
type& matrix<type>::operator ()(size_t row, size_t col)
{
	if (row >= 0 && row < numRows && col >= 0 && col < numCols)
	{
		return elem[numCols * row + col];
	}
	else
	{
		// Display error message and exit the program
		cerr	<< endl
				<< "matrix.h: In member function "
				<< "'template <class type>matrix<type>&"
				<< "matrix<type>::operator"
				<< " ()(size_t row, size_t col)'"
				<< endl
				<< "matrix.h: error: "
				<< "Index ("
				<< row << "," << col
				<< ") is out of bounds!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}


//_________________________________________________________________________
template <class type>
matrix<type> matrix<type>::operator ()(size_t firstRow,
									   size_t lastRow,
									   size_t firstCol,
									   size_t lastCol)
{
	// Check to see if indices are in range
	if (firstRow >= 0 && firstRow < numRows &&
		lastRow >= 0 && lastRow < numRows &&
		lastCol >= 0 && lastCol < numCols &&
		firstCol >= 0 && firstCol < numCols)
	{
		// Check to see that lastRow/Col >= firstRow/Col
		if (lastRow >= firstRow && lastCol >= firstCol)
		{
			// The temporary matrix to be returned
			matrix <type> tempMat((1 + lastRow - firstRow),
								  (1 + lastCol - firstCol));

			// Fill the temporary matrix
			for (size_t i = 0; i < tempMat.height(); ++i)
			{
				for (size_t j = 0; j < tempMat.width(); ++j)
				{
					tempMat(i,j) = (*this)(i + firstRow, j + firstCol);
				}
			}

			// Return the temporary matrix
			return tempMat;
		}
		else
		{
			// Display error message and exit the program
			cerr	<< endl
					<< "matrix.h: In member function "
					<< "'template <class type>matrix<type>&"
					<< "matrix<type>::operator"
					<< " ()(size_t firstRow, size_t lastRow,"
					<< " size_t firstCol, size_t lastCol)'"
					<< endl
					<< "matrix.h: error: "
					<< "In index range ("
					<< firstRow << "," << lastRow << ","
					<< firstCol << "," << lastCol
					<< ") starting index must be less than ending index!"
					<< endl << endl;

			// Exits program
			std::exit(EXIT_FAILURE);
		}
	}
	else
	{
		// Display error message and exit the program
		cerr	<< endl
				<< "matrix.h: In member function "
				<< "'template <class type>matrix<type>&"
				<< "matrix<type>::operator"
				<< " ()(size_t firstRow, size_t lastRow,"
				<< " size_t firstCol, size_t lastCol)'"
				<< endl
				<< "matrix.h: error: "
				<< "Index range ("
				<< firstRow << "," << lastRow << ","
				<< firstCol << "," << lastCol
				<< ") is out of bounds!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}
}


//_________________________________________________________________________
template <class type>
void matrix<type>::fill(const type& val)
{
	// Fills elem[] with val
	for (size_t i = 0; i < (numRows * numCols); ++i)
	{
		elem[i] = val;
	}
}


//_________________________________________________________________________
// Swaps row 'i' with row 'j'
template <class type>
void matrix<type>::rowSwap(size_t i, size_t j)
{
	// Temporary element
	type temp;
	
	// Swap the elements in row 'i' with the elements in row 'j'
	for (size_t k = 0; k < numCols; ++k)
	{
		temp = (*this)(i,k);
		(*this)(i,k) = (*this)(j,k);
		(*this)(j,k) = temp;
	}
}


//_________________________________________________________________________
// Swaps column 'i' with column 'j'
template <class type>
void matrix<type>::colSwap(size_t i, size_t j)
{
	// Temporary element
	type temp;
	
	// Swap the elements in column 'i' with the elements in column 'j'
	for (size_t k = 0; k < numRows; ++k)
	{
		temp = (*this)(k,i);
		(*this)(k,i) = (*this)(k,j);
		(*this)(k,j) = temp;
	}
}


//_________________________________________________________________________
template <class type>
void matrix<type>::crop(size_t minrow,
			size_t maxrow,
			size_t mincol,
			size_t maxcol)
{
	// Checks that max's are greater than min's
	if ((minrow > maxrow) || (mincol > maxcol))
	{
		// Display error message and exit the program
		cerr	<< endl
				<< "matrix.h: In member function "
				<< "'template <class type>void"
				<< "matrix<type>::crop(size_t minrow, "
				<< "size_t maxrow,"
				<< "size_t mincol, size_t maxcol)'"
				<< endl
				<< "matrix.h: error: "
				<< "Max must be greater than min!"
				<< endl << endl;

		// Exits program
		std::exit(EXIT_FAILURE);
	}

	// New dimensions
	size_t	height = 1 + maxrow - minrow,
			width = 1 + maxcol - mincol;

	// Temporary matrix
	type *temp;
	temp = new type [height * width];

	// Copies data from elem to temp
	for (size_t row = minrow; row <= maxrow; ++row)
	{
		for (size_t col = mincol; col <= maxcol; ++col)
		{
			temp[width * (row - minrow) + col - mincol] = (*this)(row,col);
		}
	}

	// Frees the old elem memory
	delete[] elem;

	// Sets elem to temp
	elem = temp;

	// Sets new dimensions
	numRows = height;
	numCols = width;
}


//_________________________________________________________________________
template <class type>
void matrix<type>::grow(size_t above,
						size_t below,
						size_t left,
						size_t right)
{
	// New dimensions
	size_t	height = numRows + above + below,
			width = numCols + left + right;

	// Temporary matrix
	type *temp;
	temp = new type [height * width];

	// Copies data from elem to temp
	for (size_t row = above; row < (numRows + above); ++row)
	{
		for (size_t col = left; col < (numCols + left); ++col)
		{
			temp[width * row + col] = (*this)((row - above), (col - left));
		}
	}

	// Frees the old elem memory
	delete[] elem;

	// Sets elem to temp
	elem = temp;

	// Sets new dimensions
	numRows = height;
	numCols = width;
}


//_________________________________________________________________________
template <class type>
void matrix<type>::grow(size_t above,
						size_t below,
						size_t left,
						size_t right,
						type val)
{
	// Increases size of matrix
	grow(above,below,left,right);

	// Fills extra spaces
	// (Top rows)
	for (size_t row = 0; row < above; ++row)
	{
		for (size_t col = 0; col < numCols; ++col)
		{
			(*this)(row,col) = val;
		}
	}

	// (Bottom rows)
	for (size_t row = (numRows - 1); row >= (numRows - below); --row)
	{
		for (size_t col = 0; col < numCols; ++col)
		{
			(*this)(row,col) = val;
		}
	}

	// (Left-over spaces)
	for (size_t row = above; row < (numRows - below); ++row)
	{
		// Left spaces
		for (size_t col = 0; col < left; ++col)
		{
			(*this)(row,col) = val;
		}

		// Right spaces
		for (size_t col = (numCols - 1); col >= (numCols - right); --col)
		{
			(*this)(row,col) = val;
		}
	}
}


//_________________________________________________________________________
template <class type>
matrix<type>& matrix<type>::
	operator =(const matrix<type>& other)
{
	// Check for self assignment
	if (this != &other)
	{
		// Frees elem memory
		delete[] elem;

		// Copies height and width values
		numRows = other.height();
		numCols = other.width();

		// Allocates memory
		elem = new type [numRows * numCols];

		// Copies data from other into element array
		for (size_t i = 0; i < numRows; ++i)
		{
			for (size_t j = 0; j < numCols; ++j)
			{
				(*this)(i,j) = other(i,j);
			}
		}
	}

	return *this;
}

} // namespace matrix

#endif
