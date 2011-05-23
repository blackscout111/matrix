#include <iostream>
#include <string>
#include <iomanip>
#include "matrix.h"
#include "matrixmath.h"

using namespace std;

void dispNumMat(const matrix<int>& mat)
{
	for (_DIM i = 0; i < mat.height(); i++)
	{
		for (_DIM j = 0; j < mat.width(); j++)
		{
			cout << setw(3) << mat(i,j);
		}

		cout << endl;
	}

	cout << endl;
}

void dispNumMat(const matrix<float>& mat)
{
	for (_DIM i = 0; i < mat.height(); i++)
	{
		for (_DIM j = 0; j < mat.width(); j++)
		{
			cout << setw(4) << mat(i,j);
		}

		cout << endl;
	}

	cout << endl;
}

void dispNumMat(const matrix<double>& mat)
{
	for (_DIM i = 0; i < mat.height(); i++)
	{
		for (_DIM j = 0; j < mat.width(); j++)
		{
			cout << setw(4) << mat(i,j);
		}

		cout << endl;
	}

	cout << endl;
}

void dispStrMat(const matrix<string>& strmat)
{
	for (_DIM i = 0; i < strmat.height(); i++)
	{
		for (_DIM j = 0; j < strmat.width(); j++)
		{
			cout << setw(15) << strmat(i,j);
		}

		cout << endl;
	}

	cout << endl;
}

//////////////////////////////////
// Main Function Starts			//
//////////////////////////////////
int main()
{
	matrix <int> a(4,9);

	cout << "a.width() =" << a.width() << endl;
	cout << "a.height() =" << a.height() << endl;

	// Fill a
	int count = 0;
	for (_DIM i = 0; i < a.height(); i++)
	{
		for (_DIM j = 0; j < a.width(); j++)
		{
			a(i,j) = count;
			count++;
		}
	}

	// Display a
	dispNumMat(a);

	// Display max(a)
	dispNumMat(max(a));

	// Display max2d(a)
	printf("%d \n\n",max2d(a));

	// Display min(a)
	dispNumMat(min(a));

	// Display min2d(a)
	printf("%d \n\n",min2d(a));

	// Display sum(a)
	dispNumMat(sum(a));

	// Display sum2d(a)
	printf("%d \n\n",sum2d(a));

	// Display mean(a)
	dispNumMat(mean(a));

	// Display mean2d(a)
	printf("%d \n\n",mean2d(a));

	// Display vari(a)
	dispNumMat(vari(a));

	// Display vari2d(a)
	printf("%d \n\n",vari2d(a));

	// Display stdev(a)
	dispNumMat(stdev(a));

	// Display stdev2d(a)
	printf("%d \n\n",stdev2d(a));

	// Create other matrix
	matrix <int> b(5,5);

	// Fill other matrix with 2's
	b.fill(2);

	// Display b
	dispNumMat(b);

	// Sets b equal to a
	b = a;

	// Display b
	dispNumMat(b);

	// Crops b
	b.crop(1,2,2,4);

	// Display b
	dispNumMat(b);

	// Grows b
	b.grow(1,2,3,4,88);

	// Display b
	dispNumMat(b);

	// Grows b
	b.grow(2,3,3,2,-1);

	// Display b
	dispNumMat(b);

	// Create another matrix and set equal to a portion of b
	matrix <int> c = b(1,3,4,7);

	// Display c
	dispNumMat(c);

	// Make copy of c
	matrix <int> d = c;

	// Fill d with -5
	d.fill(1);

	// Display d
	dispNumMat(d);

	// Display d + c
	dispNumMat(d + c);

	// Display d - d
	dispNumMat(d - d);

	// Display 10*d
	dispNumMat((int)10*d);

	// Display d*20
	dispNumMat(d*(int)2);

	// A new matrix
	matrix<int> e(4,3);
	e.fill(1);
	dispNumMat(e);

	// Display (d*trans(d))
	dispNumMat(d * trans(d));

	// Set d = d*(2*e)
	dispNumMat(d *= (int)2*e);

	// Displays (d*2 - 3*(ones(d.height(),d.width())))
	dispNumMat(d*(int)2 - (int)3*idnty<int>(d.height(),d.width()));

	printf("%d\n\n",trace(d*(int)2 - (int)3*idnty<int>(d.height(),d.width())));

	dispNumMat(idnty<int>(10));

	// Creates and displays string matrix
	matrix <string> strmat(3,2);

	// Fills strmat with "*"'s
	strmat.fill("*");

	// Displays strmat
	dispStrMat(strmat);

	// Fills string matrix
	strmat(0,0) = "Hi";
	strmat(0,1) = "How are you?";
	strmat(1,0) = "Hola";
	strmat(1,1) = "?Como estas?";
	strmat(2,0) = "Yo";
	strmat(2,1) = "Sup?";

	// Displays strmat
	dispStrMat(strmat);

	// Creates other string matrix
	matrix <string> otherstrmat(strmat);
	
	// Displays otherstrmat
	dispStrMat(otherstrmat);

	// Grows otherstrmat
	otherstrmat.grow(4,3,1,1,"silence");

	// Displays otherstrmat
	dispStrMat(otherstrmat);

	// Swap row 6 with row 5 of otherstrmat
	otherstrmat.rowSwap(6,5);

	// Swap col 6 with col 5 of otherstrmat
	otherstrmat.colSwap(1,2);

	// Displays otherstrmat
	dispStrMat(otherstrmat);

	////////////////////////////////////////////////////////////////////////////
	// Test the determinant function
	////////////////////////////////////////////////////////////////////////////
	matrix <float> fltMat2x2(2,2);
	fltMat2x2(0,0) = 1.3;
	fltMat2x2(0,1) = 2;
	fltMat2x2(1,0) = 3.8;
	fltMat2x2(1,1) = 4;
	dispNumMat(fltMat2x2);
	printf("det(fltMat2x2) = %f\n\n",det(fltMat2x2));


	matrix <float> fltMat3x3(3,3);
	fltMat3x3(0,0) = 1.3;
	fltMat3x3(0,1) = 2;
	fltMat3x3(0,2) = -1;
	fltMat3x3(1,0) = 5;
	fltMat3x3(1,1) = 3.8;
	fltMat3x3(1,2) = 4;
	fltMat3x3(2,0) = -2;
	fltMat3x3(2,1) = 0;
	fltMat3x3(2,2) = 1;
	dispNumMat(fltMat3x3);
	printf("det(fltMat3x3) = %f\n\n",det(fltMat3x3));
	printf("det(fltMat3x3) = %f\n\n",det(fltMat3x3,2,1));
	printf("det(fltMat3x3) = %f\n\n",det(fltMat3x3,2,2));
	////////////////////////////////////////////////////////////////////////////

	return 0;
}
