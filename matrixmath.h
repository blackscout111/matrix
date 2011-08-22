#ifndef matrixmath_h
#define matrixmath_h

/*******************************************************************************
* This file contains include statements for header files that include functions
* for performing mathematical operations with the matrix class.
* (e.g. multiplication, addition, row reduction).
*
* The functions are divided into multiple sections, with each section contained
* in its own separate header file. (Section # XX contained in matrixmathXX.h)
*
* ==============================================================================
* Section #:	Description:
* ==============================================================================
* 00			Functions for checking the types of matrices and #define values
* 01			Operator overloads(+,+=,-,-=,*)
* 02			Special matrix functions (e.g. transpose, reduction, statistics)
* 03            Functions that generate special matrices (e.g. identity matrix)
*
*******************************************************************************/

// Include statements for different sections
#include "matrixmath00.h"
#include "matrixmath01.h"
#include "matrixmath02.h"
#include "matrixmath03.h"

#endif
