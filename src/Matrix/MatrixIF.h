/* matrixIF.h ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	Matrices & Vectors 			    Full Interface 	**
**								 	**
**	Copyright (c) 1991, 1992, 1997				 	**
**	Tilo Levante, Scott A. Smith            		 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** This file contains all the include statements associated with the	**
** matrices and vectors contained in the GAMMA simulation platform.	**
** Thus, including this file in a program will includes everything	**
** needed to use 1 & 2 dimensional arrays in GAMMA.			**
**								 	**
*************************************************************************/

#ifndef GmatrixIF_H__			// Is this already included?
#define GmatrixIF_H__ 1			// If no, then include it

#include <Matrix/MxModBas.h>            // Module errors
#include <Matrix/complex.h>             // Complex numbers
#include <Matrix/_matrix.h>             // Base matrix class
#include <Matrix/i_matrix.h>		// Identity matrix
#include <Matrix/d_matrix.h>		// Diagonal matricies
#include <Matrix/n_matrix.h>		// Normal (full) matrix
#include <Matrix/h_matrix.h>		// Hermitian matrix
#include <Matrix/matrix.h>              // Matrix (collection of all types)
#include <Matrix/col_vector.h>		// Column vectors (derived from matrix)
#include <Matrix/row_vector.h>		// Row vectors    (derived from matrix)

#endif 					// MatrixIF.h 
