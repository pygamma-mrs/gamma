/* SpinOpSng.h **************************************************-*-c++-*-
**									**
**                                  G A M M A				**
**									**
**      Single Spin Operators				Interface	**
**									**
**      Copyright (c) 1990, 1991, 1992					**
**      Tilo Levante, Scott A. Smith					**
**      Eidgenoessische Technische Hochschule				**
**      Labor fuer physikalische Chemie					**
**      8092 Zuerich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
** This module contains functions which define the common single spin	**
** operators for GAMMA.  These are returned as matrices, represented in	**
** the common basis: alpha = (1,0,...,0) & beta = (0,1,0,...,0)		** 
**									**
** Currently provided are the following:				**
**									**
**                     Ix, Iy, Iz, Ie, Im, Ip, Ru(phi)                  **
**									**
*************************************************************************/

#ifndef   single_spin_op_h_		// Is file already included?
#  define single_spin_op_h_ 1		// If no, then remember that it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <Basics/Gconstants.h>		// Include constant PI 

/*	     Input          qn : The spin Hilbert space, (2I+1).
	     Output         I* : I* matrix of dimension qn
	     Note              : Included to span operator space
	     Note              : I=1/2 -> qn = 2, I=1 -> qn = 3
	     Note              : List mx_list herein maintains
				    all single spin I*'s produced                  */

MSVCDLL matrix Ie(int qn);		// Identity matrix of dimension qn
MSVCDLL matrix Ix(int qn);		// Ix matrix of dimension qn
MSVCDLL matrix Iy(int qn);		// Iy matrix of dimension qn
MSVCDLL matrix Iz(int qn);		// Iz matrix of dimension qn
MSVCDLL matrix Ip(int qn);		// I+ matrix of dimension qn
MSVCDLL matrix Im(int qn);		// I- matrix of dimension qn

/*                              [                                            ]
                       I=1/2    | cos(B/2)-iu sin(B/2)  (-iu -u )sin(B/2)    |
 R (B) = exp(-i*B*I ) ------->  |            z              x  y             |
  u                u            |                                            |
                                | (-iu +u )sin(B/2)     cos(B/2)+iu sin(B/2) |
                                [     x  y                         z         ] */

MSVCDLL matrix Raxis(int qn, double beta, char axis);
 
  //  Input          qn : The spin Hilbert space, (2I+1).
  //               beta : A rotation angle
  //               axis : A rotation axis
  //  Output    R(beta) : A rotation operator for spin
  //                      rotation about the defined axis
  //                      by angle beta.
  //  Note              : I=1/2 -> qn = 2, I=1 -> qn = 3
  //  Note              : Unlike other functions in this module,
  //                      there is no list of rotation operators
  //                      maintained herein.  That is because
  //                      these lists are ordered by qn and the
  //                      rotations must also be ordered by angle
  //                      and axis.
 
#endif						// SpinOpSng.h (Single Spin Operator)

