/* Level1IF.h ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA Level 1 Library				Interface 	**
**								 	**
**      Scott Smith                                                     **
**	Copyright (c) 1999						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** This file contains all include statements associated with the GAMMA	**
** level 1 library.  This file includes classes and functions that are	**
** used in higher level GAMMA workings.  Herein are coordinates and	**
** coordinate vectors.							**
**								 	**
*************************************************************************/

#ifndef __Level1IF_H__				// Is this already included?
#define __Level1IF_H__				// If no, then include it

#include <Level1/coord.h>			// Include coordinates
#include <Level1/coord_vec.h>			// Include coordinate vectors
#include <Level1/SphHarmic.h>			// Include spherical harmonics
#include <Level1/Wigner.h>			// Include Wigner rotations
#include <Level1/SpinT.h>			// Include spin tensors
#include <Level1/SpaceT.h>			// Include space tensors
#include <Level1/nmr_tensor.h>			// Include common MR spin tensors
#include <Level1/Exponential.h>			// Include exponential functions
#include <Level1/Lorentzian.h>			// Include Lorentzian functions
#include <Level1/WindowFct.h>			// Include Windowing functions
#include <Level1/ExProcessM.h>			// Include mutual exchange processes

#endif 						// Level1IF
