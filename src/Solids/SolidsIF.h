/* IntRank2IF.h *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	Rank 2 Interactions 				Interface 	**
**								 	**
**      Scott Smith                                                     **
**	Copyright (c) 1998						**
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
** class supporing irreducible rank 2 spin interactions.  Including 	**
** this file in a program will include everything associated with 	**
** such rank 2 interactions.						**
**								 	**
*************************************************************************/

#ifndef __IntRank2IF_H__		// Is this already included?
#define __IntRank2IF_H__		// If no, then include it

#include <IntRank2/IntRank2A.h>		// Rank 2 interation space tensors
#include <IntRank2/IntRank2T.h>		// Rank 2 interaction spin tensors
#include <IntRank2/IntRank2.h>		// Rank 2 interactions
#include <IntRank2/IntDip.h>		// Dipolar interactions
#include <IntRank2/IntQuad.h>		// Quadrupolar interactions
#include <IntRank2/IntCSA.h>		// Shift Anisotropy interactions
#include <IntRank2/IntG.h>		// Electron G interactions
#include <IntRank2/IntHF.h>		// Hyperfine interactions
#include <IntRank2/IntDipVec.h>		// Vector of Dipolar interactions
#include <IntRank2/IntCSAVec.h>		// Vector of SA interactions
#include <IntRank2/IntQuadVec.h>	// Vector of Quadrupolar interacts
#include <IntRank2/IntGVec.h>		// Vector of G interactions
#include <IntRank2/IntHFVec.h>		// Vector of HF interactions
#include <IntRank2/SolidSys.h>		// Solid state systems

#endif 					// IntRank2IF.h
