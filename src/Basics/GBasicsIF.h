/* GBasicsIF.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA Basics					Interface 	**
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
** modules setting up the platform basics.  Basics include commonly 	**
** required constants in MR calculations, often needed functions for	**
** interactive programs, spin isotope definitions, string parsing, as	**
** well as parameter and parameter set I/O to ASCII files.		**
**								 	**
*************************************************************************/

#ifndef GBasicsIF_H__				// Is this already included?
#define GBasicsIF_H__				// If no, then include it

#include <Basics/Gutils.h>			// Include query functions
#include <Basics/Gconstants.h>			// Include GAMMA constants

#include <Basics/StringCut.h>			// Include string parsing
#include <Basics/SinglePar.h>			// Include single parameters
#include <Basics/ParamSet.h>			// Include parameters sets

#include <Basics/IsotopeData.h>			// Single Isotope Definition
#include <Basics/Isotope.h>			// Single Isotope Class

#endif 						// GBasicsIF
