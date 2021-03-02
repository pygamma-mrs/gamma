/* TestingIF.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA Testing Module				Interface 	**
**								 	**
**      Scott Smith                                                     **
**	Copyright (c) 2000						**
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
** testing module.  This file includes classes and functions that are	**
** used in testing the GAMMA platform. 					**
**								 	**
*************************************************************************/

#ifndef __TestingIF_H__				// Is this already included?
#define __TestingIF_H__				// If no, then include it

#include <Testing/SingleTest.h>			// Include single tests
#include <Testing/SectTest.h>			// Include section tests
#include <Testing/ClassTest.h>			// Include class tests
#include <Testing/ModTest.h>			// Include module tests
#include <Testing/GamTest.h>			// Include GAMMA tests
#include <Testing/ConstTest.h>			// Include directory names

#endif 						// TestingIF
