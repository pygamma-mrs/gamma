/* relaxBWRIF.h *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA BWR Relaxation & Exchange			Definition	**
**								 	**
**      Scott Smith                                                     **
**      Copyright (c) 1998                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**									**
** Support For BWR Treatment of Spin Relaxation in GAMMA 		**
**									**
*************************************************************************/

#ifndef __BWR_H__			// Is this already included?
#define __BWR_H__			// If no, then include it

#include <BWRRelax/relaxNMR.h>		// Liquid nmr relaxation
#include <BWRRelax/relaxRF.h>		// Relaxation with an exteral field
//#include <BWRRelax/relaxAux.h>		// Relaxation auxiliary functions
#include <BWRRelax/relaxJ.h>		// Relaxation spectral density functions
#include <BWRRelax/relaxDip.h>		// Dipolar relaxation
#include <BWRRelax/relaxCSA.h>		// Chemical shift anisotropy relaxation
#include <BWRRelax/relaxQuad.h>		// Quadrupolar relaxation
#include <BWRRelax/relaxRand.h>		// Random field relaxation
#include <BWRRelax/relaxDCSA.h>		// Dipolar-csa cross relaxation
//#include <BWRRelax/relaxDQuad.h>	// Dipolar-quadrupolar cross relaxation
#include <BWRRelax/relaxProp.h>		// Relaxation propagators
#include <BWRRelax/relaxExch.h>		// Exchange superoperators
#include <BWRRelax/relaxBWRexch.h>	// BWR relaxation and exchange

#endif 					// relaxBWRIF.h
