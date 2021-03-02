/* GradIF.h *****************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Gradients				     Interface 		**
**                                                                      **
**      Copyright (c) 1999						**
**      Scott A. Smith							**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     **
** Resonance Experiments and Other Associated   Mathematical            **
** Capabilities.  This file contains all include statements associated	**
** with the GAMMA module supporting computation in Hilbert space.	**
**								 	**
*************************************************************************/

#ifndef __GradIF_H__			// Is this already included?
#define __GradIF_H__			// If no, then include it

#include <Gradients/sys_gradz.h>	// System in z-gradient
#include <Gradients/Gradients2.h>
#include <Gradients/GrdPulses.h>	// Pulses in gradients
#include <Gradients/GrdEvolve.h>	// Time evolutions in gradients
#include <Gradients/GrdAcquire.h>	// Acqisitions in gradients
#include <Gradients/GrdDeprec.h>	// Deprecated Functions

#endif 					// GradIF.h 
