/* HSLib2IF.h ***************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      MR Hilbert Space Library 		  Interface 		**
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

#ifndef __HSLibIF_H__			// Is this already included?
#define __HSLibIF_H__			// If no, then include it

#include <HSLib/SpinOpSng.h>		// Single spin operators
#include <HSLib/SpinSys.h>		// Base spin systems
#include <HSLib/SpinOp.h>		// Generic spin operators
#include <HSLib/SpinOpCmp.h>		// Composite spin operators
#include <HSLib/SpinOpRot.h>		// Spin rotation operators
#include <HSLib/Basis.h>		// Operator bases
#include <HSLib/GenOpRep.h>		// Operator representations
#include <HSLib/GenOp.h>		// Operators
#include <HSLib/HSacquire.h>		// Acquisitions in Hilbert space
//#include <HSLib/HSevolve.h>		// Evolution in Hilbert space
#include <HSLib/HSdetect.h>		// Detection operators
#include <HSLib/HSdecomp.h>		// Density operator decomposition
#include <HSLib/HSprop.h>		// Hilbert space propagators
#include <HSLib/HSanalyze.h>		// System Analysis
#include <HSLib/SpinSystem.h>		// Include isotropic systems
#include <HSLib/HSham.h>		// Common MR Hamiltonians
#include <HSLib/HSauxil.h>		// More generic HS stuff
#include <HSLib/PulseI.h>		// Ideal pulses
#include <HSLib/PulseS.h>		// Real (rectangular) pulses
#include <HSLib/PulseShp.h>		// Shaped pulses
 


#endif 					// HSLibIF.h 
