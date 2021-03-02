/* gamma.h ******************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA 			                      Definition	**
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
** This file contains all include statements associated with the GAMMA	**
** simulation platform.  Thus, including this file in a program will	**
** include everything associated with GAMMA.				**
**								 	**
*************************************************************************/

#ifndef GAMMA_H__			// Is this already included?
#define GAMMA_H__			// If no, then include it


/*************************************************************************
**									**
**             System Dependent And Global Definitions			**
**									**
*************************************************************************/

#include <GamGen.h>

/*************************************************************************
**									**
**                      Support For GAMMA Basics			**
**									**
** 1. Constants used in MR calculations					**
** 2. Functions for interactive I/O					**
** 3. Spin isotopes							**
** 4. String cutting utilities 						**
** 5. Single Parameter & Parameter Set I/O To ASCII			**
**									**
*************************************************************************/

#include <Basics/GBasicsIF.h>		// All GAMMA basics headers

/*************************************************************************
**									**
**                      Support For GAMMA Matrices			**
**									**
** 1. Complex numbers							**
** 2. Row and column vectors						**
** 3. Generic matrices							**
** 4. Specialized matrix structures					**
**									**
*************************************************************************/

#include <Matrix/MatrixIF.h>		// All GAMMA matrix related headers

/*************************************************************************
**									**
**                Support Of Input/Output Formats		**
**									**
**									**
**									**
*************************************************************************/

#include <GamIO/GammaIOIF.h>		// All program interfaces

/*************************************************************************
**									**
**                Support For Hilbert Space Calculations		**
**									**
** 1. Spin Operators							**
** 2. Base Spin Systems							**
** 3. General Operators							**
** 4. Specialized matrix structures					**
**									**
*************************************************************************/

#include <HSLib/HSLibIF.h>		// All H.S. Library related headers

/*************************************************************************
**									**
**                      Support For Level 1 Calculations		**
**									**
** 1. Coordinates and coordinate vectors.				**
**									**
*************************************************************************/

#include <Level1/Level1IF.h>		// All Level 1 Library headers

/*************************************************************************
**									**
**                Support For Liouville Space Calculations		**
**									**
**									**
**									**
*************************************************************************/

#include <LSLib/LSLibIF.h>		// All L.S. Library related headers

/*************************************************************************
**									**
**                      Support For Level 2 Calculations		**
**									**
** 1. Transition Tables.						**
** 2. Rapid 1D Acquisitions						**
** 3. Quaternions                                                       **
**									**
*************************************************************************/

#include <Level2/Level2IF.h>		// All Level 2 Library headers

/*************************************************************************
**									**
**             Support For Bloch Equation Simulations			**
**									**
**									**
*************************************************************************/

#include <Bloch/BlochIF.h>		// All Bloch Library headers

/*************************************************************************
**									**
**             Support For Bloch Wangness Redfield Relaxation		**
**									**
**									**
*************************************************************************/

#include <BWRRelax/relaxBWRIF.h>	// All BWR relaxation 

/*************************************************************************
**									**
**                    Support For Rank 2 Interactions			**
**									**
** 1. Rank 2 Spatial Tensors						**
** 2. Rank 2 Spin Tensors						**
** 3. Rank 2 Interactions						**
** 4. Shift Anisotropy Interactions					**
** 5. Dipolar Interactions						**
** 6. Electron G Interactions						**
** 7. Hyperfine Interactions						**
** 8. Quadrupolar Interaction						**
** 9. Solid State Spin System						**
**									**
*************************************************************************/

#include <IntRank2/IntRank2IF.h>	// All Rank 2 Interactons Stuff

#include <Floquet/FloquetIF.h>		//All Floquet operator support

#include <Pulses/PulsesIF.h>		//All Pulses operator support

#include <MultiSys/MultiSysIF.h>	// The multi-spinsystem module

#include <Deprecated/DeprecIF.h>	// Include deprecated codes

#include <ESRLib/ESRIF.h>		// Include the ESR library

#include <Gradients/GradIF.h>		// Include gradients

#include <Testing/TestingIF.h>		// Include testing

#endif 							// gamma.h
