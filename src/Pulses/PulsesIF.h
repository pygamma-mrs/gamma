/* PulsesIF.h ***************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      NMR Pulses                                   Interface 		**
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
** Resonance Experiments and Other Associated Mathematical            **
** Capabilities.  This file contains all include statements associated	**
** with the GAMMA module supporting multiple step pulses.		**
**								 	**
*************************************************************************/

#ifndef __PulPkg_H__			// Is this already included?
#define __PulPkg_H__			// If no, then include it

//#include <Pulses/PSeqAux.h>
//#include <Pulses/PTrain.h>
//#include <Pulses/PTrainAux.h>
//#include <Pulses/P_Gauss.h>
//#include <Pulses/P_Shaped.h>
#include <Pulses/PulAuxil.h>
#include <Pulses/PulCHIRP.h>
#include <Pulses/PulComposite.h>
#include <Pulses/PulCycle.h>
#include <Pulses/PulDANTE.h>
#include <Pulses/PulGARP.h>
#include <Pulses/PulGauss.h>
#include <Pulses/PulSinc.h>
#include <Pulses/PulMLEV.h>
#include <Pulses/PulSupCycle.h>
#include <Pulses/PulTrain.h>
#include <Pulses/PulTrainSCyc.h>
#include <Pulses/PulWALTZ.h>
#include <Pulses/PulWaveform.h>
#include <Pulses/Pulse.h>

#endif 					// PulsesIF.h 
