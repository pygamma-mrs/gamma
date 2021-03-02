/* MultiSysIF.h *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	GAMMA Multi-Systems                             Interface 	**
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
** modules setting up the multiple-system spin systems.                 **
**								 	**
*************************************************************************/

#ifndef __MultiSysIF_H__			// Is this already included?
#define __MultiSysIF_H__			// If no, then include it

#include <MultiSys/SpinMap.h>			// Include spin mappings
#include <MultiSys/ExProcess.h>			// Include exchange processes
#include <MultiSys/MultiSys.h>			// Include spin system
#include <MultiSys/MultiLib.h>			// Include library functions
#include <MultiSys/MultiSOp.h>			// Include spin operators
#include <MultiSys/MultiHam.h>			// Include HS Hamiltonians
#include <MultiSys/MultiIPul.h>			// Include ideal pulses
#include <MultiSys/MultiHSLib.h>		// Include HS library functions
#include <MultiSys/MultiLOp.h>			// Include super operators
#include <MultiSys/MultiWBR.h>			// Include BWR relaxation
#include <MultiSys/MultiExch.h>			// Include non-mutual exchange

#endif 						// MultiSysIF.h
