/* BlochIF.h ****************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Bloch Module 				Interface 		**
**                                                                      **
**      Copyright (c) 2001						**
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
** This module facilitates simulations that intend to utilize the	**
** phenomenological Bloch equations.					**
**								 	**
*************************************************************************/

#ifndef  BlochIF_H__			// Is this already included?
#define  BlochIF_H__			// If no, then include it

#include <Bloch/BlochMx.h>		// Bloch matrices
#include <Bloch/MagVec.h>		// Magnetization vectors
#include <Bloch/DetVec.h>		// Detection vectors
#include <Bloch/BlochSys.h>		// Bloch system of vectors
#include <Bloch/BlochAcq.h>		// Bloch acquisitions
#include <Bloch/BlochB.h>		// Bloch field & offset matrices
#include <Bloch/BlochR.h>		// Bloch relaxation matrices
#include <Bloch/BlochK.h>		// Bloch exchange matrices
#include <Bloch/BlochM.h>		// Bloch special mag. vectors
#include <Bloch/BlochTraj.h>		// Bloch trajectories
#include <Bloch/Bloch.h>		// Various Bloch equation functions

#endif 					// BlochIF.h 
