/* MxModBas.h ***************************************************-*-c++-*-
**                                                                      **
**	                               G A M M A                         **
**                                                                      **
**	Matrix Module Basics                           Interface         **
**                                                                      **
**	Copyright (c) 2000                                               **
**	Scott A. Smith                                                   **
**	Eidgenoessische Technische Hochschule                            **
**	Labor fur physikalische Chemie                                   **
**	8092 Zurich / Switzerland                                        **
**                                                                      **
**      $Header:  $
**                                                                      **
*************************************************************************/
 
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** This file contains some basic functions that are used throughout     **
** the GAMMA Matrix module. It is included in most Matrix module files  **
** and allows the module to be independent of the rest of GAMMA.        **
** Furthermore, it lists the type of GAMMA matrices available. These    **
** matrices are the types derived from the base class _matrix and which **
** fit transparently under the umbrella of class matrix.  This file is  **
** intended to be included in each and every GAMMA matrix class, i.e.   **
** all classes derived from _matrix, so that each will know about the   **
** other. A necessity when performing matrix-matrix operations. 	**
**                                                                      **
*************************************************************************/

#ifndef   MxModBas_h_			// Is this file already included?
#  define MxModBas_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include string manipulations
#include <iostream>			// If we don't already know cout/cin

// _____________________________________________________________________________
// A              Matrix Module Generic Error Handling Functions
// _____________________________________________________________________________
 
/*
           Input                hdr     : Class or Module header name
                               *msg     : The error message to output
                               *eidx    : The error index
                               *pname   : Additional error string
                                noret   : Flag whether line feed needed
           Output               none    : An error message of some sort
                                          is sent to standard output          */
 
MSVCDLL void MxModError(const std::string& hdr, const std::string& msg, int noret=0);
MSVCDLL void MxModError(const std::string& hdr, int eidx, int noret=0);
MSVCDLL void MxModError(const std::string& hdr, int eidx, 
                                      const std::string& pname, int noret=0);
MSVCDLL volatile void MxModFatal();

// ____________________________________________________________________________ 
// B             Matrix Module Generic String Formatting Functions
// ____________________________________________________________________________
 
/* The functions dec and form are C functions which used to be also defined in
   the Gnu C++ library (g++).  However, at various times in GCC evolution they
   seem to have disappeared.  They are back in now, but are still missing in
   MSCV++.  As a consequence, they are redefined here.  With the proper includes
   current GCC still knows dec and form conversions to string, but other
   compilers (e.g. MSVC) do not so it is set up here.....                    */

MSVCDLL std::string MxModdec(int i);
MSVCDLL std::string MxModdec(const std::string& fmt, int i);
MSVCDLL std::string MxModdec(int i, int digs);
MSVCDLL std::string MxModform(const std::string& fmt, double d);

#endif							// MxModBas.h
