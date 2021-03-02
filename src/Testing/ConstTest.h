/* ConstTest.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      GAMMA Testing Constants				Interface	**
**                                                                      **
**      Copyright (c) 2002 						**
**      Scott A. Smith 							**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**      email: gamma@magnet.fsu.edu					**
**      www: http://gamma.magnet.fsu.edu				**
**                                                                      **
**      $Header: 
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains constants that are of use in the GAMMA testing   **
**  framework.								**
**                                                                      **
*************************************************************************/
 
#ifndef   ConstTest_H__			// Is this module known?
#  define ConstTest_H__	1		// If not then define & include
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Know about libstdc++ strings

// ____________________________________________________________________________ 
//                    Testing FrameWork Base Directory
// ____________________________________________________________________________

/* Many of the GAMMA tests read local parameter files. When the tests are run
   from any arbitrary directory they must be told where the parameter files
   can be found. This constant sets the global base directory for the GAMMA
   tests. By default this is set to ".", that is it assumes things are run 
   locally.  However, it may be set to something like 
   /gamma/gamma-4.1.0/test.	                                             */ 

extern MSVCDLL std::string GTESTDIR;		

MSVCDLL void        SetGTestDir(const std::string& tdir); 
MSVCDLL std::string GetGTestDir(); 

// ____________________________________________________________________________ 
//                    Testing FrameWork ParameterSet Directory
// ____________________________________________________________________________

/* Many of the GAMMA tests read local parameter files. When the tests are run
   from any arbitrary directory they must be told where the parameter files
   can be found. This constant sets the global base directory for the GAMMA
   tests. By default this is set to ".", that is it assumes things are run 
   locally.  However, it may be set to something like 
   /gamma/gamma-4.1.0/test.	                                             */ 

extern MSVCDLL std::string GTESTPDIR;		

MSVCDLL void        SetGTestPDir(const std::string& tdir); 
MSVCDLL std::string GetGTestPDir(); 

#endif							// ConstTest.h
