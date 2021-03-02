/* ConstTest.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      GAMMA Constants				Implementation		**
**                                                                      **
**      Copyright (c) 2002						**
**      Scott A. Smith							**
**      National High Magnetic Field Laboratory                         ** 
**      1800 E. Paul Dirac Drive                                        ** 
**      Tallahassee Florida, 32306-4005                                 ** 
**      email: gamma@magnet.fsu.edu                                     **
**      www: http://gamma.magnet.fsu.edu                                **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains constants and generic functions that are of use  **
**  in the GAMMA testing framework.					**
**                                                                      **
*************************************************************************/

#ifndef   CONSTTEST_CC__		// Is this file known?
#  define CONSTTEST_CC__ 1		// If not, define & include
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Testing/ConstTest.h>			// Include our own header

// ____________________________________________________________________________
//                    Testing FrameWork Base Directory
// ____________________________________________________________________________

/* Many of the GAMMA tests read local parameter files. When the tests are run
   from any arbitrary directory they must be told where the parameter files
   can be found. This constant sets the global base directory for the GAMMA
   tests. By default this is set to ".", that is it assumes things are run
   locally.  However, it may be set to something like
   /gamma/gamma-4.1.0/test.                                                  */

std::string GTESTDIR = std::string("./Standards");

void        SetGTestDir(const std::string& tdir) { GTESTDIR = tdir; }
std::string GetGTestDir()                        { return GTESTDIR; }

// ____________________________________________________________________________
//                    Testing FrameWork ParameterSet Directory
// ____________________________________________________________________________

/* Many of the GAMMA tests need/read/use parameter files. When the tests are
   run from any arbitrary directory they must be told where the parameter files
   can be found. This is because the "local" directory during execution will be
   just that, the directory where the program is executed. The directory where
   the parameter files exist must be set relative to the local (execution) 
   directory.

   So, how then do we tell each test program where to find the directory that
   contains any of its parameter files? For all tests the parameter files will
   be placed in the directory "Standards" stemming off of this location (this
   directory). There are several cases that must then be considered:

      Test Type        Exec. Dir.     GTESTPDIR              Result
   ===============  ================  =========  ==============================
   Module Specific  Module Directory    Empty    Pgm sets where Standards is.
 * Module Specific   Any Directory    Relative   Pgm uses GTESTPDIR setting
** GAMMA  Summed    This Directory   ./Standards Pgm uses GTESTPDIR setting
 * GAMMA  Summed     Any Directory 
  
   *: In these cases a function call to SetGTestPDir must be made in order to
      specify where the Standards directory is relative to the test execution
      directory.
  **: In this case the flag  
   
                                                */

std::string GTESTPDIR = std::string("./Standards");

void        SetGTestPDir(const std::string& tdir) { GTESTPDIR = tdir; }
std::string GetGTestPDir()                        { return GTESTPDIR; }

  
#endif						// ConstTest.cc
