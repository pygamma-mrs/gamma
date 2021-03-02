/* MxModBas.cc **************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A				**
**                                                                      **
**      Matrix Module Basics                        Implementation	**
**                                                                      **
**      Copyright (c) 2000						**
**      Scott A. Smith                                                  **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fur physikalische Chemie                                  **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header:  $
**                                                                      **
*************************************************************************/
 
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** This file contains some basic functions that are used throughout	**
** the GAMMA Matrix module. It is included in most Matrix module files  **
** and allows the module to be independent of the rest of GAMMA.	**
**                                                                      **
*************************************************************************/
 
#ifndef   MxModBas_cc_			// Is this file already included?
#  define  MxModBas_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the interface
#  endif

#include <Matrix/MxModBas.h>		// Include the implementation
#include <stdlib.h>
#include <string>			// Know about string class
#ifndef MSC_VER				// If not using MS Visual C++
  #include <stdio.h>			// know about sprintf function
#endif

// ____________________________________________________________________________ 
// A                Matrix Module Generic Error Handling Functions
// ____________________________________________________________________________
 
/*
           Input                hdr     : Class or Module header name
                               *msg     : The error message to output
                               *eidx    : The error index
	  		       *pname   : Additional error string
	  			noret   : Flag whether line feed needed
           Output               none    : An error message of some sort
				          is sent to standard output          */

 
void MxModError(const std::string& hdr, const std::string& msg, int noret)
  {
  std::string s = "\n" + hdr + ": " + msg;		// Set error string
  std::cout << s; if(!noret) std::cout << std::endl;	// GAMMA error stuff
  }   
  
void MxModError(const std::string& hdr, int eidx, int noret)
  { 
  std::string msg;
  switch(eidx)
    {
    case 0:  msg = std::string("Program Aborting.....");            break; 	// (0)
    case 1:  msg = std::string("Problems With Input File Stream");  break; 	// (1)
    case 2:  msg = std::string("Problems With Output File Stream"); break; 	// (2)
    case 3:  msg = std::string("Can't Construct From Parameter Set"); break;	// (3)
    case 4:  msg = std::string("Cannot Construct From Input File"); break;	// (4)
    case 5:  msg = std::string("Cannot Write To Parameter File");   break;	// (5)
    case 6:  msg = std::string("Cannot Write To Output FileStream");break;	// (6)
    case 7:  msg = std::string("Cannot Produce Any Output");	       break;	// (7)
    case 8:  msg = std::string("Problems With Parameter Set");      break;	// (8)
    case 9:  msg = std::string("Problems During Construction");     break;	// (9)
    case 10: msg = std::string("Cannot Access Internal Component"); break;	// (10)
    case 11: msg = std::string("Cannot Read From Input FileStream");break;	// (11)
    case 12: msg = std::string("End of File Reached");	       break;	// (12)
    case 13: msg = std::string("Cannot Open ASCII File For Read");  break;	// (13)
    case 14: msg = std::string("Cannot Read From Input File");      break;	// (14)
    case 15: msg = std::string("Cannot Write To Output File");      break;	// (15)
    case 16: msg = std::string("Problems With Array Conversion");   break;	// (16)
    case 17: msg = std::string("Dimensioning Problems");            break;	// (17)
    case 18: msg = std::string("Division By Zero");                 break;	// (18)
    case 19: msg = std::string("Problems During Diagonalization");  break;	// (19)
    case 20: msg = std::string("Too Many Function Iterations");     break;	// (20)
    case 23: msg = std::string("Cannot Perform Division");          break;	// (23)
    case 24: msg = std::string("Function Not Yet Implemented");     break;	// (24)
    case 25: msg = std::string("Function Not Fully Implemented?");  break;	// (25)
    case 26: msg = std::string("Binary File Stream Not Open?");  break;	// (26)
    case 49: msg = std::string("Problems With Block Put");          break;   // (49)
    case 50: msg = std::string("Conversion To Unknown Type");       break;   // (50)
    case 51: msg = std::string("Cannot Continue Calculations");     break;   // (81)
    case 52: msg = std::string("Fundamental GAMMA Problems!");      break;   // (52)
    case 53: msg = std::string("Please Report This Trouble");       break;   // (53)
    case 54: msg = std::string("http://gamma.magnet.fsu.edu");      break;   // (54)
    case 55: msg = std::string("Unlikely Error! Bad Input?");       break;   // (55)
    case 56: msg = std::string("Unlikely Error! Check Indices?");   break;   // (56)
    default: msg = std::string("Unknown Error", noret);             break;	// (-1)
    }
  MxModError(hdr, msg, noret);
  }  

  
void MxModError(const std::string& hdr, int eidx, const std::string& pname, int noret)
  { 
  std::string msg;
  switch(eidx)
    {
    case 1:  msg = std::string("Problems With File ");		break;	// (1)
    case 2:  msg = std::string("Cannot Read Parameter ");		break;	// (2)
    case 3:  msg = std::string("Invalid Use Of Function ");		break;	// (3)
    case 4:  msg = std::string("Use Of Deprecated Function ");	break;	// (4)
    case 5:  msg = std::string("Please Use Class Member Function ");	break;	// (5)
    case 6:  msg = std::string("Failed To Perform ");		break;	// (6)
    case 7:  msg = std::string("Sorry, No Implementation Of ");	break;	// (7)
    default: msg = std::string("Unknown Error - ");			break;	// (-1)
    }
  MxModError(hdr, msg+pname, noret);
  }  
 
volatile void MxModFatal()
  {
   std::cout << std::endl;				//   Add a linefeed
   exit(-1);						//   Abort the program
  }

// ____________________________________________________________________________ 
// B             Matrix Module Generic String Formatting Functions
// ____________________________________________________________________________
 
/* The functions dec and form are C functions which used to be also defined in
   the Gnu C++ library (g++).  However, these don't seem to be set up for
   strings in ANSI C++ so they are redefined here. These are given the names
   MxModdec & MxModform so that the GAMMA matrix module is standalone.       */

std::string MxModdec(int i)
  { 
	char buffer[201]; 

#ifdef _MSC_VER
	sprintf_s(buffer, 200, "%d", i); 
#else
	sprintf(buffer, "%d", i); 
#endif

	return std::string(buffer); 
	}
std::string MxModdec(const std::string& fmt, int i)
  { 
	char buf[201]; 

#ifdef _MSC_VER
	sprintf_s(buf, 200, fmt.c_str(),i); 
#else
	sprintf(buf, fmt.c_str(),i); 
#endif

	return std::string(buf); 
	}
std::string MxModdec(int i, int digs)
  { 
	std::string fmt=std::string("%")+MxModdec(digs)+std::string("d"); 
	return MxModdec(fmt,i); 
	}
std::string MxModform(const std::string& fmt, double d)
  { 
	char buf[201]; 

#ifdef _MSC_VER
	sprintf_s(buf, 200, fmt.c_str(), d); 
#else
	sprintf(buf, fmt.c_str(), d);
#endif

	return std::string(buf); 
	}

   
#endif							// MxModBas.cc
