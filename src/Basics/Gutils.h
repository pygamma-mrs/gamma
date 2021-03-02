/* Gutils.h *****************************************************-*-c++-*-
**									**
**	                         G A M M A				**
**								 	**
**	GAMMA Utilities 	 	              Interface 	**
**								 	**
**	Copyright (c) 1990				 		**
**	Tilo Levante, Scott A. Smith		 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fur physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**									**
**      $Header: $
**                                                     			**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
**  This file contains some utilities which have a general use in	**
**  GAMMA programs.  There doesn't seem to be any other better place	**
**  to put them except in their own file.				**
** 									**
*************************************************************************/

#ifndef   Gutils_h_				// Is this file already included?
#  define Gutils_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// Then this is the interface
#  endif

#include <string>				// Include libstdc++ strings
#include <iostream>				// If we don't already know cout/cin
#include <GamGen.h>				// Know OS dependent stuff

// _____________________________________________________________________________
// A                   Funtions For Interactive Requests
// _____________________________________________________________________________

        // Input                argc : Integer, # of arguments in argv
        //                      argv : Character array of arguments/parameters
        //                       par : Argument/parameter index
        //                         Q : A string which asks for a parameter
        //                         V : A string/Double/Integer parameter value
        // Output               none : Value of parameter is set to the par
        //                             value in argv if it exists, or it
        //                             is set interactively after the user is
        //                             prompted by the string question.

MSVCDLL void query_parameter(int argc,char* argv[],
                                  int par,const std::string& Q, std::string& V);
MSVCDLL void query_parameter(int argc,char* argv[],
                                  int par,const std::string& Q, double&      V);
MSVCDLL void query_parameter(int argc,char* argv[],
                                  int par,const std::string& Q, int&         V);

MSVCDLL bool ask_set(int argc,char* argv[],int par,const std::string& Q, std::string& V);
MSVCDLL bool ask_set(int argc,char* argv[],int par,const std::string& Q, double&      V);
MSVCDLL bool ask_set(int argc,char* argv[],int par,const std::string& Q, int&         V);

// _____________________________________________________________________________
// B                GAMMA Generic Error Handling Functions
// _____________________________________________________________________________
 
/*         Input                hdr     : Class or Module header name
                               *msg     : The error message to output
                               *eidx    : The error index
	  		       *pname   : Additional error string
	  			noret   : Flag whether line feed needed
           Output               none    : An error message of some sort
				          is sent to standard output
					  Program terminated if fatal
				Note	: New C++ compilers replace volatile
					  --attribute--((fatal)	              */
 
MSVCDLL void GAMMAerror(const std::string& hdr, const std::string& msg,             int noret=0);
MSVCDLL void GAMMAerror(const std::string& hdr, int eidx,                           int noret=0);
MSVCDLL void GAMMAerror(const std::string& hdr, int eidx, const std::string& pname, int noret=0);
MSVCDLL volatile void GAMMAfatal();


#endif								// Gutils.h
