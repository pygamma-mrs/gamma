/* Gutils.cc ****************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**      GAMMA Utilities                            Implementation	**
**									**
**      Copyright (c) 1990						**
**      Tilo Levante, Scott A. Smith					**
**      Eidgenoessische Technische Hochschule   	                **
**      Labor fur physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/
 
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains some utilities which have a general use in       **
**  GAMMA programs.  There doesn't seem to be any other better place    **
**  to put them except in their own file.                               **
**                                                                      **
*************************************************************************/
 
#ifndef   Gutils_cc_			// Is this file already included?
#  define Gutils_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Basics/Gutils.h>		// Include the implementation
#include <GamGen.h>			// Include OS specifics
#include <string>			// Know about string class
#include <iostream>			// Know about endl
#include <stdlib.h>
#include <math.h>

using std::string;			// Using libstdc++ strings
using std::cout;			// Using libstdc++ standard output
using std::cin;				// Using libstdc++ standard input
using std::endl;			// Using libstdc++ line endings

// _____________________________________________________________________________
// A              Interactive Parameter Request Functions
// _____________________________________________________________________________
 
/* 	   Input                argc : The # of arguments/parameters in argv
                                argv : Character array of arguments/parameters
	  			 par : Argument/parameter index
	  		           Q : String which asks for a parameter
	  		           V : - string parameter value
	  		               - double parameter value
	  		               - int parameter value
           Output               none : Parameter V is set to the par
	  			       parameter in argv if it exists, or it
	  			       is set interactively after the user
	  			       is prompted by the question Q.         */

void query_parameter(int argc,char* argv[],int par,const string& Q,string& V)
  {
  if(argc>par) V = argv[par];		// If parameter input, set as string
  else { cout << Q; cin >> V; }		// else ask for it & read it
//while(cin.good()&&(cin.get()!='\n'));	// Clear out any linefeeds
  }

void query_parameter(int argc,char* argv[],int par,const string& Q,double& V)
  {
  if(argc>par) V = atof(argv[par]);	// If parameter input, set as float
  else { cout << Q; cin >> V; }		// else ask for it & read it
//while(cin.good()&&(cin.get()!='\n'));	// Clear out any linefeeds
  }

void query_parameter(int argc,char* argv[],int par,const string& Q,int& V)
  {
  if(argc>par) V = atoi(argv[par]);	// If parameter input, then set as int
  else { cout << Q; cin >> V; }		// else ask for it & read it
//while(cin.good()&&(cin.get()!='\n'));	// Clear out any linefeeds
  }

bool ask_set(int argc,char* argv[],int par,const string& Q,string& V)
  {
  if(par < 0) return false;		// Dont read if par negative index
  if(argc>par) V=string(argv[par]);	// If parameter input, set as string
  else 					// If parameter Not Input
    {
    char s[100];			//   A character array for input
    cout << Q; 				//   Ask for the value
    cin.getline(s, 100);		//   Read line containing the value
    string ss(s);			//   Convert to C++ string
    if(!ss.length()) return false;	//   Here if only a return was input
    V = s;				//   Set the variable to string
    }
  return true;
  }

bool ask_set(int argc,char* argv[], int par,const string& Q, int& V)
  {
  if(par < 0) return false;		// Dont read if par negative index
  if(argc>par) V = atoi(argv[par]);	// If parameter input, set as int
  else 					// If parameter Not Input
    {
    char s[100];			//   A character array for input
    cout << Q; 				//   Ask for the value
    cin.getline(s, 100);		//   Read line containing the value
    string ss(s);			//   Convert to C++ string
    if(!ss.length()) return false;	//   Here if only a return was input
    V = atoi(s);			//   Convert to integer point #
    }
  return true;
  }

bool ask_set(int argc,char* argv[], int par, const string& Q, double& V)
  {
  if(par < 0) return false;		// Dont read if par negative index
  if(argc>par) V = atof(argv[par]);	// If parameter input, set as float
  else 					// If parameter Not Input
    {
    char s[100];			//   A character array for input
    cout << Q; 				//   Ask for the value
    cin.getline(s, 100);		//   Read line containing the value
    string ss(s);			//   Convert to C++ string
    if(!ss.length()) return false;	//   Here if only a return was input
    V = atof(s);			//   Convert to floating point #
    }
  return true;
  }

// ____________________________________________________________________________ 
// B                GAMMA Generic Error Handling Functions
// ____________________________________________________________________________
 
/*
           Input                hdr     : Class or Module header name
                               *msg     : The error message to output
                               *eidx    : The error index
	  		       *pname   : Additional error string
	  			noret   : Flag whether line feed needed
           Output               none    : An error message of some sort
				          is sent to standard output          */

void GAMMAerror(const string& hdr, const string& msg, int noret)
  {
	  /*
string errstr;
errstr += string("\n");
errstr += hdr;
errstr += string(": ");
errstr += msg;
  cout << errstr; if(!noret) cout << endl;	// GAMMA error stuff
  */

  string s = string("\n") + hdr + string(": ") + msg;	// Set error string
  cout << s; if(!noret) cout << endl;	// GAMMA error stuff
  }  

void GAMMAerror(const string& hdr, int eidx, int noret)
  { 
  string msg;
  switch(eidx)
    {
    case 0:  msg = string("Program Aborting.....");              break; // (0)
    case 1:  msg = string("Problems With Input File Stream");    break; // (1)
    case 2:  msg = string("Problems With Output File Stream");   break; // (2)
    case 3:  msg = string("Can't Construct From Parameter Set"); break;	// (3)
    case 4:  msg = string("Cannot Construct From Input File");   break;	// (4)
    case 5:  msg = string("Cannot Write To Parameter File");     break;	// (5)
    case 6:  msg = string("Cannot Write To Output FileStream");  break;	// (6)
    case 7:  msg = string("Cannot Produce Any Output");		 break;	// (7)
    case 8:  msg = string("Problems With Parameter Set");        break;	// (8)
    case 9:  msg = string("Problems During Construction");       break;	// (9)
    case 10: msg = string("Cannot Access Internal Component");   break;	// (10)
    case 11: msg = string("Cannot Read From Input FileStream");  break;	// (11)
    case 12: msg = string("End of File Reached");		 break;	// (12)
    case 13: msg = string("Cannot Open ASCII File For Read");    break;	// (13)
    case 14: msg = string("Cannot Read From Input File");        break;	// (14)
    case 15: msg = string("Cannot Write To Output File");        break;	// (15)
    default: msg = string("Unknown Error", noret);               break;	// (-1)
    }
  GAMMAerror(hdr, msg, noret);
  }  
  
void GAMMAerror(const string& hdr, int eidx, const string& pname, int noret)
  { 
  string msg;
  switch(eidx)
    {
    case 1:  msg = string("Problems With File ");		break;	// (1)
    case 2:  msg = string("Cannot Read Parameter ");		break;	// (2)
    case 3:  msg = string("Invalid Use Of Function ");		break;	// (3)
    case 4:  msg = string("Use Of Deprecated Function ");	break;	// (4)
    case 5:  msg = string("Please Use Class Member Function ");	break;	// (5)
    case 6:  msg = string("Please Use Modern GAMMA Function ");	break;	// (6)
    case 7:  msg = string("Parameter ") + pname 
                 + string(" Is The Culprit!");			break;	// (7)
    default: msg = string("Unknown Error - ");			break;	// (-1)
    }
  GAMMAerror(hdr, msg+pname, noret);
  }  
 
volatile void GAMMAfatal()
  {						
     cout << endl;					//   Add a linefeed
     exit(-1);						//   Abort the program
  }

#endif							// Gutils.cc
