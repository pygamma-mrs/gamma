/* Gascii.cc ***************************************************-*-c++-*-*
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      GAMMA ASCII File Functions                   Implementation     **
**                                                                      **
**      Copyright (c) 1997                                              **
**      Dr. Scott A. Smith                                              **
**      1800 E. Paul Dirac Drive                                        **
**      National High Magnetic Field Laboratory                         **
**      Tallahassee FL 32306 USA                                        **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This module deals with I/O of ASCII files in GAMMA.  Although GAMMA **
**  data types often come with the ability to read themselves from set  **
**  parameter files (ASCII files in GAMMA parameter format), often one  **
**  needs to input a set of values (int, double, string, ....) that     **
**  control program function but have nothing to do with GAMMA objects. **
**  These functions simplify that process by allowing flexible parsing  **
**  of ASCII files.                                                     **
**                                                                      **
*************************************************************************/


#ifndef   Gascii_cc_			// Is file already included?
#  define Gascii_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <stdlib.h>
#include <GamIO/Gascii.h>		// Include the interface
#include <string>			// Know about libstc++ strings
#include <Basics/StringCut.h>		// Know about string cutting
#include <Basics/Gutils.h>		// Know about GAMMA error messages
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Matrix/matrix.h>		// Know about matrices
#include <iostream>			// Know about input/output file streams
#include <fstream>

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output
using std::ifstream;			// Using libstdc++ input file streams
using std::ios;				// Using libstdc++ file type settings

// ____________________________________________________________________________
// i                       ASCII File I/O Error Functions
// ____________________________________________________________________________

                                                                                
void ASCIIerr(int eidx, int noret)
  
        // Input 		eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

  {
  string hdr("ASCII File I/O");
  switch(eidx)
    {
    case 1:
    default: GAMMAerror(hdr, eidx, noret); break;
    }
  }  
 
 
void ASCIIdie(int eidx)
  
        // Input 		eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

  {
  ASCIIerr(eidx, 1);				// Output error message
  if(eidx) ASCIIerr(0);				// Write that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
 
void ASCIIerr(int eidx, const string& pname, int noret)
  
        // Input 		eidx    : Flag for error type
        //                      pname   : string included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
  {
  string hdr("ASCII File I/O");
  switch(eidx)
    {
    case 1:
    default: GAMMAerror(hdr,eidx,pname,noret); break;
    }
  }  
 


// ____________________________________________________________________________
// A                    ASCII File Reading & Parsing Of Strings
// ____________________________________________________________________________

 
string* Getstrings(const string& filename, int& N)

        // Input        filename: ASCII input file name
        // Return       Svect   : A string vector is filled with names
        //                        found in the file "filename"
        // Note                 : Assumes file contains 1 name per line
        // Note	           : Value of N set to # of strings read
        // Note                 : String array must be externally destruced

  {
  string input;                         // String for a whole line
  string datas;                         // String for read data
  ifstream inp;                         // Input stream
  int lines = LineCount(filename);      // Line count
  string* Svect;                        // String vector, 1 string per line
  Svect = new string[lines];            // Allocate string vector space
  lines = 0;                            // Reset line count
  inp.open(filename.c_str());           // Reopen the file
  while(Greadline(inp, input))          // Loop file, read lines
    Svect[lines++] = input;             // Store the value in vector
  N = lines;
  return Svect;
  }


// ----------------------------------------------------------------------------
//                    I/O To Obtain Doubles From An ASCII File
// ----------------------------------------------------------------------------

double* GetDoubles(const string& filename, int& N)

        // Input        filename: ASCII input file name
        // Return       Dvect   : A double vector is filled with values
        //                        found in the file "filename"
        // Note                 : Assumes file contains 1 value per line

  {
  int nlines = LineCount(filename);     // Line count
  double *Dvect = NULL;                 // Double vector, 1 value per line
  Dvect = new double[nlines];           // Allocate double vector space
  ifstream inp;                         // Input stream  lines = 0;                            // Reset line count
  inp.open(filename.c_str());           // Reopen the file
  int line=0;
  string input;                         // String for a whole line
  string datas;                         // String for read data
  while(Greadline(inp, input))          // Loop file, read lines
    {
    datas = cutDouble(input);           // Store the value in vector
    Dvect[line] = atof(datas.c_str());  // Convert from string to float
    line++;                             // Increment line count
    }
  N = nlines;
  return Dvect;
  }
 

row_vector GetDoubles(const string& filename)

        // Input        filename: ASCII input file name
        // Return       data	: A vector is filled with values
        //                        found in the file "filename"

  {                                                       
  int lines = LineCount(filename);             // Line count
  row_vector data(lines);                      // Vector for input data
  ifstream inp;                                // Input stream
  inp.open(filename.c_str());                  // Re-open the file
  if(!inp.good())                              // Insure file is still O.K.
    {
    ASCIIerr(1, filename, 1);                  //    Problems with file
    ASCIIdie(13);                              //    Can't Read, Quit
    }
  string input;                                // String for a whole line
  string datac;                                // String for read data
  double dataf;                                // Double for converted data
  lines = 0;
  while(Greadline(inp, input) !=0)
    {
    datac = cutDouble(input);                  // Get double value in line
    if(datac!="") dataf=atof(datac.c_str());   // Either set the value
    else          dataf = 0;                   // or make it zero
    data.put(dataf,lines);                     // Put data into vector
    lines++;                                   // Increment data count
    }
  return data;
  }  


void vxread(const string& filename, row_vector& vx)
 
     // Input           filename : ASCII input file name
     //                 vx       : A row vector
     // Return          none     : vx is filled with data in
     //                            file "filename"
 
  {
  string input;
  string datac;
  double datad=0;
  int pt = 0;
  ifstream inp;                     // Input stream
  inp.open(filename.c_str());       // Open the file
  if(!inp.good())                   // Insure file is O.K.
    {
    ASCIIerr(1, filename, 1);       //	Problems with file
    ASCIIdie(13);                   // 	Can't Read, Quit
    }
  while(Greadline(inp, input))
    {
    cutWhite(input);                // Remove blanks at line beginning
    datac = cutDouble(input);       // Cut out double precsion number
    datad = 0;
    if(datac!="")
      datad = atof(datac.c_str());
    if(pt<vx.size()) vx.put(datad, pt);
    pt++;
    }    
  }
 

// ____________________________________________________________________________
// G                       ASCII File Auxiliary Functions
// ____________________________________________________________________________
 

void ASCIItell(const string& filename, int verbose)
 
	// Input           filename	: ASCII input file name
        //                 verbose	: Flag how much info to output        
	// Return          none		: The ASCII file is probed to determine
	//                                what type of information it contains
 
  {
  ifstream inp;					// Input filestream
  string input;					// string for input line
  inp.open(filename.c_str());			// Open the file
  if(!inp.good())				// Insure file is O.K.
    {
    ASCIIerr(1, filename, 1);			//	Problems with file
    ASCIIdie(13);				// 	Can't Read, Quit
    }
  int lines = 0;				// Line counter
  while(Greadline(inp, input)!=0 )		// Count lines in input file
    lines++;
  inp.close();					// Close file for reopen
  cout << "\n\tASCII File " << filename
       << " Information\n";
  cout << "\n\t\tNumber of Lines: " << lines;
  if(verbose)
    {
    cout << "\n\t\tFile Data:       ";    
    inp.open(filename.c_str());			// Open the file
    if(!inp.good())				// Insure file is O.K.
      {
      ASCIIerr(1, filename, 1);			//	Problems with file
      ASCIIdie(13);				// 	Can't Read, Quit
      }
    int line=0;
    while(Greadline(inp, input)!=0 )		// Count lines in input file
      {
      cout << line << ". " << input
           << "\n\t\t                 ";    
      line ++;
      }
    }
  cout << "\n\n";
  }

int LineCount(const string& filename)
 
        // Input           filename	: ASCII input file name   
        // Return          int		: The nubmer of lines in the ASCII file
	//                                   named filename

  {
  ifstream inp;                         // Input stream
  inp.open(filename.c_str());           // Open the file
  if(inp.rdstate() & ios::failbit)      // Insure file is O.K.
    {
    ASCIIerr(1, filename, 1);           //	Problems with file
    ASCIIdie(13);                       // 	Can't Read, Quit
    }
  int lines = -1;                       // Line count
  char buf[200];			     // Buffer for line read
  while(!(inp.rdstate() & ios::failbit))
    {
    inp.getline(buf, 200);              // Read next line
    lines++;                            // to get a total line count
    }
  inp.close(); 
  return lines;
  }

// ____________________________________________________________________________
// G                     ASCII File Input Query Functions
// ____________________________________________________________________________
 

//string ASCask_read(int argc, char* argv[], int argn, const string& quest)
 
        // Input 		argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
	//			quest	: Question to ask interactively
        // Output               string  : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the spin system is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The spin system is modifed (filled)
 
//    {
//    string filename;                            // Name of spin system file
//    query_parameter(argc, argv, argn,           // Get filename from command
//       quest, filename); // Or ask for it
//    read(filename);                             // Read system from filename
//    return filename;
//    }

 
#endif  	                                           // Gascii.cc

