/* ML5AF.cc *****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**	MATLAB Array Flags (MAT Version 5)        Implementation	**
**                                                                      **
**      Scott Smith                                                     **
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
** This MATLAB class provides functions for the input and output of a   **
** MATLAB array flags sub-element in MAT version 5 files.  MATLAB is a 	**
** commercial product and can be purchaced from The MathWorks, Inc.  	**
** This  module  deals with binary MAT files version 5 circa 1999.      **
**                                                                      **
*************************************************************************/

#ifndef   GML5AF_CC_				// Is file already included?
#  define GML5AF_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/ML5AF.h>			// Include the header file
#include <GamIO/ML5SubE.h>			// Include base class header
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i             MATLAB MAT 5 Array Flags Sub-Element Error Handling
// ____________________________________________________________________________

/*              Input           MLAF5   : MAT version 5 sub-element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5AF::MLAF5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Array Flags Sub-Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Read Data Element Header", noret); // (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5AF::MLAF5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Array Flags Sub-Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5AF::MLAF5fatality(int eidx)
  {
  MLAF5error(eidx, 1);
  if(eidx) MLAF5error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              MATLAB MAT 5 Array Flags Sub-Element Constructors
// ____________________________________________________________________________


MatLab5AF::MatLab5AF() : MatLab5SE() { Class=0; Fs=0; }
MatLab5AF::MatLab5AF(int cmplx, int C, int global, int logical) :MatLab5SE()
  {
              Fs  = 0;			// Set integer equivalent
  if(cmplx)   Fs += 8;			// of flags byte
  if(global)  Fs += 4;
  if(logical) Fs += 2;
              Class = C;		// Store integer equal of class byte
  MLData = new char[8];			// Allocate new data array
  MLData[2] = char(Fs);			// Set flags array (3rd byte)
  MLData[3] = char(Class);		// Set array type (4th byte)
  }

MatLab5AF::MatLab5AF(const matrix& mx, int cmplx) :MatLab5SE()

        // Input                ML5AF	: MAT version 5 array flags SE (this)
	//			mx	: A GAMMA matrix
	//			cmplx 	: Flag for whether real vs complex
	// Output		this	: Constructed appropriately
	//				  for the array flags sub-element
	//				  of a GAMMA array
	// Note				: The tag is set to be of type
	//				  6 (i.e. MiUINT32) and the size
	//				  is set to eight bytes long
	// Note				  The data is set to be of type
	//				  6 (i.e. mxDOUBLE_CLASS) and the
	//				  output size is eight bytes as is
	//				  indicated in the preceding tag

  {
  MLTag = MatLab5Tag(6,8);		// Set up the tag
  if(cmplx) Fs = 8; 			// Set integer equivalent
  else	    Fs = 0; 			// of flags byte
  Class     = 6;			// Store integer equal of class byte
  MLData    = new char[8];		// Allocate new data array
  MLData[2] = char(Fs);			// Set flags array (3rd byte)
  MLData[3] = char(Class);		// Set array type (4th byte)
  }
 
// ____________________________________________________________________________
// B           MATLAB MAT 5 Array Flags Sub-Element Access Functions
// ____________________________________________________________________________
   

std::string MatLab5AF::SClass() const

        // Input                ML5AF	: MAT version 5 array flags SE (this)
        // Output               string  : String labeling the type

  {
  switch(Class)
    {
    case 0: return std::string("None"); break;
    case 1: return std::string("Cell Array"); break;
    case 2: return std::string("Structure"); break;
    case 3: return std::string("Object"); break;
    case 4: return std::string("Character Array"); break;
    case 5: return std::string("Sparse Array"); break;
    case 6: return std::string("Double Precision Array"); break;
    case 7: return std::string("Single Precision Array"); break;
    case 8: return std::string("8 Bit, Signed Integer"); break;
    case 9: return std::string("8 Bit, Unsigned Integer"); break;
    case 10: return std::string("16 Bit, Signed Integer"); break;
    case 11: return std::string("16 Bit, Unsigned Integer"); break;
    case 12: return std::string("32 Bit, Signed Integer"); break;
    case 13: return std::string("32 Bit, Unsigned Integer"); break;
    }
  return std::string("Unknown");
  }
   

std::string MatLab5AF::Symbol() const

        // Input                ML5AF	: MAT version 5 array flags SE (this)
        // Output               string  : String labeling the type

  {
  switch(Class)
    {
    case 0: return std::string("None"); break;
    case 1: return std::string("mxCELL_CLASS");    break;
    case 2: return std::string("mxSTRUCT_CLASS");  break;
    case 3: return std::string("mxOBJECT_CLASS");  break;
    case 4: return std::string("mxCHAR_CLASS");    break;
    case 5: return std::string("mxSPARSE_CLASS");  break;
    case 6: return std::string("mxDOUBLE_CLASS");  break;
    case 7: return std::string("mxSINGLE_CLASS");  break;
    case 8: return std::string("mxINT8_CLASS");    break;
    case 9: return std::string("mxUINT8_CLASS");   break;
    case 10: return std::string("mxINT16_CLASS");  break;
    case 11: return std::string("mxUINT16_CLASS"); break;
    case 12: return std::string("mxINT32_CLASS");  break;
    case 13: return std::string("mxUINT32_CLASS"); break;
    }
  return std::string("Unknown");
  }


int MatLab5AF::IsComplex() const { return IsFlag(4); } 

        // Input                ML5AF	: MAT version 5 array flags SE (this)
        // Output               TF      : True if complex flag set
	// Note				: The flags in this sub-element
	//				  are stored in the 3rd (index 2)
	//				  byte of the data.  The 5th bit
	//				  (index 4) is the complex flag

/******************************************************************************
  sosi - will want to replace this when libstdc++ bitset class works!
         or perhaps just get rid of it entirely!
******************************************************************************/

int MatLab5AF::IsFlag(int flg) const
  {
  int Flags[8];
  register int t;
  unsigned u = Fs;
  int i=0;
  for(t=128; t>0; t/=2)
    {
    if(u & t) Flags[i] = 1;
    else      Flags[i] = 0;
    i++;
    }
  return Flags[flg];
  }
 

// ____________________________________________________________________________
// C       MATLAB MAT 5 Array Flags Sub-Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output this sub-element in MATLAB ".mat"
   binary format, version 5.  This is done with a pair {tag,data}.  Both the
   tag & data will be written directly from these functions, bypassing the
   associated class structures which have no knowledge of GAMMA array types.

   Remember, the MATLAB MAT version 5 tag for this sub-element contains
   8-bytes with two important settings: 1.) data type & 2.) # of bytes.
   This will be written in compact form, i.e. the sub-element containing a
   total of eight bytes.  The first two bytes contain the # of bytes, the 2nd
   two bytes contain the data type, and  is placed in the first 4 bytes and 
   the number of bytes placed into the last 4 bytes (for non-compressed data 
   elements which we will always output...).  In this case, we output the 
   elements as doubles so that the data type should be MATLAB type miDOUBLE, 
   i.e. value 9.

                Input           ML5Re   : MAT version 5 reals SE (this)
                                fp      : Output file stream
                                mx      : GAMMA matrix
                                rv      : GAMMA row vector
                                cv      : GAMMA column vector
                Output          void    : Output file stream is modified
                                          by having the ML5Re written to it
                Note                    : No care taken herein to insure the
                                          sub-element is properly written.   */
   
int MatLab5AF::write(std::fstream& fp) const
  {
  MLTag.write(fp);				// First write the tag
  char Data[8];
  if(WeRBigEnd())				// Set array type (class) and
    {						// array flags.  These take up
    Data[3] = char(Class);			// one byte each & are witten
    Data[2] = char(Fs);			// to file together as a 4 byte
    }						// integer.  So we may have to
  else						// untangle it if byte swapped
    {
    Data[0] = char(Class);			// one byte each & are witten
    Data[1] = char(Fs);			// to file together as a 4 byte
    }
  fp.write((char*)&Data, 8*sizeof(char));
  return 1;
  } 

/* These are the functions which will return the number of bytes that are
   written upon output this sub-element in MATLAB ".mat" binary format, V5   */

int MatLab5AF::Size() const { return 2*MLTag.Size(); }


// ____________________________________________________________________________
// D               MATLAB MAT 5 Array Flags Sub-Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input an Array Flags sub-element in
   MATLAB ".mat" binary format, version 5.  Note that in a valid MATLAB MAT 
   file, a "Array Flags Sub-Element" is a combination of Tag + Data. The data
   will be 8 bytes where the 3rd byte are the Flags and the 4th byte the Class.
 
                Input           ML5AF	: MAT version 5 array flags SE (this)
                                fp      : Input file stream
				bigend  : Data storage endian type 
                                warn    : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5AF is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          tag is read from the top of a data
                                          element (as it should be in MATLAB)*/

int MatLab5AF::read(std::fstream& fp, int bigend, int warn)
  {
  if(!MatLab5SE::read(fp, bigend, (warn)?1:0))	// Read sub-element
    {                                           // If this has problems, then
    if(warn)                                    // we must deal with them
      {                                         // appropriately now
      if(warn==1) MLAF5error(10, 1);
      else        MLAF5fatality(10);
      }
    return 0;
    } 						// Now set array type (class) &
  if(bigend) 					// array flags.  These take up
    { 						// 1 byte each & may need swap
    Class = int(MLData[3]); 			// Here if big -> little/big
    Fs    = int(MLData[2]); 			// (this is tested & works)
    }
  else 						// Here if little -> little/big
    { 						// (this is tested & works)
    Fs    = int(MLData[1]);
    Class = int(MLData[0]);
    }
  return 1;
  }
 
// ____________________________________________________________________________
// E         MATLAB MAT 5 Array Flags Sub-Element ASCII Output Functions
// ____________________________________________________________________________

void MatLab5AF::print(std::ostream& ostr) const

        // Input                ML5AF	: MAT version 5 sub-element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                sub-element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each sub-element

 {
 ostr << "\n\t\tSE Array Flags  ";		// Output header
 MLTag.print(ostr, 0);				// Output tag values
 ostr << "\n\t\t  Array Type:   " << Class
      << " - " << SClass();

/******************************************************************************
  sosi - will want to replace this when libstdc++ bitset class works!
******************************************************************************/

if(Class)
{
int Flags[8];
register int t;
unsigned u = MLData[2];
int i=0;
for(t=128; t>0; t/=2)
  {
  if(u & t) Flags[i] = 1;
  else      Flags[i] = 0;
  i++;
  }

/******************************************************************************
******************************************************************************/

   ostr << "\n\t\t  Real/Cmplx:   ";
   if(Flags[4]) ostr << "Complex";
   else         ostr << "Real";
   ostr << "\n\t\t  Global/Local: ";
   if(Flags[5]) ostr << "Global";
   else         ostr << "Local";
   ostr << "\n\t\t  Logical:      ";
   if(Flags[6]) ostr << "Yes";
   else         ostr << "No";
   }
 }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab5AF& ML5AF)
 
        // Input                ostr    : An output stream
        //                      ML5AF	: MAT version 5 sub-element
        // Output               ostr    : The output stream modified by
        //                                the sub-element parameters

  { ML5AF.print(ostr); return ostr; }


#endif							// ML5AF.cc
