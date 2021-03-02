/* ML5Tag.cc ****************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**	MATLAB Tag (MAT Version 5)                 Implementation	**
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
** This MATLAB class provides functions for the input and output of     **
** MAT file version 5 tags assicated with a data element. MATLAB is a 	**
** commercial product and can be purchaced from The MathWorks, Inc.  	**
** This  module  deals with binary MAT files version 5 circa 1999.      **
**                                                                      **
*************************************************************************/

#ifndef _ML5TAG_CC_				// Is file already included?
#define _ML5TAG_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5Tag.h>			// Include the header file
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <stdio.h>
#include <time.h>				// Include time/date functions
#if !defined(_MSC_VER) && !defined(__MINGW32__) // MSVC++ & MinGW do not use this
# include <sys/utsname.h>                       // Include platform information
#endif
#include <Matrix/matrix.h>			// Include GAMMA matrices


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      MATLAB MAT 5 Tag Error Handling
// ____________________________________________________________________________
 
/*              Input           ML5T    : MAT version 5 header
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
		Output		none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void MatLab5Tag::MLT5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Tag");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, eidx, noret); break;   // Program Aborting     (0)
    case 1: GAMMAerror(hdr, eidx, noret); break;   // Input FileStream Bad (1)
    case 10: GAMMAerror(hdr, "Unknown Data Type Specified", noret);      //(10)
    case 11: GAMMAerror(hdr, "Setting Data Type To Double", noret);      //(11)
    case 12: GAMMAerror(hdr, "Non-positive Data Size Set", noret);       //(12)
    case 13: GAMMAerror(hdr, "Unreasonably Large Data Size", noret);     //(13)
    case 14: GAMMAerror(hdr, "Compressed Tag With >4 Bytes", noret);     //(14)
    case 20: GAMMAerror(hdr, "Error During Construction", noret);        //(20)

    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5Tag::MLT5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Tag");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5Tag::MLT5fatality(int eidx)
  {
  MLT5error(eidx, 1);
  if(eidx) MLT5error(0);
  GAMMAfatal();					// Clean exit from program
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                    MATLAB MAT 5 Tag Constructors
// ____________________________________________________________________________


MatLab5Tag::MatLab5Tag()
  {
  for(int i=0; i<8; i++) EightBytes[i] = '-'; 	// NULL the bytes
  DType = 0;					// Data type (unknown)
  NBytes = 0; 					// Number of bytes (none)
  Compressed = 0;				// Not compressed
  }


MatLab5Tag::MatLab5Tag(int TY, int NB, int CMP)
  {
int MLMaxBytes = 2000000000;
  longchars LC;
  LC.longval = TY;
  EightBytes[0] = LC.chars[0];
  EightBytes[1] = LC.chars[1];
  EightBytes[2] = LC.chars[2];
  EightBytes[3] = LC.chars[3];
  LC.longval = NB;
  EightBytes[4] = LC.chars[0];
  EightBytes[5] = LC.chars[1];
  EightBytes[6] = LC.chars[2];
  EightBytes[7] = LC.chars[3];
  DType = TY;				// Data type (unknown)
  if(TY<1 && TY>14)			// Odd construction as we don't know
    { 					// what this data type is
    MLT5error(10, 1);			//   Unknown data type
    MLT5error(20, 1);			//   Error during construction
    MLT5error(11);			//   Set data type default to double
    TY = 9;				// Set the type to double
    }
  NBytes = NB; 				// Number of bytes (none)
  if(NB<1)				// Trying to construct with <1 byte!
    {
    MLT5error(12, 1);			//   Cannot construct with <1 byte
    MLT5fatality(20);			//   Fatal during construction
    }
  if(NB>MLMaxBytes)			// Cannot construct a tag
    { 					// with more than MLMaxBytes
    MLT5error(13, 1); 			//   Too large of data size specified!
    MLT5fatality(20);			//   Fatal during construction
    }
  Compressed = CMP;			// Set compressed flag
  if(CMP && NB>4)			// Cannot construct compressed
    {					// data tag with more than 4
    MLT5error(14, 1);			//   bytes present in data!
    MLT5fatality(20);			//   Fatal during construction
    }
  }

MatLab5Tag::MatLab5Tag(const MatLab5Tag& T)
  {
  DType  = T.DType;				// Copy data type
  NBytes = T.NBytes;				// Copy # of bytes
  Compressed = T.Compressed;			// Copy compressed flag
  for(int i=0; i<8; i++) 			// Copy all 8 bytes
    EightBytes[i] = T.EightBytes[i];
  }

void MatLab5Tag::operator=(const MatLab5Tag& T)
  {
  DType  = T.DType;				// Copy data type
  NBytes = T.NBytes;				// Copy # of bytes
  Compressed = T.Compressed;			// Copy compressed flag
  for(int i=0; i<8; i++) 			// Copy all 8 bytes
    EightBytes[i] = T.EightBytes[i];
  }

// ____________________________________________________________________________
// B                    MATLAB MAT 5 Tag Access Functions
// ____________________________________________________________________________


int MatLab5Tag::Type(int LE) const

        // Input 		ML5T    : MAT version 5 tag (this)
	//			LE      : Little Endian flag
        // Output               type	: An integer for the data type
	// Note				: Byte swapping done if the endian
	//				  specified doesn't match here

  {
  longchars LC;
  int i;
  if(!LE &&  this->IsBigEndian())
    for(i=0; i<4; i++) LC.chars[i] =  EightBytes[3-i];
  else
    for(i=0; i<4; i++) LC.chars[i] =  EightBytes[i];
  return int(LC.longval);
  }


std::string MatLab5Tag::DataType() const

        // Input 		ML5T    : MAT version 5 tag (this)
        // Output               string  : String labeling the type

  {
  switch(DType)
    {
    case 0: return std::string("empty"); break;
    case 1: return std::string("8 bit, signed"); break;
    case 2: return std::string("8 bit, unsigned"); break;
    case 3: return std::string("16 bit, signed"); break;
    case 4: return std::string("16 bit, unsigned"); break;
    case 5: return std::string("32 bit, signed"); break;
    case 6: return std::string("32 bit, unsigned"); break;
    case 7: return std::string("IEEE 754 single format"); break;
    case 8: return std::string("Reserved"); break;
    case 9: return std::string("IEEE 754 double format"); break;
    case 10: return std::string("Reserved"); break;
    case 11: return std::string("Reserved"); break;
    case 12: return std::string("64 bit, signed"); break;
    case 13: return std::string("64 bit, unsigned"); break;
    case 14: return std::string("MATLAB array"); break;
    }
  return std::string("Unknown");
  }


std::string MatLab5Tag::DataSymbol() const

        // Input 		ML5T    : MAT version 5 tag (this)
        // Output               string  : String labeling the type
	//				  according to MATLAB symbol

  {
  switch(DType)
    {
    case 0: return std::string("none"); break;
    case 1: return std::string("miINT8"); break;
    case 2: return std::string("miUINT8"); break;
    case 3: return std::string("miINT16"); break;
    case 4: return std::string("miUINT16"); break;
    case 5: return std::string("miINT32"); break;
    case 6: return std::string("miUINT32"); break;
    case 7: return std::string("miSINGLE"); break;
    case 8: return std::string("--"); break;
    case 9: return std::string("miDOUBLE"); break;
    case 10: return std::string("--"); break;
    case 11: return std::string("--"); break;
    case 12: return std::string("miINT64"); break;
    case 13: return std::string("miUINT64"); break;
    case 14: return std::string("miMATRIX"); break;
    }
  return std::string("Unknown");
  }


int MatLab5Tag::Bytes() const { return NBytes; }

        // Input 		ML5T    : MAT version 5 tag (this)
        // Output               NBytes	: Number of bytes 


void MatLab5Tag::SetCompressed(int bigend)

        // Input 		ML5T    : MAT version 5 tag (this)
	//			bigend  : Data storage endian type
        // Output               void	: Determines if input tag has been
	//				  read in compressed format or not 
	//				  and sets flag Compressed accordingly
	// Note				: For a non-compressed tag, a 4 integer	
	//				  [0, 14] for the data type is written.
	//				  Thus, part of these bytes will be 
	//				  zero.  For big endian,bytes 3&4
	//				  are zero & for little endian bytes
	//				  1&2 are zero if tag is NOT compressed

  { 
  shortchars IC;			// We check for 0 in two of the bytes
  if(bigend)				// If data in was big endian, then the
    {					// 3rd & 4th bytes will have the data
    IC.chars[0] = EightBytes[1];	// type leaving bytes 1&2 empty (0)
    IC.chars[1] = EightBytes[0];	// the tag isn't compressed, if it
    }
  else					// If data in was little endian, then
    {					// the 1st & 2nd bytes will have the
    IC.chars[0] = EightBytes[2];	// data type in them, leaving the 3rd &
    IC.chars[1] = EightBytes[3];	// 4th bytes empty (0) if uncompressed
    }					
  (IC.shortval)?Compressed=1:Compressed=0;
  }


// ____________________________________________________________________________
// C                  MATLAB MAT 5 Tag Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a tag in MATLAB ".mat" binary
   format, version 5.  Note that in a valid MATLAB MAT file, this tag will
   be associated with a "Data Element" or a "Sub-Element".
 
                Input           ML5H    : MAT version 5 tag (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5H written to it
                Note                    : No care is taken to insure the
                                          tag is written tothe top of the
                                          file (as it should be in MATLAB)   */

int MatLab5Tag::write(std::fstream& fp) const                                       
  {
  if(!Compressed)				// Output format is simple if
    {						// the tag isn't compressed.
    int32_t LI = DType;				// Two int32_t integers (8 bytes),
    fp.write((char*)&LI, sizeof(int32_t));		// data type followed by the
    LI = NBytes;				// size of data in bytes
    fp.write((char*)&LI, sizeof(int32_t));
    }
  else						// If the tag is compressed
    {						// we just write two short
    short NB = NBytes; 				// data followed by the data
    short DT = DType; 				// 4 bytes which MUST BE
    if(WeRBigEnd()) 				// WRITTEN by whoever called us
      {
      fp.write((char*)&NB, sizeof(short));		// If we write big endian then
      fp.write((char*)&DT, sizeof(short));		// output shorts NB DT
      }
    else					// If we write small endian
      {						// then output shorts DT NB
      fp.write((char*)&DT, sizeof(short));
      fp.write((char*)&NB, sizeof(short));
      }
    }
  return 1;
  }

  
/* These are the functions which will return the number of bytes that are 
   written upon output a tag in MATLAB ".mat" binary format, Ver. 5          */ 

int MatLab5Tag::Size() const { return 8; }	// Bytes used in a file


// ____________________________________________________________________________
// D                MATLAB MAT 5 Tag Binary Input Functions
// ____________________________________________________________________________
 
/* These are the functions which will input a tag in MATLAB ".mat" binary
   format, version 5.  Note that in a valid MATLAB MAT file, this tag will
   associated with a "Data Element", a combination of Tag + Data 

		Input           ML5T    : MAT version 5 tag (this)
           			fp	: Input file stream
				bigend  : Data storage endian type
                                warn     : Warning output level
                                                0 = no warnings
                                                1 = warnings
                                               >1 = fatal warnings
		Output		void    : ML5T is set from values found
					  in the input stream
		Note			: No care is taken to insure the
					  tag is read from the top of a data
					  element (as it should be in MATLAB)*/
 
int MatLab5Tag::read(std::fstream& fp, int bigend, int warn)

  {
  for(int i=0; i<8; i++)			// Just input 8 char bytes
    { 						// as the tag text field
    if(!fp || fp.eof())				// Insure filestream is open
      {						// Here if bad filestream
      if(warn)					// We'll output a warning
        {					// if needed and even die
        MLT5error(1, 1);			// if the error is fatal
        if(warn>1) MLT5fatality(11);
        else       MLT5error(11,1);
        }
      return 0;
      }
    fp.read(&(EightBytes[i]), sizeof(char));
    }
  SetCompressed(bigend);			// See if tag input compressed
  int swap = 0;					// Assume no byte swapping,
  if(bigend != int(WeRBigEnd())) swap++;		// but check and see
  if(Compressed)				// See if tag indicates the
    {						// data element is comrpessed
    shortchars IC; 				// If so, NBytes will be in the
    IC.chars[0] = EightBytes[0];		// first two bytes and the data
    IC.chars[1] = EightBytes[1];		// type will be in the 2nd pair
    if(swap) Swap(IC.shortval);			// of bytes
    if(bigend) NBytes = IC.shortval;
    else       DType  = IC.shortval; 
    IC.chars[0] = EightBytes[2];
    IC.chars[1] = EightBytes[3];
    if(swap) Swap(IC.shortval);
    if(bigend) DType  = IC.shortval;
    else       NBytes = IC.shortval;
    }
  else
    {
    longchars LC;				// Now we clip out the data
    LC.chars[0] = EightBytes[0];		// type. This tells us how
    LC.chars[1] = EightBytes[1];		// the data is stored. Most
    LC.chars[2] = EightBytes[2];		// times this is MATLAB form,
    LC.chars[3] = EightBytes[3];		// indicated by type 14
    if(swap) Swap(LC.longval);
    DType = LC.longval;
    LC.chars[0] = EightBytes[4];		// Next we clip out the # of
    LC.chars[1] = EightBytes[5];		// bytes stored in the data
    LC.chars[2] = EightBytes[6];		// element.
    LC.chars[3] = EightBytes[7];
    if(swap) Swap(LC.longval);
    NBytes = LC.longval;
    }
  return 1;
  } 

// ____________________________________________________________________________
// E                  MATLAB MAT 5 Tag ASCII Output Functions
// ____________________________________________________________________________
 

void MatLab5Tag::print(std::ostream& ostr, int hf, int bf) const

	// Input 		ML5T    : MAT version 5 tag (this)
	// 			ostr	: An output stream
	//			hf 	: Flag to print header   
	//                      bf      : Flag to print tag chars
	// Output		ostr	: Output stream with MATLAB MAT file
	//				  tag version 5 info placed into it
	// Note				: This information exists at the
	//				  beginning of each data element

  {
  if(hf) ostr << "\n\t\tTag             ";		// 8 byte field
  if(bf)						// If desired then
    for(int i=0; i<8; i++) ostr << EightBytes[i]; 	// output all bytes
  if(Compressed)     ostr << "(Compressed)";
  ostr << "\n\t\t  Data Type:    " << DType 		// The Data type
       << " - " << DataSymbol() << " " << DataType();
  ostr << "\n\t\t  No. Bytes:    " << NBytes;		// # of Bytes
  ostr.flush();
  }


std::ostream& operator<< (std::ostream& ostr, const MatLab5Tag& MLT5)

        // Input                ostr	: An output stream
	// 			ML5T    : MAT version 5 tag
        // Output               ostr	: The output stream modified by
        //                                the tag parameters

  { MLT5.print(ostr); return ostr; }


#endif							// ML5Tag.cc
