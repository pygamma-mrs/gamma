/* ML5DElem.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**     MATLAB Data Element (MAT Version 5)        Implementation	**
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
** MATLAB data element in MAT version 5 files.  MATLAB is a commercial  **
** product and can be purchaced from The MathWorks, Inc.		**
** This module deals with binary MAT files version 5 circa 1999.      	**
**                                                                      **
*************************************************************************/

#ifndef _ML5DE_CC_				// Is file already included?
#define _ML5DE_CC_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <GamIO/ML5DElem.h>			// Include the header file
#include <GamIO/ML5SubE.h>			// Include base class header
#include <Basics/Gutils.h>                      // Include GAMMA errors
#include <string>				// Include libstdc++ strings
#include <Basics/StringCut.h>	     // Include Gdec function


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              MATLAB MAT 5 Data Element Error Handling
// ____________________________________________________________________________

/*              Input           MLDE5   : MAT version 5 data element
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void MatLab5DE::MLDE5error(int eidx, int noret)
  {                                                     
  std::string hdr("MATLAB MAT V5 Data Element");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 1: GAMMAerror(hdr, "End of File Reached Before Data Found", noret);
      break;                                     //                        (1)
    case 10: GAMMAerror(hdr, "Cannot Read Data Element Header", noret);	// (10)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

 
void MatLab5DE::MLDE5error(int eidx, const std::string& pname, int noret)
  {
  std::string hdr("MATLAB MAT V5 Data Element");
  std::string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr,   1, pname, noret); break;	// File Problems   (1)
    case 2:  GAMMAerror(hdr,   2, pname, noret); break;	// !Read SingPar   (2)
    default: GAMMAerror(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
     
 
volatile void MatLab5DE::MLDE5fatality(int eidx)
  {
  MLDE5error(eidx, 1);
  if(eidx) MLDE5error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  MATLAB MAT 5 Data Element Constructors
// ____________________________________________________________________________


MatLab5DE::MatLab5DE() { BigEndIn = 1; }
MatLab5DE::MatLab5DE(const matrix& mx, const std::string& Name, int cmplx)
          :Tag(14, DataSize(mx,Name,cmplx)), AF(mx,cmplx),
           DA(mx), AN(Name), RData(mx), IData(mx) { MX = mx; }
MatLab5DE::~MatLab5DE() { }
 

// ____________________________________________________________________________
// B                MATLAB MAT 5 Data Element Access Functions
// ____________________________________________________________________________
   

void MatLab5DE::whos(std::ostream& ostr, std::fstream& fp, int bigend)

        // Input                ML5DE	: MAT version 5 data element (this)
        //                      ostr    : An output stream
	//			fp	: A file stream
	//			bigend	: Stored data endian type	
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 5 info placed
	//				  into it in MATLAB "whos" fashion
        // Note                         : It is assumed that fp points to the
        //                                beginning of a data element and that
	//				  it is open for reading.
	// Note				: The file stream position will be
	//				  advanced to the end of the element

  {
  if(!fp) return;			// No output if no filestream
  int pos = fp.tellp();			// Get current file position
  int TF = Tag.read(fp, bigend, 1);	// Read data element tag
  TF *= AF.read(fp, bigend, 1);		// Read array flags sub-element
  TF *= DA.read(fp, bigend, 1);		// Read dim. array sub-element
  TF *= AN.read(fp, bigend, 1);		// Read array name sub-element
  ostr << "\n" << std::string(12, ' ');	// Start new line & output margin
  std::string name = AN.Name();		// Get element array name
  int len = name.length();		// Get array name length
  if(len > 10)				// Output the array name
    ostr << name.substr(0,7) << "...";	// within the 9 columns in whos
  else					// output allocated for it
    ostr << name << std::string(10-len, ' ');
  int ND = DA.Dims();			// Get number of dimensions
  std::string dimstr ;			// String for dimension output
  ostr << "  ";				// Output a spacer
  for(int i=0; i<ND; i++)		// Output all of the individual
    {					// array dimensions
    dimstr += Gdec(DA.Dim(i));
    if(i<ND-1) dimstr += std::string(" x ");
    }
  len = dimstr.length();		// Get string length
  ostr << dimstr << std::string(10-len, ' ');// Output dimensions
  TF = Tag.read(fp, bigend, 1);		// Read data element tag
  int NB = Tag.Bytes();			// Get the bytes
  if(AF.IsComplex()) NB *= 2;		// Twice a many if complex
  std::string byt = Gdec(NB);			// Set bytes in string form
  len = byt.length();			// Length of bytes string
  ostr << "  ";				// Output a spacer
  ostr << byt << std::string(10-len, ' ');	// Output the bytes
  ostr << "  ";				// Output a spacer
  ostr << AF.SClass();			// Ouput the MATLAB class type
  if(AF.IsComplex())
    ostr << " (complex)";		// Ouput the MATLAB class type
  fp.seekp(pos);			// Set original current file position
  }


int MatLab5DE::Size(std::fstream& fp, int bigend)
  {
  if(!fp) return 0;			// No output if no filestream
  int pos = fp.tellp();			// Get current file position
  Tag.read(fp, bigend, 1);		// Read data element tag
  fp.seekp(pos);			// Set original current file position
  return Tag.Bytes();			// Return specified size in tag
  }


std::string MatLab5DE::Name(std::fstream& fp, int bigend)
  {
  if(!fp) return std::string("");		// No output if no filestream
  int pos = fp.tellp();			// Get current file position
  int TF = Tag.read(fp, bigend, 1);	// Read data element tag
  TF *= AF.read(fp, bigend, 1);		// Read array flags sub-element
  TF *= DA.read(fp, bigend, 1);		// Read dim. array sub-element
  TF *= AN.read(fp, bigend, 1);		// Read array name sub-element
  std::string name = AN.Name();		// Get element array name
  fp.seekp(pos);			// Set original current file position
  return name;				// Return requested element name
  }


void MatLab5DE::skip(std::fstream& fp, int bigend)
  {
  if(!fp) return;			// No output if no filestream
  Tag.read(fp, bigend, 1);		// Read data element tag
  if(!Tag.Compressed)			// If there is data (& there 
    fp.seekp(Tag.Bytes(), std::ios::cur);	// should be in data element) then  
  }					// jump the file pointer after it



// ____________________________________________________________________________
// C             MATLAB MAT 5 Data Element Binary Output Functions
// ____________________________________________________________________________

/* These are the functions which will output a data element in MATLAB ".mat"
   binary format, version 5. This is done with a pair {tag,data}.  However,
   the data itself can be broken up into sub-elements, each of which will have
   the {tag,data} structure as well.  For matrix output the data is actually
   5 sub-elements of the following:

          Sub-Element     GAMMA Class    Data Type        Number of Bytes
     -------------------  -----------  ------------   ------------------------
     1. Array Flags        MatLab5AF   miUINT32 (6)   2*miUINT32  (8 bytes)
     2. Dimensions Array   MatLab5DA   miINT32  (5)   ND*miINT32  (ND*4 bytes)
     3. Array Name         MatLab5AN   miINT8   (1)   NC*miINT8   (NC bytes)
     4. Real Part          MatLab5Re   miDOUBLE (9)   NE*miDOUBLE (NE*8 bytes)
     5. Imaginary Part     MatLab5Im   miDOUBLE (9)   NE*miDOUBLE (NE*8 bytes)

   Thus our array output will consist of a Tag followed by the 5 Sub-Elements.
 
                Input           ML5DE	: MAT version 5 data element (this)
                                fp      : Output file stream
                Output          void    : Output file stream is modified
                                          by having the ML5DE written to it
                Note                    : No care is taken to insure the
                                          data element is written &
					  written at data element start.      */
 
int MatLab5DE::write(std::fstream& fp) const
  {
  int TF = Tag.write(fp);		// Write main data element tag
  TF *= AF.write(fp);			// The array flags sub-element
  TF *= DA.write(fp);			// The dim. array sub-element
  TF *= AN.write(fp);			// The array name sub-element
  TF *= RData.write(fp, MX);		// The reals sub-element
  if(AF.IsComplex())
    TF *= IData.write(fp, MX);		// The imags sub-element
  return TF;
  }

int MatLab5DE::write(std::fstream& fp,const matrix& mx,    const std::string& N) const
  {
  int TF = Tag.write(fp);		// Write main data element tag
  TF *= AF.write(fp);			// The array flags sub-element
  TF *= DA.write(fp, mx);		// The dim. array sub-element
  TF *= AN.write(fp, N);		// The array name sub-element
  TF *= RData.write(fp, mx);		// The reals sub-element
  if(AF.IsComplex())
    TF *= IData.write(fp, mx);		// The imags sub-element
  return TF;
  }

int MatLab5DE::write(std::fstream& fp,const row_vector& rv,const std::string& N) const
  {
  int TF = Tag.write(fp);		// Write main data element tag
  TF *= AF.write(fp);			// The array flags sub-element
  TF *= DA.write(fp, rv);		// The dim. array sub-element
  TF *= AN.write(fp, N);		// The array name sub-element
  TF *= RData.write(fp, rv);		// The reals sub-element
  if(AF.IsComplex())
    TF *= IData.write(fp, rv);		// The imags sub-element
  return TF;
  }

int MatLab5DE::write(std::fstream& fp,const col_vector& cv,const std::string& N) const
  {
  int TF = Tag.write(fp);		// Write main data element tag
  TF *= AF.write(fp);			// The array flags sub-element
  TF *= DA.write(fp, cv);		// The dim. array sub-element
  TF *= AN.write(fp, N);		// The array name sub-element
  TF *= RData.write(fp, cv);		// The reals sub-element
  if(AF.IsComplex())
    TF *= IData.write(fp, cv);		// The imags sub-element
  return TF;
  }

/* These are the functions which will return the number of bytes that are
   written upon output a data element in MATLAB ".mat" binary format, Ver. 5 
   Note that Size is the complete size of the data element (Tag + Data), but
   that DataSize is the value in bytes that will/should be present in the 
   data element tag.  This latter value will be eight bytes samaller because
   the tag size is not included in that value.                               */

int MatLab5DE::DataSize(const matrix& mx, const std::string& name, int cmplx) const
  {
  int NB  = AF.Size();			// Bytes in array flags sub-element
      NB += DA.Size(mx);		// Bytes in dim. array sub-element
      NB += AN.Size(name);		// Bytes in array name sub-element
      NB += RData.Size(mx);		// Bytes in reals sub-element
  if(cmplx)
      NB += IData.Size(mx);		// Bytes in imags sub-element
  return NB;
  }

int MatLab5DE::Size(const matrix& mx, const std::string& name, int cmplx) const
  {
  int NB  = Tag.Size();			// Bytes from data element tag
      NB += AF.Size();			// Bytes in array flags sub-element
      NB += DA.Size(mx);		// Bytes in dim. array sub-element
      NB += AN.Size(name);		// Bytes in array name sub-element
      NB += RData.Size(mx);		// Bytes in reals sub-element
  if(cmplx)
      NB += IData.Size(mx);		// Bytes in imags sub-element
  return NB;
  }

// ____________________________________________________________________________
// D             MATLAB MAT 5 Data Element Binary Input Functions
// ____________________________________________________________________________

/* These are the functions which will input a data element in MATLAB ".mat"
   binary format, version 5.  Note that in a valid MATLAB MAT file a 
   "Data Element" is a combination of Tag + Data and the Data itself may
   be broken up into Sub-Elements, each of which contains a Tag and Data
 
                Input           ML5DE   : MAT version 5 dim. array (this)
                                fp      : Input file stream
				bigend	: Stored data endian type	
                                warn	: Warning output level
						0 = no warnings
						1 = warnings
                                               >1 = fatal warnings
                Output          void    : ML5DE is set from values found
                                          in the input stream
                Note                    : No care is taken to insure the
                                          dim. array is read from the top of
 					  a data element (as it should be 
					  in MATLAB)                         */

int MatLab5DE::read(std::fstream& fp, int bigend, int warn)
  {
  int TF  = Tag.read(fp, bigend, (warn)?1:0);	// Read data element tag
      TF *= AF.read( fp, bigend, (warn)?1:0);	// Read array flags sub-element
      TF *= DA.read( fp, bigend, (warn)?1:0);	// Read dim. array sub-element
      TF *= AN.read( fp, bigend, (warn)?1:0);	// Read array name sub-element
      TF *= RData.fread(fp,bigend,(warn)?1:0);	// Pseudo read of reals SE
      if(AF.IsComplex()) 			// Pseudo read of imaginary SE
        TF *= IData.fread(fp,bigend,(warn)?1:0);
  if(!TF)
    {						// If this has problems, then
    if(warn)					// we must deal with them
      {						// appropriately now
      if(warn==1) MLDE5error(10, 1);
      else        MLDE5fatality(10);
      }
    return 0;
    }
  return 1;
  }

matrix MatLab5DE::GetMatrix(std::fstream& fp, int bigend, int warn)
  {
  int TF  = Tag.read(fp, bigend, (warn)?1:0);	// Read data element tag
      TF *= AF.read( fp, bigend, (warn)?1:0);	// Read array flags sub-element
      TF *= DA.read( fp, bigend, (warn)?1:0);	// Read dim. array sub-element
      TF *= AN.read( fp, bigend, (warn)?1:0);	// Read array name sub-element
  matrix mx = RData.read(fp,bigend,DA,warn?1:0);
  if(AF.IsComplex())
    IData.read(fp,bigend,mx,warn?1:0);
  return mx;
  }
 
// ____________________________________________________________________________
// E            MATLAB MAT 5 Data Element ASCII Output Functions
// ____________________________________________________________________________

                                                                                
void MatLab5DE::print(std::ostream& ostr) const

        // Input                ML5DE	: MAT version 5 data element (this)
        //                      ostr    : An output stream
        // Output               ostr    : Output stream with MATLAB MAT file
        //                                data element version 5 info placed
	//				  into it
        // Note                         : This information exists at the
        //                                beginning of each data element

  {
  ostr << "\n\t\tData Element";		// Output a header
  Tag.print(ostr, 0);			// Print data element tag
  AF.print(ostr); 			// Print array flags SE
  DA.print(ostr);			// Print dim. array SE
  AN.print(ostr);			// Print name array SE
  RData.print(ostr);			// Print real data SE
  if(AF.IsComplex()) IData.print(ostr);	// Print imag data SE
  ostr.flush();
  }
 
std::ostream& operator<< (std::ostream& ostr, const MatLab5DE& ML5DE)
 
        // Input                ostr    : An output stream
        //                      ML5DE	: MAT version 5 data element
        // Output               ostr    : The output stream modified by
        //                                the data element parameters

  { ML5DE.print(ostr); return ostr; }

 
// ____________________________________________________________________________
// F       MATLAB MAT 5 Data Element ASCII Output Auxiliary Functions
// ____________________________________________________________________________

/* These functions are just used for testing purposes.  The mirror the binary
   write functions, but instead of writing in binary to a filestream the write
   in ASCII to an output stream.  Thus it is easy to see most of what is
   being written, at least in principle..... The function overload with matrix
   should work with classes row_vector and col_vector also.                  */
                                                                                
void MatLab5DE::print(std::ostream& ostr,
                            const matrix& mx, const std::string& N, int cmplx) const
  {
  MatLab5Tag T(14, Size(mx,N,cmplx));	// Construct main tag
  T.print(ostr);			// Print main data element tag
  MatLab5AF A(mx, cmplx);		// Set array flags SE for matrix
  A.print(ostr);			// Print array flags sub-element
  MatLab5DA D(mx);			// Set dim. array SE for matrix
  D.print(ostr);			// Print dim. array sub-element
  MatLab5AN NAME(N);			// Set dim. array SE for matrix
  NAME.print(ostr);			// Print array name sub-element
  MatLab5Re R(mx);			// Set dim. array SE for matrix
  R.print(ostr);			// Print reals sub-element
  if(cmplx)				// See if imaginaries are desired
    {					// and deal with them if so
    MatLab5Im I(mx);			// Set dim. array SE for matrix
    I.print(ostr);			// Print imaginaries sub-element
    }
  }

#endif							// ML5DE.cc
