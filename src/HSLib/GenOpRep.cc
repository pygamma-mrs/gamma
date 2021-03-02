/* GenOpRep.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      General Operator Represenation              Implementation 	**
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zuerich / Switzerland                                      **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** The class genoprep defines a single operator representaion for a	**
** variable of type gen_op, i.e. a GAMMA general operator.  Each 	**
** representation consists of a matrix, a basis, and a priority #.	**
**                                                                      **
*************************************************************************/

#ifndef   GenOpRep_cc_ 			// Is file already included?
#  define GenOpRep_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA) 		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <HSLib/GenOpRep.h>             // Include the interface
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <Basics/StringCut.h>		// Include CenterString function

bool genoprep::BSPrnt = false;		// Flag NOT to print basis array

// ____________________________________________________________________________
// i                CLASS OPERATOR REPRESENTATION ERROR HANDLING
// ____________________________________________________________________________
 
	// Input		OpRep	: Operator representation (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               void    : An error message is output
        //                                Program execution stopped if fatal
 
/* The following error messages use the defaults set in the Gutils package
 
                Case                          Error Message
 
                (0)                     Program Aborting.....
    		(9)                     Problems During Construction
                default                 Unknown Error                        */
 
void genoprep::OpReperror(int eidx, int noret) const
  {
  std::string hdr("Operator Representation");
  std::string msg;
  switch (eidx)
    {
    case 50:GAMMAerror(hdr,"Rectangular Array Construction",noret);break;//(50)
    case 51:GAMMAerror(hdr,"Matrix-Basis Dimension Mismatch",noret);break;//(51)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

volatile void genoprep::OpRepfatal(int eidx) const
  {
  OpReperror(eidx, 1);				// Normal non-fatal error
  if(eidx) OpReperror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// A            OPERATOR REPRESENTATION CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________


genoprep::genoprep( ) { RepPty=0; }

genoprep::genoprep(const genoprep& OpRep)
  {
  RepMx  = OpRep.RepMx;				// Copy OpRep matrix 
  RepBs  = OpRep.RepBs;				// Copy OpRep basis
  RepPty = OpRep.RepPty;			// Copy OpRep priority
  }

genoprep::genoprep(const matrix& mx, const basis& bs, int pty)
  {
  RepMx  = mx;					// Set OpRep matrix
  RepBs  = bs; 					// Set OpRep default basis
  RepPty = pty;					// Set default prority
  if(!OpRepCheck(1)) OpRepfatal(9);
  }

genoprep& genoprep::operator= (const genoprep& OpRep)
  {
  if(this == &OpRep) return *this;		// Avoid self copy
  RepMx  = OpRep.RepMx;				// Copy OpRep matrix 
  RepBs  = OpRep.RepBs;				// Copy OpRep basis
  RepPty = OpRep.RepPty;			// Copy OpRep priority
  return *this;
  }

genoprep::~genoprep() {}

// ____________________________________________________________________________
// J                        CLASS OPERATOR CHECKS
// ____________________________________________________________________________

	// Input		OpRep: A operator representation
	// 			warn : A warning level
	//			         0 = no warning
	//			         1 = non-fatal warning
	//			         2 = fatal warnings
	// Output		bool : True if OprRep matrix mx is square
        //			       & its dimension matches its basis

bool genoprep::OpRepCheck(int warn) const
  {
  if(RepMx.cols() != RepMx.rows())	// Insure OpRep matrix is square
    {					// If not square, then issue some
    if(warn)				// warnings and either exit or
      {					// return that we failed.
      if(warn > 1) OpRepfatal(50);
      else         OpReperror(50,1);
      }
    return false;
    }
  if(RepMx.cols()!=RepBs.size())	// Insure OpRep matrix dimension
    {					// matches the OpRep basis dimension
    if(warn)				// If they do not match, issue some
      {					// warnings and either exit or
      if(warn > 1) OpRepfatal(51);	// return that we failed
      else         OpReperror(51,1);
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// K                  CLASS OPERATOR REPRESENTATION I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input           OpRep: Operator representation (this)
                                ostr : Output ASCII file stream
                                full : Flag for amount of output
                Return          void : OpRep is sent to the output stream    */

std::ostream& genoprep::print(std::ostream& ostr, int full) const
  {
  std::string hdr;
  if(!RepMx.rows() && full)
    {
    hdr = std::string("Empty Operator Representation");
    ostr << CenterString(hdr) << std::endl;
    return ostr;
    }
  std::string bn = RepBs.name();
  if(!BSPrnt && bn.length())
    {
    hdr = std::string("(") + bn + std::string(" Representation)");
    ostr << CenterString(hdr) << std::endl;
    }
  ostr << RepMx;
  if(BSPrnt) ostr << RepBs;
  return ostr;
  }
 
std::ostream& operator<< (std::ostream& ostr,const genoprep &OpRep)
  {return OpRep.print(ostr);}

// ------------------------ Binary Output Functions ---------------------------
 
/*              Input           OpRep: Operator representation (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at OpRep spot)
                Return          void : OpRep is written to either the
                                       specified file or filestream.
                Note                 : Output format is partially set by
                                       class matrix (matrix typing)          */
 
void genoprep::write(const std::string& fn) const
  {
  std::ofstream fp;					// Construct a file stream
  fp.open(fn.c_str(),std::ios::out|std::ios::binary);	// Open the file for output
  write(fp);                   				// Write OpRep to file stream
  fp.close();						// Close file stream
  }
 
std::ofstream& genoprep::write(std::ofstream& fp) const
  {
  fp.write((char*)&RepPty,sizeof(int));	// Write the operator rep. priority
  RepMx.write(fp);  	                // Write the operator rep. matrix
  RepBs.write(fp);			// Write the operator rep. basis
  return fp;
  }

// ------------------------ Binary Input Functions ----------------------------
 
/*              Input           OpRep: Operator representation (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at OpRep)
                Return          void : OpRep is read in from either the
                                       specified file or filestream.         */
 
void genoprep::read(const std::string& fn)
  {
  std::ifstream fp;					// Construct a file stream
  fp.open(fn.c_str(),std::ios::in|std::ios::binary);	// Open file for reading
  read(fp);                     			// Read OpRep, use overload
  fp.close();						// Close the file stream
  }
 
std::ifstream& genoprep::read(std::ifstream &fp)
  {
  fp.read((char*)&RepPty,sizeof(int));		// Read the operator rep. priority
  RepMx.read(fp);				// Read the operator rep. matrix
  RepBs.read(fp);				// Read the operator rep. basis
  return fp;
  }

// ____________________________________________________________________________
// L        OPERATOR REPRESENTATION LIST/VECTOR SUPPORT FUNCTIONS
// ____________________________________________________________________________

bool genoprep::operator==(const genoprep& OpRep) const
  {
  if(RepMx != OpRep.RepMx) return false;
  if(RepBs != OpRep.RepBs) return false;
  return true;
  }

bool genoprep::operator!=(const genoprep& OpRep) const
  { return (!((*this) == OpRep)); }
 
bool genoprep::operator<(const genoprep& OpRep) const
  { return (RepPty < OpRep.RepPty); }

bool genoprep::operator>(const genoprep& OpRep) const
  { return (RepPty > OpRep.RepPty); }

#endif							// GenOpRep.cc
