/* FrameMakerC.h ************************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
*      FrameMaker Constructs			Interface 		**
**                                                                      **
**      Copyright (c) 1997                                              **
**      Scott A. Smith                                                  **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  GAMMA's FrameMaker modules provide functions to produce FrameMaker  **
**  files from within GAMMA programs.  This particular module contains  **
**  the functions that produce stack plots from vectors and matrices.   **
**                                                                      **
**  Note that this module relies on the two other GAMMA FM modules,     **
**  FrameMakerP and FrameMakerM, they which deal with FrameMaker        **
**  plotting parameters and MIF output respectively.                    **
**                                                                      **
*************************************************************************/

#ifndef   GFrameMakerC_h_		// Is this file already included?
#  define GFrameMakerC_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <GamIO/FrameMakerP.h>		// Know FM plotting parameters
#include <GamIO/FrameMakerM.h>		// Know FM MIF functions
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Matrix/complex.h>		// Know about GAMMA complex #s

class FMPL : public row_vector
  {
  double PWidth; 			// PolyLine pen width
  int PColor;				// PolyLine pen color
  int PFill;				// PolyLine pen fill
  int PID;				// PolyLine group ID
  int Ppts;				// PolyLine points
  int PLdebug;				// PolyLine debugging flag
  int Preduce;				// PolyLine reduction flag
 
                
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
//                    CLASS FM STACK PARAMETER ERROR HANDLING
// ____________________________________________________________________________
 

void PLerror(int eidx, int noret=0) const;

        // Input                FMPL	: FrameMaker PolyLine (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
 
volatile void PLfatality(int eidx) const;

        // Input                FMSTK   : FrameMaker stack parameters (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped
 

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:

// ____________________________________________________________________________
// A                 FM STACK PARAMETERS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 

MSVCDLC FMPL();

        // Input                none    : No arguments required
        // Output               none    : FrameMaker stack parameters
        ///F_list FMP                   - Constructor


MSVCDLC FMPL(int npts);

        // Input                PL      : FM PolyLine (this)
        //                      npts    : Number of points
        // Output               none    : FrameMaker PolyLine
        ///F_list FMPL                  - Constructor

                                                      
MSVCDLC ~FMPL();

        // Input                FMSTK   : FrameMaker stack parameters (this)
        // Output               none    : FMSTK is destructed

                                                            
MSVCDLL void operator=(const FMPL& FMSTK1);

        // Input                FMSTK	: FrameMaker stack parameters (this)
        //                      FMSTK1  : FrameMaker stack parameters
        // Output               void    : FrameMaker stack parameters FMSTK1
        //                                copied into FMSTK
        ///F_list =                     - Assignment

 
// ____________________________________________________________________________
//                      POLYLINE PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// ----------------------------- Pen Width Access -----------------------------
 
MSVCDLL void SetWidth(double pw);
 
        // Input        PL      : FM PolyLine (this)
        //              pw      : Input pen width
        // Output       void    : The PL pen width is set to pw
 
         
MSVCDLL double GetWidth();
 
        // Input        PL      : FM PolyLine (this)
        // Output       pw      : The current pen width
 

// ----------------------------- Pen Color Access -----------------------------
 
 
MSVCDLL void SetColor(int pc);
 
        // Input        PL      : FM PolyLine (this)
        //              pc      : Input pen color
        // Output       void    : The PL pen color is set to pc
 
 
MSVCDLL int GetPen();
 
        // Input        PL      : FM PolyLine (this)
        // Output       pc      : The current pen color

// --------------------------- PolyLine Fill Access ---------------------------
 
MSVCDLL void SetFill(int pf);
 
        // Input        PL      : FM PolyLine (this)
        //              pf      : Input pen fill
        // Output       void    : The PL pen fill is set to pf
 
 
MSVCDLL int GetFill();
 
        // Input        PL      : FM PolyLine (this)
        // Output       pf      : The current pen fill
 
 
// ---------------------------- PolyLine ID Access ----------------------------
 
 
MSVCDLL void SetID(int pid);
 
        // Input        PL      : FM PolyLine (this)
        //              pid     : Input ID
        // Output       void    : The PL ID is set to pid
 
 
MSVCDLL int GetID();
 
        // Input        PL      : FM PolyLine (this)
        // Output       PID     : The current PL ID


// -------------------------- PolyLine Point Access ---------------------------
 
 
MSVCDLL int GetPoints();
 
        // Input        PL      : FM PolyLine (this)
        // Output       Ppts    : The current # points


// ----------------------- PolyLine Debug Flag Access -------------------------
 
 
MSVCDLL int GetDebug();
 
        // Input        PL      : FM PolyLine (this)
        // Output       PLdebug	: The current debugging level
 
 
MSVCDLL void SetDebug(int db);
 
        // Input        PL      : FM PolyLine (this)
	//		db	: A debugging level
        // Output       void 	: Debugging level set to db
 

 
// ------------------------- PolyLine Data Reduction --------------------------
 
  
MSVCDLL int GetReduce();
 
        // Input        PL      : FM PolyLine (this) 
        // Output       Preducd : The current reductoin level
 
 
MSVCDLL void SetReduce(int rf);
 
        // Input        PL      : FM PolyLine (this) 
        //              rf      : Data reduction flag
        // Output       void    : Reduction level set to rf 
 

// ____________________________________________________________________________
//                       POLYLINE AUXILARY FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL std::string Colors(int pc) const;
 
        // Input        PL      : FM PolyLine (this)
        //              pc      : Pen color 
        //              SPC     : string for Pen color
 
 
MSVCDLL std::string Filling(int pf) const;
 
        // Input        PL      : FM PolyLine (this)
        //              pf      : PolyLine fill
        //              SPC     : string for PL fill


MSVCDLL std::string Reduction() const;

        // Input        PL      : FM PolyLine (this)
        //              SRed    : string for point reduction


MSVCDLL void Set(double pw, int pc, int pf, int id);
 
        // Input        PL      : FM PolyLine (this)
        //              pw      : Pen width (pts)
        //              pc      : Pen color
        //              pf      : PolyLine fill
        //              id      : PolyLine ID
        // Output       void    : PL parameters are set
 


// ____________________________________________________________________________
// C                     ROW BLOCK DETERMINATION FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void Zero();

        // Input        PL      : FM PolyLine (this)
        //              z       : PolyLine point
        // Return       void    : The point z is added to PL
 

MSVCDLL void AddPt(const complex& z);

        // Input        PL      : FM PolyLine (this)
        // 		z	: PolyLine point
        // Return	void    : The point z is added to PL


MSVCDLL void WriteMIF(std::ostream& ostr);

        // Input        PL      : FM PolyLine (this)
	//		ostr	: An output stream
	// Output	void	: The PolyLine is written into
	//			  the output stream ostr in MIF format

// ____________________________________________________________________________
// Z                  FM POLYLINE STANDARD I/O FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL std::ostream& print(std::ostream& ostr) const;
 
        // Input                PL      : FM PolyLine (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                the FrameMaker PolyLine
 
    
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const FMPL& PL);
 
        // Input                ostr    : An output stream
        //                      PL      : FM PolyLine (this)
        // Output               ostr    : The output stream modified by
        //                                the FM PolyLine
 

                                                     
  };


#endif							// FrameMakerC.h
