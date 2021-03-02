/* FrameMakerS.h ************************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
*      FrameMaker Stack Plots                    Interface 		**
**                                                                      **
**      Copyright (c) 1997                                              **
**      Scott A. Smith                                                  **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
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
**  Note that this module relies on other GAMMA FM modules, FrameMakerC **
**  FrameMakerP and FrameMakerM, they which deal with FrameMaker        **
**  constructs, plotting parameters, and MIF output respectively.       **
**                                                                      **
*************************************************************************/

#ifndef   GFrameMakerS_h_		// Is this file already included?
#  define GFrameMakerS_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <GamIO/FrameMakerC.h>		// Know FM MIF constructs (PolyLines)
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Matrix/matrix.h>		// Know about matrices
#include <string>
//class ostream;				// Know output streams exist

//forward declarations
class FM_stack;

MSVCDLL void FM_stack(const std::string& filename, const matrix &mx, double xinc, 
          double yinc, int RI=1, double xsize=14, double ysize=14, int grid=0, int CI=1);

// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
//             CLASS RowBlock: FRAMEMAKER STACK PLOT ROW
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________

 
enum AreaType { _above, 			// Point above row level
                _baseplane,			// Point below row level
                _below				// Point below base plane
              };
 

enum PtTrend { _up, 				// Point in area up
               _same,				// Point in same area
               _down				// Point in area below
              };


class RowBlk : public row_vector
  {
  int* Bareas;				// Pointer to array of area flags
  int Bpts;				// Points in Block and Areas
  int Brow;				// Block row index
  double Bhoff;				// Block horizontal offset
  double Bvoff;				// Block vertical offset
  double Bhsf;				// Block horzontal scaling factor
  double Bvsf;				// Block vertical scaling factor
  int Bids[4];				// Block associated ID values
  double Blevs[4];			// Block associated levels
  int lastarea;				// Last point area flag
  double lasty;				// Last point vertical height
  int lastlev;				// Last point level flag
  int Bdebug;				// Debugging level
  int BaseID;				// FM Base ID value

  friend class FMStack;			// Allow FMstack FULL Access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------


void RBerror(int eidx, int noret=0) const;

        // Input                Blk     : Stack plot row block (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
        
volatile void RBfatality(int eidx) const;
 
        // Input                Blk     : Stack plot row block (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped
 

public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   ROW BLOCK CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


MSVCDLC RowBlk( );

        // Input                none    : No arguments required
        // Output               none    : FrameMaker parameters

    
MSVCDLC RowBlk(int npts);
 
        // Input                none    : No arguments required
        // Output               none    : FrameMaker parameters
 
 

MSVCDLL void operator=(const RowBlk& Blk);
 
        // Input                Blk     : Stack plot row block (this)
        //                      Blk1	: Stack plot row
        // Output               void    : FrameMaker stack plot row Blk1
        //                                is copied into Blk
        ///F_list =                     - Assignment
    

MSVCDLC ~RowBlk();
 
        // Input                Blk     : Stack plot row block (this)
        // Output               none    : Blk is destructed

    
// ____________________________________________________________________________
// B                       ROW BLOCK ACCESS FUNCTIONS
// ____________________________________________________________________________

// ---------------------------- Row Index Access ------------------------------


MSVCDLL void SetRow(int row);
 
        // Input                Blk     : Stack plot row block (this)
        //                      row     : Row index
        // Output               void    : Block row index set to row
 

MSVCDLL int GetRow() const;
 
        // Input                Blk     : Stack plot row block (this)
        // Output               Brow    : Block row index
        // Note                         : This will set the IDs also

 
// ----------------------------- Row ID Access -------------------------------


MSVCDLL void SetIDs();

        // Input                Blk     : Stack plot row block (this)
        // Return               void    : Possible block ID values are set

 
// --------------------------- Row Levels Access -----------------------------
 
MSVCDLL void SetLevels(double* Levels);
 
        // Input                Blk     : Stack plot row block (this)
        // Return               void    : Possible block ID values are set
 


// -------------------------- Row Offsets Access -----------------------------


MSVCDLL void SetOffsets(double ho, double vo);
 
        // Input                Blk     : Stack plot row block (this)
        //                      ho      : Horizontal offset (cm) 
        //                      vo      : Vertical offset (cm)
        // Output               void    : Block offsets set (cm)
 

MSVCDLL double VOffset() const;

        // Input                Blk     : Stack plot row block (this)
        // Output               Bvoff   : Block row vertical offset (cm)


MSVCDLL double HOffset() const;
 
        // Input                Blk     : Stack plot row block (this)
        // Output               Bhoff   : Block row horizontal offset (cm)
 
 
// ---------------------- Row Scaling Factor Access --------------------------
 
 
MSVCDLL void SetScaling(double HSF, double VSF);

        // Input                Blk     : Stack plot row block (this)
        //                      ho      : Horizontal offset (cm)
        //                      vo      : Vertical offset (cm)
        // Output               void    : Block offsets set (cm)

 
MSVCDLL double VScale() const;
 
        // Input                Blk     : Stack plot row block (this)
        // Output               Bvsf    : Block vertical scaling (cm/asbunits)
 
 
MSVCDLL double HScale() const;
 
        // Input                Blk     : Stack plot row block (this)
        // Output               Bhsf    : Block horiz. scaling (cm/asbunits)
 
 
// ------------------------ RowBlk Debug Flag Access --------------------------
 
 
MSVCDLL int GetDebug();  
 
        // Input                Blk     : Stack plot row block (this)
        // Output               Bdebug : The current debugging level
 
 
MSVCDLL void SetDebug(int db);
 
        // Input                Blk     : Stack plot row block (this)
        //                      db      : A debugging level 
        // Output               void    : Debugging level set to db

// ____________________________________________________________________________
// C                     ROW BLOCK DETERMINATION FUNCTIONS
// ____________________________________________________________________________

 
void Zero();
 
        // Input                Blk     : Stack plot row block (this)
        // Return               void    : The point count is set to zero
 

void AddPt(int col, const matrix& mx);

        // Input                Blk     : Stack plot row block (this)
        //                      col     : Column index
	//			mx	: Data array for stack plot
        // Return               void    : The point <row|mx|col> is added to
        //                                the block Blk in plot coordinates


MSVCDLL void AddPt(complex& pt);
 
        // Input                Blk     : Stack plot row block (this)
        //                      pt      : Point (FM units cm)
        // Return               void    : The point pt is added to
        //                                the block Blk in plot coordinates
 

MSVCDLL void AddPt(complex& pt, int area);
 
        // Input                Blk     : Stack plot row block (this)
        //                      pt      : Point (FM units cm)
        //                      area    : Point area flag
        // Return               void    : The point pt is added to
        //                                the block Blk in plot coordinates
        //                                in the area specified
 
 
MSVCDLL complex GetPt(int col, const matrix& mx) const;
 
        // Input                Blk     : Stack plot row block (this)
        //                      col     : Column index 
        //                      mx      : FM stack plot array (this)
        // Return               void    : The point <row|mx|col> is 
        //                                returned in plot coordinates
        // Note                         : The block is unaffected              
 
 
MSVCDLL int GetArea(complex& z) const;

        // Input                Blk     : Stack plot row block (this)
        //                      z       : Block point
        // Return               area    : An area flag for the point is
        //                                returned
        // Note                         : The block is unaffected
 

MSVCDLL PtTrend Trend(int area) const;

        // Input                Blk     : Stack plot row block (this)
        //                      area    : Area flag
        // Return               PtTrend : Row trend if area is associated
        //                                with the next pt in Blk
        // Note                         : The block is unaffected

 

MSVCDLL complex ExtraPt(complex& pt1, double y);

        // Input        Blk     : Stack plot row block (this)
        //              pt1     : A point (FM units, cm)
        //              y       : Desired y coordinate (cm)
        // Output       pt      : Point on the line between pt1
        //                        last point of pt in Blk having
        //                        a height of y

//                         y1-y2                x1-x2
// Line: y-y2 = m*(x-x2) = ----- (x-x2) --> x = ----- (y-y2) + x2
//                         x1-x2                y1-y2


MSVCDLL void SetPL(FMPL& PL, int area);
 
        // Input                Blk     : Stack plot row block (this)
        //                      area    : Area flag
        // Output               void    : Block is filled for mx row


 
MSVCDLL void Fill(int cols, const matrix& mx);

        // Input                Blk     : Stack plot row block (this)
	//			cols	: Number of matrix columns
	//                      mx      : Data array for stack plot 
	// Output               void    : Block is filled for mx row


MSVCDLL void Plot(const std::string& filename);

        // Input        Blk     : Stack plot row block (this)
        //              filename: Output MIF filename
        // Return       void    : Block Blk, a stack plot row,
        //                        is put into a file filename
        //                        in FrameMaker MIF format.


// ____________________________________________________________________________
// D                      ROW BLOCK STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& print(std::ostream& ostr) const;
 
        // Input                Blk     : Stack plot row block (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                the FrameMaker block parameters
 
 
MSVCDLL friend std::ostream& operator<<(std::ostream& ostr, const RowBlk& Blk);
 
        // Input                ostr    : An output stream
        //                      Blk     : Stack plot row block (this)
        // Output               none    : Modifies output stream
        ///F_list <<                    - FrameMaker row block sent to ostr

  };						// End Class RowBlk
 

// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
//             CLASS FMStack: FRAMEMAKER STACK PLOTTING PARAMETERS
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
     

class FMStack: public FMPar
  {
  int rows;                             // Row span of plotted array
  int cols;                             // Column span of plotted array
  double vinc;                          // Vertical increment (cm/row)
  double hinc;                          // Horizontal increment (cm/row)
  int rowinc;                           // Row increment in plot
  double hdelta;                        // Total horizontal offset (cm)
  double vdelta;                        // Total vertical offset (cm)
  double hwidth;                        // Actual plotted data width (cm)
  double vheight;                       // Actual plotted data height (cm)
  double HSF;                           // Horizontal scaling (cm/pt)
  double VSF;                           // Vertical scaling (cm/abs.units)
  double Bback;				// Base plane back coord. (cm)
  double Bfront;			// Base plane front coord. (cm)
  double BbackW;			// Base plane back, west (cm)
  double BbackE;			// Base plane back, east (cm)
  double BfrontW;			// Base plane front, east (cm)
  double BfrontE;			// Base plane front, east (cm)
  int BPID;				// Base plane ID value
  int FRID;				// Anchored frame ID value
  int gridding;                         // Column gridding flag
  int colinc;                           // Stack column grid increment
  matrix MX;				// Data matrix (absolute untis)
  row_vector OffSets;			// Row offsets (cm/row)
  row_vector MaxMins;			// Row maxima (absolute units)

  RowBlk RB;				// FM row block
  FMPL RPL;				// FM row polyline
  int HLF;				// Hidden Line Algorithm Flag

double bhoff;				// Block horizontal offset
double bvoff;				// Block vertical offset
double Levels[4];			// Array of row levels

  friend class RowBlk;			// Allow RowBlk FULL Access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
//                    CLASS FM STACK PARAMETER ERROR HANDLING
// ____________________________________________________________________________
 

void FMSTKerror(int eidx, int noret=0) const;

        // Input                FMSTK   : FrameMaker stack parameters (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
 
volatile void FMSTKfatality(int eidx) const;

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
 

MSVCDLC FMStack();

        // Input                none    : No arguments required
        // Output               none    : FrameMaker stack parameters
        ///F_list FMP                   - Constructor

                                                      
MSVCDLC ~FMStack();

        // Input                FMSTK   : FrameMaker stack parameters (this)
        // Output               none    : FMSTK is destructed

                                                            
MSVCDLL void operator=(const FMStack& FMSTK1);

        // Input                FMSTK	: FrameMaker stack parameters (this)
        //                      FMSTK1  : FrameMaker stack parameters
        // Output               void    : FrameMaker stack parameters FMSTK1
        //                                copied into FMSTK
        ///F_list =                     - Assignment

                                                     
// ____________________________________________________________________________
// B                    FM STACK PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// ---------------------------- Stack Plot Offsets ----------------------------

 
MSVCDLL double VInc() const;
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               vinc    : Vertical increment (cm/row)
 
 
MSVCDLL void VInc(double vi);
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      vi      : Vertical increment (cm/row)
        // Output               vinc    : Vertical increment set to vi
 
 
MSVCDLL double HInc() const;
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               hinc    : Horizontal increment (cm/row)
 
 
MSVCDLL void HInc(double hi);
               
        // Input                FMSTK   : FM stack parameters (this)
        //                      hi      : Horizontal increment (cm/row)
        // Output               hinc    : Horizontal increment set to hi
 
 
// ----------------------- Stack Plot Row Increment ---------------------------
 
MSVCDLL int RInc() const;
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               rowinc  : The stack plot row increment
 
 
MSVCDLL void RInc(int ri);
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      ri      : Row increment
        // Output               void    : The row increment is set to ri
 
 
// ---------------------- Stack Plot Column Increment -------------------------
 
 
MSVCDLL int CInc() const;
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               colinc  : The stack plot column increment
 
 
MSVCDLL void CInc(int ci);
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      ci      : Column increment
        // Output               void    : The column increment is set to ci
 
 
// -------------------------- Stack Plot HLA Usage ----------------------------
 
 
MSVCDLL int HLAF() const;
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               HLF	: The hidden line use flag
 
 
MSVCDLL void HLAF(int hf);
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      hf      : Hidden line use flag
        // Output               void    : The hidden line use flag is set

 
// ------------------------ Stack Plot Gridding Flag --------------------------
 
         
MSVCDLL int Gridding() const;

        // Input                FMSTK   : FM stack parameters (this)
        // Output               gridding: The stack plot gridding flag

                                                                       
MSVCDLL void Gridding(int gr);

        // Input                FMSTK   : FM stack parameters (this)
        //                      gr      : Gridding flag
        // Output               void    : The gridding flag is set to gr


// ____________________________________________________________________________
// C                    STACK PARAMETER SETTING FUNCTIONS
// ____________________________________________________________________________

// ------------------------- Stack Plot Data Matrix ---------------------------

MSVCDLL void SetArray(const matrix& mx);

        // Input        FMSTK   : FM stack parameters (this)
        //              mx      : Input data matrix to plot
        // Output       void    : The following values in FMSTK are
        //                        set according to the input matrix

        //                              MX   = The array to be plotted
        //                              rows = The # of rows to plot
	//                              cols = The # of cols to plot

// ------------------------- Stack Plot Row Offsets ---------------------------

MSVCDLL void SetOffsets();

        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       void    : The vector of row offsets
        //                        is filled (in cm/row)
        // Note                 : The offsets for row i will be
        //                        <offsets|i> in cm. The real part is
        //                        the horizontal offset while the imaginary
        //                        part is the vertical offset


MSVCDLL row_vector GetOffsets();

        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       OffSets	: The vector of row offsets
        //                        is returned (in cm/row)
        // Note                 : The offsets for row i will be
        //                        <offsets|i> in cm. The real part is
        //                        the horizontal offset while the imaginary


// ------------------------- Stack Plot Row Maxima ---------------------------


MSVCDLL void SetMaxima();

        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       void    : The vector of row max & min values
        // Note                 : The maxima for row i will be
        //                        <maxmins|i> in absolute units. The
        //                        real part is the row maximum and
        //                        the imaginary part is the row minimum.
        // Note                 : Some rows may be flagged so that they
        //                        are not plotted (via row increment)
        //                        Their maxima are still determined!


MSVCDLL row_vector GetMaxima();

        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       MaxMins: The vector of row maxima
        //                        is returned (in absolute units)


// -------------------------- Stack Plot Scaling -----------------------------


MSVCDLL void SetScale();

        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       void    : Scaling parameters in FMSTK are set
        //                              hwidth: Plotted points width (cm)
        //                              HSF:    Horiz. scaling (cm/pt)
        //                              VSF:    Vertical scaling (cm/absunits)

 
// ____________________________________________________________________________
// D                       FM STACK BASE PLANE FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL void SetBasePlane();
 
        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       void    : Scaling parameters in FMSTK are set
        //                              hwidth: Plotted points width (cm)
        //                              HSF:    Horiz. scaling (cm/pt)
        //                              VSF:    Vertical scaling (cm/absunits)
        // Note                 : The "skew" herein is horizontal shift (cm)
        //                        difference between the front base & back base
 

MSVCDLL void BaseStart(std::ostream& ostr) const;

 
        // Input        FMSTK   : FM stack plot parameters (this)
        //              ostr    : Open output stream
        // Return       void    : The output stream is modified
        //                        to contain the start of the base
        //                        plane in a stack plot in FM MIF.
        // Note                 : This plots two of the 4 lines that
        //                        make up the base plane borders. The
        //                        two plotted will be those NOT hidden
        //                        by the stack plot, and which two that
        //                        will be depends upon the plot skew.


MSVCDLL void BaseEnd(std::ostream& ostr) const;

        // Input        FMSTK   : FM stack plot parameters (this)
        //              ostr    : Open output stream
        // Return       void    : The output stream is modified
        //                        to contain the finish of the base
        //                        plane in a stack plot in FM MIF.
        // Note                 : This plots two of the 4 lines that
        //                        make up the base plane borders. The
        //                        two plotted will be those that ARE hidden
        //                        by the stack plot, and which two that
        //                        will be depends upon the plot skew.


MSVCDLL int BasePlaneID() const;
 
        // Input        FMSTK   : FM stack plot paramters (this)
        // Return       BPID    : The base plane ID number is returned


MSVCDLL void BasePlaneID(int bid);

        // Input        FMSTK   : FM stack plot paramters (this)
        //              bid     : FM group ID number
        // Return       BPID    : The base plane ID number is set to bid
 

// ____________________________________________________________________________
// E               FM STACK PLOT INITIALIZATON FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void maxima(double& max, double& min, int& imax, int& imin);
 
        // Input        FMSTK	: FM stack parameters (this)
        //              max	: Maximum vertical value (cm)
        //              min	: Minimum vertical value (cm)
        //              imax	: Index of row containing maximum
        //              imin	: Index of row containing minimum
        // Return		: Void, imax and imin are modified
        // Note			: Formulae for getting intensities are
        //                          I = Re<MaxMins|i>*VSF + Im<OffSets|i>
        //                          I = Im<MaxMins|i>*VSF + Im<OffSets|i>
 

MSVCDLL void CheckRI();

        // Input        FMSTK   : FM stack parameters (this)
        // Output       void    : The current row increment
        //                        is checked to insure at least
        //                        two rows are plotted


MSVCDLL void StkBegin(const matrix& mx);

        // Input        FMSTK   : FM stack parameters (this)
        //              mx      : Input data matrix to plot


MSVCDLL void Initialize(const matrix& mx);

        // Input        FMSTK   : FM stack plot parameters (this)
        //              mx      : Data matrix
        // Return               : Void, all internal FMSTK initial values
        //                        that are dependent upon the matrix mx
        //                        are set and/or checked


MSVCDLL void StartOutput(const std::string& filename, std::ostream& ostr);
 
        // Input        FMSTK   : FM stack plot parameters (this)
        //              filename: Output filename
	//		ostr	: An output stream
        // Return       void  	: The output stream is begun with
        //                        external name "filename".  The
        //                        basics of a FrameMaker MIF file
        //                        are written into the output stream


// ____________________________________________________________________________
// F                   FM STACK PLOT DEBUGGING FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void StkInfo() const;

        // Input        FMSTK   : FM stack plot parameters (this)
        // Output       void    : Some information output to screen


MSVCDLL void ScaleInfo(int iter, double max, double min, int imax, int imin);
 
        // Input        FMSTK   : FrameMaker stack plot parameters
        // Return       void    : Vertical scaling information output


// ____________________________________________________________________________
// G                      FM STACK ROW BLOCK FUNCTIONS
// ____________________________________________________________________________
 

MSVCDLL void FillBlock(int row);

        // Input        FMSTK   : FM stack plot paramters (this)
        //              row     : Matrix row index
        // Return       void	: Row block is filled for row. Additional
	//			  points may be added to input array row
        //                        in order to define the row for HLA


MSVCDLL void PlotRow(int row, std::ostream& out);

        // Input        FMP     : FrameMaker plot parameters (this)
        //              row     : Matrix row index
        //              out     : Output stream
        // Return       void    : Row "row" of stack plot is sent into the
        //                        output stream in FrameMaker MIF format.


// ____________________________________________________________________________
// H                   FM STACK PLOT CREATION FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void FM_stack(const std::string& filename, const matrix &mx);
MSVCDLL void Plot(const std::string& filename, const matrix &mx);

        // Input        FMSTK   : FM stack plot parameters (this)
        //              filename: Output filename
        //              mx      : Data matrix
        // Return               : Void, file filename created and
        //                        will contain a stack plot of the
        //                        array mx in FrameMaker MIF format


MSVCDLL friend void FM_stack(const std::string& filename, const matrix &mx, FMStack& FMSTK);
 
        // Input        filename: Output filename
        //              mx      : Data matrix
        //              FMSTK   : FM stack plot parameters
        // Return               : Void, file out is modified
 

//friend void FM_stack(const string& filename, const matrix &mx, double xinc,
//            double yinc, int RI, double xsize, double ysize, int grid, int CI);


MSVCDLL friend void FM_stack(const std::string& filename, const matrix &mx, double xinc, 
          double yinc, int RI, double xsize, double ysize, int grid, int CI);


 
        // Input        filename  : Output filename
        //              mx        : Data matrix
        //              xinc      : Delta x in cm
        //              yinc      : Delta y in cm
        //              RI        : Row increment
        //              xsize     : Plot horizontal (x) dimension in cm
        //              ysize     : Plot vertical (y) dimension in cm
        //              grid      : Flag for gridding
        //              CI        : Column increment (only with gridding!)
        // Return                 : Void, file out is modified


// ____________________________________________________________________________
// I                  FM STACK HIDDEN LINE FUNCTIONS
// ____________________________________________________________________________


MSVCDLL complex xings(complex& z1, complex& z2, const complex &z3, const complex &z4,
                         int row, double xoff, double yoff, int k, int& cross);

        // Input        z1      : First point, on line 1
        //              z2      : Second point, on line 1
        //              z3      : First point, line 2
        //              z4      : Second point, line 2
        //              xoff    : Horizontal offset for line 1 (cm)
        //              yoff    : Vertical offset for line 1 (cm)
        //              FMxscale: x-dimension scale (cm/point)
        //              FMyscale: y-dimension scale (cm/Re(mx(i,j)))
        //              mx      : Matrix of data values
        //              row     : Row number of line 1
        //              k       : Point being checked
        //              cross   : Flag for line crossing
        // Output       z       : Point at the intersection of
        //                        the line connecting z1 & z2 with
        //                        the line connecting z3 and z4
        // Note                 : There may be no intersection !
        //                        in which case cross=FALSE
        // Note                 : Used exclusively by set_hidden
 


MSVCDLL void set_hidden(int &pts, int row, int *hidden, row_vector& crossing);

        // Input        FMSTK   : FM stack parameters (this)
        //              pts     : Number of points in PolyLine
        //              row     : Current row
        //              hidden  : Flags whether points are hidden
        //              crossing: Crossing points between visible & hidden
        // Output       none    : This function provided to services -
        //                        1.) Each point in the PolyLine is flagged
        //                            whether it should be visigle or hidden
        //                        2.) Two adjancent points which have a
        //                            transition between visible and hidden
        //                            are extrapolated to find this crossing
        // Note                 : AUXILIARY FUNCTION FOR FM_Stack



MSVCDLL int continuous(int row, double xoff, double yoff);

        // Input        FMSTK   : FM stack plot parameters (this)
        //              row     : Row index
        //              xoff    : Horizontal offset for row (in cm)
        //              yoff    : Vertical offset for row (in cm)
        // Output       npts    : Total number of working points in
        //                        vector "block".
        // Note                 : This routine alters the values of
        //                        in arrays block &  areas



MSVCDLL void HLA(std::ostream& ostr);

        // Input        FMSTK   : FM stack plot parameters
        //              ostr    : Output stream
        // Output       none    : Outputs current PolyLine to ostr
        //                        in FrameMaker MIF format using a HLA
        //                        (Hidden Line Agorithm) to determine
        //                        whether each point is visible.

// The points in the input PolyLine PL all reside in "area 1", the part of
// the stack plot which is below the current rows height yet above the bottom
// of the base plane.  That means that the points will be "hidden" (drawn with
// a white pen) if they happen to be behind any of the rows that are in front
// of them.  This reasoning is NOT applicable to any points above "area 1"
// because PolyLines in front are drawn with a white fill so they will auto-
// matically hide those behind them.


 
MSVCDLL void HLA(int &pts, int row, std::ostream & out);
 
        // Input        FMSTK   : FM stack plot parameters
        //              pts     : Number of points in PolyLine
        //              row     : Current row
        //              out     : Output stream for PolyLine
        // Output       none    : Outputs PolyLine to ostream "out" in
        //                        FrameMaker MIF format using an HLA
        //                        (Hidden Line Agorithm) to determine
        //                        whether each point is visible.

// ____________________________________________________________________________
// Z                  FM STACK PARAMETERS STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& print(std::ostream& ostr) const;

        // Input                FMSTK   : FrameMaker stack parameters (this)
        //                      ostr    : An output stream
        // Output               none    : Modifies output stream
        ///F_list print                 - FrameMaker parameters sent to ostr


MSVCDLL friend std::ostream& operator<<(std::ostream& ostr, const FMStack& FSTK);

        // Input                ostr    : An output stream
        //                      FMSTK   : FrameMaker stack parameters
        // Output               none    : Modifies output stream
        ///F_list <<                    - FrameMaker parameters sent to ostr

  };


#endif							// FrameMakerS.h
