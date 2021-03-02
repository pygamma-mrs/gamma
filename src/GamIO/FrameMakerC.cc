/* FrameMakerC.cc ***********************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
*      FrameMaker Constructs 			Implementation		**
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
**  Note that this module relies on the two other GAMMA FM modules,     **
**  FrameMakerP and FrameMakerM, they which deal with FrameMaker        **
**  plotting parameters and MIF output respectively.                    **
**                                                                      **
*************************************************************************/

#ifndef   GFrameMakerC_cc_			// Is file already included?
#  define GFrameMakerC_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#if defined(_MSC_VER)				// If we areusing MSVC++
 #pragma warning (disable : 4786)		// Kill STL namelength warnings
#endif

#include <GamIO/FrameMakerC.h>			// Include the header
#include <GamIO/FrameMakerM.h>			// Include FM MIF functions
#include <GamIO/FrameMakerP.h>			// Include FM parameters
#include <GamIO/FrameMaker.h> 			// Include FM 1D plotting
#include <Basics/Gconstants.h> 			// Include PI
#include <Basics/StringCut.h>			// Gdec & Gform (dec/form) fncts
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <string>				// Include libstdc++ STL strings
#include <iostream>				// Include input output streams

using std::ostream;				// Using libstdc++ output streams
using std::string;				// Using libstdc++ strings
using std::cout;				// Using libstdc++ standard output

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  CLASS FM STACK PARAMETER ERROR HANDLING
// ____________________________________________________________________________
 

        // Input                PL	: FrameMaker PolyLine (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
        //                                Program execution stopped (fatal)

void FMPL::PLerror(int eidx, int noret) const
  {                                                     
  cout << "\nFrameMaker PolyLine: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      cout << "Program Aborting.....";
      break;
   case 1:								// (1)
      cout << "Allocated PolyLine Dimension Exceeded!";
      break;
   case 2:								// (2)
      cout << "Cannot Add Point To PolyLine.";
      break;
   case 3:								// (3)
      cout << "Less Than 2 Points In PolyLine.";
      break;
   case 4:								// (4)
      cout << "Cannot Write PolyLine To File.";
      break;
   case 5:								// (5)
      cout << "Output Stream Is Problematic?";
      break;
    default:
      cout << "Unkown Error, Number " << eidx;
    }
  if(!noret) cout << ".\n";
  }  

 
volatile void FMPL::PLfatality(int eidx) const

  {
  PLerror(eidx,1);				// Output error message
  if(eidx) PLerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  FM POLYLINE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 


        // Input                none    : No arguments required
        // Output               none    : FrameMaker PolyLine
        ///F_list FMPL                  - Constructor

FMPL::FMPL() : row_vector()		// Call base constructor
  {
  PWidth = 1.0;				// Default Pen Width
  PColor = 0;				// Default Pen Color (black)			
  PFill = 0;				// Default Line Fill (none)
  Ppts = 0;				// Default No Points 
  PID = 1;				// Default Line ID
  PLdebug = 0;				// Default debugging
  Preduce = 0;				// Default point reduction
  }
 


        // Input                PL	: FM PolyLine (this)
	//			npts	: Number of points
        // Output               none    : FrameMaker PolyLine
        ///F_list FMPL                  - Constructor

FMPL::FMPL(int npts) : row_vector(npts)	// Call base constructor
  {
  PWidth = 1.0;				// Default Pen Width
  PColor = 0;				// Default Pen Color (black)			
  PFill = 0;				// Default Line Fill (none)
  Ppts = 0;				// Default No Points 
  PID = 1;				// Default Line ID
  PLdebug = 0;				// Default debugging
  Preduce = 0;				// Default point reduction
  }
                                                      

FMPL::~FMPL() { }

        // Input                PL	: FrameMaker PolyLine (this)
        // Output               none    : FMPL is destructed

                                                            

        // Input                PL	: FrameMaker PolyLine (this)
        //                      PL1	: FrameMaker PolyLine
        // Output               void    : FrameMaker PolyLine PL1
        //                                copied into FMPL
        ///F_list =                     - Assignment

void FMPL::operator=(const FMPL& PL1)
  {
  row_vector::operator=(PL1);		// Copy the row vector
  PWidth  = PL1.PWidth;			// Copy pen width
  PColor  = PL1.PColor;			// Copy pen color
  PFill   = PL1.PFill;			// Copy PL fill
  PID     = PL1.PID;			// Copy PL group ID
  Ppts    = PL1.Ppts;			// Copy PL points
  PLdebug = PL1.PLdebug;		// Copy debugging flag
  Preduce = PL1.Preduce;		// Copy point reduction
  }

// ____________________________________________________________________________
//                      POLYLINE PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________

// ----------------------------- Pen Width Access -----------------------------

//
void FMPL::SetWidth(double pw)

        // Input	PL	: FM PolyLine (this)
        // 		pw	: Input pen width
	// Output	void	: The PL pen width is set to pw

  { PWidth = pw; }


double FMPL::GetWidth()

	// Input	PL	: FM PolyLine (this)
 	// Output	pw      : The current pen width

  { return PWidth; }


// ----------------------------- Pen Color Access -----------------------------


void FMPL::SetColor(int pc)

        // Input	PL	: FM PolyLine (this)
        // 		pc	: Input pen color
	// Output	void	: The PL pen color is set to pc

  { PColor = pc; }


int FMPL::GetPen() { return PColor; }

	// Input	PL	: FM PolyLine (this)
 	// Output	pc      : The current pen color



// --------------------------- PolyLine Fill Access ---------------------------


void FMPL::SetFill(int pf) { PFill = pf; }

        // Input	PL	: FM PolyLine (this)
        // 		pf	: Input pen fill
	// Output	void	: The PL pen fill is set to pf



int FMPL::GetFill() { return PFill; }

	// Input	PL	: FM PolyLine (this)
 	// Output	pf      : The current pen fill



// ---------------------------- PolyLine ID Access ----------------------------


void FMPL::SetID(int pid)

        // Input	PL	: FM PolyLine (this)
        // 		pid	: Input ID
	// Output	void	: The PL ID is set to pid

  { PID = pid; }


int FMPL::GetID()

	// Input	PL	: FM PolyLine (this)
 	// Output	PID     : The current PL ID

  { return PID; }


// -------------------------- PolyLine Point Access ---------------------------


int FMPL::GetPoints()

	// Input	PL	: FM PolyLine (this)
 	// Output	Ppts	: The current # points

  { return Ppts; }



// ----------------------- PolyLine Debug Flag Access -------------------------


int FMPL::GetDebug()

        // Input        PL      : FM PolyLine (this)
        // Output       PLdebug : The current debugging level
 
  { return PLdebug; }
 

void FMPL::SetDebug(int db)

        // Input        PL      : FM PolyLine (this)
        //              db      : A debugging level 
        // Output       void    : Debugging level set to db
 
  { PLdebug = db; }


// ------------------------- PolyLine Data Reduction --------------------------

 
int FMPL::GetReduce()

        // Input        PL      : FM PolyLine (this)
        // Output       Preducd : The current reductoin level

  { return Preduce; }


void FMPL::SetReduce(int rf)

        // Input        PL      : FM PolyLine (this)
        //              rf      : Data reduction flag
        // Output       void    : Reduction level set to rf

  { Preduce = rf; }



// ____________________________________________________________________________
//                       POLYLINE AUXILARY FUNCTIONS
// ____________________________________________________________________________


void FMPL::Set(double pw, int pc, int pf, int id)

	// Input	PL	: FM PolyLine (this)
	//		pw	: Pen width (pts)
	//		pc	: Pen color
	//		pf	: PolyLine fill
	//		id	: PolyLine ID
 	// Output	void	: PL parameters are set

  {
  PWidth = pw;
  PColor = pc;
  PFill = pf;
  PID = id;
  }


string FMPL::Colors(int pc) const

	// Input	PL	: FM PolyLine (this)
	//		pc	: Pen color
	//		SPC	: string for Pen color

  {
  switch(pc) 
    {
    case 0: return string("Black");   break;
    case 1: return string("White");   break;
    case 2: return string("Red");     break;
    case 3: return string("Blue");    break;
    case 4: return string("Green");   break;
    case 5: return string("Cyan");    break;
    case 6: return string("Magenta"); break;
    case 7: return string("Yellow");  break;
    }
  return string("Black");
  }


string FMPL::Filling(int pf) const

	// Input	PL	: FM PolyLine (this)
	//		pf	: PolyLine fill
	//		SPC	: string for PL fill

  {
  switch(pf) 
    {
    case 7:  return string("Solid"); break;
    case 15: return string("Clear"); break;
    default: return string(Gdec(pf)); break;
    }
  return string(Gdec(pf));
  }


string FMPL::Reduction() const

	// Input	PL	: FM PolyLine (this)
	//		SRed	: string for point reduction

  {
  switch(Preduce) 
    {
    case 0: return string("None"); break;
    case 1: return string("Horizonal"); break;
    case 2: return string("Vertical"); break;
    case 3: return string("Linear"); break;
    }
  return string("None");
  }


// ____________________________________________________________________________
//                    POLYLINE POINT MANIPULATION FUNCTIONS
// ____________________________________________________________________________


void FMPL::Zero() { Ppts = 0; }      
 
        // Input        PL      : FM PolyLine (this)
        //              z       : PolyLine point
        // Return       void    : The point z is added to PL
 


void FMPL::AddPt(const complex& z)
 
        // Input        PL      : FM PolyLine (this)
        //              z       : PolyLine point
        // Return       void    : The point z is added to PL
 
 
  {      
  if(Ppts >= size())                            // Insure we haven't added too
    {                                           // many points to the PL
    PLerror(1, 1);                              //      PL dimen. exceeded
    PLfatality(2);                              //      Can't add this point
    }
  put(z, Ppts);					// Add the point
  Ppts++;                                       // Update the point count
  }


void FMPL::WriteMIF(ostream& ostr)
 
        // Input        PL      : FM PolyLine (this)
        //              ostr    : An output stream
        // Output       void    : The PolyLine is written into 
        //                        the output stream ostr in MIF format 
 
  {
  if(Ppts < 2)					// Cannot output PL if
    {						// less than 2 points
    PLerror(3,1);				// 	Note < 2 points
    PLerror(4);					//	Note cannot output
    return;					//	Exit without output
    }
  if(!ostr.good())				// Insure ostream OK
    {
    PLerror(5,1);				// 	Bad ostream
    PLfatality(4);				//	Can't write to ostr
    }
  ostr << "  <PolyLine\n";
  ostr << Gform( "    <GroupID %i>\n", PID);
  ostr << Gform( "    <Pen %i>\n", PColor);
  ostr << string("    <ObColor Colors(PColor)>\n");
  ostr << Gform( "    <PenWidth %3.1f pt>\n", PWidth); 
  ostr << Gform( "    <Fill %i>\n", PFill); 
  int npts = Ppts;
  int npl = Ppts;
  double x,y; 
double VRES = 0.001;
  int k=0;
  switch(Preduce)
    {
    case 0:					// No data reduction
    default:
      ostr << Gform ( "    <NumPoints %i>\n", Ppts); 
      for(k=0; k<Ppts; k++) 
        { 
        x = getRe(k); 
        y = getIm(k); 
        ostr << Gform("    <Point %3.3f cm ", x); 
        ostr << Gform("%3.3f cm>\n", y); 
        } 
      break;
    case 1:					// Horizontal data reduction
      FMPL PL(Ppts);
      double lastx, lasty;			//	Track last pts
      int skip = 0;				//	Flag for point skip
      lastx = getRe(0); 			//	Get first x
      lasty = getIm(0); 			//	Get first y
      PL.AddPt(get(0));				// 	Store the first point
      for(k=1; k<Ppts; k++) 
        { 
        x = getRe(k); 
        y = getIm(k); 
        if(fabs(y-lasty) < VRES)
          {
          skip = 1;
          lastx = x;
          lasty = y;
          }
        else
          {
          if(skip)
            {
            PL.AddPt(complex(lastx,lasty));	// 	Store this point
            skip = 0;
            }
          PL.AddPt(complex(x,y));		// 	Store this point
          lastx = x;
          lasty = y;
          }
        } 
      if(skip)					// Make sure last point
        PL.AddPt(complex(lastx,lasty));		// is not skipped
      ostr << Gform ( "    <NumPoints %i>\n", PL.Ppts); 
      for(k=0; k<PL.Ppts; k++) 
        { 
        x = PL.getRe(k); 
        y = PL.getIm(k); 
        ostr << Gform("    <Point %3.3f cm ", x); 
        ostr << Gform("%3.3f cm>\n", y); 
        } 
      npl = PL.Ppts;
      break;
    }
  ostr << "    > # end of PolyLine                   \n"; 
  ostr.flush();
  if(PLdebug) 
    { 
    cout << "\n\t\t" << npl << " Point PolyLine Output";
    cout << "\n\t\t  ID = "    << PID    << ", Pen Color = "  << Colors(PColor) 
         << "\n\t\t  Width = " << PWidth << ", Fill = " << Filling(PFill) 
         << "\n\t\t  Points = " << Ppts; 
    cout << "\n\t\t  First Point " << get(0); 
    cout << "\n\t\t  Last Point  " << get(Ppts-1);
    if(Preduce)
      {
      cout << "\n\t\t  Applied " << Reduction() << " Reduction";
      cout << "\n\t\t  " << npts-npl << " Points Cropped in Reduction";
      }
    } 
  cout.flush();
  }


// ____________________________________________________________________________
// Z                  FM POLYLINE STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                PL	: FM PolyLine (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                the FrameMaker PolyLine

ostream& FMPL::print(ostream& ostr) const
  {
  ostr << "\n\tPolyLine Pen Width: " << PWidth;
  ostr << "\n\tPolyLine Pen Color: " << Colors(PColor);
  ostr << "\n\tPolyLine Pen Fill:  " << Filling(PFill);
  ostr << "\n\tPolyLine Group ID:  " << PID;
  ostr << "\n\tPolyLine Points:    " << Ppts;
  ostr << "\n\tPolyLine Size:      " << size();
  ostr << "\n\tPolyLine Reduction: " << Reduction();
  ostr << "\n";
  return ostr;
  }
 


        // Input                ostr    : An output stream
        //                      PL	: FM PolyLine (this)
        // Output               ostr    : The output stream modified by
        //                                the FM PolyLine

ostream& operator<< (ostream& ostr, const FMPL& PL)
  {                                                                 
  PL.print(ostr);
  return ostr;
  }
 

#endif 						// FrameMakerC.cc
