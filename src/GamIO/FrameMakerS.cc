/* FrameMakerS.cc ***********************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
*      FrameMaker Stack Plots                    Implementation		**
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

#ifndef   FrameMakerS_cc_			// Is file already included?
#  define FrameMakerS_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/FrameMakerS.h>			// Include the header
#include <GamIO/FrameMakerC.h>			// Include FM MIF constucts
#include <GamIO/FrameMakerM.h>			// Include FM MIF functions
#include <GamIO/FrameMakerP.h>			// Include FM parameters
#include <GamIO/FrameMaker.h>			// Include FM 1D plotting
#include <Basics/Gconstants.h>			// Include PI
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <Basics/StringCut.h>			// Include Gdec and Gform
#include <iostream>				// Include input output streams
#include <string>				// Include libstdc++ STL strings

using std::string;				// Using libstdc++ strings
using std::ostream;				// Using libstdc++ output streams
using std::ofstream;				// Using libstdc++ output file streams
using std::cout;				// Using libstdc++ standard output

// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
//             CLASS RowBlock: FRAMEMAKER STACK PLOT ROW
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
 
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS ROWBLOCK ERROR HANDLING
// ____________________________________________________________________________


void RowBlk::RBerror(int eidx, int noret) const

	// Input                Blk     : Stack plot row block (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

                                                        
  {
  cout << "\nFrameMaker Stack Row Block: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      cout << "Program Aborting.....";
      break;
    case 1:								// (1)
      cout << "Allocated Block Dimension Exceeded!";
      break;
    case 2:								// (2)
      cout << "Cannot Add Point To Row Block.";
      break;
    case 5:								// (5)
      cout << "Problems Adding Downward Trend Point From Area 0 in Block";
      break;
    case 6:								// (6)
      cout << "Problems Adding Downward Trend Point From Area 1 in Block";
      break;
    case 7:								// (7)
      cout << "Problems Adding Downward Trend Point From Area 2 in Block!?";
      break;
    case 8:								// (8)
      cout << "Problems Adding Upward Trend Point From Area 2 in Block";
      break;
    case 9:								// (9)
      cout << "Problems Adding Upward Trend Point From Area 1 in Block";
      break;
    case 10:								// (10)
      cout << "Problems Adding Upward Trend Point From Area 0 in Block!?";
      break;
    default:
      cout << "Unkown Error, Number " << eidx;
    }
  if(!noret) cout << ".\n";
  }  

 
volatile void RowBlk::RBfatality(int eidx) const

	// Input                Blk     : Stack plot row block (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped

  {
  RBerror(eidx, 1);				// Output error message
  if(eidx) RBerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                   ROW BLOCK CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
RowBlk::RowBlk( ) : row_vector()

	// Input                Blk     : Stack plot row block (this)
        // Output               none    : Stack plot row block started

  {                                                             
  Bareas = NULL;			// No array of area flags
  Bpts = 0;				// No points in block
  Bpts = 0;				// There are no points
  Brow = 0;				// There is no row designation
  Bhoff = 0;				// Nope, no horizontal offset
  Bvoff = 0;				// Nor vertical offset 
  Bhsf = 1;				// Horizontal scaling (cm/pt)
  Bvsf = 1;				// Vertical scaling (cm/absunits)
  Bdebug = 0;				// Set for no debugging
  BaseID = 21;				// Default FM ID
  }
 

RowBlk::RowBlk(int npts) : row_vector(npts, complex0)

	// Input                Blk     : Stack plot row block (this)
        // Output               none    : Stack plot row block started

  {                                                             
  Bareas = new int[npts];		// Start the block area flags
  Bpts = 0;				// There are no points
  Brow = 0;				// There is no row designation
  Bhoff = 0;				// Nope, no horizontal offset
  Bvoff = 0;				// Nor vertical offset 
  Bhsf = 1;				// Horizontal scaling (cm/pt)
  Bvsf = 1;				// Vertical scaling (cm/absunits)
  Bdebug = 0;				// Set for no debugging
  BaseID = 21;				// Default FM ID
  }

                                                            
void RowBlk::operator=(const RowBlk& Blk1)

	// Input                Blk     : Stack plot row block (this)
        //                      Blk1    : Stack plot row
        // Output               void    : FrameMaker stack plot row Blk1
        //                                is copied into Blk

  {
  row_vector::operator=(Blk1);		// Copy the row vector
  if(Bareas) delete [] Bareas;		// Remove any old areas
  Bareas = new int[size()];		// Allocate space for areas
  for(int i=0; i<Blk1.Bpts; i++)	// Copy all the area flags
    Bareas[i] = Blk1.Bareas[i];
  Bpts  = Blk1.Bpts;			// Copy the number of points
  Brow  = Blk1.Brow;			// Copy row index
  Bhoff = Blk1.Bhoff;			// Copy horizontal offset (cm)
  Bvoff = Blk1.Bvoff;			// Copy vertical offset (cm)
  for(int j=0; j<4; j++)		// Copy the possible block IDs
    Bids[j] = Blk1.Bids[j];
  Bhsf  = Blk1.Bhsf;			// Copy the horizontal scaling
  Bvsf  = Blk1.Bvsf;			// Copy the vertical scaling
  for(int k=0; k<3; k++)		// Copy the block levels
    Blevs[k] = Blk1.Blevs[k];
  Bdebug = Blk1.Bdebug;			// Copy the debugging level
  BaseID = Blk1.BaseID;			// Copy the FM base ID
  }
  
    
RowBlk::~RowBlk()
 
	// Input                Blk     : Stack plot row block (this)
        // Output               none    : Blk is destructed
 
  { if(Bareas) delete [] Bareas; }
 
 
// ____________________________________________________________________________
// B                       ROW BLOCK ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// ---------------------------- Row Index Access ------------------------------

void RowBlk::SetRow(int row)

	// Input		Blk	: Stack plot row block (this)
	//			row	: Row index
 	// Output		void	: Block row index set to row
	// Note				: This will set the IDs also

  {
  Brow = row;
  SetIDs();
  }


int RowBlk::GetRow() const

	// Input		Blk	: Stack plot row block (this)
 	// Output		Brow	: Block row index

  { return Brow; }
 
 
// ----------------------------- Row ID Access -------------------------------


void RowBlk::SetIDs()

	// Input		Blk	: Stack plot row block (this)
        // Return		void    : Possible block ID values are set
	// Note				: This must match stack plot #s

  {
  for(int area=0; area<4; area++)
    Bids[area] = 5*Brow + area + BaseID;
  }

 
// --------------------------- Row Levels Access -----------------------------


void RowBlk::SetLevels(double* Levels)

	// Input		Blk	: Stack plot row block (this)
        // Return		void    : Possible block ID values are set

  {
  for(int lev=0; lev<4; lev++)
    Blevs[lev] = Levels[lev];
  }

// -------------------------- Row Offsets Access -----------------------------


void RowBlk::SetOffsets(double ho, double vo)

	// Input		Blk	: Stack plot row block (this)
	//			ho	: Horizontal offset (cm)
	//			vo	: Vertical offset (cm)
 	// Output		void 	: Block offsets set (cm)

  {
  Bhoff = ho;                  		// Set horizontal offset (cm)
  Bvoff = vo; 				// Set vertical offset (cm)
  }
 

double RowBlk::VOffset() const

	// Input		Blk	: Stack plot row block (this)
 	// Output		Bvoff	: Block row vertical offset (cm)

  { return Bvoff; }
 

double RowBlk::HOffset() const

	// Input		Blk	: Stack plot row block (this)
 	// Output		Bhoff	: Block row horizontal offset (cm)

  { return Bhoff; }
 

// ---------------------- Row Scaling Factor Access --------------------------


void RowBlk::SetScaling(double HSF, double VSF)

	// Input		Blk	: Stack plot row block (this)
	//			ho	: Horizontal offset (cm)
	//			vo	: Vertical offset (cm)
 	// Output		void 	: Block offsets set (cm)

  {
  Bhsf = HSF;                  		// Set horizontal scaling (cm/pt)
  Bvsf = VSF; 				// Set vertical offset (cm/absunits)
  }
 

double RowBlk::VScale() const

	// Input		Blk	: Stack plot row block (this)
 	// Output		Bvsf	: Block vertical scaling (cm/asbunits)

  { return Bvsf; }
 

double RowBlk::HScale() const

	// Input		Blk	: Stack plot row block (this)
 	// Output		Bhsf	: Block horiz. scaling (cm/asbunits)

  { return Bhsf; }
 

// ------------------------ RowBlk Debug Flag Access --------------------------


int RowBlk::GetDebug()

	// Input		Blk	: Stack plot row block (this)
        // Output       	Bdebug : The current debugging level

  { return Bdebug; }


void RowBlk::SetDebug(int db)

	// Input		Blk	: Stack plot row block (this)
        //			db      : A debugging level
        // Output       	void    : Debugging level set to db

  { Bdebug = db; }


// ____________________________________________________________________________
// C                     ROW BLOCK DETERMINATION FUNCTIONS
// ____________________________________________________________________________
 

void RowBlk::Zero()

	// Input		Blk	: Stack plot row block (this)
        // Return		void    : The point count is set to zero

  { Bpts = 0; }


void RowBlk::AddPt(int col, const matrix& mx)

	// Input		Blk	: Stack plot row block (this)
	//			col	: Column index
        // 			mx	: FM stack plot array (this)
        // Return		void    : The point <row|mx|col> is added to
	//				  the block Blk in plot coordinates

  {
  if(Bpts >= size())				// Insure we haven't added too
    { 						// many points to the block
    RBerror(1, 1);				//	Block dimen. exceeded
    RBfatality(2);				//	Can't add this point
    }
  double x = Bhoff + col*Bhsf;			// Compute x coordinate (cm)
  lasty = Bvoff - mx.getRe(Brow,col)*Bvsf;	// Compute y coordinate (cm)
  put(complex(x,lasty), Bpts); 			// Add the point
  if(lasty <= Blevs[1])      Bareas[Bpts] = 0;	// Area y on or above row
  else if(lasty <= Blevs[2]) Bareas[Bpts] = 1;	// Area y in base plane
  else                       Bareas[Bpts] = 2;	// Area y below base plane 
  lastarea = Bareas[Bpts];			// Store this as last area
  Bpts++;					// Update the point count
  }
 

void RowBlk::AddPt(complex& pt)

	// Input		Blk	: Stack plot row block (this)
	//			pt	: Point (FM units cm)
        // Return		void    : The point pt is added to
	//				  the block Blk in plot coordinates

  {
  if(Bpts >= size())				// Insure we haven't added too
    { 						// many points to the block
    RBerror(1, 1);				//	Block dimen. exceeded
    RBfatality(2);				//	Can't add this point
    }
  lasty = zIm(pt);				// Set last y coordinate (cm)
  put(pt, Bpts);		 		// Add the point
  if(lasty <= Blevs[1])      Bareas[Bpts] = 0;	// Area y on or above row
  else if(lasty <= Blevs[2]) Bareas[Bpts] = 1;	// Area y in base plane
  else                       Bareas[Bpts] = 2;	// Area y below base plane 
  lastarea = Bareas[Bpts];			// Store this as last area
  Bpts++;					// Update the point count
  }


void RowBlk::AddPt(complex& pt, int area)

	// Input		Blk	: Stack plot row block (this)
	//			pt	: Point (FM units cm)
	//			area	: Point area flag
        // Return		void    : The point pt is added to
	//				  the block Blk in plot coordinates
	//				  in the area specified

  {
  if(Bpts >= size())				// Insure we haven't added too
    { 						// many points to the block
    RBerror(1, 1);				//	Block dimen. exceeded
    RBfatality(2);				//	Can't add this point
    }
  lasty = zIm(pt);				// Set last y coordinate (cm)
  put(pt, Bpts);		 		// Add the point
  Bareas[Bpts] = area;				// Set the point area
  lastarea = area;				// Store this as last area
  Bpts++;					// Update the point count
  }


complex RowBlk::GetPt(int col, const matrix& mx) const

	// Input		Blk	: Stack plot row block (this)
	//			col	: Column index
        // 			mx	: FM stack plot array (this)
        // Return		void    : The point <row|mx|col> is 
	//				  returned in plot coordinates
	// Note				: The block is unaffected

  {
  double x = Bhoff + col*Bhsf;			// Compute x coordinate (cm)
  double y = Bvoff - mx.getRe(Brow,col)*Bvsf;	// Compute y coordinate (cm)
  return complex(x,y);				// Return the point
  }


int RowBlk::GetArea(complex& z) const

	// Input		Blk	: Stack plot row block (this)
        // 			z	: Block point
        // Return		area	: An area flag for the point is 
	//				  returned
	// Note				: The block is unaffected

  {
  double y = zIm(z);				// Get y coordinate (cm)
  if(y <= Blevs[1])      return 0;		// On or above row height
  else if(y <= Blevs[2]) return 1;		// In base plane (BP)
  else                   return 2;		// Completely below BP
  }


PtTrend RowBlk::Trend(int area) const

	// Input		Blk	: Stack plot row block (this)
        // 			area	: Area flag
        // Return		PtTrend : Row trend if area is associated
	//				  with the next pt in Blk
	// Note				: The block is unaffected

  {
  if(lastarea < area)			// If lastarea less than area
    return _down;			// then the row goes down 
  if (lastarea > area)			// If laster area greater than
    return _up;				// area then the row goes up
  return _same;
  }
 

complex RowBlk::ExtraPt(complex& pt1, double y)

	// Input	Blk	: Stack plot row block (this)
 	// 		pt1	: A point (FM units, cm)
	//		y 	: Desired y coordinate (cm)
	// Output	pt	: Point on the line between pt1
	//			  last point of pt in Blk having
	//			  a height of y

//			   y1-y2                x1-x2
// Line: y-y2 = m*(x-x2) = ----- (x-x2) --> x = ----- (y-y2) + x2
//                         x1-x2                y1-y2

  {
  complex pt2 = get(Bpts-1);
  double m = (zIm(pt1)-zIm(pt2))/(zRe(pt1)-zRe(pt2));
  return complex(zRe(pt2)+(y-zIm(pt2))/m, y);
  }

void RowBlk::Fill(int cols, const matrix& mx) 

	// Input		Blk	: Stack plot row block (this)
        //                      cols    : Number of matrix columns
        //                      mx      : Data array for stack plot  
 	// Output		void	: Block is filled for mx row

  {
  AddPt(0, mx);					// Add the first point
  complex pt, ptx;				// For working data points
  int area;					// Area flag, # of added points
  double ylev, ylevup, ylevdwn;			// Working plot vertical levels
  for(int col=1; col<cols; col++)		// Loop full row (past 1st pt)
    {
    pt = GetPt(col, mx);			// Get point for this column
    area = GetArea(pt);				// Get area flag this column
    ylev = Blevs[area];				// Level line below area
    switch(Trend(area))				// Now, add this pt depending
      {						// on where it is in the plot
      case _down:				// Point is in lower area
        {
        switch(lastarea)			//	How to add pt depends
          {					//	on start & finish areas
          case 0:				// 	Start in area 0
            {
            ylevdwn = Blevs[1];			//	Height at area 0 bottom
            if(area == 1)			//	New point is in area 1
              {					//	Get level below area 0
              if(lasty != ylevdwn)		//	if last pt isn't on the
                {				//	level line, project pt
                ptx = ExtraPt(pt,ylevdwn);	//	that's on level line
                AddPt(ptx);			//	and add it (area is 0)
                }
              AddPt(pt);			//	Now add new pt (area 1)
              }
            else if(area == 2)
              {
              if(lasty != ylevdwn)		//	If last pt isnt' right
                {				//	level line, project pt
                ptx = ExtraPt(pt,ylevdwn);	//	that's on level line
                AddPt(ptx);			//	and add it (area is 0)
                }
              ylevdwn = Blevs[2];		//	Get level below area 1
              ptx = ExtraPt(pt,ylevdwn);	//	On level under area 1
              AddPt(ptx);			//	Add pt on level below
              AddPt(pt);			//	Add new pnt (area 2)
              }
            else RBerror(5);			//	Shouldn't get to here
            }
            break;
          case 1:				//	Start in area 1
            {
            ylevdwn = Blevs[2];			//	Height at area 1 bottom
            if(area == 2)			//	New point is in area 2
              {					//	Get level below area 1
              if(lasty != ylevdwn)		//	If last pt not just on
                {				//	level line, project pt
                ptx = ExtraPt(pt,ylevdwn);	//	that's on level line
                AddPt(ptx);			//	and add it (area is 1)
                }
              AddPt(pt);			//	Add new point (area 2)
              }
            else RBerror(6);			//	Shouldn't get to here
            }
            break;
          case 2:				// 	Start in area 2
          default:				//	Start in unknown area
            RBerror(7);				//	Shouldn't get to here?
            break;				//	Cannot go "down" from
          }					//	bottom area!
        }
        break;
      case _up:					// Point is in an above area
        {
        switch(lastarea)			//	How to add pt depends
          {					//	on start & finish areas
          case 2:				// 	Start in area 2
            {
            ylevup = Blevs[2];			//	Height at area 1 bottom
            if(area == 1)			//	New point is in area 1
              {					//	Get level below area 0
              if(zIm(pt) != ylevup)		//	if new pt isnt' just on
                {				//	level line, project pt
                ptx = ExtraPt(pt,ylevup);	//	that's on level line
                AddPt(ptx);			//	and add it (area is 1)
                }
              AddPt(pt);			//	Add new pt (area is 1)
              }
            else if(area == 0)			//	New point is in area 0
              {
              ptx = ExtraPt(pt,ylevup);		//	Proj. up (below area 1)
              AddPt(ptx);			//	and add it (area is 1)
              ylevup = Blevs[1];		//	Height at area 0 bottom
              if(zIm(pt) != ylevup)		//	If new pt isnt' just on
                {				//	level line, project pt
                ptx = ExtraPt(pt,ylevup);	//	that's on level line
                AddPt(ptx);			//	and add it (area is 0)
                }
              AddPt(pt);			//	Now add new pt (area 0)
              }
            else RBerror(8);			//	Shouldn't get to here
            }
            break;
          case 1:				//	Start in area 1
            {
            ylevup = Blevs[1];			//	Height at area 0 bottom
            if(area == 0)			//	New point is in area 0
              {
              if(zIm(pt) != ylevup)		//	If new pt isn't just on
                {				//	level line, project up
                ptx = ExtraPt(pt,ylevup);	//	to level below area 0
                AddPt(ptx);			//	Add proj. pt (area 0)
                }
              AddPt(pt);			//	Now add new pt (area 0)
              }
            else RBerror(9);			//	Shouldn't get to here
            }
            break;
          case 0:				// 	Start in area 0
          default:				//	Start in unknown area
            RBerror(10);			//	Shouldn't get to here?
            break;				//	Cannot go "up" from
          }					//	top area!
        }
        break;
      case _same:				// It the new point doesn't
      default:					// change its area, then just
        AddPt(pt);				// add it to to Block
        break;
      }
    }

// One Last Case To Deal With.  If The First Point Is Exactly On The  
// Level Line Below It's Area and the 2nd Point Is In The Next Lower
// Area ---> The 1st Point is Isolated.  If That Happens We Just Set
// The 1st Points Area Flag To Be In The Next Lower Level (With 2nd Pt).
// sosix

  if(getIm(0) == Blevs[Bareas[0]+1])
    if(Bareas[1] == Bareas[0]+1)
      (Bareas[0])++;
  }


void RowBlk::Plot(const string& filename) 

	// Input	Blk	: Stack plot row block (this)
        //		filename: Output MIF filename
        // Return	void 	: Block Blk, a stack plot row,
	//			  is put into a file filename
	//			  in FrameMaker MIF format.
	// Note			: In principle, the block Blk
	//			  switched area only after a
	//			  point on a level is reached

  {
  double pw = 1.0;				// For pen width 1 pt 
  int pf = 15;					// For pen fill to none
  int pc[4] = { 3, 4, 5, 6 };			// For pen colors
  int FMID = 11;				// Frame ID in plot
  double hsize = 15.0;				// Plot horizontal size
  double vsize = 10.0;				// Plot vertical size

//row_vector x = get_block(0,0,1,Bpts);
//GP_1D("crap.asc", x, -1);
//GP_1Dplot("crap.gnu", "crap.asc");
//for(int t=0; t<Bpts; t++)
//x.put(complex(Bareas[t]), t);
//GP_1D("crap.asc", x);
//GP_1Dplot("crap.gnu", "crap.asc");

  FMPL PL(Bpts);				// PolyLine for FM output
  PL.SetDebug(Bdebug);				// Set debugging level
  PL.SetReduce(1);				// Set data reduction
  ofstream ostr(filename.c_str());		// Define output stream
  FM_Begin(ostr);				// FrameMaker begin comment
  FM_AFrames_Begin(ostr);			// Anchored frame begin
  FM_AFrame_Set(ostr,hsize,vsize,FMID);		// Set anchored frame
  complex pt = get(0);				// Working data point
  PL.AddPt(pt);					// Set the 1st point
  int area, prevarea = Bareas[0];		// Set 1st pt area flag
  PL.Set(pw,pc[prevarea],pf,Bids[prevarea]);	// Set PL pen,area,fill
  for(int j=1; j<Bpts; j++)			// Loop columns of this block
    {   
    area = Bareas[j];				// Get "area" of point
    if(prevarea != area)			// Point doesn't change area
      {						// from the previous one.
      PL.WriteMIF(ostr);	 		// Output this PolyLine
      PL.Zero();				// Zero the PL
      PL.AddPt(pt);				// Set the 1st point
      prevarea = area;				// Set new area
      PL.Set(pw,pc[prevarea],pf,Bids[prevarea]);// Set PL pen,area,fill
      }
    pt = get(j);				// This is the next pt
    PL.AddPt(pt);				// Set the next point
    }   
  PL.WriteMIF(ostr);		 		// Output any end PolyLine
  int a;
  for(a=0; a<3; a++) FM_Group(ostr,Bids[a]);	// Group PolyLines by area
  for(a=0; a<3; a++)				// Group PolyLines together 
    FM_Group(ostr,Bids[a],Bids[3]);	
  FM_AFrame_End(ostr);				// Frame end
  FM_AFrames_End(ostr);				// Anchored Frame end
  FM_ParaText_End(ostr);			// Create TextFlow with Frame
  FM_End(ostr);					// Frame End
  }  


// ____________________________________________________________________________
// D                      ROW BLOCK STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


ostream& RowBlk::print(ostream& ostr) const

	// Input		Blk	: Stack plot row block (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                the FrameMaker block parameters

  {
  ostr << "\n\t\t\tStack Plot Row Block\n";
  ostr << "\n\t\tBlock Associated Row:   " << Brow;
  ostr << "\n\t\tBlock Points:           " << Bpts;
  ostr << "\n\t\tBlock Vector Size:      " << size();
  ostr << "\n\t\tBlock Horiz. Offset:    " << Bhoff << " cm";
  ostr << "\n\t\tBlock Vertical Offset:  " << Bvoff << " cm";
  ostr << "\n\t\tBlock Horiz. Scaling:   " << Bhsf << " cm/pt";
  ostr << "\n\t\tBlock Vertical Scaling: " << Bvsf << " cm/absunits";
  ostr << "\n\t\tBlock ID Area 0:        " << Bids[0];
  ostr << "\n\t\tBlock ID Area 1:        " << Bids[1];
  ostr << "\n\t\tBlock ID Area 2:        " << Bids[2];
  ostr << "\n\t\tBlock Line ID:          " << Bids[3];
  ostr << "\n\t\tBlock Level 0:          " << Blevs[0] << " cm";
  ostr << "\n\t\tBlock Level 1:          " << Blevs[1] << " cm";
  ostr << "\n\t\tBlock Level 2:          " << Blevs[2] << " cm";
  ostr << "\n\t\tBlock Level 3:          " << Blevs[3] << " cm";
  ostr << "\n";
  return ostr;
  }


ostream& operator<<(ostream& ostr, const RowBlk& Blk)

        // Input                ostr    : An output stream
	// 			Blk	: Stack plot row block (this)
        // Output               none    : Modifies output stream
        ///F_list <<                    - FrameMaker row block sent to ostr

  {                                                                 
  Blk.print(ostr);
  return ostr;
  }
                                                    

// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
//             CLASS FMStack: FRAMEMAKER STACK PLOTTING PARAMETERS
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
 

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  CLASS FM STACK PARAMETER ERROR HANDLING
// ____________________________________________________________________________
 

void FMStack::FMSTKerror(int eidx, int noret) const

        // Input                FMSTK   : FrameMaker stack parameters (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 

  {                                                     
  cout << "\nFrameMaker Stack Parameters: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      cout << "Program Aborting.....";
      break;
    case 1:								// (1)
      cout << "Stack Plot |yinc| too Large For Plot Width! Adjusting";
      break;
    case 2:								// (2)
      cout << "Stack Plot |xinc| too Large For Plot Height! Adjusting";
      break;
    case 5:								// (5)
      cout << "Data Vector Found to Contain Only 0 or 1 point.";
      break;
    case 10:								// (10)
      cout << "Having Trouble With The Output Stream!";
      break;
    case 11:								// (11)
      cout << "Cannot Write Data To Stack Plot File!";
      break;
    case 49:                                                            // (49)
      cout << "Submitted Array For Plot Contains < 2 Columns!";
      break;
    case 50:                                                            // (50)
      cout << "Submitted Array For Plot Contains < 2 Rows!";
      break;
    case 51:                                                            // (51)
      cout << "Cannot Produce A Stack Plot.";
      break;
    case 60:								// (60)
      cout << "Cannot Find Stack Plot Vertical Scaling Iteratively!";
      break;
    default:
      cout << "Unkown Error, Number " << eidx;
    }
  if(!noret) cout << ".\n";
  }  

 
volatile void FMStack::FMSTKfatality(int eidx) const

        // Input                FMSTK   : FrameMaker stack parameters (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped

  {
  FMSTKerror(eidx,1);				// Output error message
  if(eidx) FMSTKerror(0);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                 FM STACK PARAMETERS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 

FMStack::FMStack() : FMPar()			// Call base constructor

        // Input                none    : No arguments required
        // Output               none    : FrameMaker stack parameters
        ///F_list FMSTK                 - Constructor

  {
  rows = 0;                             // Row span of plotted array
  cols = 0;                             // Column span of plotted array
  vinc = 0.1;                           // Stack vertical increment (cm)
  hinc = 0.1;                           // Stack horizontal increment (cm)
  rowinc = 1;                           // Stack row increment
  hdelta = 0;				// Total horizontal offset over rows
  vdelta = 0;				// Total vertical offset over rows
  hwidth = 0;				// Actual plotted data width
  vheight = 0;				// Actual plotted data height
  HSF = 1;				// Horizontal scaling (cm/abs.units)
  VSF = 1;				// Vertical scaling (cm/abs.units)
  Bback = 0;				// Base plane back coord. (cm)
  Bfront = 0;				// Base plane front coord. (cm)
  BbackW = 0;				// Base plane back, west (cm)
  BbackE = 0;				// Base plane back, east (cm)
  BfrontW = 0;				// Base plane front, east (cm)
  BfrontE = 0;				// Base plane front, east (cm)
  BPID = 10;				// Base plane ID default is 10
  gridding = 0;                         // Column gridding flag
  colinc = 1;                           // Stack column grid increment
//  Areas = NULL;				// Insure the areas are not set				
  FRID = 11;				// Anchored frame ID value
//  bpts = 0;				// No points in Block
  bhoff = 0;				// No block horizontal offset
  bvoff = 0;				// No block vertical offset
  HLF = 1;				// Set to use HLA
  }
                                                      

FMStack::~FMStack()

        // Input                FMSTK   : FrameMaker stack parameters (this)
        // Output               none    : FMSTK is destructed

 {  }
// { if(Areas) delete [] Areas; }

                                                            
void FMStack::operator=(const FMStack& FMSTK1)

        // Input                FMSTK	: FrameMaker stack parameters (this)
        //                      FMSTK1  : FrameMaker stack parameters
        // Output               void    : FrameMaker stack parameters FMSTK1
        //                                copied into FMSTK
        ///F_list =                     - Assignment

  {
  rows     = FMSTK1.rows;	// Copy row span of plotted array
  cols     = FMSTK1.cols;	// Copy column span of plotted array
  vinc     = FMSTK1.vinc; 	// Copy vertical increment (cm/row)
  hinc     = FMSTK1.hinc;	// Copy horizontal increment (cm/row)
  rowinc   = FMSTK1.rowinc;	// Copy stack row increment
  hdelta   = FMSTK1.hdelta;	// Copy total horizontal offset over rows
  vdelta   = FMSTK1.vdelta;	// Copy total vertical offset over rows
  hwidth   = FMSTK1.hwidth;	// Copy actual plotted data width
  vheight  = FMSTK1.vheight;	// Copy actual plotted data height
  HSF      = FMSTK1.HSF;	// Copy horizontal scaling (cm/abs.units)
  VSF      = FMSTK1.VSF;	// Copy vertical scaling (cm/abs.units)
  Bback    = FMSTK1.Bback;	// Copy base plane back coord. (cm)
  Bfront   = FMSTK1.Bfront;	// Copy base plane front coord. (cm)
  BbackW   = FMSTK1.BbackW;	// Copy base plane back, west (cm)
  BbackE   = FMSTK1.BbackE;	// Copy base plane back, east (cm)
  BfrontW  = FMSTK1.BfrontW;	// Copy base plane front, east (cm)
  BfrontE  = FMSTK1.BfrontE;	// Copy base plane front, east (cm)
  gridding = FMSTK1.gridding;	// Copy column gridding flag
  colinc   = FMSTK1.colinc;	// Copy column grid increment
  HLF      = FMSTK1.HLF;	// Copy the HLA flag
  }
                                                     
// ____________________________________________________________________________
// B                    FM STACK PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________
 
 
// ---------------------------- Stack Plot Offsets ----------------------------
 
 
double FMStack::VInc() const { return vinc; }
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               vinc    : Vertical increment (cm)
 
 
void FMStack::VInc(double vi) { vinc = vi; }
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      vi      : Vertical increment (cm)
        // Output               vinc    : Vertical increment set to vi
 
 
double FMStack::HInc() const { return hinc; }
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               hinc    : Horizontal increment (cm)
 
               

 
void FMStack::HInc(double hi) { hinc = hi; }

        // Input                FMSTK   : FM stack parameters (this)
        //                      hi      : Horizontal increment (cm)
        // Output               hinc    : Horizontal increment set to hi
 
                  
// ----------------------- Stack Plot Row Increment ---------------------------
 
int FMStack::RInc() const
         
        // Input                FMSTK   : FM stack parameters (this)
        // Output               rowinc  : The stack plot row increment
 
  { return rowinc; }
 
 
void FMStack::RInc(int ri)
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      ri      : Row increment
        // Output               void    : The row increment is set to ri
 
  { rowinc = ri; }

                 
// ---------------------- Stack Plot Column Increment -------------------------
 
int FMStack::CInc() const
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               colinc  : The stack plot column increment
 
  { return colinc; }
 
 
void FMStack::CInc(int ci)
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      ci      : Column increment
        // Output               void    : The column increment is set to ci
 
  { colinc = ci; }

 
// -------------------------- Stack Plot HLA Usage ----------------------------
                                                                                
 
int FMStack::HLAF() const
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               HLF     : The hidden line use flag

  { return HLF; }
         
 
void FMStack::HLAF(int hf)
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      hf      : Hidden line use flag
        // Output               void    : The hidden line use flag is set
 
  { HLF = hf; }

                   
// ------------------------ Stack Plot Gridding Flag --------------------------
 
 
int FMStack::Gridding() const
 
        // Input                FMSTK   : FM stack parameters (this)
        // Output               gridding: The stack plot gridding flag
 
  { return gridding; }
 
 
void FMStack::Gridding(int gr)
 
        // Input                FMSTK   : FM stack parameters (this)
        //                      gr      : Gridding flag
        // Output               void    : The gridding flag is set to gr
 
  { gridding = gr; }


// ____________________________________________________________________________
// C                    STACK PARAMETER SETTING FUNCTIONS
// ____________________________________________________________________________

// ------------------------- Stack Plot Data Matrix ---------------------------

void FMStack::SetArray(const matrix& mx)

        // Input	FMSTK	: FM stack parameters (this)
        // 		mx	: Input data matrix to plot
	// Output	void	: The following values in FMSTK are
	//			  set according to the input matrix

	//				MX   = The array to be plotted
	//				rows = The # of rows to plot
	//				cols = The # of cols to plot

  {
  MX = mx;					// Store the data array
  rows = mx.rows();				// Set # of rows spanned
  cols = mx.cols();				// Set number of columns
  if(rows <= 1)					// Check matrix rows
    {						// If less than 2 cannot
    FMSTKerror(50); 				// make a stack plot
    FMSTKfatality(51);
    }
  if(cols <= 1)					// Check matrix columns
    {						// If less than 2 cannot
    FMSTKerror(49);				// make a stack plot
    FMSTKfatality(51);
    }
  }

// ------------------------ Stack Plot Row Offsets ---------------------------


void FMStack::SetOffsets()

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	void    : The vector of row offsets
	//			  is filled (in cm/row)
	// Note			: The offsets for row i will be
	//			  <offsets|i> in cm. The real part is  
	//			  the horizontal offset while the imaginary
	//			  part is the vertical offset

  {
  OffSets = row_vector(rows);			// Vector of x,y offsets (cm)
  double xoff, yoff;				// We'll keep {x,y} pairs for
  double ip1, im1;				// values (horiz., vert.)
  for(int i=0; i<rows; i++)
    {						// Note that this is done for
    ip1 = double(i+1);				// all rows whether they are
    im1 = double(i-1);				// plotted or not (gridding
    if(vinc>0) yoff =  ip1*vinc;		// 1st row down, last row up
    else       yoff = -im1*vinc;		// 1st row up, last row down
    if(hinc>0) xoff =  ip1*hinc;		// 1st row left, last row right
    else       xoff = -im1*hinc;		// 1st row right, last row left
    OffSets.put(complex(xoff, yoff), i);	// Store row i offset values
    }
  return;
  }


row_vector FMStack::GetOffsets()

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	Offsets	: The vector of row offsets
	//			  is returned (in cm/row)
	// Note			: The offsets for row i will be
	//			  <offsets|i> in cm. The real part is  
	//			  the horizontal offset while the imaginary
	//			  part is the vertical offset

  { return OffSets; }


// ------------------------- Stack Plot Row Maxima ---------------------------


void FMStack::SetMaxima()

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	void    : The vector of row max & min values
	// Note			: The maxima for row i will be
	//			  <maxmins|i> in absolute units. The
	//			  real part is the row maximum and  
	//			  the imaginary part is the row minimum.
	// Note			: Some rows may be flagged so that they
	//			  are not plotted (via row increment) 
	//			  Their maxima are still determined!

  {
  MaxMins = row_vector(rows, 0);		// Vector for max and min
  double pt, vxre, vxim;			// Values of each row
  int i, j;					// Indexing variables
  for(i=0; i<rows; i+=rowinc)			// Maxima are reals
    {						// Minima are imaginarys
    vxre = 1.e-50;				// Need these values so we
    vxim = 1.e50;				// get the true maxima
    for(j=0; j<cols; j++)
      {
      pt = MX.getRe(i,j);			//	Get point j
      if(pt > vxre) vxre = pt;			//	Check maximum
      if(pt < vxim) vxim = pt;			//	Check minimum
      }
    MaxMins.put(complex(vxre, vxim),i);		// Store row i maxima
    }
  }


row_vector FMStack::GetMaxima()

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	MaxMins	: The vector of row max & min values
	// Note			: The maxima for row i will be
	//			  <maxmins|i> in absolute units. The
	//			  real part is the row maximum and  
	//			  the imaginary part is the row minimum.
	// Note			: Some rows may be flagged so that they
	//			  are not plotted (via row increment) 
	//			  Their maxima are still determined!

  { return MaxMins; }


// -------------------------- Stack Plot Scaling -----------------------------


void FMStack::SetScale()

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	void    : Scaling parameters in FMSTK are set
	//				hwidth: Plotted points width (cm)
	//				HSF:	Horiz. scaling (cm/pt)
	//				VSF:	Vertical scaling (cm/absunits)

// The ith Point in Vector MaxMins, <MaxMins|1> Contains Maxima For Row i In A
// Stack Plot in ABSOLUTE UNITS (those input).  These, As Well As All Plotted
// Points, Must Be Converted Into Plotting Units, Namely To Centimeters (cm).
// The Conversion On the Horizontal Axis Is Simple Because We Just Use Points
// As Absolute Units and Get The Simple Conversion Factor HSF in cm/point.  The
// Total Plot Width (hwidth) Is Simply The Length Of A Row + Total Horizontal
// Offset Over All Rows (rows*hinc).  So the Horizontal Scaling is Obtained By

//	     VSF = (hwidth - rows*hinc)/cols = (hwidth-hdelta)/cols

// Conversion On the Vertical Axis Is Quite a Bit More Tricky Because We Have
// To Keep The Total Vertical Plot Span To Be Within the Specified Size:

//		VPlotSize = <mxmn|imax>*VSF + offsets(imax)
//                        - <mxmn|imin>*VSF + offsets(imin)

// Where VSF is The Vertical Scaling Factor (cm/pt), imax The Plotted Row That
// Has The Highest Point, and imin The Plotted Row With The Lowest Point.  We
// Dont Know Any Of These Three Values, Although We Only Want VSF.  But, Note
// That As VSF Changes The Rows With Absolute Maxima (imax & imin) May Change. 

// Consider A Plot Where the Lowest Row Has Noise Whereas The Third Row Has
// A Small Peak  Which Doesn't Go Down And Reach The First Row.  Then, The
// First Row Contains The Mimimum.  As The Vertical Scaling Increases, The Peak
// Will Grow Down Becomeing the Lowest Point So Then That Row Has The Mimimum.

// The End Result Is That We Have To Scale Iteratively for VSF......


  {
  hwidth = hsize-hdelta;			// Horizontal Row Width
  HSF = hwidth/(cols-1);			// Horiz. scaling cm/point
  row_vector mxmn = GetMaxima();		// For convenience
  int scale = 1;				// Iterative scaling flag
  int iter = 0;					// Interation count
  VSF = 1;					// Initilize scaling factor
  double max, min;				// Input data global maxima
  int imax, imin;				// Row indices of glob. maxima
  maxima(max, min, imax, imin);			// Get maxima for this VSF
  double span = max-min;			// The span for this VSF (abs)
  double last_span;				// To track spans vs. VSF
  if(FMdebug>1)					// Output vertical scaling
    ScaleInfo(iter,max,min,imax,imin); 		// info if debugging
  while(scale)					// Iterate for VSF
    {
    last_span = span;				// Previous iteration span
    iter++;					// Increment iteration count
    VSF = vsize - OffSets.getIm(imax)		// Guess at a new vertical
                + OffSets.getIm(imin);		// scaling factor (VSF) in
    VSF /= (mxmn.getRe(imax)-mxmn.getIm(imin));	// cm/absunits
    maxima(max, min, imax, imin); 		// Get maxima for this VSF
    span = max-min;				// Span for this scaling
    if(FMdebug>1)				// Output vertical scaling
      ScaleInfo(iter,max,min,imax,imin); 	// info if debugging
    if(fabs(span-vsize) < VRES) scale = 0;	// Stop if VSF correct (10 um!)
    else if(iter>100)				// Stop if too many iterations
      { 					// for a decent scaling factor
      scale = 0;
      FMSTKerror(60);
      }
    else if(fabs(span-last_span) < VRES)	// Stop if we aren't finding
      { 					// any differences in the span
      scale = 0; 				// compared to VSF iterations
      FMSTKerror(60);
      }
    }
  }

// ____________________________________________________________________________
// D                       FM STACK BASE PLANE FUNCTIONS
// ____________________________________________________________________________


void FMStack::SetBasePlane()

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	void    : Scaling parameters in FMSTK are set
	//				hwidth: Plotted points width (cm)
	//				HSF:	Horiz. scaling (cm/pt)
	//				VSF:	Vertical scaling (cm/absunits)
	// Note			: The "skew" herein is horizontal shift (cm)
	//			  difference between the front base & back base

  {
  double max, min, skew=0;
  int imax, imin;
  maxima(max,min,imax,imin);
  Bback = VSF*MaxMins.getRe(imax)
        - (hdelta - OffSets.getIm(imax));	// Back vert. coord (cm)
  Bfront = Bback + hdelta;			// Frnt vert.coord. (cm)
  if(vinc*hinc < 0) skew = hdelta;		// Adjust skewing (cm)
  BbackW  = hdelta - skew;			// BackWest (cm) 
  BbackE  = BbackW + hwidth;			// BackEast (cm) 
  BfrontW = skew;				// FrontWest (cm)
  BfrontE = BfrontW + hwidth;			// FrontEast (cm)
  }


void FMStack::BaseStart(ostream& ostr) const

        // Input	FMSTK   : FM stack plot parameters (this)
        // 		ostr	: Open output stream
        // Return	void    : The output stream is modified
	//			  to contain the start of the base
	//			  plane in a stack plot in FM MIF.
	// Note			: This plots two of the 4 lines that
	//			  make up the base plane borders. The
	//			  two plotted will be those NOT hidden
	//			  by the stack plot, and which two that
	//			  will be depends upon the plot skew.

  {
  double plB  = Bback;				// Base plane back coord (cm)
  double plF  = Bfront;				// Base plane front coord (cm)
  double plBW = BbackW;				// Base plane front left (cm)
  double plBE = BbackE;				// Base plane front right (cm)
  double plFW = BfrontW;			// Base plane back left (cm)
  double plFE = BfrontE;			// Base plane back right (cm)
  if(vinc >= 0)					// Draw 1 plane back (skew!)
    FM_Line(ostr,BPID,-1,2,plBW,plB,plBE,plB);
  else
    FM_Line(ostr,BPID,-1,2,plFW,plF,plFE,plF);
  if(hinc >= 0)					// Draw 1 plane side (skew!)
    FM_Line(ostr,BPID,-1,2,plFW,plF,plBW,plB);
  else
    FM_Line(ostr,BPID,-1,2,plFE,plF,plBE,plB);
  }
     

void FMStack::BaseEnd(ostream& ostr) const

        // Input	FMSTK   : FM stack plot parameters (this)
        // 		ostr	: Open output stream
        // Return	void    : The output stream is modified
	//			  to contain the finish of the base
	//			  plane in a stack plot in FM MIF.
	// Note			: This plots two of the 4 lines that
	//			  make up the base plane borders. The
	//			  two plotted will be those that ARE hidden
	//			  by the stack plot, and which two that
	//			  will be depends upon the plot skew.

  {
  double plB  = Bback;				// Base plane back coord (cm)
  double plF  = Bfront;				// Base plane front coord (cm)
  double plBW = BbackW;				// Base plane front left (cm)
  double plBE = BbackE;				// Base plane front right (cm)
  double plFW = BfrontW;			// Base plane back left (cm)
  double plFE = BfrontE;			// Base plane back right (cm)
  if(hinc >= 0)					// Draw 1 plane side (skew!)
    FM_Line(ostr,BPID,-1,2,plFE,plF,plBE,plB);
  else
    FM_Line(ostr,BPID,-1,2,plFW,plF,plBW,plB);
  if(vinc >= 0)					// Draw 1 plane top (skew!)
    FM_Line(ostr,BPID,-1,2,plFW,plF,plFE,plF);
  else
    FM_Line(ostr,BPID,-1,2,plBW,plB,plBE,plB);
  FM_Group(ostr, BPID);				// Group all 4 plane sides
  }


int FMStack::BasePlaneID() const

        // Input	FMSTK	: FM stack plot paramters (this)
        // Return	BPID	: The base plane ID number is returned

  { return BPID; }


void FMStack::BasePlaneID(int bid)

        // Input	FMSTK	: FM stack plot paramters (this)
	//		bid	: FM group ID number
        // Return	BPID	: The base plane ID number is set to bid

  { BPID = bid; }

 
// ____________________________________________________________________________
// E               FM STACK PLOT INITIALIZATON FUNCTIONS
// ____________________________________________________________________________


void FMStack::maxima(double& max, double& min, int& imax, int& imin)


        // Input        FMSTK   : FM stack parameters (this)
        //              max     : Maximum vertical value (cm)
        //              min     : Minimum vertical value (cm)
        //              imax    : Index of row containing maximum
        //              imin    : Index of row containing minimum
        // Return               : Void, imax and imin are modified
        // Note                 : Formulae for getting intensities are
        //                          I = Re<MaxMins|i>*VSF + Im<OffSets|i>
        //                          I = Im<MaxMins|i>*VSF + Im<OffSets|i>

   {
   int dim = MaxMins.elements();		// Number of maxima
   max = -1.e50;				// Very small maximum
   min = 1.e50;					// Very large miminum
   double pt;					// Scratch point
   for(int i=0; i<dim; i+=rowinc)		// Loop over maxima
     {
     pt = MaxMins.getRe(i)*VSF + OffSets.getIm(i);
     if(pt > max)
       {
       max = pt;
       imax = i;
       }
     pt = MaxMins.getIm(i)*VSF+ OffSets.getIm(i);
     if(pt < min)
       {
       min = pt;
       imin = i;
       }
     }
   }

void FMStack::CheckRI()

	// Input	FMSTK	: FM stack parameters (this)
 	// Output	void	: The current row increment
	//			  is checked to insure at least 
	//			  two rows are plotted

  {
  if(rowinc >= rows/2) rowinc = rows/10-1;	// We want at least 2 rows
  if(rowinc < 1)       rowinc = 1; 		// plotted. 
  }


void FMStack::StkBegin(const matrix& mx)

        // Input        FMSTK   : FM stack parameters (this)
        //              mx      : Input data matrix to plot

  {
  plotsize();				// Insure plot size O.K.
  SetArray(mx);				// Set data array to plot
  CheckRI();				// Insure row increment O.K.

//	   	         Determine The Last Plotted Row
//		     (Last Matrix Row If All Rows Plotted)

  int lastrc;
  if(rowinc > 1)
    {
    lastrc = 0;
    for( ; lastrc<rows; lastrc += rowinc) {}	// Get last row for plot
    rows -= rowinc;				// Adjust row dimension
    rows++;
    }
	
//        Determine Total Base Offsets From Row & Column Incremnts

  double ip1 = double(rows+1); 
  hdelta = fabs(hinc)*ip1;			// Determine total x offset
  vdelta = fabs(vinc)*ip1;			// Determine total y offset

//         Insure Stack Plot Will Stay With Choosen Plot Height

  if(vsize < vdelta)				// Check bad vertical offset
    {                                           // This occurs if summed row
    vdelta = 0.30*vsize;			// offset surpasses plot height
    vinc *= vdelta/(ip1*fabs(vinc));		// So, adjust the increment
    } 						// to span  30 % of total plot
 
//         Insure Stack Plot Will Stay With Choosen Plot Width
 
  if(hsize < hdelta)				// Check bad horizontal offset
    {                                           // This occurs if summed row
    hdelta = 0.30*hsize; 			// offset surpasses plot width
    hinc *= hdelta/(ip1*fabs(hinc));		// So, adjust the increment
    } 						// to span  30 % of total plot

//	     Insure The Specified Column Increment is Reasonable
//		  (This Is Only Active If Gridding Is On)

//  if(gridding)
//    {
//    if(cinc >= cols/2) cinc = cols/10-1;	
//    if(cinc <= 0)      cinc = 1;
//    }
//  else                 cinc = 1;

//	          Determine The Last Plotted Column
//	     (Last Matrix Column If All Columns Plotted)

//  if(cinc > 1)
//    {
//    lastrc = 0;
//    for(; lastrc<cols; lastrc+=cinc) {}	// Find last plot col
//    cols = lastrc+1;
//    }
  }


void FMStack::Initialize(const matrix& mx)

	// Input 	FMSTK   : FM stack plot parameters (this)
	// 		mx      : Data matrix
	// Return		: Void, all internal FMSTK initial values
	//			  that are dependent upon the matrix mx
	//			  are set and/or checked

  {
  StkBegin(mx);				// Set up plot parameters
  SetOffsets();				// Set up row offsets(cm)
  SetMaxima();				// Set up row maxima(cm)
  SetScale();				// Set up scaling factors
  SetBasePlane();			// Set up the base plane
  RB = RowBlk(2*cols-1);		// Set up the row block
  RPL = FMPL(2*cols-1);			// Set up the PolyLine
  RB.SetDebug(FMdebug);			// Set block debugging level
  Levels[0] = -1;			// Initialize top to above plot
  Levels[2] = Bfront;			// Initialize to base plane front
  Levels[3] = vsize+1;			// Initialize to below plot
  if(FMdebug) StkInfo();		// If debugging, output info
  }


void FMStack::StartOutput(const string& filename, ostream& ostr)

	// Input 	FMSTK   : FM stack plot parameters (this)
	//		filename: Output filename
	// Return	ostr	: An output stream is begun with
	//			  external name "filename".  The
	//			  basics of a FrameMaker MIF file
	//			  are written into the output stream

  {
  if(!ostr.good())			// Insure ostream OK
    {
    FMSTKerror(10,1);			//      Bad ostream
    FMSTKfatality(11);			//      Can't write to ostr
    }                
  FM_Begin(ostr);			// FrameMaker begin comment
  FM_AFrames_Begin(ostr);		// Anchored frame begin
  FM_AFrame_Set(ostr,hsize,vsize,FRID);	// Set anchored frame
  }


// ____________________________________________________________________________
// F                   FM STACK PLOT DEBUGGING FUNCTIONS
// ____________________________________________________________________________


void FMStack::StkInfo() const

        // Input	FMSTK   : FM stack plot parameters (this)
	// Output	void	: Some information output to screen

  {
  int rl = rows-1; 
  cout << "\n\n\tBasic Stack Plot Parameters\n";
  cout << "\n\t\tHorizontal Plot Size:       " << hsize << " cm";
  cout << "\n\t\tVertical Plot Size:         " << vsize << " cm";
  cout << "\n\t\tHorizontal Axis Start:      " << haxmin;
  cout << "\n\t\tHorizontal Axis Finish:     " << haxmax;
  cout << "\n\t\tVertical Axis Start:        " << vaxmin;
  cout << "\n\t\tVertical Axis Finish:       " << vaxmax;
  cout << "\n\t\tVertical Plot Resolution:   " << VRES << " cm";
  cout << "\n\t\tHorizontal Plot Resolution: " << HRES << " cm";
  cout << "\n\t\tData Useage In Plot:        " << int(duse);
  cout << "\n\t\tMaximum PolyLine Size:      " << PLmax << " Points";
  cout << "\n\t\tDebugging Level:            " << FMdebug;
  cout << "\n\t\tMatrix Rows To Plot:        " << rows;
  cout << "\n\t\tMatrix Columns To Plot:     " << cols;
  cout << "\n\t\tTotal Horizontal Offset:    " << hdelta << " cm";
  cout << "\n\t\tTotal Vertical Offset:      " << vdelta << " cm";
  cout << "\n\t\tFirst Row Offset Over:      " << OffSets.getRe(0) << " cm";
  cout << "\n\t\tFirst Row Offset Up:        " << OffSets.getIm(0) << " cm";
  cout << "\n\t\tLast Row Offset Over:       " << OffSets.getRe(rl) << " cm";
  cout << "\n\t\tLast Row Offset Up:         " << OffSets.getIm(rl) << " cm";
  cout << "\n\t\tVertical Scaling Factor:    " << VSF << " cm/absunits";
  cout << "\n\t\tHorizontal Scaling Factor:  " << HSF << " cm/pt";
  cout << "\n\t\tFront Left Border:         (" << BfrontW
                                       << ", " << Bfront << ") cm";
  cout << "\n\t\tFront Right Border:        (" << BfrontE
                                       << ", " << Bfront << ") cm";
  cout << "\n\t\tBack Left Border:          (" << BbackW
                                       << ", " << Bback<< ") cm";
  cout << "\n\t\tBack Right Border:         (" << BbackE
                                       << ", " << Bback<< ") cm";
  cout << "\n";  
  }


void FMStack::ScaleInfo(int iter, double max, double min, int imax, int imin) 

        // Input 	FMSTK	: FrameMaker stack plot parameters
        // Return	void    : Vertical scaling information output

  {
  cout << "\n\t\tVertical Scaling Iteration " << iter
       << ": " <<  VSF << " cm/absunits"; 
  cout << "\n\t\t\tMaximum Is:          " << max << ", Row " << imax; 
  cout << "\n\t\t\tMinimum Is:          " << min << ", Row " << imin; 
  cout << "\n\t\t\tTotal Scaled Span:   " << max-min << " cm";
  }


// ____________________________________________________________________________
// G                      FM STACK ROW/BLOCK FUNCTIONS
// ____________________________________________________________________________


void FMStack::FillBlock(int row)

        // Input	FMSTK	: FM stack plot paramters (this)
	//		row	: Matrix row index
        // Return	void 	: The row block RBL is filled from data array
	//			  in FM coordinates.

// Note that the integer array "Blevs" contains four values.  These are
// FrameMaker vertical plot coordinates. They were originally set in the
// function StkRowPlot and one appropriate for the input row will be
// set in this routine.  These levels are defined to be the following:

//		0. Above entire plot		2. Base plane front
//		1. Row baseline			3. Below entire plot

// Based on the these four values, each point is assigned an "area" flag.
// These area flags are as follows:

// 			0: Point is above row baseline
//		 	1: Point below baseline, above bottom of base plane
//			2: Below bottom of base plane

// This function takes the row "row" (of the stack plot data matrix) and
// converts its points from absoulute to FM coordinates.  These points are
// placed into the RowBlock "RB" and to each point is assigned an
// "area" flag which is stored in the integer array "Bareas". Additional
// points may be interdispersed into block (& areas) whenever the row
// points change the area in which they reside. This is in order to
// facilitate implmentation of a HLA when the block is plotted.

// Whenever the row crosses from one area into the next points are added at
// the intersection with the specified row levels.

  {
  RB.SetRow(row);				// Set the block row index
  bhoff = OffSets.getRe(row);			// Row's horiz. offset (cm)
  bvoff = Bback + vdelta - OffSets.getIm(row);	// Row vertical offset (cm)
  RB.SetOffsets(bhoff, bvoff);			// Set block row offsets
  RB.SetScaling(HSF,VSF);			// Set block scaling
  Levels[1] = bvoff;				// Store row baseline
  RB.SetLevels(Levels);				// Set block levels
  RB.Fill(cols, MX);				// Set block values 
  if(FMdebug>1) cout << RB; 			// If debug: output block 
  if(FMdebug>2)					// If more debugging, plot
    {						// the block
    string FMname="R"+string(Gdec(row))+".mif";	// This will be output file
    cout << "\n\t\tPlotting Block To File: "	// Say what we are doing
         << FMname;
    RB.Plot(FMname); 				// Plot Block to FM
    cout << "\n\t\tBlock MIF File " << FMname
         << " Complete.";
    cout.flush();
    }
  }


void FMStack::PlotRow(int row, ostream& ostr)

        // Input	FMP	: FrameMaker plot parameters (this)
	//		row	: Matrix row index
	//		ostr	: Output stream
        // Return	void 	: Row "row" of stack plot is sent into the
	//			  output stream in FrameMaker MIF format.

// Note that the integer array "limits" contains four values.  These are
// FrameMaker vertical plot coordinates. They are set in the functions
// FillBlock and StkRowPlot and defined to be the following:

//		0. Above entire plot		2. Base plane front
//		1. Row baseline			3. Below entire plot

// Based on the these for values, each point is assigned an "area" flag.
// These area flags are as follows:

// 			0: Point is above row baseline
//		 	1: Point below baseline, above bottom of base plane
//			2: Below bottom of base plane

// First we take the row and convert it into plotting coordinates.  This will
// be the vector "block".  Additional points may be interdispersed in the
// block.  These occurs whenever any two point cross one of the designated
// levels.  The added point(s) will lie on the line connecting the two points
// and reside on one of the levels.  To each of the points in "block" and area
// flag will be set so we know which of the desinated levels it is above.

// After the row is converted into a block (PLB) it is output as FrameMaker
// PolyLines in sections.  Connected points in area 0 and 2 are written out
// with a solid black pen with a solid white fill.  This partially creates
// the hidden lines because these hide whatever is behind them.  Connected
// point is area 1 are run through another HLA. These sections will be plotted
// with no fill and will either have a clear or black pen.  They will be clear
// if they lie behind any successive rows (as yet not plotted) or they will be
// black if nothing is in front of them. 

  {
  FillBlock(row);				// Set Levels, Block, & Area
  double pw = 1.0;				// Set pen width 1 pt
  int pf[3] = { 15, 15, 15 };			// Pen filling (15 = none)   
  int pc = 0;					// Set pen color (black)
  if(HLF)					// If HLA then lines areas
    {						// 0 (above row base) and
    pf[0] = 7;					// 2 (below base plane)
    pf[2] = 7;					// should be solid filled
    }
  complex pt = RB.get(0);			// Working data point
  RPL.AddPt(pt);				// Set the 1st point
  int area, parea = RB.Bareas[0];		// Set 1st pt area flag
  RPL.Set(pw,0,pf[parea],RB.Bids[parea]);	// Set PL pen,area,fill
  RPL.SetReduce(1);				// Set data reduction
  for(int j=1; j<RB.Bpts; j++)			// Loop columns of this block
    {
    area = RB.Bareas[j];			// Get "area" of point
    if(parea != area)				// Point doesn't change area
      {						// from the previous one.
      if(parea > area) RPL.AddPt(RB.get(j));	// Include next pt if going up
      if(HLF>0 && parea==1) HLA(ostr); 		// Use HLA for area 1 pts
      else RPL.WriteMIF(ostr); 			// No HLA for areaa 0 & 2
      RPL.Zero();				// Zero the PL point count
      if(parea < area) RPL.AddPt(pt);		// Set 1st pt if going down
      parea = area;				// Set new area
      RPL.Set(pw,pc,pf[parea],RB.Bids[parea]);	// Set PL pen,area,fill
      }
    pt = RB.get(j);				// This is the next pt
    RPL.AddPt(pt);				// Set the next point
    }
  if(RPL.GetPoints()>1) RPL.WriteMIF(ostr);	// Output any end PolyLine
  RPL.Zero();					// Zero the PL point count
  for(int a=0; a<3; a++)			// Group PolyLines together 
    FM_Group(ostr,RB.Bids[a],RB.Bids[3]);	// each with new ID Bids[3]
  RB.Zero();
  }


// ____________________________________________________________________________
// H                   FM STACK PLOT CREATION FUNCTIONS
// ____________________________________________________________________________


void FMStack::FM_stack(const string& fn, const matrix &mx) { Plot(fn, mx); }
void FMStack::Plot(const string& filename, const matrix &mx)

	// Input 	FMSTK   : FM stack plot parameters (this)
	//		filename: Output filename
	// 		mx      : Data matrix
	// Return		: Void, file filename created and
	//			  will contain a stack plot of the 
	//		 	  array mx in FrameMaker MIF format

  {
  Initialize(mx);			// Initialize FMSTK for mx
  ofstream ostr(filename.c_str());	// Open output stream
  StartOutput(filename, ostr);		// Begin MIF output file
  BaseStart(ostr);			// Draw base plane start
  int row, nr=rows-1;			// Row index, # of rows
  for(row=nr; row>=0; row-=rowinc)	// Loop entire matrix and
    PlotRow(row, ostr);			// plot rows back to front
  BaseEnd(ostr);			// Draw base plane start
  for(row=nr; row>=0; row-=rowinc)	// Group each row
    FM_Group(ostr, 5*row+3+RB.BaseID, FRID);
  FM_Group(ostr, FRID);			// Group all rows 
  FM_AFrame_End(ostr);			// Frame end
  FM_AFrames_End(ostr);			// Anchored Frame end
  FM_ParaText_End(ostr);		// Create TextFlow with Frame
  FM_End(ostr);				// FrameMaker end comment
  }


void FM_stack(const string& filename, const matrix &mx, FMStack& FMSTK)

	// Input	filename: Output filename
	// 		mx      : Data matrix
	//		FMSTK   : FM stack plot parameters
	// Return		: Void, file out is modified

  { FMSTK.Plot(filename, mx); }


void FM_stack(const string& filename, const matrix &mx, double xinc, 
            double yinc, int RI, double xsize, double ysize, int grid, int CI)

	// Input	filename  : Output filename
	// 		mx        : Data matrix
	// 		xinc      : Delta x in cm
	// 		yinc      : Delta y in cm
	//		RI   	  : Row increment
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	//		grid  	  : Flag for gridding (not currently active!)
	//		CI  	  : Column increment (only with gridding!)
	// Return		  : Void, file out is modified
 
  { 
  FMStack FMSTK;			// FM stack plot parmeters
  FMSTK.HInc(xinc);			// Set horizontal offset increment
  FMSTK.VInc(yinc);			// Set vertical offset increment
  FMSTK.HPlotSize(xsize);		// Set horizontal plot size
  FMSTK.VPlotSize(ysize);		// Set vertical plot size
  FMSTK.RInc(RI);			// Set the row incrementation
  FMSTK.Gridding(grid);			// Set the gridding flag
  FMSTK.CInc(CI);			// Set column increment
  FMSTK.FM_stack(filename, mx);		// Use function overload
  }
 

// ____________________________________________________________________________
// I                  FM STACK HIDDEN LINE FUNCTIONS
// ____________________________________________________________________________


complex intersect(const complex& z1, const complex& z2,
                  const complex& z3, const complex& z4, int& crx)

 	// Input	z1	: First point, line 1
 	// 		z2	: Second point, line 1
	// 		z3	: First point, line 2
	// 		z4	: Second point, line 2
	//		crx	: Flag for line crossing
	// Output	z	: Point at the intersection of
	//			  the line connecting z1 & z2 with
	//			  the line connecting z3 and z4
	// Note			: There may be no intersection!
	// 			  In that case crx=FALSE
	// Note			: Used exclusively by function xings

  {
  complex z = 0;
  double x1, y1, x2, y2;
  double u1, v1, u2, v2;
  double xint, yint;
  double mxy, muv;
  x1 = Re(z1);
  y1 = Im(z1);
  x2 = Re(z2);
  y2 = Im(z2);
  u1 = Re(z3);
  v1 = Im(z3);
  u2 = Re(z4);
  v2 = Im(z4);
  mxy = (y1-y2)/(x1-x2);
  muv = (v1-v2)/(u1-u2);
  if(mxy == muv) crx = 0;
  else
    {
    xint = ((y1-mxy*x1)-(v1-muv*u1))/(muv-mxy);
    yint = y1 + mxy*(xint-x1);
    z = complex(xint, yint);
    crx = 1;
    }
  return z;
  }


complex FMStack::xings(complex& z1, complex& z2, const complex &z3, const complex &z4,
	                  int row, double xoff, double yoff, int k, int& cross)

 	// Input	z1	: First point, on line 1
 	// 		z2	: Second point, on line 1
	// 		z3	: First point, line 2
	// 		z4	: Second point, line 2
	// 		xoff	: Horizontal offset for line 1 (cm)
	// 		yoff	: Vertical offset for line 1 (cm)
	//		FMxscale: x-dimension scale (cm/point)
	//		FMyscale: y-dimension scale (cm/Re(mx(i,j)))
	//		mx	: Matrix of data values
	//		row	: Row number of line 1
	//              k       : Point being checked
	//		cross   : Flag for line crossing
	// Output	z	: Point at the intersection of
	//			  the line connecting z1 & z2 with
	//			  the line connecting z3 and z4
	// Note			: There may be no intersection !
	// 			  in which case cross=FALSE
	// Note			: Used exclusively by set_hidden

  {
  complex z = intersect(z3,z4,z1,z2,cross);	// Try for intersection
  if(cross)					// If intersection found
    {						// then make sure the
    if(Re(z) < Re(z3))      cross=0;		// crossing point is
    else if(Re(z) > Re(z4)) cross=0;		// reasonable
    }
  if(!cross)					// If bad crossing then
    {						// look 1 pt to the left
    z2 = z1;
    double x1 = xoff + (k-2)*HSF;
    double y1 = yoff - Re(MX.get(row,k-2))*VSF;
    z1 = complex(x1, y1);
    z = intersect(z3, z4, z1, z2, cross);	// Try again for intersection
    if(Re(z)<Re(z3) || Re(z)>Re(z4))	// Check if this one is
      {  					// O.K. or not
      cross = 0;				// Not O.K., set bad crossing
      return complex0;				// flag and just zero z
      }
    }
  return z;
  }


void FMStack::set_hidden(int &pts, int row, int *hidden, row_vector& crossing)

 	// Input	FMSTK	: FM stack parameters (this)
        //              pts     : Number of points in PolyLine
	//		row	: Current row
	//		hidden  : Flags whether points are hidden
	//		crossing: Crossing points between visible & hidden
	// Output	none    : This function provided two services -
	//			  1.) Each point in the PolyLine is flagged
	//			      whether it should be visible or hidden
	//			  2.) Two adjancent points which have a
	//			      transition between visible and hidden
	//			      are extrapolated to find this crossing

// The points in the input PolyLine PL all reside in "area 1", the part of
// the stack plot which is below the current rows height yet above the bottom
// of the base plane.  That means that the points will be "hidden" (drawn with
// a white pen) if they happen to be behind any of the rows that are in front
// of them.  This reasoning is NOT applicable to any points above "area 1"
// because PolyLines in front are drawn with a white fill so they will auto-
// matically hide those behind them.

  {
complex zint;
int nxtrow;
complex zcross;
complex z1, z2;
  int lev, search=0;
  int *levels;
  levels = new int[pts];
//  int levels[pts];
  double x,y;
  double x1,y1,x2,y2;
  complex z;
  double yint;
  double xoff,yoff;
  for(int j=0; j<pts; j++)		// Loop through entire PL
    {					// determine whether hidden or not
    nxtrow = row;				// Initialize row count
    hidden[j] = 0;				// Default is visible point
    search = 1;					// Set for crossing search
    x = RPL.getRe(j);				// Point x coord (cm)
    y = RPL.getIm(j);				// point y coord (cm)
    while(search)				// Search rows in front
      {
      nxtrow -= rowinc;				// Next plotted row
      xoff = OffSets.getRe(nxtrow);		// Next row x offset (cm)
      yoff = Bback+vdelta-Im(OffSets(nxtrow));	// Next row y offset(cm)
      if(nxtrow<0)              search = 0;	// Below border -> Visible
      else if(x<xoff)           search = 0; // Left of border -> Visible
      else if (x>(xoff+hwidth)) search = 0;	// Right of border -> Visible
      else					// Within border->Invisible?
        {
        x2 = xoff;				// Find line segment of row
        int k = 0;				// beneath current section
        while(x2 < x)
          {
          k++;
          x2 = xoff + k*HSF;
          }
        x1 = xoff + (k-1)*HSF;
        y1 = yoff - MX.getRe(nxtrow,k-1)*VSF;
        y2 = yoff - MX.getRe(nxtrow,k)*VSF;
        yint = ((y1-y2)/(x1-x2))*(x-x1) + y1;
        if(yint < y)
          {
          hidden[j] = 1;			// Point should be invisible
          levels[j] = nxtrow;			// Keep level it is under
          search = 0;
          int cross = 0;
          if((j > 0) && (hidden[j-1] == 0))
{
z1 = complex(x1, y1);
z2 = complex(x2, y2);
           zcross = xings(z1, z2, RPL.get(j-1), RPL.get(j),
                              nxtrow, xoff, yoff, k, cross);
           crossing.put(zcross, j-1);
if(zcross == complex0)
  cout << "\n\tCrossing Point Is Zero! " << j-1; 
}
          }
        }
      }						// Go to next lower level

    if((j>0) && (hidden[j-1]>hidden[j])) 	// For visible pt, if last pt
      {						// hidden get crossover point
lev = levels[j-1];
      nxtrow = levels[j-1];
      xoff = OffSets.getRe(nxtrow);		// Next row x offset (cm)
      yoff = Bback+vdelta-OffSets.getIm(nxtrow);// Next row y offset(cm)
      x2 = xoff;				// Find line segment of row
      int k = 0;				// beneath current section
      while(x2 < x)
        {
        k++;
        x2 = xoff + k*HSF;
        }
      x1 = xoff + (k-1)*HSF;
      y1 = yoff - MX.getRe(nxtrow,k-1)*VSF;
      y2 = yoff - MX.getRe(nxtrow,k)*VSF;
      int cross = 0;
z1 = complex(x1, y1);
z2 = complex(x2, y2);
      z = xings(z1,z2,RPL.get(j-1),RPL.get(j),
                     nxtrow,xoff,yoff,k,cross);
      nxtrow -= rowinc;				// Check 1 level lower
      xoff = OffSets.getRe(nxtrow);		// Set row x offset (cm)
      yoff = Bback+vdelta-OffSets.getIm(nxtrow);// Set row y offset(cm)
      x1 = xoff + (k-1)*HSF;			// Previous x
      y1 = yoff - MX.getRe(nxtrow,k-1)*VSF;	// Previous y
      x2 = xoff + k*HSF;			// This point x
      y2 = yoff - MX.getRe(nxtrow,k)*VSF;	// This point y
z1 = complex(x1, y1);
z2 = complex(x2, y2);
      zint = xings(z1,z2,RPL.get(j-1),RPL.get(j),
                       nxtrow,xoff,yoff,k,cross);
      if((cross == 0) || (Re(zint) < RPL.getRe(j-1))
                      || (Re(zint) > RPL.getRe(j)))
        zint = z;
      if(zIm(z)>zIm(zint) && cross) z = zint;
      crossing.put(z, j-1);
      }
    }						// Check next PL point
  delete [] levels;
  }


void FMStack::HLA(ostream& ostr)

 	// Input	FMSTK   : FM stack plot parameters
	//		ostr	: Output stream
	// Output	none    : Outputs current PolyLine to ostr
	//			  in FrameMaker MIF format using a HLA
	//			  (Hidden Line Agorithm) to determine
	//			  whether each point is visible.

// The points in the input PolyLine PL all reside in "area 1", the part of
// the stack plot which is below the current rows height yet above the bottom
// of the base plane.  That means that the points will be "hidden" (drawn with
// a white pen) if they happen to be behind any of the rows that are in front
// of them.  This reasoning is NOT applicable to any points either above area 1
// (in area 0) or below area 1 (in area 2) because PolyLines in front are drawn
// with a white fill so they will auto matically hide those behind them.

  {
  int *hidden;
  hidden = new int[RB.Bpts];
//  int hidden[RB.Bpts];				// Flags if points hidden
  row_vector crossing(RB.Bpts,0);		// Vector of crossing points
  int npts = RPL.GetPoints();			// Point in input RPL
  set_hidden(npts,RB.Brow,hidden,crossing); 	// Set hidden flags, crossings
  int hide, prevhide=hidden[0];			// PolyLine point hide flags
  FMPL PL1(RB.Bpts);				// Vector for HLA PolyLine
  PL1.AddPt(RPL.get(0));			// Set initial PL point
  double pw = 1.0;				// Pen width (1 pt)
  int pc = 0;					// Pen color (0 = black)
  int pf = 15;					// Pen fill (15 = solid white)
  int ID = 5*RB.Brow + 23;			// PolyLine ID for row
  PL1.Set(pw,pc,pf,ID);				// Set PL pen,area,fill
  complex pt;					// For PL points
  for(int jj=1; jj<npts; jj++)			// Loop through entire PL
    {
    hide = hidden[jj];				// Get the "hidden" flag
    if(hide != prevhide)			// Draw if hidden to visible
      {						// or visible to hidden
      PL1.AddPt(crossing.get(jj-1));		// Last pt ist crossing pt
      if(prevhide) pc = 15; 			// Switch pen color either
      else         pc = 0;			// blk->whi or whi->blk
      pf = 15;					// Switch PL fill
      PL1.Set(pw,pc,pf,ID);			// Set PL pen,area,fill
      PL1.WriteMIF(ostr);			// Output any PolyLine
      PL1.Zero();				// Zero the PL point count
      PL1.AddPt(crossing.get(jj-1));		// First pt is crossing pt
      }
    pt = RPL.get(jj);				// This is the next pt
    PL1.AddPt(pt);				// Now add in the pt
    prevhide = hide;				// Update hide flag
    }
  if(prevhide) pc = 15; 			// Switch pen color either
  else         pc = 0;				// blk->whi or whi->blk
  pf = 15;					// Switch PL fill
  PL1.Set(pw,pc,pf,ID);				// Set PL pen,area,fill
  if(PL1.GetPoints()>1) PL1.WriteMIF(ostr);	// Output any end PolyLine
  PL1.Zero();					// Zero the PL point count
  PL1.AddPt(pt);				// Reset the 1st point
  delete [] hidden;
  }

// ____________________________________________________________________________
// J                  FM STACK PARAMETERS STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


ostream& FMStack::print(ostream& ostr) const

        // Input                FMSTK	: FM stack parameters (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                the FrameMaker stack parameters

  {
  FMPar::print(ostr); 
  ostr << "\n\tNumber of Rows To Plot:       " << rows;
  ostr << "\n\tNumber of Columns To Plot:    " << cols;
  ostr << "\n\tTotal Horizontal Offset:      " << hdelta << " cm";
  ostr << "\n\tTotal Vertical Offset:        " << vdelta << " cm";
  ostr << "\n";
  return ostr;
  }
 

ostream& operator<< (ostream& ostr, const FMStack& FSTK)

        // Input                ostr    : An output stream
        //                      FSTK	: FM stack parameters (this)
        // Output               ostr    : The output stream modified by
        //                                the FM stack parameters

  {                                                                 
  FSTK.print(ostr);
  return ostr;
  }

// ____________________________________________________________________________
// K                         FM STACK LEFT OVERS
// ____________________________________________________________________________


void grids(ostream &out, const matrix &mx, int RI, int CI,
	row_vector &offsets, double FMxscale, double FMyscale,
		double deltay, double plB, int row, int jdim, int ID)
	// Note			: AUXILIARY FUNCTION FOR FM_Stack
  {
  int l;
  double x, y;
  FMPL PL(RI+1);
  for(int j=0; j<jdim; j+=CI)
    {
    for(int k=0; k<=RI; k++) 
      {
      l = row + RI - k;
      x = offsets.getRe(l) + j*FMxscale; 
      y = plB + deltay - offsets.getIm(l) - mx.getRe(l,j)*FMyscale; 
      PL.AddPt(complex(x,y));
      }
//    PL1.Set(pw,pc,pf,ID);			// Set PL pen,area,fill
//    PL.WriteMIF(ostr);			// Output any PolyLine
FM_PolyLine(out, PL, ID, 15, RI+1);
    PL.Zero();				// Zero the PL point count
    }
  }


/*************************************************************************
**								 	**
**  HLA Details								**
**						 			**
**  GAMMA's FrameMaker stack function automatically uses a Hidden Line	**
**  Algorithm (HLA).  Users can turn off hidden lines with FM so there	**
**  is not loss in flexibility this way.				**
**								 	**
**  There are two means by which hidden lines are accomplished. First	**
**  is to use a "white" pen in plotting hidden points.  Second is to	**
**  use a PolyLine fill of white to hide things behind it.  These are	**
**  fundamentally different and used accordingly.			**
**								 	**
**  White Pen: The points of each row are checked to see if they lie	**
**  1.) above the row "baseline", 2.) below the baseline but in the	**
**  "baseplane", or 3.) below the baseplane.  For cases 1. and 3. the	**
**  points are assumed visible and drawn with a black pen.  For case 2.	**
**  the points are checked to see if they are behine the ensuing row.	**
**  If they are they are drawn with a white (hidden) pen and if not	**
**  they are drawn with a black pen.					**
**								 	**
**  PolyLine Fill: The second way GAMMA makes hidden lines if FM is to	**
**  fill the PolyLines with white so they "hide" whatever is behind	**
**  them.  To use this one must plot the rows from back to front so	**
**  that the ones in front will "hide" the ones in back.  An additional	**
**  trick is to teach the line how to do the fill in sections.  For	**
**  that the row is broken up into sections so that the filling is	**
**  done properly.							**
**								 	**
*************************************************************************/

#endif 						// FrameMakerS.cc
