/* GgnuplotC.cc *************************************************-*-c++-*-
**									**
**	                           G A M M A 				**
**								 	**
**	Gnuplot Controls 			Implementation   	**
**				        	 		 	**
**      Copyright (c) 2002						**
**      National High Magnetic Field Laboratory 	                **
**      1800 E. Paul Dirac Drive                        	        **
**      Tallahassee Florida, 32306-4005                       		**
**									**
**      $Header: $
**	                         		         	 	**
*************************************************************************/

/*************************************************************************
**								 	**
**  Description						 		**
**								 	**
*************************************************************************/
 
#ifndef   GgnuplotC_cc_			// Is file already included?
#  define GgnuplotC_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#if defined(_MSC_VER)			// If we areusing MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

#include <GamIO/GgnuplotC.h>		// Include gnuplot header
#include <GamIO/Ggnuplot.h>
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <iostream>			// Include libstdc++ file streams
#include <string>			// Include libstdc++ ANSI strings

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ vectors
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output
using std::ios;				// Using libstdc++ standard input
using std::ofstream;			// Using libstdc++ output file streams
using std::endl;			// Using libstdc++ line end

// ____________________________________________________________________________
// i                 GAMMA Gnuplot Controls Error Handling
// ____________________________________________________________________________

/*         Input                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void GPControls::GPCerror(int eidx, int noret)
  {
  string hdr("Gnuplot Controls");
  switch(eidx)
    {
    case 10: GAMMAerror(hdr,"Cannot Find Gnuplot Executable",  noret);  //(10)
      break;
    case 50: GAMMAerror(hdr,"Cannot Open New Load File",  noret);  //(10)
      break;
    default: GAMMAerror(hdr,eidx, noret); break;// Usually Unknown Error  (-1)
      break;
    }
  }

void GPControls::GPCerror(int eidx, const string& pname, int noret)
  {
  string hdr("GAMMA Gnuplot Routine");
  string msg;
  switch(eidx)
    {
    case 10: msg = string("Executable Set To " + pname);                // (10)
             GAMMAerror(hdr, msg, noret); break;
      break;		
    case 50: msg = string("Cannot Open New Load File " + pname);	// (50)
             GAMMAerror(hdr, msg, noret); break;
      break;
    case 51: msg = string("Cannot Open New Data File " + pname);	// (51)
             GAMMAerror(hdr, msg, noret); break;
      break;
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }
     
volatile void GPControls::GPCfatal(int eidx)
  { GPCerror(eidx, 1); if(eidx) GPCerror(0); GAMMAfatal(); }

volatile void GPControls::GPCfatal(int eidx, const string& pname)              
  { GPCerror(eidx, pname, 1); if(eidx) GPCerror(0); GAMMAfatal(); }

// ____________________________________________________________________________ 
// i                 GAMMA Gnuplot Controls Default Settings
// ____________________________________________________________________________ 
  
void GPControls::defaults()
  {
  datafile = string("gnuplot.asc");		// Default data file name
  loadfile = string("gnuplot.gnu");		// Default load file name
  delim    = string("  ");			// Default data col delimiter

  string none("");				// Empty string
  term      = none;				// Terminal is default
  basedir   = none;				// Base directory
  basename  = none;				// Base filename
  cmd       = none;				// System command
  ylabel    = none;				// Label for yaxis
  xlabel    = none;				// Label for xaxis
  zlabel    = none;
  plottitle = none;				// Label for plot

  nokey      = false;				// Flag for key 
  parametric = false;				// Flag for parametric 
  border = 1;					// Border flag (yes/no)
  join   = 1;					// Line draw flag (yes/no)
  xtics  = 1;					// X-axis tic flag (yes/no)
  ytics  = 1;					// Y-axis tic flag (yes/no)
  xaxis  = 1;					// Flag for x-axis draw (yes/no)
  yaxis  = 1;					// Flag for y-axis draw (yes/no)
 
  compress = false;				// No data compression
  xsize  = 0;					// X-axis scaling factor
  ysize  = 0;					// Y-axis scaling factor
  }

void GPControls::copy(const GPControls& GPC)
  {
//  *Lfp      = *(GPC.Lfp);			// Load file output stream
//  *Dfp      = *(GPC.Dfp);			// Data file output stream
  loadfile  = GPC.loadfile;			// Copy loadfile name
  datafile  = GPC.datafile;			// Copy datafile name
  datafiles = GPC.datafiles;			// Copy additonal names

  plottitle = GPC.plottitle;			// Label for plot
  xlabel    = GPC.xlabel;			// Label for xaxis
  ylabel    = GPC.ylabel;			// Label for yaxis
  zlabel    = GPC.zlabel;			// Label for zaxis
  delim     = GPC.delim;			// Copy column delimiter

  term      = GPC.term;				// Terminal is default
  basedir   = GPC.basedir;			// Base directory
  basename  = GPC.basename;			// Base filename
  cmd       = GPC.cmd;				// System command

  nokey      = GPC.nokey;			// Flag for nokey (yes/no)
  border     = GPC.border;			// Border flag (yes/no)
  parametric = GPC.parametric;			// Parametric flag (yes/no)
  join       = GPC.join;			// Line draw flag (yes/no)
  xtics      = GPC.xtics;			// X-axis tic flag (yes/no)
  ytics      = GPC.ytics;			// Y-axis tic flag (yes/no)
  xaxis      = GPC.xaxis;			// Flag for x-axis draw (yes/no)
  yaxis      = GPC.yaxis;			// Flag for y-axis draw (yes/no)
 
  compress = GPC.compress;			// Copy data compression flag
  xsize  = GPC.xsize;				// X-axis scaling factor
  ysize  = GPC.ysize;				// Y-axis scaling factor
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A           Gnuplot Controls Contruction, Assigment, Destruction
// ____________________________________________________________________________



// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

GPControls::GPControls()                      { defaults(); }
GPControls::GPControls(const GPControls& GPC) { copy(GPC);  }

GPControls::~GPControls() {}
GPControls& GPControls::operator= (const GPControls& GPC)
  {
  if(this == &GPC) return *this;
  copy(GPC);
  return *this;
  }


// ____________________________________________________________________________
//                        Gnuplot Load File Basics
// ____________________________________________________________________________

bool GPControls::NewLoadFile(bool warn)
  {
  if(Lfp.is_open()) Lfp.close();		// Close existing file
  Lfp.open(loadfile.c_str(), ios::out);		// Open loadfile file
  if(!Lfp.is_open())				// Insure file is OK
    {
    if(warn) GPCerror(50, loadfile, 1);		//   Problems with file
    return false;
    }
  return true;
  }

void GPControls::CloseLoadFile() { Lfp.close(); }

void GPControls::RunLoadFile()
  {
RunGnuplot(loadfile);
// sosi
  }

// ____________________________________________________________________________
//                        Gnuplot Data File Basics
// ____________________________________________________________________________

bool GPControls::NewDataFile(int idx, bool warn)
  {
  if(Dfp.is_open()) Dfp.close();		// Close existing file
  string fname = datafile;
  if(idx !=-1) fname = datafiles[idx];
// sosi check if above index doesnot exist
  Dfp.open(fname.c_str(), ios::out);		// Open loadfile file
  if(!Dfp.is_open())				// Insure file is OK
    {
    if(warn) GPCerror(51, fname, 1);		//   Problems with file
    return false;
    }
  return true;
  }

void GPControls::CloseDataFile() { Dfp.close(); }

// ____________________________________________________________________________
//            Gnuplot Controls Gnuplot Loadfile Output Commands
// ____________________________________________________________________________

void GPControls::SetLoadFile(const string& fln)  { loadfile   = fln; }
void GPControls::SetDataFile(const string& fln)  { datafile   = fln; }
void GPControls::SetDataFiles(const vector<string>& dfs) { datafiles = dfs; }
void GPControls::SetParametric(bool par)  { parametric = par?false:true; }
void GPControls::SetKey(       bool key)  { nokey      = key?false:true; }
void GPControls::SetBorder(    bool bord) { border     = bord; }
void GPControls::SetXTics(     bool xt)   { xtics      = xt; }
void GPControls::SetYTics(     bool yt)   { ytics      = yt; }
void GPControls::SetDelimiter( const string& d) {delim = d; }

// ____________________________________________________________________________
//            Gnuplot Controls Gnuplot Loadfile Output Commands
// ____________________________________________________________________________

void GPControls::WriteBorder(ofstream& ostr)
  {
  if(border)  ostr << "set border"   << endl;	// Put in a border
  else        ostr << "set noborder" << endl;	// No border
  }

void GPControls::WriteParametric(ofstream& ostr)
  {
  if(parametric) ostr << "set parametric"   << endl;
  else           ostr << "set noparametric" << endl;
  }

void GPControls::WriteKey(ofstream& ostr)
  {
  if(nokey) ostr << "set nokey" << endl;	// Don't put in any legend
  else      ostr << "set key"   << endl;	// Put in a legend
  }

void GPControls::WriteXLabel(ofstream& ostr)
  { if(xlabel.length())  ostr << "set xlabel \"" << xlabel << "\"" << endl; }

void GPControls::WriteYLabel(ofstream& ostr)
  { if(ylabel.length())  ostr << "set ylabel \"" << ylabel << "\"" << endl; }

void GPControls::WriteZLabel(ofstream& ostr)
  { if(zlabel.length())  ostr << "set zlabel \"" << zlabel << "\"" << endl; }


void GPControls::WriteXTics(ofstream& ostr)
  {
  if(xtics)  ostr << "set xtics"   << endl;	// Put in x-axis tic marks
  else       ostr << "set noxtics" << endl;	// No x-axis tic marks
  }

void GPControls::WriteYTics(ofstream& ostr)
  {
  if(xtics)  ostr << "set ytics"   << endl;	// Put in y-axis tic marks
  else       ostr << "set noytics" << endl;	// No y-axis tic marks
  }

void GPControls::WritePause(ofstream& ostr)
  {
  ostr << "pause(-1)" << endl;			// Issue pause command
  }

// ____________________________________________________________________________
// Z              Gnuplot Controls Standard Output Functions
// ____________________________________________________________________________

ostream& GPControls::print(ostream& ostr, bool phdr)
  {
  string hdr;
  int hl;
  if(phdr)
    {
    hdr = string("Gnuplot Plotting Controls");		// Output header
    hl = hdr.length();					// Header length
    ostr << "\n" << string(40-hl/2, ' ') << hdr;
    }

  ostr << "\n\tData File Name: " << datafile;
  ostr << "\n\tLoad File Name: " << loadfile;
  string On("On");
  string Off("Off");

  ostr << "\n\tDraw Key:         "; if(nokey)    ostr << Off;
                                    else         ostr << On;
  ostr << "\n\tData Compression: "; if(compress) ostr << On;
                                    else         ostr << Off;
  return ostr;
  }

ostream& operator<< (ostream& ostr, GPControls& GPC)
  { return GPC.print(ostr); }



#endif 							// GnuplotC.cc
