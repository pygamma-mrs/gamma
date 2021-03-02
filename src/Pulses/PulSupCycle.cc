/* PulSupCycle.cc ***********************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	Class Pulse Supercycle			  Implementation	**
**                                                                      **
**      Copyright (c) 1998                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
** This class sets up a pulse supercycle for GAMMA. A pulse supercycle	**
** consists of a defined pulse cycle (a cycled pulse waveform) which is	**
** repeated (i.e. again cycled) over a specified number of steps with	**
** defined phase changes.						**
**                                                                      **
*************************************************************************/

#ifndef _PulSupCycle_cc_		// Is this file already included?
#define _PulSupCycle_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulSupCycle.h>		// Include the header
#include <Pulses/PulWaveform.h>		// Include pulse waveforms
#include <Pulses/PulCycle.h>		// Include pulse cycles
#include <Matrix/row_vector.h>		// Must know about row vectors
#include <string>			// Must know stdlibc++ strings
#include <Basics/StringCut.h>           // Include Gform and Gdec functions
#include <iostream>                     // Include input output streams
#include <stdlib.h>

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                   CLASS PULSE SUPERCYCLE ERROR HANDLING
// ____________________________________________________________________________

 
void PulSupCycle::SCycerror(int eidx, int noret) const
 
        // Input                SCyc	: Pulse Supercycle (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag 
        // Output               none    : Error Message Output
 

  {
  cout << "\nClass Pulse Supercycle: ";
  switch(eidx)
    {
    case 0:								// (0)
      cout << "Program Aborting....";
      break;
    case 1:								// (1)
      cout << "Error During Construction";
      break;
    case 30:								// (30)
      cout << "Step Hamiltonian Access, Hamiltonian Does Not Exist";
      break;
    case 31:								// (31)
      cout << "Build Step Hamiltonians Before Their Access";
      break;
    case 32:								// (32)
      cout << "Step Propagator Access, Propagator Does Not Exist";
      break;
    case 52:								// (52)
      cout << "Step Propagator Access, Superop. Propagator Does Not Exist";
      break;
    default:
      cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) cout << ".\n";
  else       cout << ".";
  }
 
 
void volatile PulSupCycle::SCycfatality(int eidx) const
 
        // Input                SCyc	: Pulse Supercycle (this)
        //                      eidx    : Error flag
        // Output               none	: Stops execution & error Message
 
  {
  SCycerror(eidx);
  if(eidx) SCycerror(0);
  cout << "\n";
  exit(-1);
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            CLASS PULSE SUPERCYCLE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


PulSupCycle::PulSupCycle()

	// Input	none	:
	// Output	SCyc	: A NULL PulSupCycle (this)

  {
  SCycname    = "";			// No pulse cycle name
  SCycnosteps = 0;			// No pulse cycle steps
  }


PulSupCycle::PulSupCycle(const row_vector& pcsteps, const string& pcname)

	// Input	pcsteps	: Pulse train cycle step points (phase)
	//		pcname  : Pulse train cycle name
	// Output	SCyc	: A new pulse cycle (this)

  {
  SCycnosteps = pcsteps.size();		// Set pulse cycle # steps
  SCycsteps   = pcsteps;			// Set pulse cycle steps
  SCycname    = pcname;			// Set pulse cycle name
  }


PulSupCycle::PulSupCycle(const PulSupCycle& SCyc1)

	// Input	SCyc1	: Pulse Supercycle
	// None		SCyc	: Pulse Supercycle (this), identical
	//			  copy of SCyc1

  {
  SCycnosteps = SCyc1.SCycnosteps;	// Copy the number of steps
  SCycname    = SCyc1.SCycname;		// Copy pulse cycle name
  SCycsteps   = SCyc1.SCycsteps;	// Set pulse cycle steps
  }


PulSupCycle::PulSupCycle(const PulCycle& Cyc)

        // Input        SCyc    : Pulse supercycle (this)
        //              Cyc     : Pulse cycle
        // None         SCyc    : Pulse Supercycle (this) made from SCyc1

  {
  SCycnosteps = Cyc.steps();		// Copy the number of steps
  SCycname    = Cyc.name();		// Copy pulse cycle name
  SCycsteps   = Cyc.values();		// Copy pulse cycle steps
  }


// ------------------------------ Destruction ---------------------------------

PulSupCycle::~PulSupCycle() { }

	// Input	SCyc	: A pulse cycle (this)
	// Output	none	: SCyc is destructed


// ------------------------------- Assignment ---------------------------------


PulSupCycle& PulSupCycle::operator = (const PulSupCycle& SCyc1)

	// Input	SCyc1	: Pulse Supercycle
	// None		SCyc	: Pulse Supercycle (this) copied from SCyc1

{
  SCycnosteps = SCyc1.SCycnosteps;		// Copy the number of steps
  SCycname    = SCyc1.SCycname;		// Copy pulse cycle name
  SCycsteps   = SCyc1.SCycsteps;		// Copy pulse cycle steps

  return (*this);
}


void PulSupCycle::operator = (const PulCycle& Cyc)

	// Input	SCyc	: Pulse Supercycle (this)
	//		Cyc	: Pulse cycle
	// None		SCyc	: Pulse Supercycle copied from cycle

  {
  SCycnosteps = Cyc.steps();		// Copy the number of steps
  SCycname    = Cyc.name();		// Copy pulse cycle name
  SCycsteps   = Cyc.values();		// Copy pulse cycle steps
  }


// ____________________________________________________________________________
// B                CLASS PULSE SUPERCYCLE ACCESS FUNCTIONS
// ____________________________________________________________________________

// -------------------- Functions Over Full Supercycle ------------------------

int        PulSupCycle::steps()   const { return SCycnosteps; }
string     PulSupCycle::name()    const { return SCycname; }
row_vector PulSupCycle::values()  const { return SCycsteps; }

	// Input	SCyc	: A pulse cycle (this)
	// Output	steps	: SCyc steps
	// 		name	: SCyc name
        //              values  : Array of phi values

// ----------------- Functions For Specific Supercycle Step -------------------

complex PulSupCycle::value(int i)  const { return SCycsteps.get(i); }
double  PulSupCycle::phase(int i)  const { return SCycsteps.getRe(i); }



// ____________________________________________________________________________
// D                 CLASS PULSE SUPERCYCLE PLOTTING FUNCTIONS
// ____________________________________________________________________________


/*
void PulSupCycle::GP(int split, int ends) const

	// Input		SCyc	: Pulse Supercycle
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*step length/10
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*step length
        // Output               none	: Pulse Supercycle plot
	//				  is made interactively using
	//				  Gnuplot.  Plot is phase vs. "time"

  {
  if(!SCycnosteps) return; 			// Quit if no steps defined
  split = abs(split);				// Amount space between steps
  int spf = 1;
  if(split) spf =2;
  int endadd = 0; 
  if(ends) endadd += 4;
  double delt = double(split)/10.0;
  double endt = double(ends);
  string Aname = SCycname + "_TvsP.asc";
  string Gname = SCycname + "_TvsP.gnu";
  row_vector data(spf*2*SCycnosteps+1+endadd);
  double pttime = 0;
  double steplen = 1.0;
  int j=0;
  if(ends)					// Add horizontal line at
    {						// start of zero phase
    data.put(complex(pttime,0),j++);
    pttime += endt;
    data.put(complex(pttime,0),j++);
    }
  double ptphase = SCycsteps.getRe(0);		// Phase of 1st cycle step
  data.put(complex(pttime,ptphase),j++);	// Set start of this step
  for(int i=0; i<SCycnosteps; i++)		// Loop over rest of steps
    {
    pttime += steplen;				// Advance through this step
    data.put(complex(pttime,ptphase), j++);	// Set end of this step
    if(split)
      {
      data.put(complex(pttime,0), j++);
      pttime += delt;
      data.put(complex(pttime,0), j++);
      }
    if(i < SCycnosteps-1);
      {
      ptphase = SCycsteps.getRe(i+1);		// Set start next step
      data.put(complex(pttime,ptphase), j++);
      }
    }
  if(ends)
    {
    data.put(complex(pttime,0),j++);
    pttime += endt;
    data.put(complex(pttime,0),j++);
    }
  GP_xy(Aname, data);
  GP_xyplot(Gname, Aname, -1);
  }


void PulSupCycle::GP(const PulWaveform& PW, int split, int ends) const

	// Input		SCyc	: Pulse Supercycle
	//			PW	: Pulse Waveform 
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*step length/10
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Ends, length #*cyclelength/10
        // Output               none	: Pulse Supercycle plot
	//				  is made interactively using
	//				  Gnuplot. Plot is phase vs time

  {
  if(!SCycnosteps) return; 			// Quit if no steps defined
  if(!PW.steps)
    {
    GP(split, ends);				// If no waveform use overload
    return;
    }
  double totime = SCycnosteps*PW.length();	// Length of pulse cycle
  int endadd = 0; 				// Add no points due to ends
  if(ends) endadd += 2;				// If ends used, add 4 points
  double endt = double(ends)*totime/10.;	// Length of ends
  split = abs(split);				// Amount space between steps
  double delt = double(split)*PW.length(0)/10.0;// Length of spacers between steps
  int spf = 1;					// Assume no spacers
  if(split) spf=2;				// If spacers present double points
  string Aname = SCycname + "_TvsP.asc";		// ASCII output filename
  string Gname = SCycname + "_TvsP.gnu";		// GNUplot output filename
  int npts = spf*2*SCycnosteps*PW.steps();
  row_vector data(npts+endadd);			// Vector of plotted points
  double pttime = 0;				// Time in pulse cycle
  double Cphase, WFphase, Value;
  int k=0, j=0;
  if(ends)					// Add horizontal line at
    {						// start of zero phase
    data.put(complex(pttime,0),k++);		// Start at 0,0
    pttime += endt;				// Increment time to cycle start
    data.put(complex(pttime,0),k++);		// Set point at line end
    }
  for(int i=0; i<SCycnosteps; i++)		// Loop over cycle steps
    {
    Cphase = SCycsteps.getRe(i);			// Phase of cycle step
    for(j=0; j<PW.steps(); j++)			// Loop over waveform steps
      {
      WFphase = PW.phase(j);			// Phase of waveform step
      Value = WFphase + Cphase;
//      Value = WFphase + Cphase;
      Value = PW.strength(j);
      data.put(complex(pttime,Value),k++);	// Set start of this step
      pttime += PW.length(j);			// Advance through this step
      data.put(complex(pttime,Value), k++);	// Set end of this step
      if(split)
        {
        if((i!=SCycnosteps-1)||(j!=PW.steps()-1))// If gap desired between
          {					// steps, draw an added
          data.put(complex(pttime,0), k++);	// horizontal line
          pttime += delt;
          data.put(complex(pttime,0), k++);
          }
        }
      }
    }
  if(ends)
    {
    data.put(complex(pttime,0),k++);
    pttime += endt;
    data.put(complex(pttime,0),k++);
    }
  GP_xy(Aname, data);
  GP_xyplot(Gname, Aname, -1);
  }

*/


// ____________________________________________________________________________
// Z                 CLASS PULSE SUPERCYCLE I/O FUNCTIONS
// ____________________________________________________________________________


// --------------------------- ASCII Output -----------------------------------


ostream& PulSupCycle::printBase(ostream &ostr, double SCyclen) const

	// Input		SCyc	: Pulse Supercycle
        //                      ostr	: Output stream
	//			SCyclen	: Length of 1 cycle (sec)
	//			full	: Flag for output amount
        // Output               none	: Pulse Supercycle base info is sent
	//				  to the output stream

  {
  ostr << "\n\tSupCycle Steps:                     " << SCycnosteps;
  if(SCyclen)
    {
    ostr << "\n\tSupCycle Length:                    ";
    if(SCyclen > 0.1)         ostr << SCyclen      << " sec";
    else if(SCyclen > 0.0001) ostr << SCyclen*1.e3 << " msec";
    else                       ostr << SCyclen*1.e6 << " nsec";
    ostr << "\n\tSupCycle Spectral Width:            ";
    double SW = 1.0/SCyclen;
    if(SW < 1000)        ostr << SW << " Hz";
    else if(SW < 100000) ostr << SW*1.e-3 << " KHz";  
    else                 ostr << SW*1.e-6 << " MHz";
    ostr << "\n\tSupCycle Step Length:               ";
    double cslen = SCyclen/SCycnosteps;
    if(cslen > 0.1)         ostr << cslen      << " sec";
    else if(cslen > 0.0001) ostr << cslen*1.e3 << " msec";
    else                    ostr << cslen*1.e6 << " nsec";
    ostr << "\n\tSupCycle Step Spectral Width:       ";
    double cssw = 1.0/cslen;
    if(cssw < 1000)        ostr << cssw << " Hz";
    else if(cssw < 100000) ostr << cssw*1.e-3 << " KHz";  
    else                   ostr << cssw*1.e-6 << " MHz";
    }
  return ostr;
  }


ostream& PulSupCycle::printSteps(ostream &ostr) const

	// Input		SCyc	: Pulse Supercycle
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: SCyc step info is sent
	//				  to the output stream

  {
  string startln = "\n\t";
  string spacer = "  ";			// Column separator
  string dblsp = "    ";
  ostr << "\n\tSupCycle Step Phases:\n" << startln;
  int k, st=0, ncol=4;			// Number of columns
  for(k=0; k<ncol; k++)
    ostr << "Step" << spacer << " Phase " << dblsp;
  ostr << startln;
  for(k=0; k<ncol; k++)
    ostr << "----" << spacer << "-------" << dblsp;
  ostr << startln;
  double sphase;
  for(k=0, st=0; k<SCycnosteps; k++,st++)
    {
    if(st == ncol)
      {
      ostr << startln;
      st = 0;
      }
    sphase = SCycsteps.getRe(k);			// Step phase
    ostr << Gform("%3i.", k) << spacer;		// Output step number
    ostr << Gform("%7.2f", sphase);		// Output step phase
    ostr << dblsp;				// Write a column spacer
    }
  return ostr;
  }


ostream& PulSupCycle::print(ostream &ostr, int full) const

	// Input		SCyc	: Pulse Supercycle
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: Pulse Supercycle info is sent
	//				  to the output stream

  {
  ostr << "\n\n\t\t";
  if(!SCycnosteps)
    {
    ostr << "\tEmpty Pulse Supercycle\n\n";
    return ostr;
    }
  if(full) ostr << "\t";
  ostr << "Pulse Supercycle " << SCycname << "\n";

  printBase(ostr);			// Print out base cycle info
  if(full)				// For full output also
    printSteps(ostr);			// print cycle steps
  ostr << "\n\n"; 
  return ostr;
  }


ostream &operator << (ostream &ostr, const PulSupCycle &SCyc)

	// Input		SCyc	: Pulse Supercycle
        //                      ostr	: Output stream
        // Output               none	: Pulse Supercycle info is sent
	//				  to the output stream

  {
  SCyc.print(ostr);
  return ostr;
  }

 
#endif						// PulSupCycle.cc
