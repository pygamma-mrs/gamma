/* PulWaveform.cc ************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	Class Pulse Train Waveform             	   Implementation	**
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
** This class sets up a pulse waveform for implementation in GAMMA.	**
** The waveform in designated in steps, each step having a specified	**
** strength, rf-phase, and applied time. 				**
**                                                                      **
*************************************************************************/

#ifndef _PulWaveform_cc_			// Is this file already included?
#define _PulWaveform_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulWaveform.h>		// Include the header
#include <Pulses/PulAuxil.h>		// Include pulse auxilary
#include <Basics/Gconstants.h>		// Must know about DEG2RAD
#include <Basics/Gutils.h>		// Know about GAMMA errors
#include <Basics/StringCut.h>		// Include Gform and Gdec functions
#include <GamIO/Ggnuplot.h>		// Know aobout Gnuplot functions
#include <GamIO/FrameMaker.h>		// Know aobout FrameMaker functions
#include <stdlib.h>
#include <string>			// Must know about strings
#include <iostream>                     // Include input output streams

using std::string;			// Using libstdc++ strings	
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS PULSE WAVEFORM ERROR HANDLING
// ____________________________________________________________________________
 
/*         Input                PWF	: Pulse waveform (this)
                                eidx    : Error flag
                                noret   : Return flag 
           Output               none    : Error Message Output

  The following error messages use the defaults set in the Gutils package

               Case                          Error Message

               (0)                     Program Aborting.....
               (1)                     Problems With Input File Stream
               (2)                     Problems With Output File Stream
               (9)                     Problems During Construction
               default                 Unknown Error                        */

void PulWaveform::PWFerror(int eidx, int noret) const
  {
  string hdr("Pulse Waveform");
  string msg;
  switch(eidx)
    {
    case 20: msg = string("Error In Waveform Step Timing");								// (2)
             GAMMAerror(hdr,msg,noret);  break;				// (20)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }
 
void volatile PulWaveform::PWFfatality(int eidx) const
  {
  PWFerror(eidx, 1);				// Normal non-fatal error
  if(eidx) PWFerror(0);				// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  CLASS PULSE WAVEFORM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


// -------------------------- Empty Construction ------------------------------


PulWaveform::PulWaveform()
  {
  WFsteps   = 0;			// No pulse waveform steps
  WFname    = "";			// No pulse waveform name
  WFtp      = 0;			// No pulse waveform length
  rad       = 0;			// Pulse steps in degrees
  }

// ----------------------------------------------------------------------------
//              Constructors That Fill Up The Waveform Directly
// ----------------------------------------------------------------------------

	// Input	wfsteps	: Pulse waveform step points (gB1,phase)
	//		wftimes : Pulse waveform step lengths
	//		wfname  : Pulse waveform name
	//		wfrad   : Pulse waveform radians flag
	// Output	PWF	: A new pulse waveform (this)


PulWaveform::PulWaveform(const row_vector& wfsteps,
                    const row_vector& wftimes, const string& wfname, int wfrad)
  {
  WFsteps = wfsteps.size();		// Set pulse waveform # steps
  WFvals  = wfsteps;			// Set pulse waveform steps
  WFtimes = wftimes;			// Set pulse waveform step lengths (s)
  WFname  = wfname;			// Set pulse waveform name
  WFtp    = Re(wftimes.sum());		// Set pulse waveform total length (s)
  rad     = wfrad;			// Copy the radians flag
  }

// -------------------------- Self Construction ------------------------------

PulWaveform::PulWaveform(const PulWaveform& PWF1)
  {
  WFsteps = PWF1.WFsteps;		// Copy the number of steps
  WFname  = PWF1.WFname;		// Copy pulse waveform name
  WFtp    = PWF1.WFtp;			// Copy pulse waveform length
  WFvals  = PWF1.WFvals;		// Set pulse waveform steps
  WFtimes = PWF1.WFtimes;		// Set pulse waveform step lengths (s)
  rad     = PWF1.rad;			// Copy the radians flag
  }

// ------------------------------ Destruction ---------------------------------

PulWaveform::~PulWaveform() { }

	// Input	PWF	: A pulse waveform (this)
	// Output	none	: PWF is destructed

// ------------------------------- Assignment ---------------------------------

PulWaveform& PulWaveform::operator = (const PulWaveform& PWF1)
{
  WFsteps = PWF1.WFsteps;		// Copy the number of steps
  WFname  = PWF1.WFname;		// Copy pulse waveform name
  WFtp    = PWF1.WFtp;			// Copy pulse waveform length
  WFvals  = PWF1.WFvals;		// Set PT waveform steps
  WFtimes = PWF1.WFtimes;		// Set PT waveform step lengths (sec)
  rad     = PWF1.rad;			// Copy the radians flag

  return *this;
}


// ____________________________________________________________________________
// B               CLASS PULSE WAVEFORM ACCESS FUNCTIONS
// ____________________________________________________________________________

// --------------------- Functions Over Full Waveform -------------------------

	// Input	PWF	: A pulse waveform (this)
	// Output	steps	: PWF steps
	// 		name	: PWF name
	//		length  : PWF length (sec)
        //              values  : Array of { gamB1, phi } values
	//		lengths : Array of times

int        PulWaveform::steps()   const { return WFsteps; }
string     PulWaveform::name()    const { return WFname; }
double     PulWaveform::length()  const { return WFtp; }
row_vector PulWaveform::values()  const { return WFvals; }
row_vector PulWaveform::lengths() const { return WFtimes; }

// ----------------- Functions For Specific Waveform Step ---------------------

	// Input	PWF	: A pulse waveform (this)
	// Output	strength: Step rf-field strength (Hz)
	//		phase	: Step phase value (degrees)
        //              value   : Step { gamB1, phi } values
	//              length  : Step length (sec)

double  PulWaveform::strength(int i) const { return WFvals.getRe(i); }
complex PulWaveform::value(int    i) const { return WFvals.get(i); }
double  PulWaveform::length(int   i) const { return WFtimes.getRe(i); }
double  PulWaveform::phase(int    i) const
  { 
  double PHF = 1;
  if(rad) PHF = 1/DEG2RAD;
  return WFvals.getIm(i)*PHF;
  }

// ____________________________________________________________________________
// C               CLASS PULSE WAVEFORM AUXILIARY FUNCTIONS
// ____________________________________________________________________________

// ----------------------------- Step Lengths ---------------------------------
 
	// Input	PWF	: A pulse waveform (this)
	//		cutoff	: Minimum length to consider
        // Output (max) tsum    : Length from longest waveform step
        // Output (min) tsum    : Length from shortest waveform step

double PulWaveform::maxlength( ) const
  {
  double len = 0;
  for(int j=0; j<WFsteps; j++)
    if(length(j) > len) len = length(j);
  return len;
  }


double PulWaveform::minlength(double cutoff) const
  {
  double len=1.e50, lenj;
  for(int j=0; j<WFsteps; j++)
    {
    lenj = length(j);
    if(lenj<len && lenj>=cutoff) len = lenj;
    }
  return len;
  }

// ---------------------------- Step Strengths --------------------------------


double PulWaveform::maxgamB1( ) const
 
	// Input	PWF	: A pulse waveform (this)
        // Output       gB1     : Strength of strongest waveform step
 
  {
  double gB1 = 0;
  for(int j=0; j<WFsteps; j++)
    if(WFvals.getRe(j) > gB1) gB1 = WFvals.getRe(j);
  return gB1;
  }


double PulWaveform::mingamB1( ) const
 
	// Input	PWF	: A pulse waveform (this)
        // Output       gB1     : Strength of weakest waveform step
        // Output       tsum    : Length from shortest waveform step
 
  {
  double gB1 = 1.e30;
  for(int j=0; j<WFsteps; j++)
    if(WFvals.getRe(j) < gB1) gB1 = WFvals.getRe(j);
  return gB1;
  }


// --------------------------- Variation Checks -------------------------------

	// Input	PWF	: A pulse waveform (this)
	// Output	T/F     : Returns true if the RF field strength
	//			  is constant through out the waveform
	// Output	T/F     : Returns true if the RF field phase
	//			  is constant through out the waveform
	// Output	T/F     : Returns true if the step time
	//			  is constant through out the waveform

bool PulWaveform::gamB1const() const
  {
  double gB1 = WFvals.getRe(0);
  for(int i=1; i<WFsteps; i++)
    if(fabs(WFvals.getRe(i)-gB1) > 1.e-10) return false;
  return true;
  }

bool PulWaveform::phaseconst() const

  {
  double ph = WFvals.getIm(0);
  for(int i=1; i<WFsteps; i++)
    if(fabs(WFvals.getIm(i)-ph) > 1.e-10) return false;
  return true;
  }

bool PulWaveform::timeconst() const
  {
  double st = WFtimes.getRe(0);
  for(int i=1; i<WFsteps; i++)
    if(fabs(WFtimes.getRe(i)-st) > 1.e-10) return false;
  return true;
  }

// ---------------------------- Step Counting ---------------------------------
    
	// Input	PWF	: A pulse waveform (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of waveform steps needed
        //                        to evolve for time td

double PulWaveform::steps(double td) const
  {
  if(td<=0 || !WFtp) return 0;			// Exit if no waveform,time
  int nloops=-1;                                // First find out how many
  double tloop=0;                               // full pulses are needed
  while(tloop <= td)
    {
    nloops++;
    if(tloop == td)
      return double(nloops*WFsteps);
    tloop += length();
    }
  double lsteps = double(nloops*WFsteps);       // Steps evovled so far
  double t = tloop - length();                  // Time evolved so far
  for(int i=0; i<WFsteps; i++) 			// Now find out how many
    {                                           // waveform steps we need
    t += length(i);
    if(t == td) return double(i+1+lsteps);
    else if(t > td)
      {
      t -= length(i);
      t = td-t;
      return double(i) + t/length(i) + lsteps;
      }
    }
  PWFfatality(20); 				// Should never make it here
  return 0;
  }
 

int PulWaveform::fullsteps(double td) const

	// Input	PWF	: A pulse waveform (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full waveform steps
        //                        that can occur in the time td

 { return int(floor(steps(td))); }


// -------------------------- Waveform Counting -------------------------------


double PulWaveform::WFs(double td) const
   
	// Input	PWF	: A pulse waveform (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of waveforms needed
        //                        to evolve for time td
 
  {
  if(td<=0 || !WFtp) return 0;			// Exit if no waveform or time
  int nloops = -1;				// First find out how many
  double tloop = 0;				// full waveforms are needed
  while(tloop <= td)				// Loop full waveforms & count
    {
    nloops++;
    if(tloop == td) return double(nloops);
    tloop += WFtp;
    }
  td = tloop - WFtp;				// Time left to evolve
  return td/WFtp + nloops;			// Total wavforms needed
  } 
 

int PulWaveform::fullWFs(double td, double cut) const

	// Input	PWF	: Pulse Train Waveform
        //              td      : An evolution time (sec)
	//		cut	: Zero cutoff (sec)
        // Output       steps   : Number of full waveforms that
        //                        are within time td

  {
  if(cut < 0) cut = 0;			// Insure cutoff non-negative
  if(td<=0 || WFtp<=cut) return 0;	// Exit if no waveform,time
  int wfs=0, OK=1;			// Assume none within td
  while(OK)				// See if there are any
    {					// waveforms within td
    td -= WFtp;
//    if(td >= 0) wfs++;
    if(td >= -cut) wfs++;
    else        OK=0;
    }
  return wfs;
  }


double PulWaveform::sumlength(int i) const
 
	// Input	PWF	: A pulse waveform (this)
        //              i       : Step in waveform
        // Output       tsum    : Length from waveform start
        //                        over i steps (sec).
	// Note			: Time from start (before step i=0) to the
	//			  end of step i-1 (before step i=i)
 
  {
  if(!WFtp) return 0;			// Quit if WF has no length
  double len = 0;			// Start with length of zero 
  while(i > WFsteps)			// Loop over full waveforms
    {					// and sum up the their lengths
    len += WFtp;
    i -= WFsteps;
    }
  for(int j=0; j<i; j++)		// Loop over partial waveform steps
    len += WFtimes.getRe(j);		// and sum up their lengths
  return len;
  }


// --------------------------- Waveform Scaling -------------------------------


void PulWaveform::scalegB1(double sf)

	// Input	PWF	: A pulse waveform (this)
	// 		sf	: A scaling factor
	// Output	void 	: All step field strengths in PWF
	//			  are multiplied by sf.  The exception are
	//			  steps of zero length (ideal pulses)

  {
  for(int i=0; i<WFsteps; i++)
    {
    if(WFtimes.getRe(i))
      WFvals.put(complex(sf*WFvals.getRe(i), WFvals.getIm(i)), i);
    }
  }
 
// ____________________________________________________________________________
// D                 CLASS PULSE WAVEFORM PLOTTING FUNCTIONS
// ____________________________________________________________________________

//------------------------- Generic Plotting Functions ------------------------

void PulWaveform::getIdeal(double& gB1, double& ptt, int i) const

	// Input		PWF	: Pulse Train Waveform
	//			gB1     : Step intensity
	//			ptt     : Step length
	//			i	: Step index
	// Output		none    : gB1 and ptt are modified to
	//				  reflect how an ideal pulse should
	//				  be plotted for step i.

  {
  ptt = 1;					// Default length
  double stept = minlength();			// Shortest non-zero step
  if(stept!=0 && stept<ptt) ptt = stept;
  double stepgB1 = WFvals.getRe(i);		// Step i gamB1 -> ideal angle
  ptt *= stepgB1/90.0;				// Scale relative to 90 pulse 

  gB1 = 1;					// Default strength
  for(int j=0; j<steps(); j++)			// Set strength to strongest
    if(WFvals.getRe(j)>gB1) gB1=WFvals.getRe(j);// Non-zero strength
  }


row_vector PulWaveform::IvsT(int split, int ends, int N) const

	// Input		PWF	: Pulse Train Waveform
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*10%*biggest step
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*buggest step
	//			N       : Number of waveforms
        // Output               none	: Pulse Train Waveform plot vector
	//				  is made interactively of the
	//				  RF-intensity vs time.
	// Note				: Ideal pulses are set by a step
	//				  having zero length but having
	//				  an non-zero "gamB1" setting.  In
	//				  such cases "gamB1" is the pulse angle
	// Note				: For plotting purposes, the length of
	//				  an ideal 90 pulse will be taken to 
	//				  be the same as the shortest non-zero
	//				  waveform step.  Others scale with 
	//				  ideal pulse angle.

  {
  if(!WFsteps) return row_vector(); 		// Return empty if no points

//         Figure Out the Size Of Data Vector Output As Well As
//        Lengths Of The Ends and Interdispersed Step Splittings

  int spf = 0;					// Assume NOT splitting steps
  split = abs(split);				// Insure split length >= 0
  if(split) spf++;				// If splitting set split flag
  int endf = 0; 				// Assume NO ends plotted
  if(ends) endf++;				// If ends wanted, set end flag
  double delt, endt;				// Lengths of splits and ends
  double tstep1 = maxlength();			// Biggest WF step length
  delt = double(split)*0.1*tstep1;		// This will be split length
  endt = double(ends)*tstep1;			// This will be end length 
  int      npts  = 2*N*WFsteps;			// Plotted waveform points
           npts += 2;				// Add start & finish pts
  if(spf)  npts += 2*(N*WFsteps-1);		// More points if spaced apart
  if(endf) npts += 2;				// More points if ends added
  row_vector data(npts, complex0);		// Make array for PWF plot
  double pttime = 0;				// Start at time 0
  int j=0;					// Global point index
  if(ends)					// Begin output 
    {						// If ends to be added
    data.put(complex(pttime,0),j++);		//	Add point 0,0
    pttime += endt;				//	Next pt   endt,0
    }

//       Make One Waveform, Including Breaks Between Steps If Desired  
//	    Time On The Horizontal Axis, gamB1 On the Vertical Axis
//  (1*endf + 1@start + 2*N*WFsteps + 2*spf*[N*WFsteps-1] + 1@end + 1*endf)

  data.put(complex(pttime,0), j++);		// Start at gamB1=0
  for(int k=0; k<N; k++)			// Loop over desired WFs
    {
    double ptgB1 = WFvals.getRe(0);		// Step 0 gamB1
    double ptt = WFtimes.getRe(0);		// Step 0 length
    for(int i=0; i<WFsteps; i++)
      {
      ptgB1 = WFvals.getRe(i);			// Step i gamB1
      ptt = WFtimes.getRe(i);			// Step i length
      if(!ptt) getIdeal(ptgB1, ptt, i);		// Adjust if ideal pulse
      data.put(complex(pttime,ptgB1), j++);	// Add step i start point
      pttime += ptt;				// Time at step i end
      data.put(complex(pttime,ptgB1), j++);	// Add step i finish point
      if(split) 				// If spacer desired plotted
        {					// then add some horizontal
        if(k!=N-1 || (k==N-1 && i!=WFsteps-1))	// between steps and between 
          {					// WFs
          data.put(complex(pttime,0), j++);	//	Add point (time,0)
          pttime += delt;			//	Delay between pulses
          data.put(complex(pttime,0), j++);	//	Add point (time,0)
          }
        }
      }
    }
  data.put(complex(pttime,0),j++);		// Finish at gamB1=0

  if(ends)					// If ends to be added
    {						// draw horizontal line
    pttime += endt;				//	Next pt   endt,0
    data.put(complex(pttime,0),j++);		//	Draw to   endt,0
    }
  return data;
  }


row_vector PulWaveform::PvsT(int split, int ends, int N, double ph) const

	// Input		PWF	: Pulse Train Waveform
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*.1*1st pulse length
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*1st pulse
	//			N       : Number of waveforms
	//			ph	: Added phase
        // Output               none	: Pulse Train Waveform plot vector
	//				  is made interactively of the
	//				  RF-intensity vs time.

  {
  if(!WFsteps) return row_vector(); 		// Return empty if no points

//         Figure Out the Size Of Data Vector Output As Well As
//        Lengths Of The Ends and Interdispersed Step Splittings

  int spf = 0;					// Assume NOT splitting steps
  split = abs(split);				// Insure split length >= 0
  if(split) spf++;				// If splitting set split flag
  int endf = 0; 				// Assume NO ends plotted
  if(ends) endf++;				// If ends wanted, set end flag
  double delt, endt;				// Lengths of splits and ends
  double tstep1 = maxlength();			// Biggest WF step length
  delt = double(split)*0.1*tstep1;		// This will be split length
  endt = double(ends)*tstep1;			// This will be end length 
  int      npts  = 2*N*WFsteps;			// Plotted waveform points
           npts += 2;				// Add start & finish pts
  if(spf)  npts += 2*(N*WFsteps-1);		// More points if spaced apart
  if(endf) npts += 2;				// More points if ends added
  row_vector data(npts, complex0);		// Make array for PWF plot
  double pttime = 0;				// Start at time 0
  int j=0;					// Global point index
  if(ends)					// Begin output 
    {						// If ends to be added
    data.put(complex(pttime,0),j++);		//	Add point 0,0
    pttime += endt;				//	Next pt   endt,0
    }

//       Make One Waveform, Including Breaks Between Steps If Desired  
//	    Time On The Horizontal Axis, gamB1 On the Vertical Axis
//  (1*endf + 1@start + 2*N*WFsteps + 2*spf*[N*WFsteps-1] + 1@end + 1*endf)

  double tmp;
  data.put(complex(pttime,0), j++);		// Start at phase=0
  for(int k=0; k<N; k++)			// Loop over desired WFs
    {
    double ptphi = WFvals.getIm(0) + ph;	// Step 0 phase
    ptphi = 360.0*fmod(ptphi, 360.0);		// Keep phase in [0, 360)
    if(fabs(ptphi-360.0) < 1.e-10) ptphi=0.0;
    double ptt = WFtimes.getRe(0);		// Step 0 length
    for(int i=0; i<WFsteps; i++)
      {
      ptphi = WFvals.getIm(i) + ph;		// Step i phase
      ptt = WFtimes.getRe(i);			// Step i length
      if(!ptt) getIdeal(tmp, ptt, i);		// Adjust if ideal pulse
      data.put(complex(pttime,ptphi), j++);	// Add step i start point
      pttime += ptt;				// Time at step i end
      data.put(complex(pttime,ptphi), j++);	// Add step i finish point
      if(split) 				// If spacer desired plotted
        {					// then add some horizontal
        if(k!=N-1 || (k==N-1 && i!=WFsteps-1))	// between steps and between 
          {					// WFs
          data.put(complex(pttime,0), j++);	//	Add point (time,0)
          pttime += delt;			//	Delay between pulses
          data.put(complex(pttime,0), j++);	//	Add point (time,0)
          }
        }
      }
    }
  data.put(complex(pttime,0),j++);		// Finish at phase=0

  if(ends)					// If ends to be added
    {						// draw horizontal line
    pttime += endt;				//	Next pt   endt,0
    data.put(complex(pttime,0),j++);		//	Draw to   endt,0
    }
  return data;
  }


//------------------------- Gnuplot Plotting Functions ------------------------


void PulWaveform::GP(int type, int split, int ends, int N) const

	// Input		PWF	: Pulse Train Waveform
        //                      type    : Type of plot to create
	//				   0: Waveform time vs phase
	//				   1: Waveform time vs gamB1
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*.1*1st pulse length
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*1st pulse
	//			N       : Number of waveforms
        // Output               none	: Pulse Train Waveform plot
	//				  is made interactively using
	//				  Gnuplot.

  {
  switch(type)
    {
    case 0:					// Time vs. WF Phase
      {
      string Aname = WFname + "_TvsP.asc";
      string Gname = WFname + "_TvsP.gnu";
      row_vector data = PvsT(split, ends, N);
      GP_xy(Aname, data);
      GP_xyplot(Gname, Aname, -1);
      }
      break;
    case 1:					// Time vs. WF gamB1
      {
      string Aname = WFname + "_TvsB1.asc";
      string Gname = WFname + "_TvsB1.gnu";
      row_vector data = IvsT(split, ends, N);
      GP_xy(Aname, data);
      GP_xyplot(Gname, Aname, -1);
      }
    break;
    default:
      {
      break;
      }
    }
  }


//----------------------- FrameMaker Plotting Functions -----------------------


void PulWaveform::FM(int type, int split, int ends, int N) const

	// Input		PWF	: Pulse Train Waveform
        //                      type    : Type of plot to create
	//				   0: Waveform time vs phase
	//				   1: Waveform time vs gamB1
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*.1*1st pulse length
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*1st pulse
	//			N       : Number of waveforms
        // Output               none	: Pulse Train Waveform plot
	//				  is output in Framemaker MIF format

  {
  switch(type)
    {
    case 0:					// Time vs. WF Phase
      {
      string Aname = WFname + "_TvsP.mif";
      row_vector data = PvsT(split, ends, N);
      FM_xyPlot(Aname, data);
      }
      break;
    case 1:					// Time vs. WF gamB1
      {
      string Aname = WFname + "_TvsB1.mif";
      row_vector data = IvsT(split, ends, N);
      FM_xyPlot(Aname, data);
      }
    break;
    default:
      {
      break;
      }
    }
  }




// ____________________________________________________________________________
// Z                 CLASS PULSE WAVEFORM I/O FUNCTIONS
// ____________________________________________________________________________


ostream& PulWaveform::printBase(ostream &ostr) const

	// Input		PWF	: Pulse Train Waveform
        //                      ostr	: Output stream
        // Output               none	: Pulse Train base parameters are
	//				  sent into the output stream

  {
  string SL = "\n\t";
  string WF = "Waveform ";
  string ST = "Step";
  string SP = "Spectral ";
  string FL = "Field ";
  string STR = "Strength";
  ostr << SL << WF << ST << "s:                   " << Gdec(WFsteps,4);
  ostr << SL << WF << "Length:                  ";
  printTime(ostr, WFtp);
  if(WFtp)
    {
    double SW = 1.0/WFtp;
    ostr << SL << WF << SP << "Width:          ";
    printHz(ostr, SW);
    }
  if(timeconst())
    {
    ostr << SL << WF << ST << " Length:             ";
    printTime(ostr, WFtp/WFsteps);
    if(WFtp)
      {
      ostr << SL << "Corresponding Spectral Width:     ";
      printHz(ostr, WFsteps/WFtp);
      }
    }
  else
    {
    ostr << SL << WF << ST << " Length:             Variable";
    double maxt = WFtimes.getRe(WFtimes.max(1));
    ostr << " (max. ";
    printTime(ostr, maxt);
    ostr << ")";
    }
  if(gamB1const())
    {
    ostr << SL << WF << FL << STR << ":          ";
    printHz(ostr, WFvals.getRe(0));
    }
  else
    ostr << SL << WF << FL << STR << ":          Variable";
  return ostr;
  }


ostream& PulWaveform::printSteps(ostream &ostr, int full) const

	// Input		PWF	: Pulse Train Waveform
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: PWF step info is sent
	//				  to the output stream

  {
  string startln = "\n\t";
  string spacer = "  ";			// Column separator
  string dblsp = "    ";
  ostr << "\n\tWaveform Step Details:\n" << startln;
  int k, st=0, ncol=2;			// Number of columns
  for(k=0; k<ncol; k++)
    {
    ostr << "Step" << spacer;
    if(!gamB1const()) ostr << "  GamB1  ";
    ostr << " Phase " << spacer << " Angle ";
    if(!timeconst())
      {
      ostr << spacer << "    Length   ";
      if(full)
        ostr << spacer << "Summed Length";
      }
    ostr << dblsp;
    }
  ostr << startln;
  for(k=0; k<ncol; k++)
    {
    ostr << "----" << spacer;
    if(!gamB1const()) ostr << "-------" << spacer;
    ostr << "-------" << spacer << "-------";
    if(!timeconst())
      {
      ostr << spacer << "-------------";
      if(full)
        ostr << spacer << "-------------";
      }
    ostr << dblsp;
    }
  ostr << startln;
  double gB1, sphase, stime, sangle;
  double sumt;
  stime  = WFtimes.getRe(0);
  sumt = 0;
  for(k=0, st=0; k<WFsteps; k++,st++)
    {
    if(st == 2)
      {
      ostr << startln;
      st = 0;
      }
    gB1    = WFvals.getRe(k);
    sphase = WFvals.getIm(k);
    stime  = WFtimes.getRe(k);
    sangle = gB1*stime*360.;
    ostr << Gform("%3i.", k+1) << spacer;
    if(!gamB1const()) 
      {
      if(stime)
        ostr << Gform("%7.2f", gB1) << spacer;
      else
        ostr << Gform("%7.2f", 0.0) << spacer;
      }
    ostr << Gform("%7.2f", sphase);
    ostr << spacer;
    if(stime)
      ostr << Gform("%7.2f", sangle);
    else
      ostr << Gform("%7.2f", gB1);
    if(!timeconst())
      {
      ostr << spacer;
      printTime(ostr, stime);
      if(full)
        {
        ostr << spacer;
        sumt += stime;
        printTime(ostr, sumt);
        }
      }
    ostr << dblsp;
    }
  return ostr;
  }

	// Input		PWF	: Pulse Train Waveform
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: Pulse Train info is sent
	//				  to the output stream

std::ostream& PulWaveform::print(std::ostream& ostr, int full) const
  {
  if(!WFsteps)
    {
    ostr << "\n\n\t\t\tEmpty Pulse Waveform\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t\t  Pulse Waveform " << WFname << "\n";
  printBase(ostr);
  if(full) printSteps(ostr);
  ostr << "\n\n"; 
  return ostr;
  }
            
std::ostream &operator << (std::ostream& ostr, const PulWaveform& PWF)
  { return PWF.print(ostr); }

#endif						// PulWaveform.cc
