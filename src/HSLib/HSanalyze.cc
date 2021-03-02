/* HSanalyze.cc ************************************************-*-c++-*-*
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      MR Hilbert Space Analyze                     Implementation 	**
**                                                                      **
**      Copyright (c) 1993                                              **
**      Scott Smith                                                     **
**      University of Utah                                              **
**      Department of Chemistry                                         **
**      Salt Lake City, UT, 84112, USA                                  **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     **
** Resonance Experiments and Other Associated Mathematical              **
** Capabilities.  The Set of Functions Herein Allows for the Analysis   **
** of Spin Hilbert Spaces Associated with Composite Spin Systems.       **
** This module of the GAMMA MR Library provides functions designed to	**
** facilitate the analysis of the more common aspects of MR with	**
** respect to a quantum mechanical viewpoint.  As GAMMA provides full	**
** access to all its defined quantities, these routines serve just to 	**
** sort information and present it in an orderly fashion.		**
**                                                                      **
*************************************************************************/


#ifndef   HSanalyze_cc_			// Is this file already included?
#  define HSanalyze_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <HSLib/HSanalyze.h>		// Include the header
#include <HSLib/SpinOpCmp.h>		// Include composite spin ops
#include <Basics/StringCut.h>		// Include Gform function
#include <Basics/Gconstants.h>		// Include HUGE_VAL value
#include <stdlib.h>
#include <cmath>			// Inlcude HUGE_VAL_VAL

// ____________________________________________________________________________
// A                     WAVEFUNCTIONS / BASIS FUNCTIONS
// ____________________________________________________________________________
 
 
/* The wavefunctions will span the Hilbert space of the operator Op. Each
   wavefunction will be some linear combination of the default basis functions
   GAMMA uses that are inherent to spin system sys.  The linear combination
   will depend upon the current Op basis.  If Op is in its eigenbasis, wave-
   functions are eigenfunctions.  Each wavefunction corresponds to a state
   of the system.                                                            */

//std::string alphabeta[7] = {"a","b","g","d","e","w","x"};	// Spin state labels

 std::string qStatel(const spin_sys& sys, int bf)

	// Input  	     	sys	   : A basic spin system
	//			bf         : A state (basis function) index
        // Output               bflabel    : A label for the input state in terms
	//				     of individual spin states.  For example
	//				     ab for up-down nn {I1=I2=1/2}, or gb
	//				     for down-middle in {I1=1,I2=1/2}
	// Note				   : The string alphabeta is set to
	//				     handle spins with maximum I=7/2
	// Note				   : Individual spin labels are as follows:
	//					I=1/2: 1/2 - a; -1/2 - b;
	//					I=1:    1  - a;   0  - b;  -1  - g;
	//					I=3/2: 3/2 - a;  1/2 - b; -1/2 - g; -3/2 - d;

  {
  matrix PBF = sys.qStates();			// Retrieve the product basis functions
  double qnum, qstat;				// Spin I value, Iz value in bf
  int lev;
   std::string bflabel;				// This will be the basis function label
  for(int j=0; j<sys.spins(); j++)		// Loop over all the spins (in function bf)
    {
    qstat = Re(PBF(bf,j));			//	Get the Iz value of spin j in bf
    qnum = sys.qn(j);				//	This is the I value of spin j
    lev = int(qnum - qstat);			//	Spin state index, e.g. 1/2-1/2=0 -> "a"
    bflabel += alphabeta[lev];			//	Add this spin state to the bf label
    }
  return bflabel;
  }


void wf_labels( std::string* wflabels, const spin_sys& sys, gen_op &Op,
                                                  double cut, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//       	sys	: A basic spin system
	// 		Op	: General operator (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
        //		pfz	: Flag for inclusion of total Fz   (DEFAULT
	// Return	void    : The string array wflabels is filled with system
	//			  wavefunctions in the current working basis
	//			  of operator Op in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on Op input basis!

  {
  int hs = Op.dim();				// Hilbert space dimension
  if(!hs) return;				// Exit if an empty operator
  matrix B = Op.get_basis().U();		// Retrieve the basis array
  int* filter;				// Set up a working filter
  filter = new int[hs];			// Set up a working filter
  for(int i=0; i<hs; i++) filter[i] = 1; 	// But take everything (no filter)
  wf_labels(wflabels,filter,sys,B,cut,pbf,pfz);	// Use this function overload
  delete [] filter;
  return;
  }


void wf_labels( std::string* wflabels, int* filter, const spin_sys& sys,
                                const gen_op &Op, double cut, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//		filter  : An array of integers, output filter
	//       	sys	: A basic spin system
	// 		Op	: General operator (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
        //		pfz	: Flag for inclusion of total Fz   (DEFAULT
	// Return	void    : The string array wflabels is filled with system
	//			  wavefunctions in the current working basis
	//			  of operator Op in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on Op input basis!

  {
  int hs = Op.dim();				// Hilbert space dimension
  if(!hs) return;				// Exit if an empty operator
  matrix B = Op.get_basis().U();		// Retrieve the basis array
  wf_labels(wflabels,filter,sys,B,cut,pbf,pfz);	// Use this function overload
  return;
  }


void wf_labels( std::string* wflabels, const spin_sys& sys, const matrix &B,
                                                  double cut, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//       	sys	: A basic spin system
	// 		B	: A basis (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
        //		pfz	: Flag for inclusion of total Fz   (DEFAULT
	// Return	void    : The string array wflabels is filled with system
	//			  wavefunctions in the current working basis
	//			  of operator Op in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on Op input basis!

  {
  int hs = B.rows();				// Hilbert space dimension
  if(!hs) return;				// Exit if the basis is empty
  int* filter;				// Set up a working filter
  filter = new int[hs];				// Set up a working filter
  for(int i=0; i<hs; i++)			// Set the filter to take everything
    filter[i] = 1;
  wf_labels(wflabels,filter,sys,B,cut,pbf,pfz);	// Use this function overload
  delete [] filter;
  return;
  }


void wf_labels( std::string* wflabels, int* filter, const spin_sys& sys,
                              const matrix &B, double cutoff, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//		filter  : An array of integers, output filter
	//       	sys	: A basic spin system
	// 		B	: A basis (associated to sys)
	//		cutoff  : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
        //		pfz	: Flag for inclusion of total Fz   (DEFAULT
	// Return	void    : The string array wflabels is filled with system
	//			  wavefunctions in the current working basis
	//			  of operator Op in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on Op input basis!

  {
  int hs = B.rows();				// Basis dimension, spin Hilbert space
  if(!hs) return;				// Exit if the basis is empty
   std::string wflabel;
  int keep;

   std::string rlab, ilab;				// For basis function label coefficient(s)
   std::string startbf = "|";				// For basis function label start
   std::string bfind;					// For basis function label index
   std::string bfabg;					// For basis function label spin states
   std::string fzlabel;				// For basis function label Fz
   std::string endbf = ">";				// For basis funciton label end
   std::string bflabel = "";				// For basis funciton label (total)
   std::string blanks = "               ";		// For any spaces needed in labeling 

  complex zcoeff;				// For basis function coefficient value
  double rcoeff, icoeff;			// For basis function coefficient real & imag
  double fzval;					// For basis function coefficient Fz value
  int bfindl;					// For basis function coefficient index length

  int idig = 1;					// Basis function index printed digits
  while(pow(10.0, idig) < hs)		// Get the index printed digits from 
    idig++;					// the maximum possible index, hs
  int bfout = 0;				// How many basis functions in wf

//	   Looping Over All WaveFunctions, Then Loop Over All Basis
//                Functions and Build Each Wavefunction Label

  for(int wf=0; wf<hs; wf++)
    {
    bfout = 0;					// Set the output basis functions to 0
    wflabel = "";				// Start with no label here
    keep = filter[wf];				// Flag whether to make this one or not
    for(int bf=0; keep && bf<hs; bf++)		// Loop over each bf of wavefunction wf
      {
      zcoeff = B.getRe(bf,wf);			// The coefficient of basis function bf
      if(norm(zcoeff) >= cutoff) 		// Don't do small coefficients, usually [-1,1]
        {
//			Set   Up Real Coefficient of Basis Function bf

        rcoeff = B.getRe(bf,wf);		// Real coefficient of basis function bf
        rlab = "";				// Start with empty bf real coefficient label
        if(fabs(rcoeff) >= cutoff)		// Don't do small real coefficients, usually [-1,1]
          {						
          if(bfout == 0)			//	Set the sign of the real coeff
            if(rcoeff >= 0) rlab = " ";
            else            rlab = "-";
          else
            if(rcoeff >=0)  rlab = " + ";
            else            rlab = " - ";
          rlab += Gform("%5.3f", fabs(rcoeff));	//	Now add the value of the real coeff
          }
//		       Set Up Imaginary Coefficient of Basis Function bf

        icoeff = B.getIm(bf,wf);		// Imag coefficient of basis function bf
        ilab = "";				// Start with empty bf imaginary coefficient label
        if(fabs(icoeff) >= cutoff)		// Don't do small coefficients, usually [-1,1]
          {						
          if(bfout==0 && (fabs(rcoeff)<cutoff))	//	Set the sign of the imag coeff
            if(icoeff >= 0) ilab = " ";
            else            ilab = "-";
          else
            if(icoeff >=0)  ilab = "+  ";
            else            ilab = " - ";
          ilab += Gform("%5.3f", fabs(icoeff));	//	Now add the value of the imag coeff
          ilab += "i";				//	Indicate this is imaginary
          }
//		                Set Up Label of Basis Function bf

        bfind = Gdec(bf);			// Start bf index label
        bfindl = bfind.length();		// Get the current length
        if(pfz || pbf>=0) bfind += ", ";	// Add comma if labeling with more info
        bfind +=  std::string(idig-bfindl, ' ');	// Fill out the index to idig length

        bfabg = qStatel(sys, bf);		// Get bf spin states label , e.g. abaa
        if(pfz) bfabg += ", ";			// Add comma if labeling with Fz

        fzval = sys.qnState(bf);		// Fz of basis function bf (1/2, -3/2, ...)
        if(fzval < 0) fzlabel = "-";		// Start of bf Fz label, this sign
        else if(fzval > 0) fzlabel = "+";
        else fzlabel = " ";
        fzlabel += Gform("%g", fabs(fzval));		// Complex bf Fz label

        bflabel = startbf;			// Start the basis function bf label with |
        if(pbf <=0) bflabel += bfind;		// Label bf with an index number
        if(pbf >= 0) bflabel += bfabg;		// Label with spin states, e.g. abaa
        if(pfz) bflabel += fzlabel;		// Label bf with an Fz value
        bflabel += endbf;			// End the basis function bf label with >

//            Add the Contribution of Basis Function bf to Wavefunction wf's Label

        wflabel += rlab;			// Add the real coefficient of bf
        wflabel += ilab;			// Add the imag coefficient of bf
        wflabel += bflabel;			// Add the basis funciton label of bf
        bfout++;				// Track how may basis functions included
        }
      } 					// Do the next bf of this wavefunction
    wflabels[wf] = wflabel;
    } 						// Do the next wavefunction
  return;
  }


// ____________________________________________________________________________
// B                                EIGENVALUES
// ____________________________________________________________________________

/* The eigenvalues will span the Hilbert space of the operator Op. Each eigen-
   value will be a diagonal element of Op after Op has been placed into an 
   eigenbasis representation.  Each eigenvalue will then correspond to a
   particular eigenfunction of the operator. If the operator is a system
   Hamiltonian, then the eigenvalues are the energies of the system.         */


void ev_labels( std::string* evlabels, gen_op& Op, double cutoff)

	// Input 	evlabels: An array of strings
	//       	Op	: General operator
	//		cutoff  : Eigenvalues of magnitude below this are
	//			  assumed to be zero
	// Return	evlabels: An array of strings containg all the
	//			  operator eigenvalues (in the eigenbasis)
	//			  of operator Op.
	// Note			: Handles complex eigenvalues!
	// Note			: Eigenvalue ordering depends on Op eigenbasis.
	//			  Either order the eigenbasis prior to entering
	//			  function or sort the values after output

  {
  int hs = Op.dim();				// Operator dimension size
  if(!hs) return;				// Exit if NULL operator
  int* filter;				// Set up a working filter
  filter = new int[hs];				// Set up a working filter
  for(int i=0; i<hs; i++)			// Set the filter to take everything
    filter[i] = 1;
  ev_labels(evlabels,filter,Op,cutoff);		// Use this function overload
  delete [] filter;
  }


void ev_labels( std::string* evlabels, int* filter, gen_op& Op, double cutoff)

	// Input 	evlabels: An array of strings
	//		filter  : Array of integers to filter values
	//       	Op	: General operator
	//		cutoff  : Eigenvalues of magnitude below this are
	//			  assumed to be zero
	// Return	evlabels: An array of strings containg all the
	//			  operator eigenvalues (in the eigenbasis)
	//			  of operator Op.
	// Note			: Handles complex eigenvalues!
	// Note			: Eigenvalue ordering depends on Op eigenbasis.
	//			  Either order the eigenbasis prior to entering
	//			  function or sort the values after output

  {
  int hs = Op.dim();				// Operator dimension size
  if(!hs) return;				// Exit if NULL operator
   std::string evlabel;				// Label for one eigenvalue

//          Start By Figuring Out Largest, Smallest Eigenvalues, Etc.

  Op.set_EBR();					// Place Op into its eigenbasis
  complex ev;					// This will be an eigenvalue
  double *rev, *iev;				// For re,im eigenvalue components
  rev = new double[hs];
  iev = new double[hs];				// For re,im eigenvalue components
  double rmax=-HUGE_VAL, imax=-HUGE_VAL;	// For maximum (real,mag) eigenvalue cmpts. 
  double rmin=HUGE_VAL, imin=HUGE_VAL, x;	// For minimum (real,mag) eigenvalue cmpts. 
  int ef = 0;
  for(ef=0; ef<hs; ef++)			// Loop over all the eigenvalues,
    if(filter[ef])				// store components and get maximum values
      {
      ev = Op.get(ef,ef);			// 	Get (perhaps complex) eigenvalue of ef 
      x = Re(ev);				// 	Get the real part
      rev[ef] = x;				// 	Store it for later
      x = fabs(x);				// 	Get it's magnitude
      if(x > rmax) rmax = x;			// 	See what the largest real component is
      if(x < rmin) rmin = x;			// 	See what the smallest real component is
      x = Im(ev);				// 	Get the imaginary part
      iev[ef] = x;				// 	Store it for later
      x = fabs(x);				// 	Get its magnitude
      if(x > imax) imax = x;			// 	See what the largest imag component is
      if(x < imin) imin = x;			// 	See what the smallest imag component is
      }

  int prre = 1;					// Print flag for real components
  if(rmax < cutoff) prre = 0;			// See if real eigenvalue components should print
  if(prre)
    {
    while(pow(10.0, prre) < rmax)		// Get the # of printed real digits from 
      prre++;					// the maximum possible real value
    }
//          Generate strings for the Real Components of the Eigenvalues
  
   std::string *revlabels;				// Set the system real eigenvalue components
  revlabels = new  std::string[hs];				// Set the system real eigenvalue components
   std::string frm;
  for(ef=0; ef<hs; ef++)			// Loop over all the eigenvalues,
    {
    if(!prre || !filter[ef])			//	Eigenvalue real components not output
      revlabels[ef] = "";
    else
      {
      x = rev[ef];				// 	Real component printed, this is it
      if(rmin < cutoff) rmin = cutoff;		//	Insure minimum not below cutoff
      if(x >= 0)				//	Start with the sign of the component
        revlabels[ef] = " ";
      else
        revlabels[ef] = "-";
      x = fabs(x);				//	Now just use the absolute value
      if((rmax > 1.e4) || (rmin < 1.e-3))	// 	If allowed numbers can be real big,
        revlabels[ef] += Gform("%8.3e", x); 	// 	or real small then use 'e' format
      else 					//	If they are all in a good range
        {
        frm =  std::string("%")+Gdec(prre)+ std::string(".3f");
        revlabels[ef] += Gform(frm, x); 	// 	then use 'f' format
        }
      }
    }

  int prim = 1;					// Print flag for imag components
  if(imax < cutoff) prim = 0;			// See if imag eigenvalue components should print
  if(prim)
    {
    while(pow(10.0, prim) < imax)		// Get the # of printed imaginary digits from 
      prim++;					// the maximum possible imag value
    }
//          Generate  std::strings for the Imaginary Components of the Eigenvalues

   std::string *ievlabels;				// Set the system imag eigenvalue components
  ievlabels = new  std::string[hs];				// Set the system imag eigenvalue components
  for(ef=0; ef<hs; ef++)			// Loop over all the eigenvalues,
    {
    if(!prim || !filter[ef])			//	Eigenvalue imag components not output
      ievlabels[ef] = "";
    else
      {
      x = iev[ef];				// 	Imag component printed, this is it
      if(imin < cutoff) imin = cutoff;		//	Insure minimum not below cutoff
      if(x >= 0)				//	Start with the sign of the component
        ievlabels[ef] = " ";
      else
        ievlabels[ef] = "-";
      x = fabs(x);				//	Now just use the absolute value
      if((imax > 1.e4) || (imin < 1.e-3))	// 	If allowed numbers can be real big,
        ievlabels[ef] += Gform("%8.3e", x); 	// 	or real small then use 'e' format
      else 					//	If they are all in a good range
        {
        frm =  std::string("%")+Gdec(prre+4)+ std::string(".3f");
        ievlabels[ef] += Gform(frm, x);	 	// 	then use 'f' format
        }
      ievlabels[ef] += "i";			//	Make sure this shows up imaginary
      }
    }

//      Now Assemble The Eigenvalue strings Real & Imaginary Components

  for(ef=0; ef<hs; ef++)			// Loop through all the eigenfunctions
    {
    evlabel = "";				// Start with nothing as the label
    if(filter[ef])
      {
      evlabel += revlabels[ef];			// Start with the real eigenvalue component
      if(prre && prim) evlabel += ", ";		// Add a comma if the eigenvalue is complex
      evlabel += ievlabels[ef];			// Add in the imaginary eigenvalue component
      }
    evlabels[ef] = evlabel;
    }
  delete [] iev;
  delete [] rev;
  delete [] revlabels;
  delete [] ievlabels;
  return;
  }

 
// ____________________________________________________________________________
// C                             TRANSITIONS
// ____________________________________________________________________________
 
/* Transitions between eigenstates are associated with two different viewpoints
 
   1.) It is taken as a jump from one eigenfunction to another eigenfunction,
       so there are two eigenfunctions associated with each transition.
 
   2)  The transition has a value assocated with it, the value being the
       difference between two eigenvalues.  There is only one transition value,
       the net of two eigenvalues.  If operator Op is the system Hamiltonian
       then the transition value is the transition frequency.  If operator Op
       is Fz, then the value is the net Fz associated with the transtion.
 
   Note that transition values in Hilbert space (above) are equivalent to
   single eigenvalues of superoperators in Liouville space.  Similarly,
   transitions between two eigenstates in Hilbert space are associated with
   a single wavefunction in Liouville space.                                 */


void tref_labels( std::string* trlabels, const spin_sys& sys, gen_op &Op,
                                        int type, double cut, int pbf, int pfz)

	// Input 	trlabels: An array of strings
	//       	sys	: A basic spin system
	// 		Op	: General operator (associated to sys)
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//			        INACTIVE IF TYPE = 0
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE = 0
        //		pfz	: Flag for inclusion of total Fz 
	//			        INACTIVE IF TYPE = 0
	// Return	void    : The string array trlabels is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!

  {
  int hs = Op.dim();				// Operator size, spin Hilbert space
  if(!hs) return;				// Exit if Null Operator
  Op.set_EBR();					// Insure Op is in its eigenbasis
  basis B = Op.get_basis();			// Get the eigenbasis matrix
  matrix Bmx = B.U();
  tref_labels(trlabels,sys,Bmx,type,cut,pbf,pfz);	// Use function overload
  return;
  }


void tref_labels( std::string* trlabels, const spin_sys& sys, const matrix &B,
                                     int type, double cutoff, int pbf, int pfz)

	// Input 	trlabels: An array of strings
	//       	sys	: A basic spin system
	// 		B 	: Eigenbasis matrix (associated to sys)
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
	//		cutoff  : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//			        INACTIVE IF TYPE = 0
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE = 0
        //		pfz	: Flag for inclusion of total Fz 
	//			        INACTIVE IF TYPE = 0
	// Return	void    : The string array trlabels is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!
	// Note			: There is no check, but B should be an eigenbasis
	//			  of some operator associated with sys

  {
  int hs = B.rows();				// Basis size, spin Hilbert space
  if(!hs) return;				// Exit if Null Operator
   std::string *wflabels;				// string for wavefunction labels
  wflabels = new  std::string[hs];				// string for wavefunction labels
  int len, lenmax = 0;				// Maximum length of the state labels
  if(type)					// Get the wavefuntion labels if
    { 						// labeling states by wavefunction
    wf_labels(wflabels,sys,B,cutoff,pbf,pfz);	// and find the longest one
    for(int wf=0; wf<hs; wf++)
      {
      len = (wflabels[wf]).length();
      if(len > lenmax) lenmax = len;
      }
    }
  else						// Get the maximum index length if
    { 						// labeling states by number
    lenmax = 1;
    while(pow(10.0, lenmax) < hs)
    lenmax++;
    }

   std::string trlabel;				// string for a transition
   std::string arrow = " --> ";			// Arrow connecting transitions
   std::string istate, fstate;			// Initial and final state strings

//	   Looping Over All Initial States , Then Loop Over All Final States
//                    Build An Expression for Each Transition

  for(int in=0, tr=0; in<hs; in++)		// Loop over the initial states
    {
    if(type)					// Set label for the initial level
      istate = wflabels[in];			// Eithe a wavefunction or number
    else
      istate = Gdec(in);
    len = istate.length();			// Current length of initial level
    while(lenmax-len > 0) 			// Fill with blanks to length lenmax
      {
      istate += " ";
      len++;
      }
    for(int fi=0; fi<hs; fi++, tr++)		// Loop over the final states
      {
      if(type)					// Set label for the final level
        fstate = wflabels[fi];			// Eithe a wavefunction or number
      else
        fstate = Gdec(fi);
      len = fstate.length();			// Current length of final level

//          Now Put Together The string For the Transition Label

      trlabel = "";				// Begin with no label
      trlabel += istate;			// Add the initial state
      trlabel += arrow;				// Add an arrow connecting the states
      trlabel += fstate;			// Add the final state
      trlabels[tr] = trlabel;
      } 					// Do the next final state
    }						// Do the next initial state
  delete [] wflabels;
  return;
  }

// ------------------------ Hilbert Transition Types --------------------------


void tran_types( std::string* trtypes, const spin_sys& sys, gen_op &Op,
                                                          int type, double cut)

	// Input 	trtypes : An array of strings
	//       	sys	: A basic spin system
	// 		Op	: General operator (associated to sys)
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//			        INACTIVE IF TYPE = 0
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE = 0
        //		pfz	: Flag for inclusion of total Fz 
	//			        INACTIVE IF TYPE = 0
	// Return	void    : The string array trtypes is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!

  {
  int hs = Op.dim();				// Operator size, spin Hilbert space
  if(!hs) return;				// Exit if Null Operator
  Op.set_EBR();					// Insure Op is in its eigenbasis
  basis B = Op.get_basis();			// Get the eigenbasis matrix
  matrix Bmx = B.U();				// Get basis in matrix form
  tran_types(trtypes,sys,Bmx,type,cut);		// Use function overload
  return;
  }


void tran_types( std::string* trtypes, const spin_sys& sys, const matrix &B,
                                                       int type, double cutoff)

	// Input 	trtypes: An array of strings
	//       	sys	: A basic spin system
	// 		B 	: Eigenbasis matrix (associated to sys)
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
	// Return	void    : The string array trtypes is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!
	// Note			: There is no check, but B should be an eigenbasis
	//			  of some operator associated with sys

  {
  int hs = B.rows();				// Basis size, spin Hilbert space
  if(!hs) return;				// Exit if Null Operator
  matrix sps = sys.qStates();			// System pdt basis function spin states
  int *pbf;					// Find the largest product basis function
  pbf = new int[hs];					// Find the largest product basis function
  double maxcoeff; 				// contributor to each of the wave functions
  int wf, bf;					// This is equivalent to finding the largest
  for(wf=0; wf<hs; wf++)			// number in each column of the basis B
    {
    maxcoeff = 0;
    pbf[wf] = 0;
    for(bf=0; bf<hs; bf++)
    if(norm(B.get(bf,wf)) > maxcoeff)
      {
      maxcoeff = norm(B.get(bf,wf));
      pbf[wf] = bf;				// pbf[wf] is the index of the product basis
      } 					// function that dominates wavefunction wf
    }

  int ns = sys.spins();				// Get the number of spins in the system
  int* flip;					// This will store indices of flipped spins
  double* trqn;				// This will store delta Iz of flipped spins
  flip = new int[ns];
  trqn = new double[ns];
  int flips, same;
  int cnt, tr=0;
  int i, j, len;
   std::string spin0, spin1;
   std::string* spins;
  spins = new  std::string[ns];
   std::string label;					// Temporary label for the transition
  int wfi, wfj;					// Transiiton Initial & final wavefunctions
  int pbi, pbj;					// Corresponding dominant basis functions
  for(wfi=0; wfi<hs; wfi++)
    {
    pbi = pbf[wfi];				// Largest pdt basis function of wfi
    for(wfj=0; wfj<hs; wfj++)
      {
      pbj = pbf[wfj];				// Largest pdt basis function of wfj
      flips = 0;				// Start with no spin flips in transition
      for(i=0; i<ns; i++)			// Loop through spins and find the 
        if(sps(pbi,i) != sps(pbj,i)) 		// ones that have flipped by comparing
          { 					// the spin states.  If different, store
          flip[flips] = i;			// the spins index as well as the del Iz
          trqn[flips]=Re(sps(pbi,i)-sps(pbj,i));// involved in the flip.  Then keep count
          flips++;				// of the number of spins flipped
          }
      if(!flips)				// If no spins flip, set the label to 0
        trtypes[tr] =  std::string("0");
      else if(flips == 1)			// If one just 1 spin flips, its an SQT
        {
        spin0 = sys.symbol(flip[0]);		//	Get the isotope of the spin flipped
        len = spin0.length();			//	Use only the last letter for the label
        spins[0] = spin0[len-1];
        label = "";				//	Start with no label
        if(trqn[0] < 0) label += "-";		//	Explicitly set the sign of the transition
        if(abs(int(trqn[0])) != 1)		//	If the del Iz is 1, don't write it
          label += Gdec(abs(int(trqn[0])));	//	i.e. write 3WH but not 1WH, thats just WH
        label +=  std::string("W") + spins[0];	//	Now add in the frequency part WX
        trtypes[tr] = label;			//	Set this temporary label to the stored label
        }
      else if(flips == 2)			// If two spins flipped, its either a DQT or an SQT
        {
        spin0 = sys.symbol(flip[0]);		//	Get the isotope of the first spin flipped
        len = spin0.length();			//	Use only the last letter for the label
        spins[0] = spin0[len-1];
        spin1 = sys.symbol(flip[1]);		//	Get the isotope of the second spin flipped
        len = spin1.length();			//	Again use only the last letter in the label
        spins[1] = spin1[len-1];
        label = "";				//	Start with no label
        if(spin0 != spin1)			//	If different spin types, must use both
          {					//	spins in the label.
          if(trqn[0] < 0) label += "-";		//	Have explicit sign on 1st only if negative
          if(abs(int(trqn[0])) != 1)		//	If the del Iz is 1, don't write the value
            label += Gdec(abs(int(trqn[0]))); 	//	i.e. write 3WH but not 1WH, thats just WH
          label +=  std::string("W") + spins[0]; 	//	Now add in the frequency part WX of 1st spin
          if(trqn[1] < 0) label += " - ";	//	For the 2nd spin, use the sign as a separator
          else            label += " + ";
          if(abs(int(trqn[1])) != 1)		//	Still, if del Iz is 1 for this flip, dont
            label += Gdec(abs(int(trqn[1])));	//	bother to write it out.
          label +=  std::string("W") + spins[1];	//	Now add in the frequency part WX of 2nd spin
          trtypes[tr] = label;
          }
        else					//	If two spins of the same type flipped,
          {					//	just add them together.  They may cancel
          trqn[0] += trqn[1];			//	if homonuclear ZQC.  If so just set the
          if(int(trqn[0]) == 0)			//	label to 0
            trtypes[tr] =  std::string("0");
          else 					// 	Or they may add if homonuclear DQC.
            {					//	If so then start the label with a - if
            if(trqn[0] < 0) label += "-";	//	negative.  
            if(abs(int(trqn[0])) != 1) 		//	If the del Iz is 1, don't write it
              label += Gdec(abs(int(trqn[0]))); 	//	i.e. write 3WH but not 1WH, thats just WH
            label +=  std::string("W") + spins[0]; 	//	Now add in the frequency part WX
            trtypes[tr] = label;		//	Set this temporary label to the stored label
            }
          }
        }
      else					// If more that two spins have flipped, this
        { 					// is either a MQT with |M|>2 or a lower QT
        label = "";				// that has combination flips.  Start with no label.
        cnt = 0;				// Keep count of unique spin types involved
        for(i=0; i<flips; i++)			// Loop over all spins that are flipped
          {
          spin0 = sys.symbol(flip[i]);		// This is the isotope type of flipped spin i 
          same = 0;				// Assume its type is different than all other spins
          for(j=0; j<i && !same; j++)		// Loop over all previous flipped spins
            {					// and insure that this one hasn't been
            spin1 = sys.symbol(flip[j]);	// accounted for yet
            if(spin0 == spin1)			// If it has been already accounted for, flag it
              same = 1;
            }
          if(!same)				// This will occur only if spin i is the first
            {					// of its type to be added to the label
            for(j=i+1; j<flips; j++)		// Loop over all next flipped spins
              {					// If they are the same type then just
              spin1 = sys.symbol(flip[j]);	// accounted for them here
              if(spin0 == spin1)
                trqn[i] += trqn[j];
              }
            if(trqn[i] != 0)			// Now add this spin type contribution to
              {					// the label (if it is not zero)
              if(!cnt)
                {
                if(trqn[i] < 0) label += "-";	// For the first type in the label, write
                }				// sign only if negative.  For subsequent
              else				// spin types, use the sign as a spacer
                {
                if(trqn[i] < 0) label += " - ";
                else            label += " + ";
                }
              if(abs(int(trqn[i])) != 1)	// If the spin type has a del Iz of 1
                label += Gdec(abs(int(trqn[i])));// don't bother putting that in the label
              len = spin0.length();		// Now add the WX part, wher X is spin label
              label +=  std::string("W") + spin0[len-1];
              cnt++;
              }
            }
          }
        trtypes[tr] = label;			// Now set the temp label to the stored one
        }
      tr++;
      }
    }
  delete [] pbf;
  delete [] flip;
  delete [] trqn;
  delete [] spins;
  return;
  type = 0;
  cutoff = 0;
  }


// ----------------------- Hilbert Transition Energies ------------------------

 void trev_labels( std::string* trlabels, gen_op& Op, double cutoff)

	// Input 	trlabels: An array of strings
	//       	Op	: General operator
	//		cutoff  : Transitions of magnitude below this are
	//			  assumed to be zero
	// Return	trlabels: An array of strings containg all the
	//			  operator transition values (in the eigenbasis)
	//			  of operator Op.
	// Note			: Handles complex transitions!
	// Note			: Transition ordering depends on Op eigenbasis.
	//			  Either order the Op eigenbasis prior to entering
	//			  function or sort the values after output

  {
  int hs = Op.dim();				// Operator dimension size
  if(!hs) return;				// Exit if NULL operator
  Op.set_EBR();					// Place Op into its eigenbasis

//            Start By Figuring Out All Possible Transitions

  int ls = hs*hs;				// This many possible transitions
  complex *trvals;				// For storing all transitons
  trvals = new complex[ls];
  complex inev, fiev;				// For initial & final eigenvalues
  int tri=0;
  for(int in=0; in<hs; in++)			// Loop over all initial eigenvalues
    {						//
    inev = Op.get(in,in);			// Get (perhaps complex) eigenvalue
    for(int fi=0; fi<hs; fi++, tri++)		// Loop over all final eigenvalues
      {
      fiev = Op.get(fi,fi);			// Get (perhaps complex) eigenvalue
      trvals[tri] = fiev - inev;		// Store the transition value
      }
    }

//     Now Proceed By Figuring Out Largest, Smallest Transitions, Etc.

   std::string trlabel;				// Label for one transiton value
  complex tr;					// For a single transition value
  double *rtr, *itr;			// For re,im transiton components
  rtr = new double[ls];
  itr = new double[ls];
  double rmax=-HUGE_VAL, imax=-HUGE_VAL;		// For maximum (real,mag) transition cmpts. 
  double rmin=HUGE_VAL, imin=HUGE_VAL;			// For minimum (real,mag) transition cmpts. 
  double x,y;					// For temporary storage of components
  for(tri=0; tri<ls; tri++)			// Loop over all the transitions,
    {						// store components and get maximum values
    tr = trvals[tri];				// 	Get (perhaps complex) transition value 
    x = Re(tr);					// 	Get the real part
    y = fabs(x);				//	Get the absolute value
    if(y < cutoff)				//	If magnitude below cutoff then assume
      {						//	it is zero.
      y=0;					//	Set the magnitude to zero
      rtr[tri] = y;				// 	Store the value as zero
      }
    else rtr[tri] = x;				//	Else store the true value
    if(y > rmax) rmax = y;			// 	Get the largest real component is
    if(y<rmin && y!=0) rmin = y;		// 	Get the smallest nonzero real component is

    x = Im(tr);					// 	Get the imaginary part
    y = fabs(x);				//	Get the absolute value
    if(y < cutoff)				//	If magnitude below cutoff then assume
      {						//	it is zero.
      y=0;					//	Set the magnitude to zero
      itr[tri] = y;				// 	Store the value as zero
      }
    else itr[tri] = x;				//	Else store the true value
    if(y > imax) imax = y;			// 	Get largest imag component is
    if(y<imin && y!=0) imin = y;		// 	Get smallest nonzero imag component is
    }

  int prre = 1;					// Print flag for real components
  if(rmax < cutoff) prre = 0;			// See if ther are any non-zero real
  if(prre) 					// transition components.  If so, set flag
    { 						// for printint and then get the # of  
    while(pow(10.0, prre) < rmax)		// printed real digits (above the decimal)
      prre++;					// from the maximum possible real value
    }
//          Generate strings for the Real Components of the Transitions
  
   std::string *rtrlabels;				// Set the system real transition components
  rtrlabels = new  std::string[ls];
   std::string frm = "%";				// Format string
   std::string pt = "f";				// This sets the output format to 'f'
  int len = prre + 4;				// This will be the length of the printed value
  int rlen = 3;					// This many digits past decimal
  if((rmax >= 1.e4) || (rmin < 1.e-3))		// If allowed numbers can be real big,
    { 						// or real small then use 'e' format
    pt = "e";					// for all the values instead of 'f' format
    len = 8;
    }
  for(tri=0; tri<ls; tri++)			// Loop over all the transitions
    {
    if(!prre)					//	Transition real components not output
      rtrlabels[tri] = "";			//	so leave blank
    else
      {
      x = rtr[tri];				// 	Real component printed, this is it
      if(rmin < cutoff) rmin = cutoff;		//	Insure minimum not below cutoff
      if(x >= 0)				//	Start with the sign of the component
        rtrlabels[tri] = " ";
      else
        rtrlabels[tri] = "-";
      x = fabs(x);				//	Now just use the absolute value
      frm =  std::string("%") + Gdec(len)
          +  std::string(".") + Gdec(rlen) + pt;
      rtrlabels[tri] += Gform(frm, x);		// 	Store value in consistent string format
      }
    }

  int prim = 1;					// Print flag for imag components
  if(imax < cutoff) prim = 0;			// See if imag transition components should print
  if(prim)
    {
    while(pow(10.0, prim) < imax)		// Get the # of printed imaginary digits from 
      prim++;					// the maximum possible imag value
    }
//          Generate strings for the Imaginary Components of the Transitions 

   std::string *itrlabels;				// Set the system imag transition components
  itrlabels = new  std::string[ls];
  frm = "%";					// string for print format
  pt = "f";					// Again set the output format to 'f'
  len = prim + 4;				// This will be the length of the printed value
  rlen = 3;					// This many digits past decimal
  if((imax >= 1.e4) || (imin < 1.e-3))		// If allowed numbers can be real big,
    { 						// or real small then use 'e' format
    pt = "e";					// for all the values instead of 'f' format
    len = 8;
    }
  for(tri=0; tri<ls; tri++)			// Loop over all the transitions
    {
    if(!prim)					//	Transition imag components not output
      itrlabels[tri] = "";			//	so leave blank
    else
      {
      x = itr[tri];				// 	Imag component printed, this is it
      if(imin < cutoff) imin = cutoff;		//	Insure minimum not below cutoff
      if(x >= 0)				//	Start with the sign of the component
        itrlabels[tri] = " ";
      else
        itrlabels[tri] = "-";
      x = fabs(x);				//	Now just use the absolute value
      frm =  std::string("%") + Gdec(len)
          +  std::string(".") + Gdec(rlen) + pt;
      itrlabels[tri] += Gform(frm, x);		// 	Store value in consistent string format
      itrlabels[tri] += "i";			//	Make sure this shows up as imaginary
      }
    }

//      Now Assemble The Transition  std::strings Real & Imaginary Components

  for(tri=0; tri<ls; tri++)			// Loop through all the transitions
    {
    trlabel = "";				// Start with nothing as the label
    trlabel += rtrlabels[tri];			// Start with the real transition component
    if(prre && prim) trlabel += ", ";		// Add a comma if the transition is complex
    trlabel += itrlabels[tri];			// Add in the imaginary transition component
    trlabels[tri] = trlabel;
    }
  delete [] rtr;
  delete [] itr;
  delete [] trvals;
  delete [] rtrlabels;
  delete [] itrlabels;
  return;
  }

 
// ____________________________________________________________________________
// D                       EIGENSYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
 
// ------------------ System Wavefunctions in Hilbert Space -------------------
 

void wavefunctions(std::ostream& ostr, const spin_sys& sys, gen_op &Op,
                                    double cutoff, int pbf, int pfz, int title)
 
	// Input	ostr	: An output stream
	//      	sys	: A spin system
	// 		Op	: General operator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
	//              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
	//			: supplied output stream ostr.

    {
    int hs = Op.dim();				// Operator dimension size
    if(!hs) return; 				// Exit if NULL operator
    int *filter;				// Set up a working filter
    filter = new int[hs];
    for(int i=0; i<hs; i++)			// Set the filter to take everything
      filter[i] = 1;
    wavefunctions(ostr, filter, sys, Op,	// Use this function overload
                      cutoff, pbf, pfz, title);
    delete [] filter;
    }


void wavefunctions(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op &Op,
                                    double cutoff, int pbf, int pfz, int title)
 
	// Input	ostr	: An output stream
	//		filt    : An array of integers, output filter
	//      	sys	: A spin system
	// 		Op	: General operator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
	//              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
	//			: supplied output stream ostr.

    {
    int hs = Op.dim();				// Operator dimension size
    if(!hs) return; 				// Exit if NULL operator
     std::string *wflabels;			// Get the wavefunctions
    wflabels = new  std::string[hs];
    wf_labels(wflabels,filt,sys,Op,cutoff,pbf,pfz);
    int len, idig=1, maxlen = 0;
    int wf = 0;
    for(wf=0; wf<hs; wf++) 			// Determine longest wavefunction
      {
      len = (wflabels[wf]).length();
      if(len > maxlen) maxlen = len;
      }
    while(pow(10.0, idig) < hs)		// Get wavefunction printed index size
      idig++;

    len = maxlen;				// Each wavefunction is this long
    len += (idig + 2);				// This long for "#. wavefunction"
    int ncols = 1;				// Can we put multiple wavefunctions on
    int lentot = len;	 			// the same line & span < 80 columns
    while(lentot+4+len < 80)			// using 4 spaces to separate columns?
      {
      lentot += (4+len);
      ncols ++;
      }

    len = (80 - lentot)/2;			// Amount to center at 80 columns
    if(len < 0)					// If the printed functions are
      len = 0;					// then don't bother centering
     std::string tenblanks = "          ";		// Use this to add spaces in output
     std::string center = "";				// Construct blanks to center wf's
     std::string index;				// Use this for wavefunction index
    while(len > 10)
      {
      center += tenblanks;
      len -= 10;
      }
    if(len > 0) center += std::string(len, ' ');
    if(title)					// Now write a title out if desired
      {
      ostr << "\n"; 				//	Add an additional blank line
      if(Op.in_EBR())				//	If eigenbasis, these will
        {					//	be eigenfunctions
        ostr << tenblanks << tenblanks << tenblanks
             << "System Eigenfunctions";
        }
      else if((Op.get_basis()).isDefaultBasis())//	If the default basis, these will
        ostr << tenblanks << tenblanks << "  "	//	be the GAMMA default basis functions
             << "System GAMMA Default"
             << " Basis Functions";
      else					//	Else, these are some other basis
        ostr << tenblanks << tenblanks		//	wavefunctions
             << tenblanks
             << "System WaveFunctions";
      }
    ostr << "\n";				// Output a preceding blank line
    int col = 0;				// Count the columns as we print
    for(wf=0; wf<hs; wf++)			// Now write out the wavefunctions
      {
      index = Gdec(wf); 				//	Get the wavefunction index
      len = index.length();			//	This is the index length
      if(filt[wf])				// 	Output only if not filtered
        {
        if(col == 0)				//	For 1st column, new line &
          ostr << "\n" << center;		// 	space over so centered
        else
          ostr << "    ";			//	Else 4 spaces to next column
        ostr << index << ". "			//	Write the wavefunciton index
             << std::string(idig-len, ' ')	//	Fill so all indices same width
             << wflabels[wf];			//	Write the wavefunction
        col++;					//	The next column has this index
        if(col == ncols)			//	If we've surpassed the # of cols
          col = 0;				//	we want, then just start over
        }
      }
    ostr << "\n";				// Finish with a blank line
    delete [] wflabels;
    return;
    }
 

// ---------------------- Eigensystem in Hilbert Space ------------------------


void eigensystem(std::ostream& ostr, const spin_sys& sys, gen_op& Op,
                         double cute, double cutc, int pbf, int pfz, int title)

	// Input	ostr	: An output stream
	//      	sys	: A spin system
	//       	Op	: A general operator
        //              cute 	: Eigenvalues below this are assumed 0
        //              cutc 	: Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
	//              title   : Flag for inclusion of a title
	// Return	none 	: All eigenvalues and eigenfunctions
	//			  of operator Op are sent into ostream ostr

    {
    int hs = Op.dim();				// Operator dimension size
    if(!hs) return; 				// Exit if NULL operator
    int *filter;				// Set up a working filter
    filter = new int[hs];
    for(int i=0; i<hs; i++)			// Set the filter to take everything
    filter[i] = 1;
    eigensystem(ostr,filter,sys,Op,
                 cute, cutc, pbf, pfz,title);
    delete [] filter;
    return;
    }


void eigensystem(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op& Op,
                       double cute, double cutc, int pbf, int pfz, int title)

	// Input	ostr	: An output stream
	//		filt    : An array of integers, output filter
	//      	sys	: A spin system
	//       	Op	: A general operator
        //              cute 	: Eigenvalues below this are assumed 0
        //              cutc 	: Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
	//              title   : Flag for inclusion of a title
	// Return	none 	: All eigenvalues and eigenfunctions
	//			  of operator Op are sent into ostream ostr

    {
    int hs = Op.dim();				// Operator dimension size
    if(!hs) return; 				// Exit if NULL operator
    Op.set_EBR();				// Insure Op into its eigenbasis

    std::string *wflabels;			// Get the wavefunctions
    wflabels = new std::string[hs];
    wf_labels(wflabels,filt,sys,Op,cutc,pbf,pfz);
    int len, maxlen = 0;
    for(int wf=0; wf<hs; wf++) 			// Determine longest wavefunction
      {
      len = (wflabels[wf]).length();
      if(len > maxlen) maxlen = len;
      }

    int idig=1;					// Use this for eigenstate indexing
    while(pow(10.0, idig) < hs)		// Get eigenstate printed index size
      idig++;

    std::string *evlabels;			// Get the system eigenvalues
    evlabels = new std::string[hs];
    ev_labels(evlabels, Op, cute);

    if(title)					// Now write a title out if desired
      ostr << "\n               " 		// Will be centered on column 80
           << "Current Eigensystem: "
           << " Eigenvalues and Eigenfunctions"
           << "\n";

    len = (evlabels[0]).length();		// Each eigenvalue is this long
    len += (2 + maxlen);			// This long for "eigenvalue  eigenfunction"
    len += (idig + 2);				// This long for "#. eigenvalue  eigenfunction"
    int ncols = 1;				// Can we put multiple eigenstates on
    int lentot = len;	 			// the same line & span < 80 columns
    while(lentot+4+len < 80)			// using 4 spaces to separate columns?
      {
      lentot += (4+len);
      ncols ++;
      }

    len = (80 - lentot)/2;			// Amount to center each line at 80 columns
    if(len < 0)					// If the printed eigenstates are long
      len = 0;					// then don't bother centering
    std::string tenblanks = "          ";		// Use this to add spaces in output
    std::string center = "";				// Construct blanks to center wf's
    while(len > 10)				// At the loop end, string center will
      {						// be the spaces to move each eigenstate
      center += tenblanks;			// line over so it looks centered
      len -= 10;
      }
    if(len > 0) center += std::string(len, ' ');

    std::string index;				// Use this for wavefunction index
    int col = 0;				// Count the columns as we print
    for(int es=0; es<hs; es++)			// Loop through all the eigenstates
      {
      if(filt[es])
        {
        index = Gdec(es); 			//	Get the eigenstate index
        len = index.length();			//	This is the index length
        if(col == 0)				//	For 1st column, new line &
          ostr << "\n" << center;		// 	space over so centered
        else
          ostr << "    ";			//	Else 4 spaces to next column
        ostr << index << ". "			//	Write the eigenstate index
             << std::string(idig-len, ' ')		//	Fill so all indices same width
             << evlabels[es] << "  "		//	Write the eigenvalue
             << wflabels[es];			//	Write the wavefunction
        col++;					//	The next column has this index
        if(col == ncols)			//	If we've surpassed the # of cols
          col = 0;				//	we want, then just start over
        }
      }
    ostr << "\n";				//	End with a blank line
    delete [] wflabels;
    delete [] evlabels;
    return;
    }

 
// --------------------- System Transitions in Hilbert Space ------------------


void transitions(std::ostream& ostr, const spin_sys& sys, gen_op& Op, int type,
                         double cutt, double cutc, int pbf, int pfz, int title)

	// Input	ostr	: An output stream
	//      	sys	: A spin system
	//       	Op	: A general operator
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
        //              cutt 	: Transitions below this are assumed 0
        //              cutc 	: Basis functions with coefficients below
        //                        this value are not taken as part of a wf
	//			        INACTIVE IF TYPE = 0
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE = 0
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
	//			        INACTIVE IF TYPE = 0
	//              title   : Flag for inclusion of a title
	// Return	none 	: All transitions and transition values
	//			  of operator Op are sent into ostream ostr

  {
  int hs = Op.dim();				// Operator dimension size
  if(!hs) return; 				// Exit if NULL operator
  const int ls = hs*hs;				// This is the Liouville space
  int *filter;				// Set up a working filter
  filter = new int[ls];
  for(int i=0; i<ls; i++)			// Set the filter to take everything
  filter[i] = 1;
  transitions(ostr, filter, sys, Op, type,
                   cutt, cutc, pbf, pfz, title);
  delete [] filter;
  }


void transitions(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op& Op,
               int type, double cutt, double cutc, int pbf, int pfz, int title)

	// Input	ostr	: An output stream
	//		filt    : An array of integers, output filter
	//      	sys	: A spin system
	//       	Op	: A general operator
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
        //              cutt 	: Transitions below this are assumed 0
        //              cutc 	: Basis functions with coefficients below
        //                        this value are not taken as part of a wf
	//			        INACTIVE IF TYPE = 0
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE = 0
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
	//			        INACTIVE IF TYPE = 0
	//              title   : Flag for inclusion of a title
	// Return	none 	: All transitions and transition values
	//			  of operator Op are sent into ostream ostr

    {
    int hs = Op.dim();				// Operator dimension size
    if(!hs) return; 				// Exit if NULL operator
    Op.set_EBR();				// Insure Op into its eigenbasis
    int ls = hs*hs;				// Get the Liouville space

    std::string *treflabels;			// Get the transitions
    treflabels = new std::string[ls];			// Get the transitions
    tref_labels(treflabels, sys, Op,
                        type, cutc, pbf, pfz);
    int len, maxlen = 0;
    int tr = 0;
    for(tr=0; tr<ls; tr++) 			// Determine longest transition
      { 					// wavefunction labeling
      if(filt[tr])
        {
        len = (treflabels[tr]).length();
        if(len > maxlen) maxlen = len;
        }
      }

    int idig=1;					// Use this for eigenstate indexing
    while(pow(10.0, idig) < ls)		// Get eigenstate printed index size
      idig++;

    std::string *trevlabels;			// Get the system transition values
    trevlabels = new std::string[ls];	// Get the system transition values
    trev_labels(trevlabels, Op, cutt);

    if(title)					// Now write a title out if desired
      ostr << "\n           "	 		// Will be centered on column 80
           << "Current Eigensystem: "
           << " Transition Values and Eigenstate Jumps"
           << "\n";

// sosi - Alternatively, if each transition label is too long for one line, may want
//        to parse it before the -->.  Then center on the first part + 8, insert
//        \n+center+8spaces before the --> so it prints on two lines nicely?
//	  I'll deal with that when 4 spin systems are run through this I'm sure

    len = (trevlabels[0]).length();		// Each transition value is this long
    len += (2 + maxlen);			// This long for "value  label"
    len += (idig + 2);				// This long for "#. value  label"
    int ncols = 1;				// Can we put multiple transitions on
    int lentot = len;	 			// the same line & span < 80 columns
    while(lentot+4+len < 80)			// using 4 spaces to separate columns?
      {
      lentot += (4+len);
      ncols ++;
      }
    len = (80 - lentot)/2;			// Amount to center each line at 80 columns
    if(len < 0)					// If the printed transitions are long
      len = 0;					// then don't bother centering
    std::string tenblanks = "          ";		// Use this to add spaces in output
    std::string center = "";				// Construct blanks to center transitions
    while(len > 10)				// At the loop end, string center will
      {						// be the spaces to move each transition
      center += tenblanks;			// line over so it looks centered
      len -= 10;
      }
    if(len > 0) center += std::string(len, ' ');

    std::string index;				// Use this for transition index
    int col = 0;				// Count the columns as we print
    for(tr=0; tr<ls; tr++)			// Loop through all the transitiions
      {
      if(filt[tr])
        {
        index = Gdec(tr); 			//	Get the transition index
        len = index.length();			//	This is the index length
        if(col == 0)				//	For 1st column, new line &
          ostr << "\n" << center;		// 	space over so centered
        else
          ostr << "    ";			//	Else 4 spaces to next column
        ostr << index << ". "			//	Write the transition index
             << std::string(idig-len, ' ')	//	Fill so all indices same width
             << trevlabels[tr] << "  "		//	Write the transition value 
             << treflabels[tr];			//	Write the transition eigenstates
        col++;					//	The next column has this index
        if(col == ncols)			//	If we've surpassed the # of cols
          col = 0;				//	we want, then just start over
        }
      }
    ostr << "\n";				//	End with a blank line
    delete [] trevlabels;
    return;
    }
 

// ____________________________________________________________________________
// E                             SELECTION RULES
// ____________________________________________________________________________
 
// -------------------- Selection Rules in Hilbert Space ----------------------

 
void ev_select(int* select, gen_op &Op, double val1, int type,
                                          double val2, double cutoff, int reim)

	// Input 	select  : An array of integers
	// 		Op	: General operator (associated to sys)
	//		val1    : Value used in determining the selection
        //              type    : Type of selection applied
        //                         0 - Only levels having ev == val1 (DEFAULT)
        //                         1 - All levels having ev != val1
        //                         2 - All levels having ev >= val1
        //                         3 - All levels having ev <= val1
        //                         4 - All levels having val2 <= ev <= val1
	//		val2    : Value used in determining the selection,
	//			  unused unless type=4, (0 DEFAULT) 
	//		cutoff  : If fabs(ev-val1) < cutoff, assume equal
        //              reim    : Type of ev used
        //                         0 - Use the norm of the complex ev
        //                        >0 - Use the real component of the ev (DEFAULT)
        //                        <0 - Use the imaginary component of ev
	// Return	void    : The integer array is filled with 0's and 1's
	//			  depending on whether the system sys wavefunctions
	//			  satisfy the input selection rules.
	// Note			: Selection ordering depends on Op input basis!
	// Note			: Selection is based on eigenvalues here

    {
    int hs = Op.dim();				// Operator dimension
    if(!hs) return;				// Exit if the basis is empty
    complex zev;				// For complex eigenvalue
    double rev, valtmp;				// For real eigenvalue measure
    int TF;					// True - False for selection

//	                   Looping Over All Eigenvalues of Op

    for(int ev=0; ev<hs; ev++)
      {
      TF = 0;					// Assume this value is no good
      zev = Op.get(ev,ev);			// Get the, perhaps complex, eigenvalue
      if(reim == 0)				// See if the complex magnitude is used
        rev = norm(zev);
      else if(reim > 0)				// Else use the real value
        rev = Re(zev);
      else					// Else use the imaginary value
        rev = Im(zev);
      switch(type)				// See if rev satisfies selection rules
        {
        case 0:					// The value must be val1
        default:
          if(fabs(rev-val1) < cutoff)
            TF = 1;
          break;
        case 1:					// The value must not be val1
          if(fabs(rev-val1) > cutoff)
            TF = 1;
          break;
        case 2:					// The value must >= val1
          if((rev >= val1)
                  || (fabs(rev-val1) <= cutoff))
            TF = 1;
          break;
        case 3:					// The value must <= val1
          if((rev <= val1)
                  || (fabs(rev-val1) <= cutoff))
            TF = 1;
          break;
        case 4:					// The value must <= val1 & >= val2
          if(val1 < val2)
            {
            valtmp = val2;
            val2 = val1;
            val1 = valtmp;
            }
          if((rev <= val1)
                  || (fabs(rev-val1) <= cutoff))
            if((rev >= val2)
                  || (fabs(rev-val2) <= cutoff))
              TF = 1;
          break;
        }
      select[ev] = TF;				// Set the selection for this state
      } 					// Do the next eigenvalue
    return;
    }

 
void tr_select(int* select, gen_op &Op, double val1, int type,
                                          double val2, double cutoff, int reim)

	// Input 	select  : An array of integers
	// 		Op	: General operator (associated to sys)
	//		val1    : Value used in determining the selection
        //              type    : Type of selection applied
        //                         0 - Only transitions having dev == val1 (DEFAULT)
        //                         1 - All transitions having dev != val1
        //                         2 - All transitions having dev >= val1
        //                         3 - All transitions having dev <= val1
        //                         4 - All transitions having val2 <= dev <= val1
	//		val2    : Value used in determining the selection,
	//			  unused unless type=4, (0 DEFAULT) 
	//		cutoff  : If fabs(dev-val1) < cutoff, assume equal
        //              reim    : Type of dev used
        //                         0 - Use the norm of the complex dev
        //                        >0 - Use the real component of the dev (DEFAULT)
        //                        <0 - Use the imaginary component of dev
	// Return	void    : The integer array is filled with 0's and 1's
	//			  depending on whether the system sys wavefunctions
	//			  satisfy the input selection rules.
	// Note			: Selection ordering depends on Op input basis!

    {
    int hs = Op.dim();				// Operator dimension
    if(!hs) return;				// Exit if the basis is empty
    complex zdev, zevi;				// For complex transition value
    double rdev, valtmp;			// For real transition value measure
    int TF;					// True - False for selection

//	                   Looping Over All Eigenvalues of Op

    for(int evi=0,tr=0; evi<hs; evi++)		// Loop over the initial states
      {
      zevi = Op.get(evi,evi);			// Get the, perhaps complex, eigenvalue
      for(int evf=0; evf<hs; evf++,tr++)	// Loop over the final states
        {
        TF = 0;					// Assume this value is no good
        zdev = Op.get(evf,evf) - zevi;		// Get the, perhaps complex, transition value
        if(reim == 0)				// See if the complex magnitude is used
          rdev = norm(zdev);
        else if(reim > 0)			// Else use the real value
          rdev = Re(zdev);
        else					// Else use the imaginary value
          rdev = Im(zdev);
        switch(type)				// See if rdev satisfies selection rules
          {
          case 0:				// The value must be val1
          default:
            if(fabs(rdev-val1) < cutoff)
              TF = 1;
            break;
          case 1:				// The value must not be val1
            if(fabs(rdev-val1) > cutoff)
              TF = 1;
            break;
          case 2:				// The value must >= val1
            if((rdev >= val1)
                    || (fabs(rdev-val1) <= cutoff))
              TF = 1;
            break;
          case 3:				// The value must <= val1
            if((rdev <= val1)
                    || (fabs(rdev-val1) <= cutoff))
              TF = 1;
            break;
          case 4:				// The value must <= val1 & >= val2
            if(val1 < val2)
              {
              valtmp = val2;
              val2 = val1;
              val1 = valtmp;
              }
            if((rdev <= val1)
                    || (fabs(rdev-val1) <= cutoff))
              if((rdev >= val2)
                    || (fabs(rdev-val2) <= cutoff))
                TF = 1;
            break;
          }
        select[tr] = TF;			// Set the selection for this transition
        } 					// Do the next final state
      } 					// Do the next initial state
    return;
    }

 
#endif 							// HSanalyze.cc
