/* LSanalyze.cc ************************************************-*-c++-*-*
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      MR Liouville Space Analyze                   Implementation     **
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
** of Spin Liouville Spaces Associated with Composite Spin Systems.	**
** This module of the GAMMA MR Library provides functions designed to   **
** facilitate the analysis of the more common aspects of MR with        **
** respect to a quantum mechanical viewpoint.  As GAMMA provides full   **
** access to all its defined quantities, these routines serve just to   **
** sort information and present it in an orderly fashion.               **
**                                                                      **
*************************************************************************/

#ifndef   LSanalyze_cc_			// Is this file already included?
#  define LSanalyze_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <LSLib/LSanalyze.h>		// Include the header
#include <HSLib/HSanalyze.h>		// Include Hilbert space functions
#include <HSLib/SpinOpCmp.h>		// Include composite spin ops
#include <string>			// Include libstdc++ strings
#include <Basics/StringCut.h>
#include <cmath>			// Include HUGE_VAL

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output

// ____________________________________________________________________________
// A                     WAVEFUNCTIONS / BASIS FUNCTIONS
// ____________________________________________________________________________

// The spin labels below are defined in HSLib/HSanalyze.cc
//string alphabeta[7] = {"a","b","g","d","e","w","x"};	// Spin state labels

/* The wavefunctions will span the Liouville space of the superoperator LOp.
   Each wavefunction will be some linear combination of the default super-
   operator basis functions.  The default superoperator basis functions are
   inherent to each superoperator.  They are determined directly from the 
   differences of Hilbert space wavefunctions, and these Hilbert space wave-
   functions are contained in the Hilbert space basis internally present in 
   LOp.  Each of these Hilbert space wavefunctions are related to the GAMMA
   default basis functions (inherent to spin system sys).  If LOp is in its 
   eigenbasis, its wavefunctions are eigenfunctions.  If the LOp Hilbert 
   space basis an eigenbasis, then the LOp wavefunctions are linear com-
   binations of Hilbert space eigenfunction differences.  Each wavefunction 
   corresponds to a coherence between two state of the system.               */


void wf_labels(string* wflabels, const spin_sys& sys, super_op &LOp,
                                      double cut, int type, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//       	sys	: A basic spin system
	// 		LOp	: A superoperator (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
	//			   <0 - by Liouville number, [0, ls) (DEFAULT)
	//			    0 - Hilbert transitions, {1->2, 2->3,..}
	//			   >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
	//		pbf	: Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE <= 0
        //		pfz	: Flag for inclusion of total Fz in Hilbert basis
	//			        INACTIVE IF TYPE <= 0
	// Return	wflabels: An array of strings containg all system
	//			  wavefunctions in the current working basis
	//			  of operator LOp in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on LOp input basis!

    {
    int ls = LOp.size();			// Superperator dimension size
    if(!ls) return;				// Exit if the superop is empty
    basis LB = LOp.get_Lbasis();		// Get the Liouville basis
    basis HB = LOp.get_basis();			// Get the Hilbert basis
    int *filter;				// Set up a working filter
    filter = new int[ls];
    for(int i=0; i<ls; i++)			// Set the filter to take everything
      filter[i] = 1;
    matrix LBmx = LB.U();
    matrix HBmx = HB.U();
    wf_labels(wflabels,filter,sys,LBmx,HBmx,	// Use function overload
                             cut,type,pbf,pfz);
    delete [] filter;
    return;
    }


void wf_labels(string* wflabels, int* filter, const spin_sys& sys, super_op &LOp,
                                            double cut, int type, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//		filter  : An array of integers, output filter
	//       	sys	: A basic spin system
	// 		LOp	: A superoperator (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
	//			   <0 - by Liouville number, [0, ls) (DEFAULT)
	//			    0 - Hilbert transitions, {1->2, 2->3,..}
	//			   >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
	//		pbf	: Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE <= 0
        //		pfz	: Flag for inclusion of total Fz in Hilbert basis
	//			        INACTIVE IF TYPE <= 0
	// Return	wflabels: An array of strings containg all system
	//			  wavefunctions in the current working basis
	//			  of operator LOp in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on LOp input basis!

    {
    int ls = LOp.size();			// Superperator dimension size
    if(!ls) return;				// Exit if the superop is empty
    basis LB = LOp.get_Lbasis();		// Get the LOp Liouville basis
    basis HB = LOp.get_basis();			// Get the LOp Hilbert basis
    matrix LBmx = LB.U();
    matrix HBmx = HB.U();
    wf_labels(wflabels,filter,sys,LBmx,HBmx,	// Use function overload
                             cut,type,pbf,pfz);
    return;
    }


void wf_labels(string* wflabels, const spin_sys& sys, const matrix &B,
                    const matrix& HB, double cut, int type, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//       	sys	: A basic spin system
	// 		LB	: A superoperator basis
	//		HB	: A Hilbert basis (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
	//			   <0 - by Liouville number, [0, ls) (DEFAULT)
	//			    0 - Hilbert transitions, {1->2, 2->3,..}
	//			   >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
	//		pbf	: Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE <= 0
        //		pfz	: Flag for inclusion of total Fz in Hilbert basis
	//			        INACTIVE IF TYPE <= 0
	// Return	wflabels: An array of strings containg all system
	//			  wavefunctions in the current working basis
	//			  of operator LOp in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on LOp input basis!

  {
  int ls = B.cols();				// Superperator dimension size
  if(!ls) return;
  int *filter;				// Set up a working filter
  filter = new int[ls];
  for(int i=0; i<ls; i++)			// Set the filter to take everything
    filter[i] = 1;
  wf_labels(wflabels,filter,sys,B,HB,		// Use function overload
                             cut,type,pbf,pfz);
  delete [] filter;
  }


void wf_labels(string* wflabels, int* filter, const spin_sys& sys, const matrix &B,
                        const matrix& HB, double cut, int type, int pbf, int pfz)

	// Input 	wflabels: An array of strings
	//		filter  : An array of integers, output filter
	//       	sys	: A basic spin system
	// 		LB	: A superoperator basis
	//		HB	: A Hilbert basis (associated to sys)
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
	//			   <0 - by Liouville number, [0, ls) (DEFAULT)
	//			    0 - Hilbert transitions, {1->2, 2->3,..}
	//			   >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
	//		pbf	: Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE <= 0
        //		pfz	: Flag for inclusion of total Fz in Hilbert basis
	//			        INACTIVE IF TYPE <= 0
	// Return	wflabels: An array of strings containg all system
	//			  wavefunctions in the current working basis
	//			  of operator LOp in the format specified
	// Note			: Handles complex basis function coefficients!
	// Note			: Wavefunction ordering depends on LOp input basis!

  {
  int ls = B.cols();				// Superperator dimension size
  if(!ls) return;
  string wflabel;

  string *trlabels;				// Allocate array for transitions
  trlabels = new string[ls];
  if(type >=0) 					// Get the transitions if they are
    tref_labels(trlabels, sys, HB, type,	// to be used in labeling the wfs
                                cut, pbf, pfz);

  string rlab, ilab;				// For basis function label coefficient(s)
  string startbf = "|";				// For basis function label start
  string bfind;					// For basis function label index
  string endbf = ">";				// For basis funciton label end
  string bflabel = "";				// For basis funciton label (total)
  string blanks = "               ";		// For any spaces needed in labeling 

  complex zcoeff;				// For basis function coefficient value
  double rcoeff, icoeff;			// For basis function coefficient real & imag
  int bfindl;					// For basis function coefficient index length

  int idig = 1;					// Basis function index printed digits
  while(pow(10.0,double(idig)) < ls)		// Get the index printed digits from 
    idig++;					// the maximum possible index, ls
  int bfout = 0;				// How many basis functions in wf

//	   Looping Over All WaveFunctions, Then Loop Over All Basis
//                Functions and Build Each Wavefunction Label

  int keep;
  for(int wf=0; wf<ls; wf++)			// Loop over all the wavefunctions
    {
    wflabel = "";				// Start with an empty label
    keep = filter[wf];				// Flag whether to make this one or not
    bfout = 0;					// Set the output basis functions to 0
    for(int bf=0; keep && bf<ls; bf++)		// Loop over each bf of wavefunction wf
      {
      zcoeff = Re(B.get(bf,wf));		// The coefficient of basis function bf
      if(norm(zcoeff) >= cut) 			// Don't do small coefficients, usually [-1,1]
        {
//			Set Up Real Coefficient of Basis Function bf

        rcoeff = Re(B.get(bf,wf));		// Real coefficient of basis function bf
        rlab = "";				// Start with empty bf real coefficient label
        if(fabs(rcoeff) >= cut)			// Don't do small real coefficients, usually [-1,1]
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

        icoeff = Im(B.get(bf,wf));		// Imag coefficient of basis function bf
        ilab = "";				// Start with empty bf imaginary coefficient label
        if(fabs(icoeff) >= cut)			// Don't do small coefficients, usually [-1,1]
          {						
          if(bfout==0 && (fabs(rcoeff)<cut))	//	Set the sign of the imag coeff
            if(icoeff >= 0) ilab = " ";
            else            ilab = "-";
          else
            if(icoeff >=0)  ilab = " + ";
            else            ilab = " - ";
          ilab += Gform("%5.3f", fabs(icoeff));	//	Now add the value of the imag coeff
          ilab += "i";				//	Indicate this is imaginary
          }
//		                Set Up Label of Basis Function bf

        bflabel = startbf;			// Start the basis function bf label with |
        if(type <0) 				// Label bf with an index number only
          {
          bfind = Gdec(bf);			// Start bf index label
          bfindl = bfind.length();		// Get the current length
          bfind += string(idig-bfindl, ' ');	// Fill out the index to idig length
          bflabel += bfind;
          }
        else 					// Label bf with transition
          {
          if(type > 0)				// If long version, enclose with
            { 					// brackets
            bflabel += "[";
            bflabel += trlabels[bf];
            bflabel += "]";
            }
          else					// If short version (e.g. tr1 --> tr2)
            bflabel += trlabels[bf];		// just output it
          }
        bflabel += endbf;			// End the basis function bf label with >

//            Add the Contribution of Basis Function bf to Wavefuntion wf's Label

        wflabel += rlab;			// Add the real coefficient of bf
        wflabel += ilab;			// Add the imag coefficient of bf
        wflabel += bflabel;			// Add the basis funciton label of bf
        bfout++;				// Track how may basis functions included
        }
      } 					// Do the next bf of this wavefunction
    wflabels[wf] = wflabel;
    } 						// Do the next wavefunction
  delete [] trlabels;
  return;
  }


// ____________________________________________________________________________
// B                                EIGENVALUES
// ____________________________________________________________________________


/* The eigenvalues will span the Liouville space of the superoperator LOp.
   Each eigenvalue will correspond to an element of LOp after LOp has been
   placed into it representation in an eigenbasis.  Each eigenvalue will then
   correspond to a particular eigenfunction of the superoperator.  If the 
   superoperator is a commutation superoperator of the Hamiltonian, then the
   eigenvalues will be transition energies of the system.                    */


void ev_labels(string* evlabels, super_op& LOp, double cutoff)

	// Input 	evlabels: An array of strings
	//       	LOp	: A superoperator
	//		cutoff  : Eigenvalues of magnitude below this are
	//			  assumed to be zero
	// Return	evlabels: An array of strings containg all the
	//			  operator eigenvalues (in the eigenbasis)
	//			  of operator LOp.
	// Note			: Handles complex eigenvalues!
	// Note			: Eigenvalue ordering depends on LOp eigenbasis.
	//			  Either order the eigenbasis prior to entering
	//			  function or sort the values after output
  {
  int ls = LOp.size();				// Superoperator dimension size
  if(!ls) return; 				// Exit if NULL superoperator
  int *filter;				// Set up a working filter
  filter = new int[ls];
  for(int i=0; i<ls; i++)			// Set the filter to take everything
    filter[i] = 1;
  ev_labels(evlabels, filter, LOp, cutoff);
  delete [] filter;
  }


void ev_labels(string* evlabels, int* filter, super_op& LOp, double cutoff)

	// Input 	evlabels: An array of strings
	//		filter  : Array of integers to filter values
	//       	LOp	: A superoperator
	//		cutoff  : Eigenvalues of magnitude below this are
	//			  assumed to be zero
	// Return	evlabels: An array of strings containg all the
	//			  operator eigenvalues (in the eigenbasis)
	//			  of operator LOp.
	// Note			: Handles complex eigenvalues!
	// Note			: Eigenvalue ordering depends on LOp eigenbasis.
	//			  Either order the eigenbasis prior to entering
	//			  function or sort the values after output

  {
  int ls = LOp.size();				// Superoperator dimension size
  if(!ls) return; 				// Exit if NULL superoperator
  string evlabel;				// string for one eigenvalue

//          Start By Figuring Out Largest, Smallest Eigenvalues, Etc.

  LOp.set_EBR();				// Place LOp into its eigenbasis
  complex ev;					// This will be an eigenvalue
  double *rev, *iev;			// For re,im eigenvalue components
  rev = new double[ls];
  iev = new double[ls];
  double rmax=-HUGE;
  double imax=-HUGE;		// For maximum (real,mag) eigenvalue cmpts. 
  double rmin=HUGE;
  double imin=HUGE;			// For minimum (real,mag) eigenvalue cmpts. 
  double x,y;					// For temporary eigenvalue storage
  int ef = 0;
  for(ef=0; ef<ls; ef++)			// Loop over all the eigenvalues,
    if(filter[ef]) 				// store components and get maximum values
      {
      ev = LOp.get(ef,ef);			// 	Get (perhaps complex) eigenvalue of ef 
      x = Re(ev);				// 	Get the real part
      y = fabs(x);
      if(y < cutoff)                            //      If magnitude below cutoff then assume
        {                                       //      it is zero.
        y=0;                                    //      Set the magnitude to zero
        rev[ef] = y;				//      Store the value as zero
        }
      else rev[ef] = x;				// 	Else store the true value
      if(y > rmax) rmax = y;			// 	See what the largest real component is
      if(y < rmin && y!=0) rmin = y;		// 	See what the smallest real component is

      x = Im(ev);				// 	Get the imaginary part
      y = fabs(x);
      if(y < cutoff)                            //      If magnitude below cutoff then assume
        {                                       //      it is zero.
        y=0;                                    //      Set the magnitude to zero
        iev[ef] = y;				//      Store the value as zero
        }
      else iev[ef] = x;				// 	Else store the true value
      if(y > imax) imax = y;			// 	See what the largest imag component is
      if(y < imin && y!=0) imin = y;		// 	See what the smallest imag component is
      }

  int prre = 1;					// Print flag for real components
  if(rmax < cutoff) prre = 0;			// See if real eigenvalue components should print
  if(prre)
    {
    while(pow(10.0,double(prre)) < rmax)	// Get the # of printed real digits from 
      prre++;					// the maximum possible real value
    }
//          Generate strings for the Real Components of the Eigenvalues
  
  string *revlabels;				// Set the system real eigenvalue components
  revlabels = new string[ls];
  string frm;					// Printing format
  string pt = "f";				// This sets the output format to 'f'
  int len = prre + 4;                           // This will be the length of the printed value
  int rlen = 3;                                 // This many digits past decimal
  if((rmax >= 1.e4) || (rmin < 1.e-3))          // If allowed numbers can be real big,
    {                                           // or real small then use 'e' format
    pt = string("e");				// for all the values instead of 'f' format
    len = 8;
    }
  for(ef=0; ef<ls; ef++)			// Loop over all the eigenvalues,
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
      frm = string("%") + Gdec(len)
          + string(".") + Gdec(rlen) + pt;
      revlabels[ef] += Gform(frm.c_str(), x);	// 	This is the real cmpt as a string
      }
    }

  int prim = 1;					// Print flag for imag components
  if(imax < cutoff) prim = 0;			// See if imag eigenvalue components should print
  if(prim)
    {
    while(pow(10.0,double(prim)) < imax)	// Get the # of printed imaginary digits from 
      prim++;					// the maximum possible imag value
    }
//          Generate strings for the Imaginary Components of the Eigenvalues

  string *ievlabels;				// Set the system imag eigenvalue components
  ievlabels = new string[ls];
  pt = string("f");				// Again set the output format to 'f'
  len = prim + 4;                               // This will be the length of the printed value
  rlen = 3;                                     // This many digits past decimal
  if((imax >= 1.e4) || (imin < 1.e-3))          // If allowed numbers can be real big,
    {                                           // or real small then use 'e' format
    pt = string("e");				// for all the values instead of 'f' format
    len = 8;
    }
  for(ef=0; ef<ls; ef++)			// Loop over all the eigenvalues,
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
      frm = string("%") + string(Gdec(prre+4))
          + string(".3f");
      ievlabels[ef] += Gform(frm.c_str(), x);	// 	This is the imag component as a string
      ievlabels[ef] += "i";			//	Make sure this shows up imaginary
      }
    }

//      Now Assemble The Eigenvalue strings Real & Imaginary Components

  for(ef=0; ef<ls; ef++)			// Loop through all the eigenfunctions
    {
    evlabel = "";				// Start with no label
    if(filter[ef])
      {
      evlabel += revlabels[ef];			// Add in the real eigenvalue component
      if(prre && prim) evlabel += ", ";		// Add a comma if the eigenvalue is complex
      evlabel += ievlabels[ef];			// Add in the imaginary eigenvalue component
      }
    evlabels[ef] = evlabel;			// Store the eigenvalue as a string
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

// ____________________________________________________________________________
// D                       EIGENSYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________

// ------------------ System Wavefunctions in Liouville Space -----------------


void wavefunctions(ostream& ostr, const spin_sys& sys, super_op &LOp,
                          double cutoff, int type, int pbf, int pfz, int title)
 
	// Input	ostr	: An output stream
	//      	sys	: A spin system
	// 		LOp	: A superoperator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
	//              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
	//			: supplied output stream ostr.

  {
  int ls = LOp.size();				// Superoperator dimension size
  if(!ls) return; 				// Exit if NULL superoperator
  int *filter;				// Set up a working filter
  filter = new int[ls];
  for(int i=0; i<ls; i++)			// Set the filter to take everything
    filter[i] = 1;
  wavefunctions(ostr, filter, sys, LOp,
  cutoff, type, pbf, pfz, title);
  delete [] filter;
  return;
  }


void wavefunctions(ostream& ostr, int* filter, const spin_sys& sys,
           super_op &LOp, double cutoff, int type, int pbf, int pfz, int title)
 
	// Input	ostr	: An output stream
	//		filter  : An array of integers, output filter
	//      	sys	: A spin system
	// 		LOp	: A superoperator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
	//              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
	//			: supplied output stream ostr.

    {
    int ls = LOp.size();			// Superoperator dimension size
    if(!ls) return; 				// Exit if NULL superoperator
    string *wflabels;			// Get the wavefunctions
    wflabels = new string[ls];
    wf_labels(wflabels,sys,LOp,cutoff,type,pbf,pfz);
    int len, idig=1, maxlen = 0;
    int wf = 0;
    for(wf=0; wf<ls; wf++) 			// Determine longest wavefunction
      {
      if(filter[wf])				// Consider only those not filtered
        {
        len = (wflabels[wf]).length();
        if(len > maxlen) maxlen = len;
        }
      }
    while(pow(10.0,double(idig)) < ls)		// Get wavefunction printed index size
      idig++;

    len = maxlen;				// Each wavefunction will be this long
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
    string tenblanks = "          ";		// Use this to add spaces in output
    string center = "";				// Construct blanks to center wf's
    string index;				// Use this for wavefunction index
    while(len > 10)
      {
      center += tenblanks;
      len -= 10;
      }
    if(len > 0) center += string(len, ' ');
    if(title)					// Now write a title out if desired
      {
      ostr << "\n"; 				//	Add an additional blank line
      if((LOp.get_mx()).stored_type()		//	If eigenbasis, these will
                             == d_matrix_type)	//	be eigenfunctions
        {
        ostr << tenblanks << tenblanks << tenblanks
             << "System Eigenfunctions";
        }
      else if((LOp.get_Lbasis()).isDefaultBasis())//	If the default basis, these will
        ostr << tenblanks << tenblanks << "  "	//	be the GAMMA default basis functions
             << "System GAMMA Default"
             << " Basis Functions";
      else					//	Else, these are some other basis
        ostr << tenblanks << tenblanks		//	wavefunctions
             << tenblanks
             << "System WaveFunctions";
      ostr << "\n";
      }
    ostr << "\n";				// Output a preceding blank line
    int col = 0;				// Count the columns as we print
    for(wf=0; wf<ls; wf++)			// Now write out the wavefunctions
      {
      if(filter[wf])				// Consider only those not filtered
        {
        index = Gdec(wf); 			//	Get the wavefunction index
        len = index.length();			//	This is the index length
        if(col == 0)				//	For 1st column, new line &
          ostr << "\n" << center;		// 	space over so centered
        else
          ostr << "    ";			//	Else 4 spaces to next column
        ostr << index << "."			//	Write the wavefunciton index
             << string(idig-len, ' ') 		//	Fill so all indices same width
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
 

// -------------------- Eigensystem in Liouville Space ------------------------

void eigensystem(ostream& ostr, const spin_sys& sys, super_op& LOp,
               double cute, double cutc, int type, int pbf, int pfz, int title)

	// Input	ostr	: An output stream
	//      	sys	: A spin system
	//       	LOp	: A superoperator
        //              cute 	: Eigenvalues below this are assumed 0
        //              cutc 	: Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
	//              title   : Flag for inclusion of a title
	// Return	none 	: All eigenvalues and eigenfunctions
	//			  of operator LOp are sent into ostream ostr

    {
    int ls = LOp.size();			// Superoperator dimension size
    if(!ls) return; 				// Exit if NULL superoperator
    LOp.set_EBR();				// Insure LOp into its eigenbasis

    string *wflabels;			// Get the wavefunctions
    wflabels = new string[ls];
    wf_labels(wflabels,sys,LOp,cutc,type,pbf,pfz);
    int len, maxlen = 0;
    for(int wf=0; wf<ls; wf++) 			// Determine longest wavefunction
      {
      len = (wflabels[wf]).length();
      if(len > maxlen) maxlen = len;
      }

    int idig=1;					// Use this for eigenstate indexing
    while(pow(10.0,double(idig)) < ls)		// Get eigenstate printed index size
      idig++;

    string *evlabels;			// Get the system eigenvalues
    evlabels = new string[ls];
    ev_labels(evlabels, LOp, cute);

    if(title)					// Now write a title out if desired
      ostr << "\n               " 		// Will be centered on column 80
           << "Current Eigensystem: "
           << " Eigenvalues and Eigenfunctions"
           << "\n";

    len = (evlabels[0]).length();		// Each eigenvalue is this long
    len = (idig + 2 + len + 2 + maxlen);	// Each eigenstate line is this long
    len = (80 - len)/2;				// Amount to center each line at 80 columns
    if(len < 0)					// If the printed eigenstates are long
      len = 0;					// then don't bother centering
    string tenblanks = "          ";		// Use this to add spaces in output
    string center = "";				// Construct blanks to center wf's
    while(len > 10)				// At the loop end, string center will
      {						// be the spaces to move each eigenstate
      center += tenblanks;			// line over so it looks centered
      len -= 10;
      }
    if(len > 0) center += string(len, ' ');

    string index;				// Use this for wavefunction index
    for(int es=0; es<ls; es++)			// Loop through all the eigenstates
      {
      index = Gdec(es); 				//	Get the eigenstate index
      len = index.length();			//	This is the index length
      ostr << "\n" << center			// 	New line, space so centered
           << index << "."			//	Write the eigenstate index
           << string(idig-len, ' ')		//	Fill so all indices same width
           << evlabels[es] << "  "		//	Write the eigenvalue
           << wflabels[es];			//	Write the wavefunction
      }
    ostr << "\n";				//	End with a blank line
    delete [] evlabels;
    delete [] wflabels;
    return;
    }
 

// ____________________________________________________________________________
// E                             SELECTION RULES
// ____________________________________________________________________________

//  void ev_select(int* select, gen_op &Op, double val1, int type=0,
//                             double val2=0, double cutoff=1.e-4 , int reim=1)

	// Input 	select  : An array of integers
	// 		LOp	: Superoperator (associated to sys)
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

//    {
//    int hs = Op.size();				// Operator dimension
//    if(!hs) return;				// Exit if the basis is empty
//    complex zev;				// For complex eigenvalue

// ______________________________________________________________________
//                    NOT READY FOR PRIME TIME
// ______________________________________________________________________
 

void trans_relax(spin_sys& sys, gen_op &H, super_op& R)

	// Input	sys	: A spin system
	// 		H	: General operator (associated to sys)
	// 		R	: A superoperator (associated to sys)
	// Return	none 	: Indicated wavefunction wf for the working basis
	//			  of operator Op is printed to standard output

  {
  int hs = H.dim();				// Operator dimension hs
  if(!hs) return; 				// Exit if NULL operator
  H.set_EBR();					// Insure H is in its eigenbasis
  R.set_EBR();					// Insure R is in its eigenbasis
  gen_op FzOp = Fz(sys);				// Get the Fz spin operator
//  gen_op FzOp(FZ.matrix(), sys.get_basis());	// Put Fz as a general operator
//  gen_op FzOp (FZ.matrix());			// Put Fz as a general operator
  FzOp.Op_base(H);				// Put Fz in eigebasis of H
  super_op HL = commutator(H);			// Commutation superop. of H
//  super_op FL = commutator(FzOp);		// Commutation superop. of Fz
  int ls = hs*hs;				// Liouville space
  matrix B = R.get_Lbasis().U();		// Retrieve the basis array
  double rate;					// These will be the relaxation rate
  complex zcoeff;
//  int bf = 0;
//  int bout = 0;
  int idig = 1;					// Transition index printed size
  while(pow(10.0,double(idig)) < ls)		// Get the transition index printed size
    idig++;
  string Sidig = string("%i") + Gdec(idig);
  string mode;
  string blanks = "                       ";
  int rows;
  int nrows = 3;
  double maxrate = 0;
  for(int k=0; k<ls; k++)
    {
    rate = Re(R.get(k,k));
    if(fabs(rate) > maxrate)
      maxrate = rate;
    }
  int edig = 1;
  string frm;					// Printing format
  while(pow(10.0,double(edig)) < fabs(maxrate))	// Get maximum rate printed size
    edig++;
  edig += 4;
  int tind;

  cout << "\n\n\t\tSystem Relaxation Analysis\n";
  for(int m=0; m<ls; m++)			// Go through each relaxation mode
    {
    mode = Gdec(m);
    cout << "\n" << mode << "."			// Output the mode index
          << string(idig-mode.length(), ' ');
    rate = Re(R.get(m,m));			// Output the mode relaxation rate
    frm = string("%") + Gdec(edig) + string(".3f");
    cout << Gform(frm.c_str(), rate) << ": ";
    rows = 0;
    tind = -1;
    for(int t=0; t<ls; t++)			// Go through each of the transitions
      {						// First get the contribution of
      zcoeff = B(t,m);				// transition m to mode t
      if(norm(zcoeff) > 1.e-3)
        {
        cout << "  " << Gform("%6.3f", Re(zcoeff))
             << "|" << Gdec(Sidig, t) << ", "
             << Gform("%10.3f", Re(HL.get(t,t))) << ">";
//             << Re(FL.get(t,t)) << ">";
        rows++;
        tind = t;
        if(rows >= nrows)
          {
          cout << "\n" << string(idig+edig+3, ' ');
          rows = 0;
          }
        }
//      if(t == ls-1)
//        cout << ": ~freq = " << Re(HL.get(tind,tind))
//             << "\n"; 
      }
    }
  cout << "\n";
  return;
  }


#endif 							// LSanalyze.cc
