/* SpinSystem.cc ************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	Spin System                                 Implementation	**
**                                                                      **
**	Copyright (c) 1990, 1991, 1992					**
**	Scott Smith							**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description							 	**
**                                                                      **
** Class spin system contains keeps track of a grouping of spins, their	**
** basic properties, and provides their functions & I/O routines.	**
**									**
** This class is derived from a base class "spin_sys" which follows a	**
** grouping spin isotopes, each of which has the following properties -	**
**									**
**	- an isotope type	- arbitrary t/f flag			**
**	- spin quantum number	- 					**
**									**
*************************************************************************/

#ifndef   Spin_system_cc_		// Is file already included?
#  define Spin_system_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <HSLib/SpinSystem.h>		// Includes the interface 
#include <Basics/Gconstants.h>		// Include GFREE & GAMMA1H
#include <Basics/Gutils.h>		// Include GAMMA errors & queries
#include <Basics/StringCut.h>
#include <string>			// Include libstdc++ strings
#include <cmath>			// Inlcude HUGE_VAL_VAL

double HZ2GAUSSx = HZ2GAUSS/2.0;	// For cycles/sec -> Gauss
double spin_system::DefOm=500;		// Set default spec. freq. (MHz)


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
	// Output		void	: An error message is output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(0) 			Program Aborting.....
    		(3) 			Can't Construct From Parameter Set
		(4)			Cannot Construct From Input File
    		(5)			Cannot Write To Parameter File
    		(6)			Cannot Write To Output FileStream
    		(7)			Cannot Produce Any Output
    		default			Unknown Error                        */

void spin_system::SYSTerror(int eidx, int noret) const
  {
  std::string hdr("Isotropic Spin System");
  std::string msg;
  switch (eidx)
    {
    case 8:  msg = std::string("Unspecified Shifts Set to 0 Hz");
             GAMMAerror(hdr,msg,noret);  break; // (8)
    case 9:  msg = std::string("Unspecified Hyperfine Couplings Set to 0 Gauss");
             GAMMAerror(hdr,msg,noret);  break; // (9)
    case 10: GAMMAerror(hdr,"Accessed Spin Index Out of Range",noret);
             break;                              //                        (10)
    case 11: GAMMAerror(hdr,"Accessed Spin Pair Indices Out of Range",noret);
             break;                             //                         (11)
    case 12: GAMMAerror(hdr,"Requested Isotope NULL System", noret);
             break;                             //                         (12)
    case 13: GAMMAerror(hdr,"Cannot Determine # Spins In The System", noret);
             break;                             //                         (13)
    case 14: msg = std::string("Hyperfine Couplings Between Same Spin Forbidden!");
             GAMMAerror(hdr,msg,noret);  break; //                         (14)
    case 15: msg = std::string("Scalar Coupling Between Same Spin Forbidden!");
             GAMMAerror(hdr,msg,noret);  break; //                         (15)
    case 16: msg = std::string("J Coupling Between Nucleus & e- Forbidden! (Use A)");
             GAMMAerror(hdr,msg,noret);  break; //                         (16)
    case 17: msg = std::string("Electron Shift Are NOT Defiled For Nuclei?!!");
             GAMMAerror(hdr,msg,noret);  break; //                         (17)
    case 18: msg = std::string("Setting Electron Shift On a Nucleus Is Forbidden!");
             GAMMAerror(hdr,msg,noret);  break; //                         (18)
    case 19: msg = std::string("Electron G Factors Are NOT Defined For Nuclei!??");
             GAMMAerror(hdr,msg,noret);  break; //                         (19)
    case 20: msg = std::string("Electron Chemical Shifts Are Not Defined In GAMMA");
             GAMMAerror(hdr,msg,noret);  break;	//		           (20)
    case 21: msg = std::string("Unspecified G Factors Set @ Free Electron, 2.00234");
             GAMMAerror(hdr,msg,noret);  break; //                         (21)
    case 30: msg = std::string("Nyquist Frequencies Are Not Defined For Electrons");
             GAMMAerror(hdr,msg,noret);  break; //                         (30)
    case 33: msg = std::string("PPM Shifts Forbidden When No Field Strength Set");
             GAMMAerror(hdr,msg,noret);  break;		//		   (33)
    case 34: msg = std::string("Field Values Forbidden, No Base Frequency Set");
             GAMMAerror(hdr,msg,noret);  break;		//		   (34)
    case 35: msg = std::string("Warning - Unspecified J Couplings Set to 0 Hz");
             GAMMAerror(hdr,msg,noret);  break; //                         (35)
    case 41: msg = std::string("Nuclei-Electron J Couplings Forbidden (Use A)");
             GAMMAerror(hdr,msg,noret);  break; //                         (41)
    case 42: msg = std::string("No Hyperfine Coupling To Electron Pair! (Use J)");
             GAMMAerror(hdr,msg,noret);  break; //                         (42)
    case 43: msg = std::string("No Hyperfine Coupling To Nucleus Pair! (Use J)");
             GAMMAerror(hdr,msg,noret);  break; //                         (43)
    case 44: msg = std::string("Sorry, Cannot Offset Electron Shifts");
             GAMMAerror(hdr,msg,noret);  break; //                         (44)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               none    : Error message output

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

		(1) 			Problems With File PNAME
		(2)			Cannot Read Parameter PNAME
    		default			Unknown Error - PNAME                */

void spin_system::SYSTerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Isotropic Spin System");
  std::string msg;
  switch(eidx)
    {
    case 3:  msg = std::string("Warning - Unspecified Isotope Type(s) Set To ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (3)
    case 4:  msg = std::string("Warning - No Isotope Type Information: ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (4)
    case 10: msg = std::string("Cannot Set Electron Chemical Shift. Use gfactor!");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (10)
    case 102: msg = std::string("No Scalar Coupling Information ");
              GAMMAerror(hdr,msg+pname,noret);  break;	//		   (102)
    case 103: msg = std::string("1H Based Spectrometer Frequency Set To ");
              GAMMAerror(hdr,msg+pname+" MHz",noret);  break;//		   (103)
    case 104: msg = std::string("No Isotropic Chemical Shift Information: ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (104)
    case 106: msg = std::string("No Hyperfine Coupling Information: ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (106)
    case 107: msg = std::string("No G-Factor Information Avalilable: ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (107)
    case 108: msg = std::string("Problems Accessing Chemical Shift, Spin ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (108)
    case 109: msg = std::string("Problems Accessing G-Factor, Spin ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (109)
    case 130: msg = std::string("Parameter ") + pname + std::string(" Is The Culprit");
             GAMMAerror(hdr,msg,noret);  break;		//		   (130)
    case 131: msg = std::string("Requested Parameter ") + pname +
              std::string(" Return Value Is Zero");
             GAMMAerror(hdr,msg,noret);  break;		//		   (131)
    case 200: msg = std::string("Cannot Obtain Larmor Frequency of Spin ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (200)
    case 201: msg = std::string("If Setting Value, Insure Input is Double ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (201)
    case 300: msg = std::string("Request Not Applicable On Electon, Spin ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (300)
    case 301: msg = std::string("Request Not Applicable On Nucleus, Spin ");
             GAMMAerror(hdr,msg+pname,noret);  break;	//		   (301)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error   (-1)
    }
  }

        // Input                sys     : Spin system (this)
        //                      eidx    : Error index
	//			pname	: Additional error message
        // Output               none    : Error message output
        //                                Program execution stopped

volatile void spin_system::SYSTfatality(int eidx) const
  {  
  SYSTerror(eidx, 1);				// Normal non-fatal error
  if(eidx) SYSTerror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

volatile void spin_system::SYSTfatality(int eidx, const std::string& pname) const
  {  
  SYSTerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) SYSTerror(0);			// Program aborting error
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                 CLASS SPIN SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   spin system to be set up without worrying about system consistency!       */


	// Input		sys      : Spin system (this)
	// 			pset     : A parameter set
	//			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		TF	 : Spin system is filled with
	//				   parameters in pset
	// Note				 : Three things are gleaned from
	//				   the parameter set from spin_sys
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	// 				   Three more parameters are taken
	//				   from pset to complete spin_system
	//				   4.) The chemical shifts
	//				   5.) The coupling constants
	//				   6.) A spectrometer frequency
	// Note				 : Functions which place a spin_system
	//				   into a parameter set must contain
	//				   add the information read here

int spin_system::setSsys(const ParameterSet& pset, int idx, int warn)
  {
//	Use Only Parameters With Prefix [#], So Clip [#] From Names First

  ParameterSet  subpset;			// Working parameter set
  if(idx != -1) subpset=pset.strip(idx);	// Get only params with [#]
  else          subpset=pset;			// Or use full pset

//	          Now Set Up The Spin System Directly

  int TF=1;					// Track if we read OK
  int ns = getSpins(subpset, warn?1:0);		// Get the number of spins
  if(ns<=0)
    {
    if(warn) SYSTerror(13,1);
    else     SYSTerror(13);       
    return 0;
    }

  *this = spin_system(ns);			// Set the system for ns spins
  setIs(subpset);				// Set the isotope types
  setName(subpset);				// Read in the system name
  setBasis(matrix(HS(), HS(), i_matrix_type));	// Set up default basis (matrix)
  setOm(subpset);				// Set spectrometer frequency
  setShifts(subpset);				// Set the chemical shifts
  setJs(subpset); 				// Read in scalar couplings
  if(electrons())				// Read G's & J's if electrons
    {						// present in the spin system
    setGs(subpset);				// 	Set g-factors
    setAs(subpset);				// 	Set hyperfine couplings
    } 
  return TF;
  }

// ----------------------------------------------------------------------------
//                       Set External Bo Field Strength
// ----------------------------------------------------------------------------

/* This function attempts to determine the stataic external field strength, Bo.
   The following parameters are looked for in the input parameter set:

         Omega       - Proton Larmor frequency (MHz)
       * GOmega      - EPR base frequency      (GHz)
         Field       - Bo field strength       (Gauss)
         FieldT      - Bo Field strength       (Tesla)

   The * parameter is allowed only if there are electrons present in the spin
   system. If multiple spin systems are read in a program, i.e. the above
   parameters have a prefix [#], the prefix should be stripped from the pset
   prior to calling this function.


        // Input                sys	: A spin system (this)
        //                      pset    : A parameter set
        // Output               TF      : Spin system field strength
        //                                is set from parameter in pset
        //                                TRUE if set from pset, FALSE
        //                                if set to a default value
	// Note				: Assumes that the base class spin
	//				  system has already been read
	// Note				: Parameter Omega sets the frequency
	//				  based on 1H resonance (in MHz)
	// Note				: If electrons are present in the
	//				  system, GOmega (in GHz) is allowed 
	//				  for the e- base Larmor frequency  */
 
bool spin_system::setOm(const ParameterSet& pset)
  {
  double dval = 0;				// Default spec. freq.
  std::string pstate, pname = "Omega";               // Retrieve spect. frequency
  bool Omegaf = false;				// Flag if read properly
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Pix in param. list for Omega
  if(item != pset.end()) 			// If Omega has been defined 
    {						// we'll parse the information
    (*item).parse(pname,dval,pstate);		// & set 1H based spectromter
    spectrometer_frequency(fabs(dval));		// frequency
    return true;				// Exit, we found Omega
    }
  if(!Omegaf && electrons()) 			// If electrons present, also
    {						// allow Omega to be set in GHz
    pname = "GOmega";				// for e-.  Parameter GOmega
    item = pset.seek(pname);			// Pix in param. list for GOmega
    if(item != pset.end())			// If GOmega has been defined 
      {						// we'll parse the info and
      (*item).parse(pname,dval,pstate);		// set e- based spectrometer
      Omega(fabs(dval)*1.e3, "e-");		// Convert from GHz to MHz too
      return true;				// Exit, we found GOmega
      }
    }
  if(!Omegaf) 					// We can try looking for a
    {						// field strength directly
    pname = "Field";				// This will be assumed in Gauss
    item = pset.seek(pname);			// Pix in param. list for Field
    if(item != pset.end())			// If Field has been defined 
      {						// we'll parse the info and
      (*item).parse(pname,dval,pstate);		// set field appropriately
      Omega(fabs(dval)*GAMMA1H*1.e-10/HZ2RAD);	// using proton Larmor frequency
      return true;				// Exit, we found Field
      }
    }
  if(!Omegaf) 					// Last, we'll try looking for
    {						// a field strength in Tesla
    pname = "FieldT";				// This will be assumed in Tesla
    item = pset.seek(pname);			// Pix in param. list for FieldT
    if(item != pset.end())			// If Field has been defined 
      {						// we'll parse the info and
      (*item).parse(pname,dval,pstate);		// set field appropriately
      Omega(fabs(dval)*GAMMA1H*1.e-6/HZ2RAD);	// using proton Larmor frequency
      return true;				// Exit, we found FieldT
      }
    }
  if(!Omegaf)					// If we still haven't found a
    {						// defined field strength then
    spectrometer_frequency(DefOm);		// just set a default value
    if(warnings())				// Output these warnings if we
      {						// want them......
      SYSTerror(2, "Omega", 1); 		// Can't read in Omega
      if(electrons()) SYSTerror(2, pname, 1);	// Can't read in GOmega
      SYSTerror(103, Gform("%7.3f", DefOm));	// Omega at 500 MHz (FALSE)
      } 
    }
  return Omegaf;
  }

// ----------------------------------------------------------------------------
//                       Set Nuclear Spin Chemical Shifts
// ----------------------------------------------------------------------------

/* This function attempts to determine the chemical shift for each nucleus in
   the system. The shift may be set in Hz or in PPM, the latter being valid
   only if a field strength has been set. Electrons do not have shifts assigned.
   Instead, electrons will have an istropic g factor. Prior to calling this
   function the number of spins should be known. So should the field strength
   and any isotope labels be set.

         Omega       - Proton Larmor frequency (MHz)
       * GOmega      - EPR base frequency      (GHz)
         Field       - Bo field strength       (Gauss)
         FieldT      - Bo Field strength       (Tesla)

   All parameters are allowed only if there spin is a nucleus.
   * parameters are allowed only if the Bo field has been set
 
        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system isotropic shifts
        //                                are set from parameters in pset
	// Note				: If the system spectrometer
	//				  frequency is not present, PPM
	//				  shifts are disabled
	// Note				: Isotope labels should be set
	//				  prior to calling this function!
	//				: Electrons will not have any shifts
	//				  set from this function             */

void spin_system::setShifts(const ParameterSet& pset)
  {
  int ns = spins();			// Number of spins
  std::vector<int> inPPM(ns, 0);	// Flag if shifts in PPM

  std::string Pst("PPM(");		// Start of PPM parameter name
  std::string Pfi(")");			// Finish of PPM param. names
  ParameterSet::const_iterator item;	// A pix into parameter list
  int i, notfound = ns;			// No. with no shift assinged
  std::string pname, pst;			// For parameter name & statement
  double dval;				// For parameter value

//	     First See If Shift Information Is Specified in PPM
//             (All That Are Have PPM Values Set & Flagged)

  if(Omega())				// Don't do if no Omega set
    {
    for(i=0; i<ns; i++)			// Attempt to retrieve the
      { 				// chemical shifts in PPM
      if(!electron(i))			// Look for shift in PPM (if not e-)
        {
        pname = Pst + Gdec(i) + Pfi; 	// All shifts need not be present
        item = pset.seek(pname);	// Pix in parameter list for PPM(i)
        if(item != pset.end())		// Retrieve the number of spins
          {
          (*item).parse(pname,dval,pst);// Parse out the value as double
          PPM(i, dval);			// Set shift in PPM
          inPPM[i] = 1;			// Flag that PPM units used
          notfound--;			// Decrement count of shifts
          }				// we don't know yet
        }
      }
    }

//	       See If Shift Information Is Specified in Hz 

  std::string Pstv = "v("; 		// Start/finish Hz param. names
  for(i=0; i<ns; i++)			// Attempt to retrieve chemical
    { 					// shifts in Hz. Only do this
    if(!inPPM[i] && !electron(i))	// if not input in PPM
      {
      pname = Pstv + Gdec(i) + Pfi; 	// All shifts need not be present
      item = pset.seek(pname);		// Pix in parameter list for v(i)
      if(item != pset.end()) 		// Retrieve the number of spins
        {
        (*item).parse(pname,dval,pst);
        shift(i, dval);			// Set the shift in Hz
        inPPM[i] = -1;			// Flag that this is found
        notfound--;			// Decrement not found count
        }
      }
    }

//	       See If Spins With No Shift Assigned Are Electrons

  for(i=0; i<ns; i++)			// Attempt to retrieve chemical
    if(!inPPM[i] && electron(i))	// if not input in PPM
      {
      notfound--;			// Decrement not found count
      inPPM[i] = 2;			// Flag that (no) value is OK
      cshifts[i] = 0;			// Set electron shift to 0
      }

//	       Now Output Warnings If Any Spins Not Shift Assigned 
//		       (Only Done if Warnings Are Active)

  std::string serr;
  if(notfound <= 1) return;		// Exit if all are found
  if(!warnings()) return;		// Exit if warnings not wanted
  int j=0;			
  for(i=0; i<ns; i++)			// Loop over shifts not set
    {
    j++;
    if(j==1 && !inPPM[i])		// 	Start a new line of
      { 				//	parameter warnings
      if(Omega())			// 	Output PPM(#) if Omega
        {				//	has been set and v(#)
        serr = Pst + Gdec(i) + Pfi; 	//	in all cases
        serr += std::string(" nor ");
        serr += pname = Pstv + Gdec(i) + Pfi;
        SYSTerror(104, pname, 1);
        }
      else
        {
        pname = Pstv + Gdec(i) + Pfi;
        SYSTerror(104, pname, 1);
        }
      }
    else if(!inPPM[i])			//	Continuation of line
      {					//	Output PPM(#) if Omega
      if(Omega())			//	has been set and v(#)
        {				//	if Omega has been set
        pname = Pst + Gdec(i) + Pfi; 	//	in all cases
        std::cout << "; " << pname;
        pname = Pstv + Gdec(i) + Pfi;
        std::cout << " nor " << pname;
        }
      else
        {
        pname = Pstv + Gdec(i) + Pfi;
        std::cout << "; " << pname;
        }
      if(j>=2 && Omega()) j=0;	// 	Flag that new line of warnings
      else if(j>5) j=0; 		//	should start
      }
    }
  SYSTerror(8); 			// Warning, unspecified shifts at 0 Hz
  } 

 
  
// ____________________________________________________________________________
// iii              CLASS SPIN SYSTEM SPIN CHECKING FUNCTIONS
// ____________________________________________________________________________

        // Input                sys     : Spin system (this)
        // 			spin	: Spin index
	//			warn    : Flag if fatal error can occur
	//				      0 = no warnings
	//				      1 = non-fatal warnings
	//				     >1 = fatal error
        // Output               T/F	: TRUE if spin isn't electon
	//				  FALSE otherwise

bool spin_system::checkNotE(int spin, int warn) const 
  {
  if(!electron(spin)) return true;		// All is O.K. if not e-
  if(warn)					// If warnings desired
    {						// lets output 'em
    if(warn > 1) SYSTfatality(300,Gdec(spin));	//  Here if fatal error
    else	 SYSTerror(300,Gdec(spin),1);	//  Here if only a warning 
    }
  return false;					// This spin was e-
  }


bool spin_system::checkNotN(int spin, int warn) const 

        // Input                sys     : Spin system (this)
        // 			spin	: Spin index
	//			warn    : Flag if fatal error can occur
	//				      0 = no warnings
	//				      1 = non-fatal warnings
	//				     >1 = fatal error
        // Output               T/F	: TRUE if spin isn't a nucleus
	//				  FALSE otherwise
  {
  if(electron(spin)) return true;		// All is O.K. it is e-
  if(warn)					// If warnings desired
    {						// lets output 'em
    if(warn > 1) SYSTfatality(301,Gdec(spin));	//  Here if fatal error
    else	 SYSTerror(301,Gdec(spin),1);	//  Here if only a warning 
    }
  return false;					// This spin was not e-
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

spin_system::spin_system(int spins) : spin_sys(spins)
  {
  Omega1H = 0;					// No spectrometer freq
  cshifts = std::vector<double> (spins, 0.0);	// Array of shifts
  gfacts  = std::vector<double> (spins, 0.0);	// Array of g-factors
  int nsp = (spins*spins-spins)/2;		// Number of spin pairs
  Jcouplings = std::vector<double> (nsp,0.0);	// Array of J couplings
  Acouplings = std::vector<double> (nsp,0.0);	// Array of A couplings
  _spflags   = std::vector<int>    (nsp,  1); 	// Array of spin pair flags
  }

spin_system::spin_system(const spin_system& sys) : spin_sys(sys)
  {
  Omega1H = sys.Omega1H;			// Copy spectrometer freq.
  cshifts = sys.cshifts;			// Copy array of shifts
  gfacts  = sys.gfacts;				// Copy array of g-factors
  Jcouplings = sys.Jcouplings;			// Copy array of J couplings
  Acouplings = sys.Acouplings;			// Copy array of A couplings
  _spflags   = sys._spflags;			// Copy array of spin pair flags
  }

spin_system& spin_system::operator= (const spin_system &sys)
  {
  if(this == &sys) return *this;
  spin_sys::operator=(sys);			// Copy spin sys stuff
  Omega1H = sys.Omega1H;			// Copy spectrometer freq.
  cshifts = sys.cshifts;			// Copy array of shifts
  gfacts  = sys.gfacts;				// Copy array of g-factors
  Jcouplings = sys.Jcouplings;			// Copy array of J couplings
  Acouplings = sys.Acouplings;			// Copy array of A couplings
  _spflags   = sys._spflags;			// Copy array of spin pair flags
  return *this;
  }

spin_system::~spin_system () {}

// ____________________________________________________________________________
// B                CHEMICAL SHIFT AND G FACTOR MANIPULATIONS 
// ____________________________________________________________________________

// ************* Chemical Shift Manipulations in Hertz & PPM ******************
 
/* These functions allow users to either get or set the isotropic chemical
   shift of a spin or set of spins. There are a few important items worth
   taking note of.

   1.) Electrons CANNOT have their "shift" set with these functions.  Users
       must set the electron g factor to get the ESR semi-equivalent for e-.
   2.) GAMMA assumes system shifts are stored relative to a base Larmor
       frequency particular for the associated spin's isotope type.  For a
       positive gyromagnetic ratio, a spin with a positive shift is deshielded 
       relative to what the base Larmor frequency is.
   3.) Shifts in PPM may not be set/obtained unless a spectrometer field
       strength has been specified (Omega).
  
       Function    Return   Arguments	           Results
      ----------   ------   ---------   --------------------------------------
       shift       void      spin,w	Set spin shift to w spin = [0,ns)
       shift       double     spin      Retrieve shift of spin in Hz
       lab_shift   double     spin      Retrieve shift of spin in Hz, lab frame
       PPM         spin      spin,PPM   Set spin shift to PPM
       PPM         void       spin      Retrieve shift to spin in PPM
       shifts      double       w       Set all nuclear spin shifts to w in Hz
       maxShift    double               Retrieve largest shift in Hz
       maxShift    double      iso      Get largest shift in Hz, iso spins only
       minShift    double               Retrieve smallest shift in Hz
       minShift    double      iso      Get smallest shift Hz, iso spins only
       medianShift double               Get middle of shifts in Hz
                                                                             */
void spin_system::shift(int spin, double shift)
  {
  if(!check_spin(spin,1)) SYSTfatality(108,Gdec(spin));	// Insure spin exists
  if(!checkNotE(spin,1))  SYSTfatality(10);             // Insure spin not e-
  cshifts[spin]=shift;                                  // Set chemical shift
  }
 
double spin_system::shift(int spin) const
  { check_spin(spin, 2); return cshifts[spin]; }
 
double spin_system::lab_shift(int i) const
  {
  if(!check_spin(i, 1)) SYSTfatality(108,Gdec(i));	// Insure spin exists
  if(!checkNotE(i,1))  SYSTfatality(10);                // Insure spin not e-
  return cshifts[i]+((gamma(i)/GAMMA1H)*(Omega1H*1.e6));// Lab frame shift out
  }
 
void spin_system::shifts(double shift)
  { for(int i=0; i<spins(); i++) if(!electron(i)) cshifts[i]=shift; }

double spin_system::maxShift() const
  {
  double maxf = -HUGE_VAL;				// Start with a big negative
  for(int i=0; i<spins(); i++)			// Check all spins' shifts
    if(!electron(i))				// unless they are electrons
      if(maxf < cshifts[i]) maxf=cshifts[i];
  return maxf;
  }

double spin_system::maxShift(const std::string& Iso) const
  {
  Isotope I(Iso);				// Make isotope from Iso
  if(I.electron()) SYSTfatality(20);		// Electrom shifts disallowed
  double maxf = -HUGE_VAL;				// Start with a big negative
  for(int i=0; i<spins(); i++)
    if(symbol(i)==Iso && maxf<cshifts[i])
      maxf = cshifts[i];
  return maxf;
  }

double spin_system::minShift() const
  {
  double maxf=HUGE_VAL;
  for(int i=0; i<spins(); i++)			// Check all spins' shifts
    if(!electron(i))				// unless they are electrons
      if (maxf>cshifts[i]) maxf=cshifts[i];
  return maxf;
  }

double spin_system::minShift(const std::string& Iso) const
  {
  Isotope I(Iso);				// Make isotope from Iso
  if(I.electron()) SYSTfatality(20);		// Electrom shifts disallowed
  double minf = HUGE_VAL;				// Start with a big number
  for(int i=0; i<spins(); i++)
    {
    if(symbol(i)==Iso && minf>cshifts[i])
      minf = cshifts[i];
    }
  return minf;
  }

double spin_system::medianShift() const { return (maxShift()+minShift())/2; }

void spin_system::offsetShifts(double shift, int spin)
  {
  std::string Iso = symbol(spin);			// Get isotope type of spin
  offsetShifts(shift, Iso);			// Use overload function
  }

void spin_system::offsetShifts(double shift, const std::string& Iso)
  {
  Isotope I(Iso);				// Make isotope from Iso
  if(I.electron()) SYSTfatality(20);		// Electrom shifts disallowed
  for(int i=0; i<spins(); i++)
    if(Iso == symbol(i)) cshifts[i]-=shift;
  }

void spin_system::PPM(int spin, double PPMval)
  {
  if(!check_spin(spin,1)) SYSTerror(108,Gdec(spin),1);  // Insure spin exists
  if(!checkNotE(spin,1))  SYSTfatality(10);             // Insure spin not e-
  if(Omega1H == 0.0)      SYSTfatality(33);             // Fatal if no Bo Set
  cshifts[spin] = PPMval*gamma(spin)*Omega1H/GAMMA1H;   // Set shift of spin
  }
 
double spin_system::PPM(int spin) const
  {
  if(!check_spin(spin,1)) SYSTfatality(108,Gdec(spin)); // Insure spin exists
  if(!checkNotE(spin,1))  SYSTfatality(10);             // Insure spin not e-
  if(Omega1H == 0.0)      SYSTfatality(33);             // Fatal if no Bo Set
  return cshifts[spin]*GAMMA1H/(Omega1H*gamma(spin));   // Return shift in PPM
  }

// ************************** G Factor Manipulations **************************

/* These functions allow users to obtain or set electon g-factors.  These are
   NOT applicable to nuclear spin! The g-factor is a unitless quantity that
   will typically have a value near that of g for a free electron: 2.00232.
   Frequency <-> Field conversions can be made based on this value using the
   conversion factor Gauss2Hz which is planks constant (6.6262e-27 erg-s/cycle)
   over the Bohr magneton (9.2471e-21 erg/G). Gauss2Hz = 0.714474e-6 G/Hz.
   GAMMA also can track an "electron shift" which is the electon resonance at
   a particular filed strength relative to a base Larmor frequency of a free
   electron.  Note that e- has a negative gyromagnetic ratio, thus electrons
   with a positive shift are shielded!                                       */

double spin_system::gfactor(int spin) const
  {
  if(!check_spin(spin,1)) SYSTerror(109,Gdec(spin),1);	// Insure spin exists
  if(!checkNotN(spin,1)) { SYSTerror(19); return 0; } 	// Insure not nucleus
  return gfacts[spin];					// Return g-factor
  }

void spin_system::gfactor(int spin, double g)
  {
  if(!check_spin(spin,1)) SYSTfatality(109,Gdec(spin));	// Insure spin exists
  if(!checkNotN(spin,1)) SYSTfatality(19); 		// Insure not nucleus
  else gfacts[spin] = g;				// Set g-factor
  }

double spin_system::eshift(int spin) const
  {
  check_spin(spin);					// Insure spin exists
  if(!checkNotN(spin,1)) { SYSTerror(17); return 0; } 	// Insure not nucleus
  return (gfacts[spin]/GFREE - 1.00)*Omega("e-")*1.e6;	// Return e- "shift"
  }							// (typically -, MHz!)

double spin_system::lab_eshift(int spin) const
  {
  if(!check_spin(spin,1)) SYSTerror(109,Gdec(spin),1);	// Insure spin exists
  if(!checkNotN(spin,1)) { SYSTerror(17); return 0; }	// Insure not nucleus
  return (gfacts[spin]/GFREE)*Omega("e-")*1.e6;		// Return e- shift
  }							// (typically -, GHz!)
 
double spin_system::efield(int spin) const
  {
  if(!check_spin(spin,1)) SYSTfatality(108,Gdec(spin)); // Insure spin exists
  if(!checkNotN(spin,1)) { SYSTerror(17); return 0; }	// Insure not nucleus
  if(Omega1H == 0.0)      SYSTfatality(34);             // Fatal if no GOM Set
  return Bo()*(GFREE/gfacts[spin] - 1.0);		// Resonance field
  }							// (typically in G)
 
double spin_system::efield_lab(int spin) const
  {
  if(!check_spin(spin,1)) SYSTfatality(108,Gdec(spin)); // Insure spin exists
  if(!checkNotN(spin,1)) { SYSTerror(17); return 0; }	// Insure not nucleus
  if(Omega1H == 0.0)      SYSTfatality(34);             // Fatal if no GOM Set
  return Bo()*GFREE/gfacts[spin];			// Resonance field
  }							// (typically in kG)

// ____________________________________________________________________________
// C       SCALAR & HYPERFINE COUPLING CONSTANT MANIPULATION FUNCTIONS
// ____________________________________________________________________________

// ----------------------- Scalar Coupling Functions --------------------------

        // Input                sys     : A spin system (this)
        //			Jval	: Scalar coupling value (Hz)
	// 			int	: Spin index [0, nspins)
	//       		int	: Spin index [0, nspins)
	// Output Js		none	: Sets all J couplings to Jval
	//        J(i,j,d)      none	: Sets (i,j) J coupling to Jval
	//        J(d,i,j)      none	: Sets (i,j) J coupling to Jval
	//        J(i,j)	double	: Returns (i,j) J coupling in Hz

void spin_system::Js(double Jval)
  {
  int i,j,k, ns = spins();
  for(i=0,k=0; i<ns-1; i++)
    for(j=i+1; j<ns; j++, k++)
      if(!enpair(i,j)) Jcouplings[k] = Jval;
      else             Jcouplings[k] = 0;
  }

void spin_system::J(int spin1, int spin2, double coupling)
  {
  check_spin(spin1);                            // Insure 1st spin exists
  check_spin(spin2);                            // Insure 2nd spin exists
  if(enpair(spin1, spin2))
    {
    std::string p = std::string("J(")
                  + Gdec(spin1) + std::string(",")
                  + Gdec(spin2) + std::string(")");
    SYSTerror(16,1);				// Invalid for e-/nucleus pair
    SYSTfatality(130,p);			// Abort with J value being set
    }
  if(spin1 == spin2)        SYSTfatality(15);	// Can't be the same spin
  Jcouplings[pairidx(spin1, spin2)] = coupling;	// Set hyperfine coupling
  }

void spin_system::J(double Jval, int spin1, int spin2)
  { J(spin1, spin2, Jval); }
 
double spin_system::J(int spin1, int spin2) const
  {
  check_spin(spin1);                            // Insure 1st spin exists
  check_spin(spin2);                            // Insure 2nd spin exists
  if(enpair(spin1, spin2) && warnings())
    {
    std::string p = std::string("J(") + Gdec(spin1) + std::string(",")
             + Gdec(spin2) + std::string(")");
    SYSTerror(16,1);				// Invalid for e-/nucleus pair
    SYSTerror(131,p);				// Requested value return is 0
    }
  if(spin1 == spin2)       SYSTfatality(15);	// Can't be the same spin
  return Jcouplings[pairidx(spin1,spin2)];	// Return scalar coupling
  }


// --------------------- Hyperfine Coupling Functions -------------------------

        // Input                sys     : A spin system (this)
	// 			Aval	: Hyperfine coupling in Gauss
	// 			int	: Spin index [0, nspins)
	//       		int	: Spin index [0, nspins)
	// Output As		none	: Sets all HF couplings to Aval
	//        A(i,j,d)      none	: Sets (i,j) HF coupling to Aval
	//        A(d,i,j)      none	: Sets (i,j) HF coupling to Aval
	//        A(i,j)	double	: Returns (i,j) HF coupling in G
	//        AHz(i,j)	double	: Returns (i,j) HF coupling in Hz

void spin_system::As(double Aval)
  {
  int i,j,k, ns = spins();			// Number of spins
  for(i=0,k=0; i<ns-1; i++)			// Loop spin pairs
    for(j=i+1; j<ns; j++, k++)			// Set all hyperfine couplings
      if(enpair(i,j)) Acouplings[k] = Aval;
      else            Acouplings[k] = 0;
  }

void spin_system::A(int spin1, int spin2, double Aval)
  {
  check_spin(spin1);                            // Insure 1st spin exists
  check_spin(spin2);                            // Insure 2nd spin exists
  if(!enpair(spin1, spin2))
  { if(electron(spin1)) SYSTfatality(42);	// Can't be e-/e- pair
    else                SYSTfatality(43);	// Can't be nucleus/nucleus pair
  }
  if(spin1 == spin2)    SYSTfatality(14);	// Can't be the same spin
  Acouplings[pairidx(spin1, spin2)] = Aval;	// Set hyperfine coupling
  } 

void spin_system::A(double Aval, int spin1, int spin2)
  { A(spin1, spin2, Aval); }
 
double spin_system::A(int spin1, int spin2) const
  {
  check_spin(spin1);                            // Insure 1st spin exists
  check_spin(spin2);                            // Insure 2nd spin exists
  if(!enpair(spin1, spin2))
  { if(electron(spin1)) SYSTfatality(42);	// Can't be e-/e- pair
    else                SYSTfatality(43);	// Can't be nucleus/nucleus pair
  }
  if(spin1 == spin2)    SYSTfatality(14);	// Can't be the same spin
  return Acouplings[pairidx(spin1,spin2)];	// Return hyperfine coupling
  }

double spin_system::AHz(int spin1, int spin2) const
 { return A(spin1, spin2)/HZ2GAUSSx; }


// ____________________________________________________________________________
// D              SPECTROMETER FREQUENCY AND FIELD MANIPULATIONS
// ____________________________________________________________________________

/* These functions allow the user to set or obtain the spectrometer field
   strength.  Typically this is done based on the proton Larmor frequency,
   using the function Omega.  To simulate working on a 500 MHz NMR spectrometer
   you simply set the system Omega to 500.  Additional functions allow users
   to set the field strength based on a resonance frequency of ANY isotope
   (if 1H doesn't suit your fancy).  Last but not least you can access, but not
   set, the field strength in Tesla using the function Bo. Function Bo will of
   course always return a non-negative number.

   --> IMPORTANT NOTE: (Re-)Setting the field strength using an Omega function
       DOES NOT CHANGE EXISTING SYSTEM SHIFTS IN HZ!!!  That means that it does
       change the system shifts in PPM.  In contrast, the function OmegaAdjust
       KEEPS THE SYSTEM SYSTEM SHIFTS IN PPM CONSTANT. For example, say you
       which to do a field study by looping through a simulation at different
       field strengths ---> use OmegaAdjust to change the field strength so
       so that you system shifts will be constant on a PPM scale but properly
       change on a Hz scale to reflect altered field.

       Function      Return   Arguments	           Results
      ----------     ------   ---------   -------------------------------------
       Omega         void        Om	  Set 1H Larmor frequency to Om (MHz)
       Omega         void      Om, iso	  Set iso Larmor frequency to Om (MHz)
       Omega         double     spin	  Get Larmor frequency to spin
       Omega         double     iso	  Get Larmor frequency for isotope iso
       Bo	     double        	  Get Field Strength (Tesla)
       OmegaAdjust   double      Om       Set 1H Larmor freq. to Om (MHz)    */

void   spin_system::Omega(double freq) { Omega1H=freq; }
void   spin_system::Omega(double FR, const std::string& I)
                                       { Omega1H = GAMMA1H*FR/fabs(gamma(I)); }
double spin_system::Omega(int i) const
  {
  if(i<0) return Omega1H;			// Omega for 1H if i negative
  if(!check_spin(i,1))				// Insure spin accessed O.K.
    {
    SYSTerror(200, std::string(Gdec(i)), 1);		//   Can't set i based Omega
    SYSTfatality(201,std::string(Gdec(i))+std::string(".0"));// Warn if setting with int
    }						//   because Omega overload
  return gamma(i)*Omega1H/GAMMA1H;		//   can cause trouble here
  }
double spin_system::Omega(const std::string& iso) const
                                       { return gamma(iso)*Omega1H/GAMMA1H; }
double spin_system::Bo() const         { return Omega1H*HZ2RAD*(1.e10/GAMMA1H); }
void   spin_system::OmegaAdjust(double Om)
  {
  int i, ns=spins();				// Number of spins
  double* ppmvals;
  ppmvals = new double[ns];
  for(i=0; i<ns; i++) 				// Store PPM shift values
    if(!electron(i)) ppmvals[i] = PPM(i);	// for all nuclei
  Omega(Om);					// Reset field strength
  for(i=0; i<ns; i++) 				// Restore PPM shift values
    if(!electron(i)) PPM(i, ppmvals[i]);	// for all nuclei
  delete [] ppmvals;
  }
void   spin_system::FieldAdjust(double B)
  {
  double Om1H = Omega();
  double B0 = Bo();
  OmegaAdjust(Om1H*fabs(B)/B0);
  }
 
/* These functions work but are now deprecated, they're identical to Omega. */

void   spin_system::spectrometer_frequency(double freq) { Omega1H=freq; }
double spin_system::spectrometer_frequency() const      { return Omega1H; }

// ____________________________________________________________________________
// E                         SPIN PAIR FLAG FUNCTIONS
// ____________________________________________________________________________
 
// --------- Functions Which Get/Set Spin Pair Flags Within The System --------
 
void spin_system::spflags(int TF)
 
        // Input                   TF : TRUE/FALSE status
        // Output                none : All spin pairs  have their flags set
        //                              to TRUE/FALSE
  {
  int s=spins();				// Number of spins
  int nsp = (s*s-s)/2;				// Number of spin pairs
  for(int i=0; i<nsp; i++) _spflags[i] = TF;
  }
 
 
void spin_system::spflag(int spin1, int spin2, int TF)
 
        // Input               spin1 : Spin in the system
        //                     spin2 : Another spin in the system
        //                      int  : TRUE/FALSE status
        // Output               none : Flag for spin pair spin1-spin2 has
        //                             its spin pair flag set to TF
        ///F_list flag               - Set/Get spin pair's T/F flag status.
	// Note			     : This allows the index to span beyond
	//			       the number of spin pairs available

  {
  check_spin(spin1);			// Insure spin 1 exists
  check_spin(spin2);			// Insure spin 2 exists
  _spflags[pairidx(spin1, spin2)] = TF;	// Toggle spin pair flag
  }


int spin_system::spflag(int spin1, int spin2) const

        // Input               spin1 : Spin in the system
        //                     spin2 : Another spin in the system
        // Output               int  : Current T/F status of spin
        //                             pair spin1-spin2

  {
  check_spin(spin1);			// Insure spin 1 exists
  check_spin(spin2);			// Insure spin 2 exists
  return _spflags[pairidx(spin1, spin2)];
  }


// ____________________________________________________________________________
// F                      SPIN SYSTEM ASSOCIATED FUNCTIONS
// ____________________________________________________________________________

	// Input	sys	: Spin system(this)
	// 		spin	: Spin number
	// Return	centerf : Spin system frequency center based
	//			  on the shifts of all spins of the
	//			  same isotope type as input spin
	// Note		     	: Returned value is in Hertz
	// Note			: This is not allowed for electrons
	//			  since we do not store any relative
	//			  shifts

double spin_system::center(int spin)
  {
  if(electron(spin))
    {
    SYSTerror(20, 1);			// No electron shifts defined
    SYSTerror(44);			// Cant offset electron shifts
    }
  double centerf, maxf, minf;
  Isotope iso = isotope(spin);
  minf = HUGE_VAL;
  maxf = -HUGE_VAL;
  for(int i=0; i<spins(); i++)		// Get max & min shifts
    {
    if(iso == isotope(i))
      {
      if(shift(i) > maxf) maxf = shift(i);
      if(shift(i) < minf) minf = shift(i);
      }
    }
  centerf = (maxf - minf)/2.0 + minf;
  return centerf;
  }

// ____________________________________________________________________________
// G	                     Nyquist Frequency Functions
// ____________________________________________________________________________

	// Input	sys	: Spin system (this)
	//		<***>   : Selectivity
	//		fact    : Extension Factor (%)
	//		lwhh    : Anticipated half-height linewidth
	// Return	Nyqf    : Approximate Nyquist frequency needed to
	//			  insure proper quadrature acquisitions
	//			  without foldover.  This is based on the
	//			  spin system shifts and coupling constants
	//			  for the spins of the same isotope type as
	//			  the input spin
	// Note		     	: Returned value is in Hertz
	// Note		     	: For broad peaks due to relaxation or
	//			  severe apodization this can return value 
	//			  too small if no lwhh is set.  Compensation
	//			  can be made externally by adjusting the
	//			  returned value.

	// Overloads	<***>	: Selectivity 
	// 1. int	 spin   : Based on this spin's isotope type
	// 2. string	 iso    : Based on this isotope designation
	// 3. Isotope	 iso    : Based on this isotope


double spin_system::Nyquist(int spin, double fact, double lwhh) const
  { return Nyquist(isotope(spin),fact,lwhh); }

double spin_system::Nyquist(const std::string& iso, double fact, double lwhh) const
  {
  Isotope I;
  if(!I.exists(iso))
    {
    double Nyqf = 5.0*fabs(lwhh)*fact;
    if(Nyqf == 0) Nyqf = 100.;
    return Nyqf;
    }
  return Nyquist(Isotope(iso),fact,lwhh);
  }

double spin_system::Nyquist(const Isotope& iso, double fact, double lwhh) const
  {
  if(iso.electron()) SYSTfatality(30);		// No Nyquist for e-
  double Nyqf = 0;				// Initialize frequency to 0
  double maxf = -HUGE_VAL;				// Set max. frequency way high
  int maxi = 0;
  for(int i=0; i<spins(); i++)			// Get max shift magnitude
    {						// over the isotope type
    if(iso == isotope(i))			// desired.  This will produce
      {						// the largest frequency that
      if(fabs(shift(i)) > maxf)			// one has to deal with due to
        {					// chemical shifts
        maxf = fabs(shift(i));
        maxi = i;
        }
      }
    }
  for(int j=0; j<spins(); j++)			// Account for Js too
    if(!electron(j) && maxi!=j)			// which will cause frequencies
      maxf += fabs(J(maxi,j)/2.0);		// beyond the farthest shift
  maxf += 5.0*fabs(lwhh);			// Account for broad lines 
  Nyqf = maxf*fact;				// Add a bit extra (10% default)
  if(Nyqf == 0) Nyqf = 100.;			// If Nyquist = 0 set to default
  return Nyqf;
  }

// ____________________________________________________________________________
// H                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Functions To Make A Parameter Set From A Spin System
// ----------------------------------------------------------------------------

	// Input		sys	: Spin system
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with
	//				  only spin system parameters
	// Input		sys	: Spin system
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with spin system
	//				  parameters added to it


spin_system::operator ParameterSet( ) const
  {
  ParameterSet pset;			// A Null parameter set
  pset += *this; 		 	// Add in spin system parameters
  return pset;
  }

void operator+= (ParameterSet& pset, const spin_system &sys){sys.PSetAdd(pset);}

void spin_system::PSetAdd(ParameterSet& pset, int idx) const

        // Input                ss      : A spin system
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : Spin system parameters are
        //                                are added to the parameter set
        //                                with interaction prefix idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting spin systems
        //                                from parameters sets

// sosi needs e- adjustment
  {
  spin_sys::PSetAdd(pset, idx);			// Add base system parameters
  std::string prefx;                                 // Parameter prefix
  if(idx != -1)                                 // Only use suffix if idx
    prefx = std::string("[")+Gdec(idx)+std::string("]");   // is NOT -1
  std::string pname;
  std::string pstate;
  SinglePar par;

  int ns = spins();
  pstate = std::string("Chemical Shift in Hz");	// Add chemical shifts
  int i=0;					// Only done for nuclei
  for(i=0; i<ns; i++)
    {
    if(!electron(i))
      {
      pname = prefx + std::string("v(");
      pname += Gdec(i);
      pname += std::string(")");
      par = SinglePar(pname,shift(i),pstate);
      pset.push_back(par);
      }
    }

  pstate = std::string("G Factor (unitless)");	// Add g factors
  for(i=0; i<ns; i++)				// Only done for electrons
    {
    if(electron(i))
      {
      pname = prefx + std::string("g(");
      pname += Gdec(i);
      pname += std::string(")");
// sosi - need access to gfactor still
//      par = SinglePar(pname,shift(i),pstate);
      }
    }
 
  pstate = std::string("Coupling Constants in Hz");	// Add coupling constants
  for(i=0; i < ns-1; i++)			// Only done for nuclear
    {						// or electron spin pairs
    for(int j=i+1; j<ns; j++)
      {
      if(!enpair(i,j))
        {
        pname = prefx + "J(";
        pname += Gdec(i);
        pname += ",";
        pname += Gdec(j);
        pname += ")";
        par = SinglePar(pname,J(i,j),pstate);
        pset.push_back(par);
        }
      }
    }
 
  pstate = std::string("Hyperfine Coupling in Gauss");// Add HF coupling constants
  for(i=0; i < ns-1; i++)			// Only done for e-/nucleus
    {						// spin pairs
    for(int j=i+1; j<ns; j++)
      {
      if(enpair(i,j))
        {
        pname = prefx + "A(";
        pname += Gdec(i);
        pname += ",";
        pname += Gdec(j);
        pname += ")";
        pname += ")";
        par = SinglePar(pname,A(i,j),pstate);
        pset.push_back(par);
        }
      }
    }

// sosi - need to adjust for electron omega output
  pstate = "Spectrometer Frequency in MHz (1H based)";
  pname = prefx + "Omega";
  par = SinglePar(pname,spectrometer_frequency(),pstate);
  pset.push_back(par);
  return;
  } 

 
// ----------------------------------------------------------------------------
//            Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------
 
void spin_system::setJs(const ParameterSet& pset)
 
        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system scalar couplings
        //                                are set from parameters in pset
	// Note				: This restricts the setting of J(i,j)
	//				  to have j>i for all spin pairs
	// Note				: Hyperfine couplings are NOT read
	//				  in with this function
	// Note				: Error output messages are herein
	//				  controled by warnings() as
	//				     0: no warnings whatsoever
	//				    !0: warn if J involves e-
	//				    >1: also  warn if J absent

  {
  std::string Jst("J(");				// Start of J parameter names
  std::string Jmi(","), Jfi(")");			// Finish of J parameter names
  int notfound = 0;				// Count Js found
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname;					// Use for the parameter name 
  std::string pstate;				// Use for parameter statment
  double Jval = 0;				// User for parameter value
  int ns = spins();				// Number of spins in system
//  std::string NoJSet[ns*ns];				// string of Js not found
  std::string* NoJSet;
  NoJSet = new std::string[ns*ns];
  int i, j=0;
  for(i=0; i<ns-1; i++)				// Loop throught the spin pairs
    { 						// and retrieve the couplings
    for(j=i+1; j<ns; j++)			// All Js need not be present
      {
      pname = Jst+Gdec(i)+Jmi+Gdec(j)+Jfi;	// Parameter name
      item = pset.seek(pname);			// Pix in list for J(i,j)
      if(item != pset.end())
        {
        if(enpair(i,j))				// Insure not e- to nucleus
          {					// as these take hyperfine 
          if(warnings())			// couplings not J couplings
            {
            SYSTerror(41, 1);
            SYSTerror(130, pname);
            }
          Jcouplings[pairidx(i, j)]= 0;		// Set hyperfine coupling
          }					// to 0 for e-/nucleus
        else
          {
          (*item).parse(pname,Jval,pstate);
          J(i, j, Jval);
          }
        }
      else
        {
        if(!enpair(i,j))			// Insure either e- to e-
          { 					// or nucleus to nucleus
          NoJSet[notfound] = pname; 		// Store that this J not set
          J(i, j, 0.0);				// Set to default 0
          notfound++;				// Add to count of no J spins
          }
        }
      }
    }

//	    Now Output Warnings If Any Spin Pairs Not J Assigned 
//	           (Only Done if Warnings Are Active > 1)

  if(!notfound) { delete [] NoJSet; return; }	// Exit if all J's set
  if(warnings() < 2) { delete [] NoJSet; return; } 	// Exit if no warnings desired
  for(i=0,j=0; i<notfound; i++,j++)		// Loop over Js not set
    {						// and issue warnings
    if(!j) SYSTerror(102, NoJSet[i], 1); 	// (five J labels per line)
    else
      {
      std::cout << "; " << NoJSet[i];
      if(j>=2) j=-1;
      }
    }
  SYSTerror(35); 				// Warning, Unspecified Js
  delete [] NoJSet;
  } 


// sosi - needs cleanup 
void spin_system::setAs(const ParameterSet& pset)
 
        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system hyperfine couplings
        //                                are set from parameters in pset
	// Note				: This restricts the setting of A(i,j)
	//				  to have j>i for all spin pairs
	// Note				: Scalar couplings are NOT read
	//				  in with this function

  {
  std::string Ast = "A(";				// Start of A parameter names
  std::string AMst = "AMHz(";			// Alternate start of A p names
  std::string Ami = ",", Afi = ")";                  // Finish of A parameter names
  int notfound = 0;				// Count As found
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname;					// Use for the parameter name 
  std::string pstate;				// Use for parameter statment
  double Aval;					// User for parameter value
  int ns = spins();				// Number of spins in system
  std::string* NoASet;
  NoASet = new std::string[ns*ns];			// string of As not found
  int i, j, onee;
  for(i=0; i<ns-1; i++)				// Loop through the spin pairs
    { 						// and retrieve the couplings
    if(electron(i)) onee=1;			// Flag if first spin is an
    else            onee=0;			// electron
    for(j=i+1; j<ns; j++)			// All Js need not be present
      {
      if(onee && electron(j));			// Insure not both electrons
      else if(!onee && !electron(j));		// Insure not both nuclei
      else
        {
        pname = Ast+Gdec(i)+Ami+Gdec(j)+Afi;	// Parameter name (Gauss)
        item = pset.seek(pname);		// Pix in list for A(i,j)
        if(item != pset.end())
          {
          (*item).parse(pname,Aval,pstate);
          A(i, j, Aval);
          }
        else
          {
          pname = AMst+Gdec(i)+Ami+Gdec(j)+Afi;	// Parameter name (MHz)
          item = pset.seek(pname);		// Pix in list AMHz(i,j)
          if(item != pset.end())
            {
            (*item).parse(pname,Aval,pstate);
            A(i, j, Aval*1.e6*HZ2GAUSSx);
            }
          else
            {
            NoASet[notfound] = pname; 		// Store that this A not set
            A(i, j, 0.0);			// Set to default 0
            notfound++;				// Add to notfound of no A spins
            }
          }
        }
      }
    }

  if(!notfound) { delete [] NoASet; return; }	// Exit if all have been set
  if(!warnings()) { delete [] NoASet; return; }	// Exit if warnings off
  for(i=0, j=0; i<notfound; i++,j++)		// Loop over As not set
    { 						// and issue warnings about
    if(!j) SYSTerror(106, NoASet[i], 1); 	// all A's not set (5 per line)
    else
      {
      std::cout << "; " << NoASet[i];
      if(j>=2) j=0;
      }
    }
  SYSTerror(9); 				// Warning, Unspecified As
  delete [] NoASet;
  return;
  } 


        // Input                ssys    : A spin system (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system isotropic g-factors
        //                                are set from parameters in pset
	//				: A field strength must have been
	//				  specified before this function

void spin_system::setGs(const ParameterSet& pset)
  {
  int ns = spins();			// Number of spins
  int notfound = 0;			// Number of g-factors not found
  std::string* NoGSet;
  NoGSet = new std::string[ns];			// string of g factors not found
  std::string pname, pstate;			// For parameter name and statement
  double dval;
  std::string Gst = "g(", Gfi = ")";		// Start/finish g-factor param. names
  SinglePar par;			// Working single parameter
  ParameterSet::const_iterator item;	// A pix into parameter list
  int i,j;
  for(i=0; i<ns; i++)			// Attempt to retrieve the
    { 					// g-factors (unitless)
    if(electron(i))			// Look for g only if e-
      {
      pname = Gst + Gdec(i) + Gfi; 	// Parameter name for g of spin i
      item = pset.seek(pname);		// Pix in parameter list for g(i)
      if(item != pset.end())		// Retrieve the number of spins
        {
        (*item).parse(pname,dval,pstate);
        gfacts[i] = dval;		// Set g-factor (unitless)
        }
      else
        {
        NoGSet[i] = pname;		// Store parameter not found
        notfound++;			// Count parameter not found
        gfacts[i] = GFREE;		// Set default g value
        }
      }
    else gfacts[i] = 0;			// Set nuclear g-factors 0
    }
  if(notfound)
    {
    for(i=0,j=0; i<notfound; i++,j++)	// Loop over g factors not set
      {					// and issue warnings
      if(!j) SYSTerror(107,NoGSet[i],1);// (five G labels per line)
      else
        {
        std::cout << "; " << NoGSet[i];
        if(j>=4) j=0;
        }
      }
    SYSTerror(21);			// Warning, Unspecified Gs set
    }
  delete [] NoGSet;
  }


void spin_system::operator= (const ParameterSet& pset) { setSsys(pset); } 

	// Input		sys      : Spin system (this)
	// 			pset     : A parameter set
	// Output		none	 : Spin system filled with
	//				   parameters n pset
	// Note				 : Three things are gleaned from
	//				   the parameter set from spin_sys
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	// 				   Three more parameters are taken
	//				   from pset to complete spin_system
	//				   4.) The chemical shifts
	//				   5.) The coupling constants
	//				   6.) A spectrometer frequency
	// Note				 : Functions which place a spin_system
	//				   into a parameter set must contain
	//				   add the information read here

// ----------------------------------------------------------------------------
//    Functions To Output Isotropic System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

int spin_system::write(const std::string &filename, int idx, int warn) const 
 
        // Input                ss      : Spin system (this) 
        //                      filename: Output file name
        //                      idx     : Parameter index value used for 
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Spin system is written as a 
        //                                parameter set to file filename 
        ///F_list write                 - Write system to file (as  pset). 

  {
  if(!spins()) return 1;                // Nothing if no spins
  std::ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))            // If file bad then exit
    {
    if(warn)
      {
      SYSTerror(101, filename, 1);	// Problems with file
      if(warn>1) SYSTfatality(5); 	// !Write to parameter file, fatal
      }  
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }
 
        // Input                ss      : Spin system (base)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Spin system is written as a
        //                                parameter set to output filestream
        // Note                         : This depends on function PSetAdd!

int spin_system::write(std::ofstream& ofstr, int idx, int warn) const
  {
  if(!spins()) return 1;                // Nothing if no spins
  if(!ofstr.good())                     // If file bad then exit
    {
    if(warn) SYSTerror(6);		//	Problems with file stream
    if(warn>1) SYSTfatality(7);		//      Can't write anything, fatal
    }
  ParameterSet pset;			// Declare a parameter set
  PSetAdd(pset, idx);                   // Add system to parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      SYSTerror(6,1);			//	Problems with file stream
      SYSTerror(6, 1);		// 	Problems writing to filestream
      if(warn>1) SYSTfatality(7);	//      Can't write anything, fatal
      }
    return 0;
    }  
  return 1;
  }
 
// ____________________________________________________________________________
//                          SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________

	// Input		sys      : Spin system
	// 			filename : Input filename
	// 			pset	 : A parameter set
	// 			idx	 : Parameter index value used for
	//				   prefix [#] in input names
	//			warn	 : Warning output level
	//					0 = no warnings
	//					1 = warnings
	//				       >1 = fatal warnings
	// Output		none	 : Spin system filled with
	//				   parameters read from file or
	//				   from parameters in pset
	// Output		TF	 : Spin system is filled with
	//				   parameters read from file
	//				   TRUE if read is successful
	// Note			 	 : The file should be an ASCII file
	//				   containing recognized sys parameters

int spin_system::read(const std::string& filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Try and read in pset
    {					// If we cannot read the file then
    if(warn)				// we'll issue warnings as desired
      {
      SYSTerror(1,filename,1);		//	Problems with file
      if(warn>1) SYSTfatality(4);	//	This is a fatal problem
      else       SYSTerror(4,1);	//	Or maybe it aint so bad
      }
    return 0;
    }
  return read(pset, idx, warn);
  }

int spin_system::read(const ParameterSet& pset, int idx, int warn)
  {
  int TF = setSsys(pset, idx, warn?1:0);	// Use overload to read
  if(!TF)					// If setSsys didn't handle
    {						// the system read from pset
    if(warn)					// then we'll issue some
      {
      SYSTerror(8, 1);				//    Problems with pset
      if(warn>1) SYSTfatality(4);		//    This is a fatal problem
      if(warn>1) SYSTerror(4);			//    Or maybe it isn't so bad..
      }
    return 0;
    }
  return TF;
  }

	// Input		sys     : Spin system (this)
	//			argc	: Number of arguments
	//			argv    : Vecotr of argc arguments
	//			argn    : Argument index
	// Output		string  : The parameter argn of array argc
	//				  is used to supply a filename
	//				  from which the spin system is read
	//				  If the argument argn is not in argv,
	//				  the user is asked to supply a filename
	//				  The set filename is returned
	// Note			 	: The file should be an ASCII file
	//				  containing recognized sys parameters
	// Note			 	: The spin system is modifed (filled)

std::string spin_system::ask_read(int argc, char* argv[], int argn)
  {
  std::string filename;				// Name of spin system file  
  query_parameter(argc, argv, argn,		// Get filename from command
       "\n\tSpin system filename? ", filename);	// Or ask for it
  read(filename);		           	// Read system from filename
  return filename;
  }

std::string spin_system::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string msg = "\n\tSpin system filename ["     // Query we will ask if
             + def + "]? ";                     // it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename);		           	// Read system from filename
  return filename;                              // Return filename
  }

// ____________________________________________________________________________
// I                        STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                sys	: A spin system (this)
        //                      ostr    : Output stream
	// Note (printvs)		: Printed in both Hz and PPM if
	//				  the field strength is set
	// Note (printO)		: All base freqencies are written
	//				  as positive values even though
	//				  spins having a negative gamma
	//				  will have a negative Omega

/*  Function  Arguments                           Action
    --------  ---------  -------------------------------------------------------
    printvs   sys, ostr  Write system chemical shifts (PPM,Hz) to ostream ostr
    printgs   sys, ostr  Write system g-factors (unitless) to ostream ostr
    printJs   sys, ostr  Write system scalar couplings (Hz) to ostream ostr
    printAs   sys, ostr  Write system hyperfine couplings (G) to ostream ostr
    printO    sys, ostr  Write system field strength to ostream ostr
    print     sys, ostr  Write spin system to ostream ostr
    <<        sys, ostr  Write spin system to ostream ostr                     */

std::ostream& spin_system::printvs(std::ostream& ostr) const
  {
  std::vector<std::string> Strs = VStrings(8, 2);

  int i, ns=spins();
  ostr << "\nShifts   :  ";
  for(i=0; i<ns; i++) ostr << Strs[i];
  std::string NA = " -------  ";			// Print if Not Applicable
  if(spectrometer_frequency())
    {
    ostr << "\nPPM      :  ";
    for(i=0; i<ns; i++)
      {
      if(electron(i)) ostr << NA;
      else            ostr << Gform("%8.2f  ", PPM(i));
      }
    }
  return ostr;
  }

std::ostream& spin_system::printGs(std::ostream& ostr) const
  {
  std::string NA = " -------  ";			// Print if Not Applicable
  ostr << "\ng Factors:  ";
  int i, ns=spins();
  for(i=0; i<ns; i++)
    {
    if(!electron(i)) ostr << NA;
    else             ostr << Gform("%8.5f  ", gfacts[i]);
    }
  if(spectrometer_frequency())
    {
    ostr << "\nB (Gauss):  ";
    for(i=0; i<ns; i++)
      {
      if(!electron(i)) ostr << NA;
      else             ostr << Gform("%8.2f  ", efield_lab(i));
      }
    }
  return ostr;
  }

std::ostream& spin_system::printJs(std::ostream& ostr) const
  {
  std::string NA = " -------  ";			// Print if Not Applicable
  std::string EMP = "          ";			// Print if EMPTY
  ostr << "\nJ (Hz)";					// nucleus nucleus spin pairs
  int i, j, ns=spins();
  for(i=0; i<ns-1; i++)
    {
    ostr << "\nSpin " << Gdec("%2d",i)
         << "  :  ";
    for(j=0; j<ns; j++)
      {
      if(i<j)
        {
        if(enpair(i,j))      ostr << NA;
        else if(eepair(i,j)) ostr << NA;
        else                 ostr << Gform("%8.2f  ", J(i,j));
        }
      else   ostr << EMP;
      }
    }
  return ostr;
  }

std::ostream& spin_system::printAs(std::ostream& ostr) const
  {
  std::string NA = " -------  ";			// Print if Not Applicable
  std::string EMP = "          ";			// Print if EMPTY
  ostr << "\nA (Gauss)";
  int i, j, ns=spins();
  for(i=0; i<ns-1; i++)
    {
    ostr << "\nSpin " << Gdec("%2d",i)
         << "  :  ";
    for(j=0; j<ns; j++)
      {
      if(i<j)
        {
        if(enpair(i,j)) ostr << Gform("%8.2f  ", A(i,j));
        else            ostr << NA;
        }
      else   ostr << EMP;
      }
    }
  return ostr;
  }

std::ostream& spin_system::printO(std::ostream& ostr) const
  {
  ostr << "\nOmega    :  ";
  int i, ns=spins();
  double abss;
  for(i=0; i<ns; i++)
    {
    abss = fabs(Omega(i));
    if(abss > 1.e3)        ostr << Gform("%8.2f G", abss*1.e-3);// GHz
    else if(abss >= 1.0)   ostr << Gform("%8.2f M", abss);	// MHz
    else if(abss >= 1.e-3) ostr << Gform("%8.2f K", abss*1.e3);	// KHz
    else                   ostr << Gform("%8.2f H", abss*1.e6);	// Hz
    }
  if(electrons())
    {
    ostr << "\nField    :  ";
    double B = Bo();
    if(B >= 1.e4)     ostr << Gform("%8.2f T", B*1.e-4);
    else if(B < 1.e2) ostr << Gform("%8.2f mG", B*1.e3);
    else              ostr << Gform("%8.2f G", B);
    }
  return ostr;
  }
 
std::ostream& spin_system::print(std::ostream& ostr, bool hdr) const
  {
  if(hdr)					// Write out a header if
    {						// desired
    std::string sn = name();
    std::string hdr("Isotropic Spin System");
    if(sn.length())
       hdr += std::string(" ") + sn;
    ostr << CenterString(hdr);
    if(Omega())
      {
      int ne = electrons();
      int nn = nucleons();
      std::string SBo, SVo;
      if(nn)
        {
        double B = fabs(Bo());					// Bo in Gauss
        if(B >= 1.e7)       SBo = Gform("%8.4f kT", B*1.e-7);
        else if(B >= 1.e6)  SBo = Gform("%8.4f T",  B*1.e-4);
        else if(B >= 1.e5)  SBo = Gform("%8.5f T",  B*1.e-4);
        else if(B >= 1.e4)  SBo = Gform("%8.6f T",  B*1.e-4);
        else if(B >= 1.e3)
          if(!ne)  SBo = Gform("%8.7f T",  B*1.e-4);
          else     SBo = Gform("%8.3f G",  B);
        else if(B >= 1.e2)  SBo = Gform("%8.4f G",  B);
        else if(B >= 1.e1)  SBo = Gform("%8.5f G",  B);
        else if(B >= 1.0)   SBo = Gform("%8.6f G",  B);
        else if(B >= 1.e-1) SBo = Gform("%8.7f G",  B);
        else                SBo = Gform("%8.2f mG", B*1.e3);
        }
      if(ne)
        {
        double F = fabs(Omega("e-"));				// Omega in MHz
        if(F >= 1.e7)       SVo = Gform("%7.2f GHz", F*1.e-3);
        else if(F >= 1.e6)  SVo = Gform("%6.2f GHz", F*1.e-3);
        else if(F >= 1.e5)  SVo = Gform("%5.2f GHz", F*1.e-3);
        else if(F >= 1.e4)  SVo = Gform("%4.2f GHz", F*1.e-3);
        else if(F >= 1.e3)  SVo = Gform("%4.3f GHz", F*1.e-3);
        else if(F >= 1.e2)  SVo = Gform("%4.4f GHz", F*1.e-3);
        else if(F >= 1.e1)  SVo = Gform("%6.3f MHz", F);
        else if(F >= 1.0)   SVo = Gform("%6.4f MHz", F);
        else if(F >= 1.e-1) SVo = Gform("%6.5f MHz", F);
        else                SVo = Gform("%6.3f Hz",  F*1.e3);
        }
      if(nn && !ne)
        {
        hdr = std::string("(Static Bo Field Of ")
            + SBo + std::string(")");
        ostr << "\n" << CenterString(hdr) << "\n";
        }
      else if(!nn && ne)
        {
        hdr = std::string("(Spectrometer Frequency Of ")
            + SVo + std::string(")");
        ostr << "\n" << CenterString(hdr) << "\n";
        }
      else
        {
        hdr = std::string("(Bo Field Of ") + SBo
            + std::string("  Base Frequency Of ") + SVo
            + std::string(")");
        ostr << "\n" << CenterString(hdr) << "\n";
        }
      }
    }

  int ns = spins();				// Get total number of spins
  if(!ns)					// If the system has no spins
    {						// we'll just say so and exit
    ostr << CenterString("Spin System Is Empty ")
         << std::endl;
    return ostr;
    }

  std::vector<std::string> Strs;		// Vector of print strings
  Strs = SYSStrings(10, 12, 1);			// System print strings
  int nr = Strs.size();				// Number of rows to print
  int rl = (Strs[0]).length();			// Length or printed row
  int rs = int((80-rl)/2.0);			// Centering space
  std::string lst("\n");			// Row start
  if(rs>0) lst += std::string(rs, ' ');		// Adjust for centering
  for(int i=0; i<nr; i++)			// Print rows. Rows each have
    ostr << lst << Strs[i];			// different data, cols are
  ostr.flush();					// Flush the stream clean
  return ostr;
  }
 
std::ostream& operator<< (std::ostream& ostr, const spin_system& sys)
  { return sys.print(ostr); }

//-----------------------------------------------------------------------------
//                    Strings Used For Generic Output Functions
//-----------------------------------------------------------------------------

/* These functons return vectors of strings that can be used in functions that
   print out spin system information. The strings are of a specified width so
   that they can easily form nice columns when printed. The value of colwd sets
   the width of the srings returned. Function SYSStrings will return strings
   for printing the entire spin system. Each string in that case will appear
   as
                |<--w1-->|_:w3|<--w2-->|w3|<--w2-->|w3|<--w2-->|......
   e.g.         Isotope    :      1H          13C         2H

   where the column widths have default and minimal values built in. Note 
   this function coincides with the base class function of the same type.
   The difference is that it uses a wider width for the successive columns,
   setting those values in spin_sys as well.                                */

std::vector<std::string> spin_system::SYSStrings(int cw1, int cw2, int csp) const
  {
  std::vector<std::string> StrArray;		// Make a string vector
  std::vector<std::string> StrTemp;		// Temporary string array
  std::string semicol(" :");			// Need semicolon
  std::string colspcr(csp, ' ');		// Space between columns
  int i, sbuff;					// Array index, buffer length
  std::string Stemp;				// Temporary string

  if(cw1 < 10) cw1 = 10;			// Keep 1st column width OK
  if(cw2 < 10) cw2 = 10;			// Keep later col widths OK
  if(csp < 0)  csp = 1;				// Keep col spacer OK

//                    Fill Array With Base Class Strings

  StrArray = spin_sys::SYSStrings(cw1,cw2,csp);	// Get strings for spin sys
  int ns = spins();				// Spins in the system
  int dadp = 2;					// Digits after decimal point

//                       Add Chemical Shift Strings

  if(nucleons())
    {
    StrTemp = VStrings(cw2, dadp);		// Get strings for shifts
    Stemp = std::string("Shifts");		//  Set up the first column
    sbuff = cw1 - Stemp.length();		//  to width cw1
    if(sbuff > 0) Stemp += std::string(sbuff, ' ');
    Stemp += semicol + colspcr;
     for(i=0; i<ns; i++)			//  Set up the next columns
      {						//  using cw2 column widths
      Stemp += StrTemp[i];
      if(i+1 < ns) Stemp += colspcr;
      }
    StrArray.push_back(Stemp);			// Store shift values
    }

//                     Add Chemical Shift in PPM Strings

  if(nucleons() && spectrometer_frequency())
    {
    StrTemp = PPMStrings(cw2);			// Get strings for ppm
    Stemp = std::string(cw1, ' ');		//  Set up the first column
    Stemp += semicol + colspcr;
    for(i=0; i<ns; i++)				//  Set up the next columns
      {						//  using cw2 column widths
      Stemp += StrTemp[i];
      if(i+1 < ns) Stemp += colspcr;
      }
    StrArray.push_back(Stemp);			// Store ppm values
    }

//                     Add Electron G Factor Strings

  if(electrons())
    {
    StrTemp = GFStrings(cw2);			// Get strings for ppm
    Stemp = std::string("g Factors");		//  Set up the first column
    sbuff = cw1 - Stemp.length();		//  to width cw1
    if(sbuff > 0) Stemp += std::string(sbuff, ' ');
    Stemp += semicol + colspcr;
    for(i=0; i<ns; i++)				//  Set up the next columns
      {						//  using cw2 column widths
      Stemp += StrTemp[i];
      if(i+1 < ns) Stemp += colspcr;
      }
    StrArray.push_back(Stemp);			// Store g factor values
    }

//                     Add Electron Resonance Field Strings

  if(electrons() && spectrometer_frequency())
    {
    StrTemp = BeStrings(cw2);			// Get strings for field
    Stemp = std::string("B (Gauss)");		//  Set up the first column
    sbuff = cw1 - Stemp.length();		//  to width cw1
    if(sbuff > 0) Stemp += std::string(sbuff, ' ');
    Stemp += semicol; //+ colspcr;
    for(i=0; i<ns; i++)				//  Set up the next columns
      {						//  using cw2 column widths
      Stemp += StrTemp[i];
      if(i+1 < ns) Stemp += colspcr;
      }
    StrArray.push_back(Stemp);			// Store field values
    }

//                     Add Scalar Coupling Strings

  int nsl = Gdec(ns).length();
  std::string sbase, sfill;
  if(nucleons() > 1)
    {
    StrTemp = JStrings(cw2);			// Get strings for J couplings
    sbase = std::string("Js Spin ");		// Set up the first column
    sbuff = cw1 - sbase.length() - nsl;		// to width cw1
    if(sbuff > 0) sfill += std::string(sbuff, ' ');
    int j, k;
    for(i=0, k=0; i<ns-1; i++)			// Set up the first column
      {						// using cw2 column widths
      Stemp = sbase + Gdec(i, nsl) + sfill
            + semicol + colspcr;
      for(j=0; j<ns; j++, k++)			// Set up the next columns
        {					// using cw2 column widths
        Stemp += StrTemp[k];
        if(j+1 < ns) Stemp += colspcr;
        }
      StrArray.push_back(Stemp);		// Store J values
      }
    }

//                     Add Hyperfine Coupling Strings

  if(nucleons() && electrons())
    {
    StrTemp = AStrings(cw2);			// Get strings for A couplings
    sbase = std::string("As Spin ");		//  Set up the first column
    sbuff = cw1 - sbase.length() - nsl;		//  to width cw1
    if(sbuff > 0) sfill += std::string(sbuff, ' ');
    int j, k;
    for(i=0, k=0; i<ns-1; i++)			//  Set up the first column
      {						//  using cw2 column widths
      Stemp = sbase + Gdec(i, nsl) + sfill
            + semicol + colspcr;
      for(j=0; j<ns; j++, k++)			//  Set up the next columns
        {					//  using cw2 column widths
        Stemp += StrTemp[k];
        if(j+1 < ns) Stemp += colspcr;
        }
      StrArray.push_back(Stemp);		// Store A values
      }
    }

//                     Add Field Strength Strings

  if(Omega())					// Print base spectrometer 	
    {						// frequencies per each isotope
    StrTemp = OmStrings(cw2, 2);		// Get strings for fields
    Stemp = std::string("Omega");		//  Set up the first column
    sbuff = cw1 - Stemp.length();		//  to width cw1
    if(sbuff > 0) Stemp += std::string(sbuff, ' ');
    Stemp += semicol + colspcr;
     for(i=0; i<ns; i++)				//  Set up the next columns
      {						//  using cw2 column widths
      Stemp += StrTemp[i];
      if(i+1 < ns) Stemp += colspcr;
      }
    StrArray.push_back(Stemp);			// Store shift values
    }
  return StrArray;				// Return the string array
  }

std::vector<std::string> spin_system::VStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int ns = spins();					// Spins in the system
  std::string NA(colwd, '-');				// Use this if no shift
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".")  + Gdec(digs)
                  + std::string("f");
  double abss, shft;
  std::string pfmt;         
  for(int i=0; i<ns; i++)				// Loop spins & store v
    {							// strings (centered)
    if(electron(i)) StrArray.push_back(NA);
    else
      {
      abss = fabs(shift(i));
      if(abss > 1.e9)       { shft = shift(i)*1.e-9; pfmt = fmt + " GHz"; }
      else if(abss > 1.e6)  { shft = shift(i)*1.e-6; pfmt = fmt + " MHz"; }
      else if(abss > 1.e3)  { shft = shift(i)*1.e-3; pfmt = fmt + " KHz"; }
      else                  { shft = shift(i)*1.e-9; pfmt = fmt + " Hz "; }
      StrArray.push_back(Gform(pfmt.c_str(),shft));
      }
    }
  return StrArray;					// Return the string array
  }

std::vector<std::string> spin_system::PPMStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int ns = spins();					// Spins in the system
  std::string NA(colwd, '-');				// Use this if no shift
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".") + Gdec(digs)
                  + std::string("f ppm");       
  for(int i=0; i<ns; i++)				// Loop spins & store PPM
    {							// strings (centered)
    if(electron(i))					// These is no PPM for e-
      StrArray.push_back(NA);
    else if(!spectrometer_frequency())			// There is no PPM if no
      StrArray.push_back(NA);				// no Bo field set    
    else
      StrArray.push_back(Gform(fmt.c_str(),PPM(i)));
    }
  return StrArray;					// Return the string array
  }

std::vector<std::string> spin_system::GFStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int ns = spins();					// Spins in the system
  std::string NA(colwd, '-');				// Use this if no g factor
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".") + Gdec(digs)
                  + std::string("f    ");       
  for(int i=0; i<ns; i++)				// Loop spins & store g
    {							// strings (centered)
    if(!electron(i))					// These is no g for nuclei
      StrArray.push_back(NA);
    else
      StrArray.push_back(Gform(fmt.c_str(), gfacts[i]));
    }
  return StrArray;					// Return the string array
  }

std::vector<std::string> spin_system::BeStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int ns = spins();					// Spins in the system
  std::string NA(colwd, ' ');				// Use this if no g factor
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".")  + Gdec(digs)
                  + std::string("f    ");       
  for(int i=0; i<ns; i++)				// Loop spins & store field
    {							// strings (centered)
    if(!electron(i) || !spectrometer_frequency())	// These is no g for nuclei
      StrArray.push_back(NA);
    else
      StrArray.push_back(Gform(fmt.c_str(),efield_lab(i)));
    }
  return StrArray;					// Return the string array
  }

std::vector<std::string> spin_system::JStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int i,j, ns = spins();				// Spins in the system
  std::string NA(colwd, '-');				// Use this if no J coupling
  std::string EMP(colwd, ' ');				// Uset this if i>=j
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".") + Gdec(digs)    
                  + std::string("f Hz "); 
  for(i=0; i<ns-1; i++)					// Loop spins & store Js
    {
    for(j=0; j<ns; j++)
      {
      if(i<j)
        {
        if(enpair(i,j))      StrArray.push_back(NA);
        else if(eepair(i,j)) StrArray.push_back(NA);
        else                 StrArray.push_back(Gform(fmt.c_str(),J(i,j)));
        }
      else                   StrArray.push_back(EMP);
      }
    }
  return StrArray;					// Return the string array
  }

std::vector<std::string> spin_system::AStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int i,j, ns = spins();				// Spins in the system
  std::string NA(colwd, '-');				// Use this if no J coupling
  std::string EMP(colwd, ' ');				// Uset this if i>=j
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".") + Gdec(digs)    
                  + std::string("f G  "); 
  for(i=0; i<ns-1; i++)					// Loop spins & store Js
    {
    for(j=0; j<ns; j++)
      {
      if(i<j)
        {
        if(nnpair(i,j))      StrArray.push_back(NA);
        else if(eepair(i,j)) StrArray.push_back(NA);
        else                 StrArray.push_back(Gform(fmt.c_str(),A(i,j)));
        }
      else                   StrArray.push_back(EMP);
      }
    }
  return StrArray;					// Return the string array
  }

std::vector<std::string> spin_system::OmStrings(int colwd, int digs) const
  {
  std::vector<std::string> StrArray;			// Make a string vector
  int ns = spins();					// Spins in the system
  std::string NA(colwd, '-');				// Use this if no shift
  std::string fmt = std::string("%") + Gdec(colwd-4)	// Set # format string
                  + std::string(".") + Gdec(digs)
                  + std::string("f");
  double abss, Om;
  std::string pfmt;         
  for(int i=0; i<ns; i++)				// Loop spins & store v
    {							// strings (centered)
    abss = fabs(Omega(i));
    if(abss > 1.e3)       { Om = abss*1.e-3; pfmt = fmt + " GHz"; }
    else if(abss > 1.0)   { Om = abss;       pfmt = fmt + " MHz"; }
    else if(abss > 1.e-3) { Om = abss*1.e3;  pfmt = fmt + " KHz"; }
    else                  { Om = abss*1.e6;  pfmt = fmt + " Hz "; }
    StrArray.push_back(Gform(pfmt.c_str(),Om));
    }
  return StrArray;	
  }

#endif						// SpinSystem.cc
