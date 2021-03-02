/* SysGradz.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	Spin System in a z-Gradient                     Interface	**
**                                                                      **
**      Scott Smith                                                     **
**	Copyright (c) 1997						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $                                                      **
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description							 	**
**                                                                      **
** Class sys_gradz embodies an isotropic system which is placed in a	**
** magnetic field having a Bo gradient along the z-axis.		**
**									**
** This class is derived from a base class "spin_system".		**
**									**
*************************************************************************/

#ifndef _sys_gradz_cc_				// Is file already included?
#define _sys_gradz_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)					// Using the GNU compiler?
#    pragma implementation				// This is the implementation
#endif

#include <Gradients/sys_gradz.h>		// Includes the interface 
#include <Basics/Gutils.h>			// Include GAMMA std errors
#include <Basics/Gconstants.h>			// Inlcude GAMMA1H
#include <Basics/StringCut.h>

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS SYS Z-GRADIENT ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   sysg    : System in z-gradient (this)
                                eidx    : Error index 
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message 
        Output                  none    : Error message output 
                                          Program execution stopped if fatal

  The following error messages use the defaults set in the Gutils package

                Case                          Error Message

                (0)                     Program Aborting.....
                (3)                     Can't Construct from Parameter Set
                (4)                     Can't Construct from Input File
                (5)                     Can't Write To Parameter File
                default                 Unknown Error 

                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                (3)                     Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */  


void sys_gradz::SysGZerr(int eidx, int noret) const
  {
  std::string hdr("Z-Gradient System");
  std::string msg;
  switch(eidx)
    {
    case 7: msg = std::string("Requested Sub-System Out of Range!!");
             GAMMAerror(hdr, msg, noret);  break;			// (7)
    case 11: msg = std::string("Non-Positive # of Spin Sub-Systems!");
             GAMMAerror(hdr, msg, noret);  break;			// (11)
    case 12: GAMMAerror(hdr,"No Field Gradient Specified",noret); break;// (12)
    case 15: msg = std::string("Non-Positive Effective Sample Length!");
             GAMMAerror(hdr, msg, noret);  break;			// (15)
    case 19: msg = std::string("Cannot Perform Requested Alteration");
             GAMMAerror(hdr, msg, noret);  break;			// (19)
    case 25: msg = std::string("Effective Sample Length Not Specified");
             GAMMAerror(hdr, msg, noret);  break;			// (25)
    case 26: msg = std::string("Cannot Provide Sub-System Effective Distance");
             GAMMAerror(hdr, msg, noret);  break;			// (26)
    case 31: msg = std::string("One Sub-System? Use Class spin_system....");
             GAMMAerror(hdr, msg, noret);  break;			// (31)
    case 51: msg = std::string("Number of Sub-Systems Not Specified");
             GAMMAerror(hdr, msg, noret);  break;			// (51)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  if(!noret) std::cout << "\n";
  }  
 
void sys_gradz::SysGZerr(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Z-Gradient System");
  std::string msg;
  switch(eidx)
    {
    case 5: msg = std::string("Bad Use Of ") + pname + std::string(" Function ");
             GAMMAerror(hdr,msg+pname,noret);  break;                  // (5)
    default: GAMMAerror(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  if(!noret) std::cout << "\n";
  }

volatile void sys_gradz::SysGZfatal(int eidx) const
  {
  SysGZerr(eidx, 1);
  if(eidx) SysGZerr(0);
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                 SYS Z-GRADIENT SETUP FUNCTIONS
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   spin system to be set up without worrying about system consistency!       */  

int sys_gradz::setSsys(const ParameterSet& pset, int idx, int warn)
 
	// Input		sysg	: System in z-gradient (this)
        //                      pset     : A parameter set
        //                      idx      : Parameter index value used for
        //                                 prefix [#] in input names
        //                      warn     : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF       : Spin system is filled with
        //                                 parameters in pset
        // Note                          : Three things are gleaned from
        //                                 the parameter set from spin_sys
        //                                 1.) The number of spins
        //                                 2.) An isotope type for each spin
        //                                 3.) An optional spin system name
        //                                 Three more parameters are taken

  {
//      Use Only Parameters With Prefix [#], So Clip [#] From Names First
 
  ParameterSet  subpset;                        // Working parameter set
  if(idx != -1) subpset=pset.strip(idx);        // Get only params with [#]
  else          subpset=pset;                   // Or use full pset
 
//                Now Set Up The Spin System Directly
 
  int TF=1;                                     // Track if we read OK
  int ns = getSpins(subpset, warn?1:0);         // Get the number of spins
  if(ns<=0)
    {
    if(warn) SysGZerr(13,1);
    else     SysGZerr(13);
    return 0;
    }
  *this = sys_gradz(ns);			// Set the system for ns spins
  setIs(subpset);                               // Set the isotope types
  setName(subpset);                             // Read in the system name
  setBasis(matrix(HS(), HS(), i_matrix_type));  // Set up default basis (matrix)
  setOm(subpset);                               // Set spectrometer frequency
  setShifts(subpset);                           // Set the chemical shifts
  setJs(subpset);                               // Read in scalar couplings
  if(electrons())                               // Read G's & J's if electrons
    {                                           // present in the spin system
    setGs(subpset);                             //      Set g-factors
    setAs(subpset);                             //      Set hyperfine couplings
    }
  setSubSys(pset);				// Set the # of sub-systems
  setBoGrad(pset);				// Set the field gradient
  setLength(pset);				// Set the effective length 
  return TF;
  }  

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  SYS Z-GRADIENT CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/*  Arguments                         Constructed Spin System
    ---------           -------------------------------------------------------
      spins              An empty system with "spins" spins in a z-gradient
       sys               A duplicate of spin system sys
 
   Note that these constructors don't do much at all.  Systems in gradients
   need not only the parameters for a base spin system to describe them, they
   also need parameters for the gradient strength, number of sub-systems, etc.
   This being too much information to feed into the constructor, we do mimimal
   things in them. Rather, spin systems are usually read in from an external
   file where all their parameters are clearly written out.                  */


sys_gradz::sys_gradz(int spins) : spin_system(spins)
  {
  _NSS   = 0;					// No subsystems
  dBodm  = 0;					// No field gradient
  efflen = 0;					// No effective sample length
  }

sys_gradz::sys_gradz(const sys_gradz &sys) : spin_system(sys)
  {
  _NSS   = sys._NSS;				// Copy # of subsystems
  dBodm  = sys.dBodm;				// Copy the field gradient
  efflen = sys.efflen;				// Copy sample length
  Ps     = sys.Ps;  				// Copy subsystem populations
  }

sys_gradz& sys_gradz::operator= (const sys_gradz &sys)
  {
  spin_system::operator=(sys);			// Copy spin system stuff
  _NSS   = sys._NSS;				// Copy # of subsystems
  dBodm  = sys.dBodm;				// Copy the field gradient
  efflen = sys.efflen;				// Copy detected sample length
  Ps     = sys.Ps;  				// Copy subsystem populations
  return *this;
  }

sys_gradz::~sys_gradz () { }

// ____________________________________________________________________________
// B               CHEMICAL SHIFT MANIPULATION FUNCTIONS
// ____________________________________________________________________________

//	     ALL Handled By The Base Class, Class spin_system

// ______________________________________________________________________
// C               COUPLING CONSTANT MANIPULATION FUNCTIONS
// ______________________________________________________________________

//	     ALL Handled By The Base Class, Class spin_system

// ____________________________________________________________________________
// D                  SPECTROMETER FREQUENCY MANIPULATIONS
// ____________________________________________________________________________

//	     ALL Handled By The Base Class, Class spin_system

// ____________________________________________________________________________
// E                         SPIN PAIR FLAG FUNCTIONS
// ____________________________________________________________________________

//	     ALL Handled By The Base Class, Class spin_system

// ____________________________________________________________________________
// F                      SYS Z-GRADIENT ASSOCIATED FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                Functions To Access Number of Spin Sub-Systems
// ----------------------------------------------------------------------------

/*              Function     Arguments                   Result
                --------   -------------     ----------------------------------
                  NSS          int           Sets the number of sub-systems
                  NSS          ---           Returns the # of sub-systems    */

int  sys_gradz::NSS() const { return _NSS; }
void sys_gradz::NSS(int nss)
  {
  if(nss<=0) 
    {
    SysGZerr(10, 1);				// Non-positive NSS
    SysGZfatal(19);	 			// Can't alter system
    }
  if(nss==1) SysGZerr(31);			// Non-positive NSS
  _NSS = nss;					// Set new # of spins
  Ps.empty();					// Populations now empty
  }


// ----------------------------------------------------------------------------
//                 Functions To Access The Bo Field Gradient
// ----------------------------------------------------------------------------

/*              Function     Arguments                   Result
                --------   -------------     -------------------------------
                 BoGrad    double (T/m)      Gradient is set in T/m
                 BoGrad         ---          Gradient is returned in T/m
                 GradVal   double (m)        Gradient at specfified distance */
 
void   sys_gradz::BoGrad(double bgrad)       { dBodm = bgrad; }
double sys_gradz::BoGrad()             const { return dBodm; }
double sys_gradz::GradVal(double dist) const { return dBodm*dist*1.e4; }


// ----------------------------------------------------------------------------
//              Functions To Access The Effective Sample Length
// ----------------------------------------------------------------------------
         
/*              Function     Arguments                   Result
                --------   -------------     ------------------------------- 
                 SysLen    double (T/m)      Gradient is set in T/m
                 SysLen        ---           Effective sample length in m   
                 SysDist       int           Distance of sub-system in m 

  Note that in GAMMA the 1st sub-system (index 0) resides at -efflen/2
  whereas the last sub-system (index _NSS-1) resides at efflen/2             */
 
void sys_gradz::SysLen(double len)
  {
  if(len<=0) 
    {
    SysGZerr(15, 1);			// Non-positive sample length
    SysGZfatal(19);	 		// Cannot alter system
    }
  efflen = len;
  }

double sys_gradz::SysLen()         const { return efflen; }
double sys_gradz::SysDist(int nss) const
  {
  if(nss<0 || nss>=_NSS)
    {
    SysGZerr(7, 1);			// Requested sub-system out of range
    SysGZfatal(26);			// Cannot provide distance info
    }
  if(_NSS == 1) return 0.0;
  return(efflen/double(_NSS-1))*double(nss) - efflen/2.0;
  }


// ----------------------------------------------------------------------------
//                Functions To Access Particular Spin Sub-Systems
// ----------------------------------------------------------------------------
 
/*   Function       Arguments                      Result
    -----------   -------------     -------------------------------------------
    SubSys             int          Returns the requested sub-system
    SubSysShift     int, int        Shift (Hz) of spin in specified sub-system
    SubSysShift    double, int      Shift (Hz) of spin @ specified distance (m)
    SubSysPPM       int, int        Shift(PPM) of spin in specified sub-system
    SubSysPPM      double, int      Shift(PPM) of spin @ specified distance (m)
 
  Note that in GAMMA the ith sub-system is viewed as a base system having
  modified shifts.  From this perspective it is assumed that the user wants
  all sub-system referenced to the same rotating frame, namely the base
  Larmor frequency of each isotope type.                                    */
 
spin_system sys_gradz::SubSys(int nss) const
  {
  spin_system sys(*this);			// Copy sysg into spin-system
  double Bmod = dBodm*SysDist(nss)*RAD2HZ;	// Gradient at nss (T*Hz/rad)
  double newshift;
  for(int i=0; i<spins(); i++)			// Loop over the spins
    {
    newshift = sys.shift(i);			//	Current shift
    newshift += Bmod*gamma(i);			//	Gradient modified
    sys.shift(i, newshift);			//	Set modified shift
    }
  return sys;
  }
 
double sys_gradz::SubSysShift(int nss, int spin) const
  { return SubSysShift(SysDist(nss), spin); }	// User overload

 
double sys_gradz::SubSysShift(double dist, int spin) const
  {
  double Bmod = dBodm*dist*RAD2HZ;		// Gradient at dist (T*Hz/rad)
  return shift(spin)+ Bmod*gamma(spin);		// Gradient modified shift (Hz)
  }

double sys_gradz::SubSysPPM(int nss, int spin) const
  { return SubSysPPM(SysDist(nss), spin); }	// User overload

 
double sys_gradz::SubSysPPM(double dist, int spin) const
  {
  check_spin(spin);                             // Insure spin exists
  if(Omega("1H") == 0.0) SysGZerr(1);		// Fatal if no Bo Set
  double sHz = SubSysShift(dist, spin);		// Shift in Hz
  return sHz*GAMMA1H/(Omega("1H")*gamma(spin));	// Return shift in PPM
  }


// ____________________________________________________________________________
// G	                     Nyquist Frequency Functions
// ____________________________________________________________________________


// ____________________________________________________________________________
// G	                     Nyquist Frequency Functions
// ____________________________________________________________________________

//	     ALL Handled By The Base Class, Class spin_system

// ____________________________________________________________________________
// H                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

//	   Mostly Handled By The Base Class, Class spin_system


// ----------------------------------------------------------------------------
//            Functions To Make A Parameter Set From A Spin System
// ----------------------------------------------------------------------------

sys_gradz::operator ParameterSet( ) const

	// Input		sysg	: System in z-gradient (this)
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with
	//				  only spin system parameters

  {
  ParameterSet pset;                    // A Null parameter set
  pset += *this;                        // Add in spin system parameters
  return pset;
  }

void operator+= (ParameterSet& pset, const sys_gradz &sysg) 

	// Input		pset	: A parameter set
	// 			sysg	: System in z-gradient (this)
	// Output		pset	: The parameter set with
	//			          only sysg parameters

{ sysg.PSetAdd(pset); }

 
void sys_gradz::PSetAdd(ParameterSet& pset, int idx) const
 
	// Input		sysg	: System in z-gradient (this)
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
 
  {
  spin_system::PSetAdd(pset, idx);		// Add base system parameters

  std::string prefx;                                 // Parameter prefix
  if(idx != -1)                                 // Only use suffix if idx
    prefx = std::string("[")+Gdec(idx)+std::string("]");   // is NOT -1
  std::string pname;
  std::string pstate;
  SinglePar par;

  pstate = std::string("Number of Subsystems");	// Add # of subsystems
  pname = std::string("NSubSys");
  par = SinglePar(pname,_NSS,pstate);
  pset.push_back(par);

  pstate = std::string("Field Gradient (T/m)");	// Add field gradient
  pname = std::string("BoGrad");
  par = SinglePar(pname, dBodm, pstate);
  pset.push_back(par);

  pstate = std::string("Effective Sample Size (m)");	// Add sample length
  pname = std::string("SysLen");
  par = SinglePar(pname, efflen, pstate);
  pset.push_back(par);
  return;
  } 

 
// ----------------------------------------------------------------------------
//            Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------
 
//          Many of these handled by the base class spin_system
 
// int sys_gradz::setOm(ParameterSet& pset);
// void sys_gradz::setJs(ParameterSet& pset);
// void sys_gradz::setAs(ParameterSet& pset);
// void sys_gradz::setShifts(ParameterSet& pset);


 
	// Input		sysg	: System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system # of sub-systems 
        //                                is set from parameter in pset

void sys_gradz::setSubSys(const ParameterSet& pset)
  {
  int ival;				// For # of sub-systems
  std::string pstate, pname="NSubSys";	// For parameter name and statement
  ParameterSet::const_iterator item;	// A pix into parameter list
  item = pset.seek(pname);
  if(item == pset.end())
    {
    SysGZerr(2, pname, 1);		// Can't read NSubSys
    SysGZerr(51, 1);			// No sub-systems defined
    SysGZfatal(3);	 		// Can't make system from pset
    }
  (*item).parse(pname,ival,pstate);	// Get # of sub-systems
  _NSS = ival;				// Set # of sub_systems
  } 


 
	// Input		sysg	: System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system field gradient
        //                                is set from parameter in pset
	// Note				: Input units expected T/m

void sys_gradz::setBoGrad(const ParameterSet& pset)
  {
  double dval;				// For # of sub-systems
  std::string pstate, pname="BoGrad";	// For parameter name and statement
  ParameterSet::const_iterator item;	// A pix into parameter list
  item = pset.seek(pname);		// Try and find parameter
  if(item == pset.end())
    {
    SysGZerr(2, pname, 1);		// Can't read BoGrad
    SysGZerr(12, 1);	 		// No field gradient defined
    SysGZfatal(3);	 		// Can't make system from pset
    }
  (*item).parse(pname,dval,pstate);	// Get # of field gradient
  dBodm = dval;				// Set # of field gradient
  } 


void sys_gradz::setLength(const ParameterSet& pset)
 
	// Input		sysg	: System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system effective length
        //                                is set from parameter in pset
	// Note				: Input units normally um

  {
  double dval;					// For # of sub-systems
  std::string pstate, pname0 = "SysLen";		// For name and statement
  std::string pname = pname0;			// SysLen assumes um
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Try and find parameter
  if(item != pset.end())
    {
    (*item).parse(pname,dval,pstate);		// Get # of effective length(um)
    efflen = dval*1.e-6;			// Set # of effective length(m)
    return;
    }
  else
    {
    pname += "m";				// SysLenm assumes mm
    item = pset.seek(pname);			// Try and find parameter
    if(item != pset.end())
      {
      (*item).parse(pname,dval,pstate);		// Get # of effective length(mm)
      efflen = dval*1.e-3;			// Set # of effective length(m)
      return;
      }
    }
  SysGZerr(2, pname0, 1);			// Can't read SysLen
  SysGZerr(25, 1);	 			// No effective length defined
  SysGZfatal(3);	 			// Can't make system from pset
  return;
  } 


void sys_gradz::operator= (const ParameterSet& pset)

	// Input		sysg	: System in z-gradient (this)
	// 			pset	: A parameter set
	// Output		none	: Spin system filled with
	//				  parameters in pset
	// Note			 	: Six things are gleaned from
	//			  	  the parameter set from spin_system
	//				   1.) The number of spins
	//				   2.) An isotope type for each spin
	//				   3.) An optional spin system name
	//				   4.) The chemical shifts
	//				   5.) The coupling constants
	// 				  Three more parameters are taken
	//				  from pset to complete sys_gradz
	//				   6.) The number of sub-systems
	//				   7.) The field gradient
	//				   8.) An effective sample length
	// Note				 : Functions which place a sys_gradz
	//				   into a parameter set must contain
	//				   the information read here

  {
  int ns = spin_sys::getSpins(pset);	// Get the number of spins (spin_sys)
  *this = sys_gradz(ns);		// Set system for ns spins
  spin_sys::setIs(pset);		// Set the isotope types
  spin_system::setOm(pset);		// Set the spectrometer frequency
  spin_system::setShifts(pset);		// Set the chemical shifts
  setJs(pset);                          // Read in scalar couplings
//  if(electrons) setAs(pset);            // Set hyperfine couplings (if e- in sys)
  setSubSys(pset);			// Set the # of sub-systems
  setBoGrad(pset);			// Set the field gradient
  setLength(pset);			// Set the effective sample length 
  } 


// ----------------------------------------------------------------------------
//     Functions To Output Gradient System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------


	// Input		sysg	: System in z-gradient (this)
        //                      filename : Output file name
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	: Spin system is written as a 
	//				  parameter set to file filename

int sys_gradz::write(const std::string &filename, int idx, int warn) const
  {
  if(!spins()) return 1;                // Nothing if no spins
  std::ofstream ofstr(filename.c_str());     // Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))            // If file bad then exit
    {
    if(warn)
      {
      SysGZerr(1, filename, 1);		// Problems with file
      SysGZfatal(5); 			// !Write to parameter file, fatal
      }
    return 0;
    }  
  ofstr.close();                        // Close it now
  return 1;
  }

   

	// Input		sysg	: System in z-gradient (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Spin system is written as a
        //                                parameter set to output filestream
        // Note                         : This depends on function PSetAdd!

int sys_gradz::write(std::ofstream& ofstr, int idx, int warn) const
  {
  if(!spins()) return 1;                // Nothing if no spins
  if(!ofstr.good())			// If file bad then exit
    {
    if(warn) SysGZerr(5,1);		// Problems with file stream
    if(warn) SysGZfatal(5);		// !Write to parameter file, fatal
    }
  ParameterSet pset; 			// Declare a parameter set
  PSetAdd(pset, idx);                   // Add system to parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      SysGZerr(5,1);			//      Problems with file stream
      SysGZerr(5, 1);		//      Problems writing to filestream
      if(warn>1) SysGZfatal(7); 	//      Can't write anything, fatal
      }
    return 0;
    }  
  return 1;
  }  

// sosi - OK below here

// ____________________________________________________________________________
//                          SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________



	// Input		sysg	: System in z-gradient (this)
	// 			filename: Input filename
        //                      idx      : Parameter index value used for
        //                                 prefix [#] in input names
        //                      warn     : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
	// Output		none	: Spin system filled with
	//				  parameters read from file
        //                                 TRUE if read is successful
        // Note                          : The file should be an ASCII file
        //                                 containing recognized sys parameters

int sys_gradz::read(const std::string &filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))    // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      SysGZerr(1, filename, 1);		// 	Problems with file
      if(warn>1) SysGZfatal(4);		//      This is a fatal problem
      else       SysGZerr(4,1);		//      Or maybe it aint so bad
      }
    return 0;
    }  
  return read(pset, idx, warn);
  }


 
	// Input		sysg	: System in z-gradient (this)
        //                      pset     : A parameter set
        //                      idx      : Parameter index value used for
        //                                 prefix [#] in input names
        //                      warn     : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF       : Spin system is filled with
        //                                 parameters in pset
        //                                 TRUE if system filled properly

int sys_gradz::read(const ParameterSet& pset, int idx, int warn)
  {
  int TF = setSsys(pset, idx, warn?1:0);        // Use overload to read
  if(!TF)                                       // If setSsys didn't handle
    {                                           // the system read from pset
    if(warn)                                    // then we'll issue some
      {
      SysGZerr(8, 1);				//    Problems with pset
      if(warn>1) SysGZfatal(4);			//    This is a fatal problem
      if(warn>1) SysGZerr(4);			//    Or maybe it isn't so bad..
      }
    return 0;
    }  
  return TF;
  }  




	// Input		sysg	: System in z-gradient (this)
	//			argc	: Number of arguments
	//			argv    : Vecotr of argc arguments
	//			argn    : Argument index
	// Output		string  : The parameter argn of array argc
	//				  is used to supply a filename
	//				  from which the spin system is read
	//				  If the argument argn is not in argv,
	//				  the user is asked to supply a filename
	//				  The set filename is returned
        //
	// Note			 	: The file should be an ASCII file
	//				  containing recognized sys parameters
	// Note			 	: The spin system is modifed (filled)

std::string sys_gradz::ask_read(int argc, char* argv[], int argn)
  {
  std::string filename;				// Name of spin system file  
  query_parameter(argc, argv, argn,		// Get filename from command
    "\n\tSpin system filename? ", filename);	// Or ask for it
  read(filename);		           	// Read system from filename
  return filename;
  }

std::string sys_gradz::ask_read(int argc, char* argv[], int argn, const std::string& def)
  {
  std::string msg = "\n\tSpin system filename [" + def + "]? ";	// Query we will ask if
																// it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename);                               // Read system from filename
  return filename;                              // Return filename
  }


// ____________________________________________________________________________
// I                        STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

	// Input		sysg	: System in z-gradient (this)
        //                      ostr    : Output stream
        // Output               none    : Spin system info is sent
        //                                into the output stream

std::ostream& sys_gradz::print(std::ostream& out) const
  {
  out << "\n";					// Start with a line feed 
  spin_system::print(out);			// The use base class to print
  out << "\nMax Shift:";			// Now output shift values in
  double ssh, abssh;				// at the maximum and minimum
  int i, ns=spins();				// gradient values
  for(i=0; i<ns; i++)
    {  
    ssh = SubSysShift(_NSS-1, i);
    abssh = fabs(ssh);
    if(symbol(i) == "e-")  out << "NA";
    else if(abssh > 1.e9)  out << Gform("%10.2f G", ssh*1.e-9);	// GHz
    else if(abssh > 1.e6)  out << Gform("%10.2f M", ssh*1.e-6);	// MHz
    else if(abssh > 1.e3)  out << Gform("%10.2f K", ssh*1.e-3);	// KHz
    else                   out << Gform("%10.2f", ssh);		// Hz
    }  
  out << "\nMin Shift:";
  for(i=0; i<ns; i++)
    {  
    ssh = SubSysShift(0, i);
    abssh = fabs(ssh);
    if(symbol(i) == "e-")  out << "NA";
    else if(abssh > 1.e9)  out << Gform("%10.2f G", ssh*1.e-9);     // GHz
    else if(abssh > 1.e6)  out << Gform("%10.2f M", ssh*1.e-6);     // MHz
    else if(abssh > 1.e3)  out << Gform("%10.2f K", ssh*1.e-3);     // KHz
    else                   out << Gform("%10.2f", ssh);		    // Hz
    }  
  if(spectrometer_frequency())
    {  
    out << "\nMax PPM  :";
    for(i=0; i<ns; i++)
      out << Gform("%10.2f", SubSysPPM(_NSS-1, i));
    out << "\nMin PPM  :";
    for(i=0; i<ns; i++)
      out << Gform("%10.2f", SubSysPPM(0, i));
    }
  out << "\nMax Dist.:";
    {
    double dist = SysDist(_NSS-1);
    if(dist < 1.e-3)   out << Gform("%10.2f um", dist*1.e6);	// um
    else if(dist < 1)  out << Gform("%10.2f mm", dist*1.e3);	// mm
    else               out << Gform("%10.2f m",  dist);		// m
    }
  out << "\nField    :"
      << Gform("%10.2f T", Bo()*1.e-4);
  out << "\nSubSys # :" << Gform("%7i", _NSS);
  out << "\nGradient :"
      << Gform("%10.2f T/m", dBodm);
  out << "\nMax Grad :"
      << Gform("%10.2f G", GradVal(SysDist(_NSS-1)));
  out << "\nMin Grad :"
      << Gform("%10.2f G", GradVal(SysDist(0)));
  if(_NSS > 1)
    {
    out << "\nLength   :";
    if(efflen < 1.e-3)
      out << Gform("%10.2f um", efflen*1.e6);
    else if(efflen < 1.0)
      out << Gform("%10.2f mm", efflen*1.e3);
    else
      out << Gform("%10.2f m", efflen);
    }
  out << "\n";
  return out;
  }

std::ostream& operator<< (std::ostream& ostr, const sys_gradz& sys)
  { return sys.print(ostr); }


    
#endif							// SysGradz.cc
