/* sys_gradz.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Spin System in a z-Gradient                Implementation       **
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 1997                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** Class sys_gradz embodies an isotropic system which is placed in a    **
** magnetic field having a Bo gradient along the z-axis.  The system	**
** is assumed to represent N subsystems each of which exist in a set	**
** Bo field.  Each subsystem will be an isotropic system of the type	**
** "spin_sytem", the class from which this one is derived.              **
**                                                                      **
** Each subsystem contains the same sys_gradz information: # spins,     **
** their chemical shifts, scalar coupling constants, isotope types,     **
** spin flags, a field strength, spin pair flags, a default isotope     **
** type and a spectrometer frequency.                                   **
**                                                                      **
** In addition, the sys_gradz contains the following values which are	**
** used to distinguish systems experiencing Bo fields adjusted (from	**
** the standard sys_gradz Bo) by the z-gradient:                        **
**                                                                      **
** int N        - Number of spin systems with distinguishable Bo's      **
** double dBodm	- Gradient strength in T/m                              **
** double len	- Length of sample being detected.                      **
** double* Ps   - Populations of individual spin systems                **
**                                                                      **
*************************************************************************/

///Chapter Class System in z-Gradient
///Section Overview
///Body    None
///Section Available Spin System Functions

#ifndef   SysGradz_h_				// Is file already included?
#  define SysGradz_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// this is the interface
#  endif

#include <GamGen.h>					// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Includes parameter sets
#include <HSLib/SpinSystem.h>		// Includes spin system base class
#include <vector>					// Includes libstdc++ vectors

class sys_gradz: public spin_system
  {
  int                 _NSS;			// Number of subsystems
  double              dBodm;		// Field gradient (T/m)
  double              efflen;		// Effective sample length (m)
  std::vector<double> Ps;			// Sub-system populations

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                   CLASS SYS Z-GRADIENT ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   sysg    : System in z-gradient (this)
                                eidx    : Error index
                                nr      : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void SysGZerr (int  eidx,                        int nr=0) const;
         void SysGZerr(int   eidx, const std::string& st, int nr=0) const;
volatile void SysGZfatal(int eidx) const;

// ____________________________________________________________________________
// ii                     SYS Z-GRADIENT SETUP FUNCTIONS
// ____________________________________________________________________________

/* These are protected functions because they allow specific aspects of the
   spin system to be set up without worrying about system consistency!       */

virtual int setSsys(const ParameterSet& pset,int idx=-1,int wrn=2);

        // Input                sysg    : System in z-gradient (this)
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

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A          SYS Z-GRADIENT CONSTRUCTION, DESTRUCTION, ASSIGNMENT
// ____________________________________________________________________________

///Center Spin System Algebraic
///F_list sys_gradz    - Constructor
///F_list =            - Assignment

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

MSVCDLC      sys_gradz(int spins=0);
MSVCDLC      sys_gradz(const sys_gradz& sys);
MSVCDLC      ~sys_gradz();
MSVCDLL sys_gradz& operator= (const sys_gradz& sys);

// ____________________________________________________________________________
// B                       CHEMICAL SHIFT MANIPULATIONS
// ____________________________________________________________________________

//                 ALL handled by the base class spin_system

// ____________________________________________________________________________
// C           SCALAR & HYPERFINE COUPLING CONSTANT MANIPULATIONS
// ____________________________________________________________________________

//                 ALL handled by the base class spin_system

// ____________________________________________________________________________
// D                    SPECTROMETER FREQUENCY MANIPULATIONS
// ____________________________________________________________________________

//                 ALL handled by the base class spin_system

// ____________________________________________________________________________
// E                      SPIN PAIR FLAG FUNCTIONS
// ____________________________________________________________________________

//                 ALL handled by the base class spin_system

// ____________________________________________________________________________
// F                   SYS Z-GRADIENT ASSOCIATED FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                Functions To Access Number of Spin Sub-Systems
// ----------------------------------------------------------------------------

/*              Function     Arguments                   Result
                --------   -------------     ----------------------------------
                  NSS          int           Sets the number of sub-systems
                  NSS          ---           Returns the # of sub-systems    */

MSVCDLL void NSS(int nss);
MSVCDLL int  NSS() const;

// ----------------------------------------------------------------------------
//                 Functions To Access The Bo Field Gradient
// ----------------------------------------------------------------------------

/*              Function     Arguments                   Result
                --------   -------------     -------------------------------
                 BoGrad    double (T/m)      Gradient is set in T/m
                 BoGrad         ---          Gradient is returned in T/m
                 GradVal   double (m)        Gradient at specfified distance */

MSVCDLL void   BoGrad(double bgrad);
MSVCDLL double BoGrad()             const;
MSVCDLL double GradVal(double dist) const;


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

MSVCDLL void   SysLen(double len);
MSVCDLL double SysLen()         const;
MSVCDLL double SysDist(int nss) const;

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

MSVCDLL spin_system SubSys(int nss)                    const;
MSVCDLL double      SubSysShift(int nss, int spin)     const;
MSVCDLL double      SubSysShift(double dist, int spin) const;
MSVCDLL double      SubSysPPM(int nss, int spin)       const;
MSVCDLL double      SubSysPPM(double dist, int spin)   const;

// ____________________________________________________________________________
// G			  Nyquist Frequency Functions
// ____________________________________________________________________________

//                 ALL handled by the base class spin_system

// ____________________________________________________________________________
// H                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

//                 Mostly handled by the base class spin_system

///Center Spin System Parameter Set Functions

// ----------------------------------------------------------------------------
//            Functions To Make A Parameter Set From A Spin System
// ----------------------------------------------------------------------------

MSVCDLL operator ParameterSet( ) const;


        // Input                sysg    : System in z-gradient (this)
        // Output               pset    : A parameter set with
        //                                only spin system parameters
        // F_list =  - Conversion

MSVCDLL friend void operator+= (ParameterSet& pset, const sys_gradz &sysg);

        // Input                pset    : A parameter set
        //                      sysg    : System in z-gradient (this)
        // Output               pset    : The parameter set with
        //                                only sysg parameters
        ///F_list +=  - Unary Addition


MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;

        // Input                sysg    : System in z-gradient (this)
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



// ----------------------------------------------------------------------------
//            Functions To Make A Spin System From A Parameter Set
// ----------------------------------------------------------------------------

//          Many of these handled by the base class spin_system

// int setOm(const ParameterSet& pset);
// void setJs(const ParameterSet& pset);
// void setAs(const ParameterSet& pset);
// void setShifts(const ParameterSet& pset);


MSVCDLL void setSubSys(const ParameterSet& pset);

        // Input                sysg    : System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system # of sub-systems
        //                                is set from parameter in pset


MSVCDLL void setBoGrad(const ParameterSet& pset);

        // Input                sysg    : System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system field gradient
        //                                is set from parameter in pset
        // Note                         : Input units expected T/m


MSVCDLL void setLength(const ParameterSet& pset);

        // Input                sysg    : System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system effective length
        //                                is set from parameter in pset
        // Note                         : Input units expected T/m


MSVCDLL void operator= (const ParameterSet& pset);

        // Input                sysg    : System in z-gradient (this)
        //                      pset    : A parameter set
        // Output               none    : Spin system filled with
        //                                parameters in pset
        // Note                         : Functions which place a sysg
        //                                into a parameter set must contain
        //                                the information read here

// ----------------------------------------------------------------------------
//     Functions To Output Isotropic System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

        // F_list write  - Write to a file
        // Input                sysg    : System in z-gradient (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none	: Spin system is written as a
        //                                parameter set to file filename
        //                                or to output filestream
        // Note                         : This depends on function PSetAdd!

MSVCDLL virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
MSVCDLL virtual int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const;

// ____________________________________________________________________________
// I                        SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________

///F_list read          - Read spin system from disk file
///F_list ask_read      - Ask for file, read spin system from file

MSVCDLL virtual int read(const std::string& filename, int idx=-1, int warn=2);

        // Input                sysg    : System in z-gradient (this)
        //                      filename : Input filename
        //                      idx      : Parameter index value used for
        //                                 prefix [#] in input names
        //                      warn     : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               none     : Spin system is filled with
        //                                 parameters read from file
        // Note                          : The file should be an ASCII file
        //                                 containing recognized sys parameters


MSVCDLL virtual int read(const ParameterSet& pset, int idx=-1, int warn=2);

        // Input                sysg    : System in z-gradient (this)
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


MSVCDLL std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,
                                                     const std::string& def);

        // Input                sysg    : System in z-gradient (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               string  : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the spin system is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The spin system is modifed (filled)


// ____________________________________________________________________________
// I                            STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

        // Input                out      : output stream
        // Output               none	 : modifies output stream
        //
        // F_list   print  - Print to an output stream

MSVCDLL        std::ostream& print      (std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const sys_gradz& sys);

};

#endif						// sys_zgrad.h
