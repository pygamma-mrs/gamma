/* MultiSys.h ***************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      Class Multiple System                       Interface		**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      1800 E. Paul Dirac Drive                                        **
**      National High Magnetic Field Laboratory                         **
**      Tallahassee FL 32306 USA                                        **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This class embraces mulitple spin systems that have no scalar        **
** coupling between them. Each spin system is of type sys_dynamic,      **
** containing info on shifts, couplings, relaxtion, and mutual exchange **
** (dynamics). This class allows one to account for non-mutual exchange **
** between molecular groups that are not coupled to each other.  Note   **
** that intra-molecular exchange may be accounted for by sys_dynamic    **
** alone.                                                               **
**                                                                      **
*************************************************************************/

#ifndef   multi_sys_h_			// Is this file already included?
#  define multi_sys_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <MultiSys/ExProcess.h>		// Knowledge of exchange processes
#include <Basics/ParamSet.h>		// Know about GAMMA parameter sets
#include <string>			// Knowledge of libstdc++ strings
#include <vector>			// Know libstdc++ STL vectors
#include <LSLib/sys_dynamic.h>		// Knowledge of dynamic systems

class multi_sys

{
  std::string              _SysName;		// The spin system name
  std::vector<double>      _Pops;		// Population of components
  std::vector<sys_dynamic> _Comps;		// System components
  std::vector<ExchProc>	   _Exs;		// Exchange processes

private:

 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                    CLASS MULTI_SYS ERROR HANDLING
// ____________________________________________________________________________

         void MSYSerror(int eidx,                           int noret=0) const;
         void MSYSerror(int eidx, const std::string& pname, int noret=0) const;
volatile void MSYSfatal(int eidx)                                        const;
volatile void MSYSfatal(int eidx, const std::string& pname)              const;

// ____________________________________________________________________________
// ii         BUILD FUNCTIONS FOR MULTI_SYS FROM A PARAMETER SET
// ____________________________________________________________________________

/* The purpose of these functions is to set up a multiple spin system from
   parameters in a specified parameter set. Many of these must be kept private
   because their misuse could lead to big problems as they totally tweak with
   the class internal  structure.  For example, the number of components can be
   set without the actual number of allocated components changing. So, their
   misuse would be a very bad thing. If they are private they won't be abused
   (except herein). In general, the "get" functions do NOT alter the spin
   system, they only glean information from the parameter set for the system. 
   In contrast the "set" functions will alter the spin system using the
   parameters (gotten from "set" functions).

  Function        Type                            Purpose
 ==========  ==============  ==================================================
 getNComps   int             Get the number of components in MultiSys: NComp
 getMSName   string          Get system name:                          MSysName
 getFName    string          Get component (sub-system) file name:     Fname(#)
 getFNames   vector<string>  Get component (sub-system) file names:    Fname(#)
 getComp     sys_dynamic     Get component of MultiSys:                Various
 getComps    vector<sys_dyn> Get all components of MultiSys:           Various
 getPop      double          Get component population:                 Popul(#)
 getPops     vector<double>  Get populations of all components         Popul(#)
 getNex      int             Get the number of exchange processes:     Exch(#)
*/

// ----------------------------------------------------------------------------
//                 Reading Of Basics MultiSystem Parameters
// ----------------------------------------------------------------------------

/* These read very simple parameters associated with the muliple spin system.
   One reads the spin system name and this typically is set to not be of any
   consequence at all.  Another reads the number of components, i.e. the number
   of spin systems within the multiple spin system. This is usually critical
   and warning should be on if it is called.                                 */

bool getNComps(const ParameterSet& pset,int& ncmps,bool warn=true);
bool getMSName(const ParameterSet& pset, 
                                            std::string& name,bool warn=true);

// ----------------------------------------------------------------------------
//             Reading Of MultiSystem Components From Parameters
// ----------------------------------------------------------------------------

/* Each MultiSystem contains any number of spin systems of type sys_dynamic.
   When reading in the MultiSystem from a parameter file, all parameters for
   each component can either reside in the file or another filename can be
   specified that contains the component. If the latter is used then we must
   read in additional external files, one for each component.  If the former
   is used then we do additional parse on the input file. Of course the files
   are immediately read into GAMMA parameter sets, so a file is equivalent
   to a parameter set.                                                        */

bool getFName(const  ParameterSet& pset,
                                   std::string& name, int idx, bool warn=true);

bool getFNames(const ParameterSet& pset,
                              std::vector<std::string>& names, bool warn=true);

bool getComp(const ParameterSet& pset, int idx,
                                             sys_dynamic& cmp, bool warn=true);

bool getComps(const ParameterSet& pset, int ncmps,
                               std::vector<sys_dynamic>& cmps, bool warn=true);

// ----------------------------------------------------------------------------
//           Reading Of Component Populations And Exchange Processes
// ----------------------------------------------------------------------------

bool getPop(const ParameterSet& pset, int idx,
                                            double& pop, bool warn=true) const;

bool getPops(const ParameterSet& pset, int ncmps,
                              std::vector<double>& pops, bool warn=true) const;

int  getNex(const ParameterSet& pset) const;

bool getProcesses(const ParameterSet& pset, 
                            std::vector<ExchProc>& procs, bool warn=true) const;

// ----------------------------------------------------------------------------
//                    Reading Of Enire Multipule Spin System
// ----------------------------------------------------------------------------

bool getMsys(const ParameterSet& pset,
            std::string& name, std::vector<sys_dynamic>& cmps,
       std::vector<double>& pops, std::vector<ExchProc>& procs, bool warn=true);

bool setMsys(const ParameterSet& pset,              bool warn=true);


// ____________________________________________________________________________
// iii                  MULTISYSTEM CHECKING FUNCTIONS
// ____________________________________________________________________________

bool CheckNComps(int                nc, bool warn=true) const;
bool CheckRange(unsigned             n, bool warn=true) const;
bool CheckProc(int                  ip, bool warn=true) const;
bool CheckProcs(                        bool warn=true) const;
bool CheckProc(const ExchProc& XP,      bool warn=true) const;
bool CheckField(const spin_system& sys, bool warn=true) const;
 
// ---------------------------------------------------------------------------- 
// ---------------------------- PUBLIC FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
  
public:

// ____________________________________________________________________________ 
// A              MULTI-SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 
///Center Spin System Algebraic
///F_list =                   - Assignment

// ---------------------------------------------------------------------------- 
//                           Simple Constructors
// ---------------------------------------------------------------------------- 

MSVCDLC multi_sys();
//multi_sys(int Ncomp);
MSVCDLC multi_sys(const multi_sys& msys);


// ---------------------------------------------------------------------------- 
//        Construction From Systems, Populations, And Exchange Rates
// ---------------------------------------------------------------------------- 

MSVCDLC multi_sys(double pop1, sys_dynamic &sys1, 
                               double pop2, sys_dynamic &sys2, double krate=0);
 
        // Input                pop1  : Population system 1
	//			sys1  : Spin system 1
        //                      pop2  : Population system 2
	//			sys2  : Spin system 2
	//			krate : Exchange rate between systems 
        // Output               none  : spin system constructor


// ---------------------------------------------------------------------------- 
//                       Destruction, Assignment
// ---------------------------------------------------------------------------- 

MSVCDLL multi_sys& operator= (const multi_sys& msys);
//MSVCDLL void operator= (const multi_sys &msys);
MSVCDLC      ~multi_sys();

// ____________________________________________________________________________
// B                       MULTI-SPIN SYSTEM NAME
// ____________________________________________________________________________

        // Input                msys 	: Multi-spin system
        // Output               string  : Multi-spin system name
        ///F_list name                  - Set or retrieve spin system name.

MSVCDLL void               name(const std::string& sysname);
MSVCDLL const std::string& name() const;

// ____________________________________________________________________________
// C                     MULTI-SPIN POPULATION ACCESS
// ____________________________________________________________________________

        // Input		msys  : Multi system
	//			icomp : Component index
        //                      npop  : Population of component icomp
        // Output               none  : Popullation of component i set
	//                   or double: Component i population returned

MSVCDLL void   pop(int icomp, double npop);
MSVCDLL double pop(int icomp) const;
MSVCDLL double popmin() const;
MSVCDLL double popmax() const;
  
// ____________________________________________________________________________
// D                       MULTI-SPIN COMPONENT ACCESS
// ____________________________________________________________________________
 
MSVCDLL int                NComps() const;
MSVCDLL void               Comp(int icomp, const sys_dynamic& sys);
MSVCDLL const sys_dynamic& Comp(int icomp) const;
MSVCDLL void               AddComp(const sys_dynamic& sys, double pop=0);

MSVCDLL void CheckComp(unsigned n) const;

// ____________________________________________________________________________
// E                    MULTI-SPIN EXCHANGE PROCESS ACCESS
// ____________________________________________________________________________

MSVCDLL       int       NExProcs() const;
MSVCDLL const ExchProc& ExProc(int iex) const;
MSVCDLL       void      ExProc(const ExchProc& pr, int iex);
MSVCDLL       double    Kex(int iex) const;
MSVCDLL       void      Kex(double K, int iex);
MSVCDLL       int       NCompsLHS(int iex) const;
MSVCDLL       int       NCompsRHS(int iex) const;

// ____________________________________________________________________________
// F            GLOBAL & INDIVIDUAL COMPONENT SPIN SYSTEM FUNCTIONS
// ____________________________________________________________________________

/* This section contains typical GAMMA spin system functions but they are
   of two types: global and specific.  The Global functions are those which
   apply to all components (spin systems). That is, use of the function on an
   individual component will produce the same result as use on the multi_sys.
   The specific functions act on a single component within the system, but not
   on all components.

   Individual components are spin systems. Since class multi_sys is NOT derived
   from any base spin system class, access functions to individual components
   are NOT automatically inherited. Thus, they must either be explicitly
   provided for here or users must obtain the system (component), alter it,
   then put it back into the multiple spin system. 

   Note that for each of the functions, Fct(cmp,a,b,...), users could also just
   directly use Comp(cmp).Fct(a,b,....) to attain the same results. There are
   of course limited functions that are set here. Only those that do not ruin
   the multisys structure are allowed, & only those often used are set up.   */

// ----------------------------------------------------------------------------
//                      Homonuclear Versus. Heteronuclear
// ----------------------------------------------------------------------------

MSVCDLL bool homonuclear(int   cmp=-1) const;
MSVCDLL bool heteronuclear(int cmp=-1) const;

MSVCDLL int HS(int comp=-1) const;
MSVCDLL int LS(int comp=-1) const;
MSVCDLL std::vector<int> HSs()   const;
MSVCDLL std::vector<int> LSs()   const;

MSVCDLL const std::string symbol(int comp, int spin) const;
 
// ----------------------------------------------------------------------------
//                      SPECTROMETER FREQUENCY MANIPULATIONS
// ----------------------------------------------------------------------------

///Center Multi Spin System Spectrometer Frequency Functions

        // Input                freq : spectrometer frequency (MHz)
        // Output               none : stored in MHz
        ///F_list Omega              - Set/Retrieve spectrometer frequency
        // Input                freq : Larmor frequency (MHz)
        //                      iso  : spin isotope type
        // Output               none : 1H Larmor stored in MHz
        // Input                none   :
        // Output               double : spectrometer frequency in MHz
        // Input                none   :
        // Output               iso    : isotope label
 
//------------------------------ Set Omega ------------------------------------

MSVCDLL void Omega(double freq);
MSVCDLL void Omega(double freq, const std::string& iso);

//------------------------------ Get Omega ------------------------------------
 
MSVCDLL double Omega( ) const;
MSVCDLL double Omega(const std::string& iso) const;                                 

// ____________________________________________________________________________
// G                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//        Functions To Make A Parameter Set From A Multi Spin System
// ----------------------------------------------------------------------------

/* This class contains a name, a vector of components (spin systems), a vector
   of component populations, and a vector of component exchange processes. We
   will write our own name as well as our output our own component populations,
   however we will allow the components (sys_dynamic) and exchange processes
   (ExchProc) to write themselves since they are in a class structure.  For
   them we need only loop over all components and have themselves perform this
   task with the proper indexing. Note that the multi_sys is placed into only
   one parameter set (unlike input which may occur from multiple ASCII files).
   Also note that parameter prefixing thorough use of [#] is NOT allowed here
   Only one multi_sys may be put into the same parameter set.

     Function                                 Purpose
   ------------         -------------------------------------------------------
   ParameterSet         Convert system into a parameter set
   operator +=          Adds system to existing parameter set                */

MSVCDLL operator ParameterSet() const;
MSVCDLL friend void operator+= (ParameterSet& pset, const multi_sys &msys);

// ----------------------------------------------------------------------------
//    Functions To Output A Multi Spin System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

/* These functions write the multiple spin system into an ASCII file in GAMMA
   Parameter Set format.  That is, the resulting ASCII file may be read into
   any GAMMA program to create a multiple spin system identical to that written.

           Input                msys    : Multi_sys spin system (this)
                                filename: Output file name
                            or  ofstr   : Output file stream
                                warn    : Warning level
           Output               none    : System is written as a
                                          parameter set to file filename
                                          or into output file stream ofstr 

	//			basename: Component output filename base
        // Output               none    : Multi_sys spin system is written as a 
        //                                parameter set to file filename & each
	//				  component written to a file with name
	//				  basename+i.dsys - basename is given
	//				  as an arguemnt and i is an integer */

MSVCDLL void write(std::string& filename, std::string basename = "comp");

// ____________________________________________________________________________
// H                        SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                msys    : Multi_sys spin system (this)
        //                      filename: Input filename
        //                      pset    : Input parameter set
        //                      warn	: Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF	: Spin system is filled with
        //                                parameters read from file
        //                                TRUE if read is successful
        // Note				: The file should be an ASCII file
        //                                containing recognized sys parameters

MSVCDLL bool read(const std::string&  filename, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,     int warn=2);

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn, const std::string& def);
 
        // Input                msys    : Multi_sys spin system (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the spin system is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The spin system is modifed (filled)
 
// ____________________________________________________________________________
// I                  MULTIPLE SPIN SYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* The next function builds an array of spin mappings for a specified exchange
   process. The vector will contain a list of lines similar to

                            A 0 <---> B 2  13C

   Here the lone letter indicate components (spin systems) & the lone integers
   spins.  At the end the type of spin in the exchange is listed. Two added
   details. First, the A 0 <---> B 2 part of each line is supplied in the base
   class ExchProc so we need only ask that class to provide them. Second we
   will format the output strings so that they will precisely fit under the
   header used in the print function of this class.  Thus, the printing will
   appear akin to
                           Spin Mappings   Type
                          ---------------  ----
                           A 0 <---> B 2   13C
                           A 1 <---> B 3    1H

   Each string will be of length 22 (we use 1 column before the header too.) */

MSVCDLL std::vector<std::string> SpinMapStrs(int exp) const;

//vector<string> ExchMapStrs(int exp) const;
//vector<string> ExchProcStrs(int exp) const

// sosik
MSVCDLL std::vector<std::string> LHSStrs() const;
MSVCDLL std::vector<std::string> RHSStrs() const;
MSVCDLL std::vector<std::string> EXPStrs() const;

//-----------------------------------------------------------------------------
//              Functions To Modularize Exchange Process Output
//-----------------------------------------------------------------------------

/* The next function builds an array of spin mappings for a specified exchange
   process. Each map will appear as

                   SystemName Spin <----> SystemName Spin

   where each component (system) name will be limited to 15 characters max.  */

        // Input                msys  : Multi system
        //                      ostr  : Output stream
        //                      full  : Print amount flag
        //                                !0 = print individual systems (def)
        //                                 0 = don't print individual systems
        // Output               ostr  : The output stream  is returned
        //                              with the spin system added 

        // Input                out      : output stream;
        //                      sys      : Spin system to write
        // Output               none     : modifies output stream
        ///F_list <<                     - Send spin system to an output stream
 

MSVCDLL        std::ostream& print(std::ostream& ostr, int full=1) const;
MSVCDLL friend std::ostream& operator<<(std::ostream& out, const multi_sys& sys);

};
 
#endif								// MultiSys.h
