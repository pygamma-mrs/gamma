/* Powder.h *****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Powder						Interface	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2002                                              **
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
** A variable of type Powder defines an array of Euler angles and a	**
** corresponding array of scaling factors. Each set of Euler angles	**
** defines a rotation that takes the desired axes (LAB) into initial 	**
** axes (PAS). A corrsponding scaling factor defines the probability	**
** of finding that particular orientation relative to the other 	**
** orientations.							**
**                                                                      **
** Although the Euler angles used herein rotate (zum beispiel) LAB axes	**
** into PAS axes by definition, GAMMA uses these Euler angles to take	**
** interactions	from the PAS into the LAB frame. As is mentioned in	**
** countless places in GAMMA, this "inverse" rotation usage allows one	**
** to relate spherical coordinate angles (theta,phi) to the Euler	**
** angles (beta, alpha). 						**
**                                                                      **
** Thus, assume that some intial axes are the PAS of some interaction.	**
** Then the Euler angles {alpha,beta,0} would be the set GAMMA uses to 	**
** rotate into the laboratory axes so that the PAS z-axis will be at	**
** angle beta down form LAB +z axis and at angle alpha over from LAB	**
** +x axis. This is used for convenience even though {alpha, beta, 0}	**
** actually define the inverse rotation...... AOK? Ganz klar.		**
**                                                                      **
*************************************************************************/

#ifndef _Powder_h_			// Is file already included?
#define _Powder_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface 			// This is the interface
#endif

#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Include Euler angles
#include <string>			// Include stdlibc++ strings

class Powder
  {
  std::vector<EAngles> _EAs;		// Vector of Euler angles
  std::vector<double>  _Ps;		// Vector of probabilities
  std::string          _PN;		// Powder name
  bool                 _A3;		// Flag if three angles

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                         POWDER ERROR HANDLING
// ____________________________________________________________________________

/*       Input 		      CFrm    : Coordinate Frame (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void POWerror(int eidx,                           int noret=0) const;
         void POWerror(int eidx, const std::string& pname, int noret=0) const;
volatile void POWfatal(int eidx)                                        const;
volatile void POWfatal(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii                POWDER CHECKING FUNCTIONS
// ____________________________________________________________________________

//bool Powder::ChkIdx(int i, bool warn=true) const;	// Check index i OK

// ____________________________________________________________________________
// iii         POWDER PARAMETER SET SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions try and set a coordinate frame from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 CoordFrmi  (2) : AxesName - Initial Coordinate Axes
                 CoordFrmf  (2) : AxesName - Final Coordinate Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this coordinate frame) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Axesf Into Axesi

   where the # is used to indicate the object that uses the rotation. 

           Input        CFrm    : Coordinate Frame (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Coordinate Frame is set
                                  from parameters in pset		     */

//bool Powder::SetCoordFrm(const    ParameterSet& pset,int pfx=-1,  int warn=2);
//bool Powder::SetRotation(const    ParameterSet& pset,             int idx=-1);
//bool Powder::GetEulerAngles(const ParameterSet& pset,EAngles& EA, int idx=-1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                    POWDER CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Each powder contains an array of Euler angles and a cooresponding array
   of probabilities.                                                         */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

Powder::Powder();

// ----------------------------------------------------------------------------
//                       Constructors Using Two Vectors
// ----------------------------------------------------------------------------

Powder::Powder(const vector<EAngles>& EAs, const vector<double>& Ps,
                                                      const std::string& Name);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

//Powder::Powder(const ParameterSet& pset, int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

Powder& Powder::operator= (const Powder& Pow);
     Powder::~Powder();

// ____________________________________________________________________________
// B                        POWDER ACCESS FUNCTIONS
// ____________________________________________________________________________

EAngles Powder::EA(int i) const;
void    Powder::EA(const EAngles& EA, int i);

double  Powder::Pop(int i) const;
void    Powder::Pop(double& P, int i);

std::string  Powder::Name() const;
void         Powder::Name(const std::string& PN);

int     Powder:: size() const;

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Coordinate Frame
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a coordinate frame.
   Individual rotations will be placed into the parameter set as Euler
   angle in degrees.                                                         */

/*
            operator ParameterSet( ) const;
friend void operator+= (ParameterSet& pset, const Powder& CFrm);
       void Powder::PSetAdd(ParameterSet& pset,   int pfx=-1) const;
*/

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Coordinate Frame
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

//void Powder::write(const std::string &filename, int pfx=-1) const;

// ____________________________________________________________________________
// D                      POWDER INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the coordinate frame from either a
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple coordinate frames may be defined in the
   same file by switching the prefix indices.

           Input                CFrm    : Coordinate Frame
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                pfx     : Coordinate frame index (default -1->none)
           Output               none    : Coordinate frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

/*
bool Powder::read(const std::string& file,  int idx=-1, int warn=2);
bool Powder::read(const ParameterSet& pset, int idx=-1, int warn=2);
std::string Powder::ask_read(int argc, char* argv[], int argn,
                                                                   int idx=-1);
std::string Powder::ask_read(int argc, char* argv[], int argn,
                                                const string& def, int idx=-1);
*/
 
// ____________________________________________________________________________
// E                       POWDER OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame 
   to any specified output stream.

           Input                CFrm	: Coordinate Frame (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Coordinate Frame information
                                          placed into the output stream      */

       std::ostream& Powder::print(std::ostream& out, int fflag=100) const;
friend std::ostream& operator<<  (std::ostream& out, const Powder& CFrm);

// ____________________________________________________________________________
// F                             KNOWN POWDERS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                            Spherical Powders
//-----------------------------------------------------------------------------

/* These functions generate the simplest powder, just even increments over
   the entire sphere. Since the probablity of orientation with beta=0,PI is
   zero, these angles will NOT be included in the Euler angle set. If one has
   a system with no asymmetry then it is best to here set Nalpha to zero and
   limit the range to the upper hemisphere. This is done by default.         */

void Powder::Spherical(int  Nbeta=180, int Nalpha=1,   bool uhonly=true);
void Powder::Spherical2(int Nbeta=180, int Nalpha=720, bool uhonly=false);


//-----------------------------------------------------------------------------
//                           Cheng Two Angle Powder
//-----------------------------------------------------------------------------

/* These functions set the powder according the Cheng et. al..

   The array pa contains values related to the powder orientation angle beta,
   or equivalently the angle theta. Given a powder quality value, PAQ, there
   will be N = pa[PAQ] crystallite orientations set in powder. For the kth
   orientation, the angle beta will be given by PI*(k/N)..... this is
   just evenly incremented down from the z-axis within [0, PI) as in a
   spherical powder. The orientation probablity will be sin(beta), again as
   in a spherical powder. Of course the first value, k=0, will never be used
   since sin(theta) = 0.  The 2nd orientation angle, phi (or equivalently
   alpha), for the kth orientation is given by a somewhat more cryptic formula.
   For this powder, alpha = 2*PI*{(pb[PAQ]*k) % N}/N. This odd incrementation
   is of course the reason the Cheng powder differs from a spherical one.    */

void Powder::Cheng(int PAQ);

  };

#endif							// Powder.h
