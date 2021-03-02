/* PulseI.h *********************************************-*-c++-*-
**                                                              **
**                          G A M M A                           **
**                                                              **
**      Ideal Pulses                       Interface		**
**                                                              **
**      Copyright (c) 1990, 1991, 1992, 1993                    **
**      Scott Smith                                             **
**      Eidgenoessische Technische Hochschule                   **
**      Labor fuer physikalische Chemie                         **
**      8092 Zurich / Switzerland                               **
**                                                              **
**      $Header: $
**                                                              **
*****************************************************************/

/*****************************************************************
**                                                              **
** Description                                                  **
**                                                              **
** The functions in this module are the realization of ideal    **
** pulses.  Such pulses produce pure spin rotations with no     **
** offset effects.  They are infinitesimally short and taken    **
** as infinitely strong.  They are not affected by any relative **
** rotating frame.                                              **
**                                                              **
*****************************************************************/

#ifndef   PulseI_h_			// Is file already included?
#  define PulseI_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSys.h>		// Knowledge of spin systems (base)
#include <HSLib/GenOp.h>		// Knowledge of operators

// ____________________________________________________________________________
// A               Pulses Acting Directly On Density Operators
// ____________________________________________________________________________
 
/* The functions below apply an ideal pulse of specified rotation to an input
   density operator. Since the pulses are ideal they are applied irrespective
   of rotating frames and do not advance the density operator in time.
 
   Functions are overloaded to produces pulses which operate on a single spin,
   all spins of a specified isotope type, all spins in the system, or any
   combination of spins in the system as specified by spin flag settings.
 
           Input                sys   : A (base) spin system
                                sigma : Current density operator
                                SEL   : Pulse selectivity
                                beta  : Pulse rotation angle (degrees)
                                phi   : Pulse phase angle (degrees)
           Output               sigma : Density operator after application
                                        of ideal pulse on selected spin(s)
 
           OverLoads of Function where Type of Second Argument Changes

           OL Function   SEL    Label   Selectivity Result             
           -- --------  -----   -----   ------------------------------
           A. Iupuls     int    spin    Spin index, only spin affected
           B. Iupuls    string  iso     Isotope label, type iso spins affected
           C. Iupuls     ---    ---     All spins in sys affected
           D. Iupuls     int*   flags   Spins with flags TRUE affected
           E. Iupuls_sp  ---            Spins with sys flags TRUE affected   */
 
// ************************* Pulse Along The X Axis ***************************

MSVCDLL gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, int spin, double beta);
MSVCDLL gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, const std::string& iso,
                                                                  double beta);
MSVCDLL gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma,           double beta);
MSVCDLL gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, 
                                       const flagvec& flags, double beta);
MSVCDLL gen_op Ixpuls_sp(const spin_sys& sys, const gen_op& sigma,        double beta);

// ************************* Pulse Along the Y Axis ***************************

MSVCDLL gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, int spin, double beta);
MSVCDLL gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, const std::string& iso,
                                                                  double beta);
MSVCDLL gen_op Iypuls(const spin_sys& sys, const gen_op& sigma,           double beta);
MSVCDLL gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, 
                                       const flagvec& flags, double beta);
MSVCDLL gen_op Iypuls_sp(const spin_sys& sys, const gen_op& sigma,        double beta);

// ************* Pulse In the XY Plane On Axis Phi Degrees From X *************

MSVCDLL gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma,
                                            int spin, double phi, double beta);
MSVCDLL gen_op Ixypuls(const spin_sys& sys,const gen_op& sigma,
                                   const std::string& iso, double phi, double beta);
MSVCDLL gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma,
                                                      double phi, double beta);
MSVCDLL gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma, 
                           const flagvec& flags, double phi, double beta);
MSVCDLL gen_op Ixypuls_sp(const spin_sys& sys, const gen_op& sigma,
                                                      double phi, double beta);
 
// ____________________________________________________________________________
// B                        Ideal Pulse Propagators
// ____________________________________________________________________________
 
/* The functions below return an ideal pulse Hilbert space propagator of the
   specified rotation applicable to the to an input spin system. Since the
   pulses are ideal they are applied irrespective of rotating frames and do not
   advance the density operator in time.
 
   Functions are overloaded to produces pulse propagators which will operate on
   a single spin, all spins of a specified isotope type, all spins in the
   system, or any combination of spins in the system as specified by spin flag
   settings.
 
           Input                sys   : A (base) spin system
                                sigma : Current density operator
                                SEL   : Pulse selectivity
                                beta  : Pulse rotation angle (degrees)
                                phi   : Pulse phase angle (degrees)
           Output               sigma : Density operator after application
                                        of ideal pulse on selected spin(s)
 
           OverLoads of Function where Type of Second Argument Changes

        OL Function     SEL    Label   Selectivity Result              
        -- ----------  -----   -----   ------------------------------
        A. Iupuls_U     int    spin     Spin index, only spin affected
        B. Iupuls_U    string  iso      Isotope label, type iso spins affected
        C. Iupuls_U     ---    ---      All spins in sys affected
        D. Iupuls_U     int*   flags    Spins with flags TRUE affected
        E. Iupuls_sp_U  ---             Spins with sys flags TRUE affected   */
 
// ******************* Pulse Propagators Along The X Axis *********************

MSVCDLL gen_op Ixpuls_U(const spin_sys& sys, int spin,                  double beta);
MSVCDLL gen_op Ixpuls_U(const spin_sys& sys, const std::string& iso,         double beta);
MSVCDLL gen_op Ixpuls_U(const spin_sys& sys,                            double beta);
MSVCDLL gen_op Ixpuls_U(const spin_sys& sys, const flagvec& flags, double beta);
MSVCDLL gen_op Ixpuls_sp_U(const spin_sys& sys,                         double beta);
 
// ******************* Pulse Propagators Along The Y Axis *********************
 
MSVCDLL gen_op Iypuls_U(const spin_sys& sys, int spin,                  double beta);
MSVCDLL gen_op Iypuls_U(const spin_sys& sys, const std::string& iso,         double beta);
MSVCDLL gen_op Iypuls_U(const spin_sys& sys,                            double beta);
MSVCDLL gen_op Iypuls_U(const spin_sys& sys, const flagvec& flags, double beta);
MSVCDLL gen_op Iypuls_sp_U(const spin_sys& sys,                         double beta);
 
// ******** Pulse Propagators In the XY Plane On Axis Phi Degrees From X ******
 
MSVCDLL gen_op Ixypuls_U(const spin_sys& sys, int spin,       double phi, double beta);
MSVCDLL gen_op Ixypuls_U(const spin_sys& sys, const std::string& I,double phi, double beta);
MSVCDLL gen_op Ixypuls_U(const spin_sys& sys,                 double phi, double beta);
MSVCDLL gen_op Ixypuls_U(const spin_sys& sys, const flagvec& flags,
                                                      double phi, double beta);
MSVCDLL gen_op Ixypuls_U_sp(const spin_sys& sys,              double phi, double beta);
 
#endif							 // PulseI.h
