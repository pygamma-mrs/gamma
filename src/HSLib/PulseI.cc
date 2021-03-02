/* PulseI.cc ****************************************************-*-c++-*-
**									**
**	                           G A M M A				**
**							 		**
**	Ideal Pulses			           Implementation	**
**						        	 	**
**	Copyright (c) 1990, 1991, 1992, 1993	         		**
**	Scott Smith				 	        	**
**	Eidgenoessische Technische Hochschule	 	        	**
**	Labor fuer physikalische Chemie		         		**
**	8092 Zurich / Switzerland		         		**
**						         		**
**      $Header: $
**							         	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description						 		**
**							 		**
** The functions in this module are the realization of ideal pulses.	**
** Such pulses produce pure spin rotations with no offset effects.	**
** They are infinitesimally short and taken as infinitely strong.	**
**  They are not affected by any relative rotating frame.		**
**							         	**
*************************************************************************/

#ifndef _PulseI_cc_			// Is file already included?
#define _PulseI_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <HSLib/PulseI.h>		// Include the header
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of composite spin ops
#include <HSLib/SpinOpRot.h>		// Knowledge of spin rotation ops
#include <HSLib/GenOp.h>		// Knowledge of general operators
#include <HSLib/HSprop.h>		// Knowledge of propagators
#include <HSLib/SpinSys.h>
#include <HSLib/GenOp.h>
#include <vector>			// Include libstdc++ stl vectors

// ____________________________________________________________________________
// A               Pulses Acting Directly On Density Operators
// ____________________________________________________________________________

/* The functions below apply an ideal pulse of specified rotation to an input
   density operator. Since the pulses are ideal they are applied irrespective
   of rotating frames and do not advance the density operator in time.

   Functions are overloaded to produces pulses which operate on a single spin,
   all spins of a specified isotope type, all spins in the system, or any
   combination of spins in the system as specified by spin flag settings.

	   Input		sys   : A (base) spin system
	   			sigma : Current density operator
	   			SEL   : Pulse selectivity
	   			beta  : Pulse rotation angle (degrees)
				phi   : Pulse phase angle (degrees)
	   Output		sigma : Density operator after application
	  				of ideal pulse on selected spin(s)

           OverLoads of Function where Type of Second Argument Changes
 
           OL Function   SEL    Label   Selectivity Result
           -- --------  -----   -----   ------------------------------
           A. Iupuls     int    spin	Spin index, only spin affected
           B. Iupuls    string  iso	Isotope label, type iso spins affected
           C. Iupuls     ---    ---	All spins in sys affected
           D. Iupuls     int*   flags	Spins with flags TRUE affected
           E. Iupuls_sp  ---            Spins with sys flags TRUE affected   */

// ************************* Pulse Along The X Axis ***************************

gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, int spin,  double beta)
  {  return evolve(sigma, Rx(sys,spin,beta)); }

gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma,
                                                const std::string& Iso, double beta)
  { return evolve(sigma, Rx(sys,Iso,beta)); }

gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma,            double beta)
  { return evolve(sigma, Rx(sys,beta)); }

gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma,
                                        const flagvec& flags, double beta)
  { return evolve(sigma, Rx(sys,flags,beta)); }

gen_op Ixpuls_sp(const spin_sys& sys, const gen_op& sigma,         double beta)
  { return evolve(sigma, Rx_sp(sys,beta)); }

// ************************* Pulse Along the Y Axis ***************************

gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, int spin,  double beta)
  { return evolve(sigma, Ry(sys,spin,beta)); }

gen_op Iypuls(const spin_sys& sys, const gen_op& sigma,
                                                const std::string& Iso, double beta)
  { return evolve(sigma, Ry(sys,Iso,beta)); }

gen_op Iypuls(const spin_sys& sys, const gen_op& sigma,            double beta)
  { return evolve(sigma, Ry(sys,beta)); }

gen_op Iypuls(const spin_sys& sys, const gen_op& sigma,
                                        const flagvec& flags, double beta)
  { return evolve(sigma, Ry(sys,flags,beta)); }

gen_op Iypuls_sp(const spin_sys& sys, const gen_op& sigma,         double beta)
  { return evolve(sigma, Ry_sp(sys,beta)); }

// ************* Pulse In the XY Plane On Axis Phi Degrees From X *************

gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma,
                                             int spin, double phi, double beta)
  { return evolve(sigma, Rxy(sys,spin,phi,beta)); }

gen_op Ixypuls(const spin_sys& sys, const gen_op &sigma,
                                    const std::string& iso, double phi, double beta)
  { return evolve(sigma, Rxy(sys,iso,phi,beta)); }

gen_op Ixypuls(const spin_sys& sys, const gen_op &sigma,
                                                       double phi, double beta)
  { return evolve(sigma, Rxy(sys,phi,beta)); }

gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma,
                            const flagvec& flags, double phi, double beta)
  { return evolve(sigma, Rxy(sys,flags,phi,beta)); }

gen_op Ixypuls_sp(const spin_sys& sys, const gen_op &sigma,
                                                       double phi, double beta)
  { return evolve(sigma, Rxy_sp(sys,phi,beta)); }

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

	   Input		sys   : A (base) spin system
	   			sigma : Current density operator
	   			SEL   : Pulse selectivity
	   			beta  : Pulse rotation angle (degrees)
				phi   : Pulse phase angle (degrees)
	   Output		sigma : Density operator after application
	  				of ideal pulse on selected spin(s)

           OverLoads of Function where Type of Second Argument Changes
 
        OL Function     SEL    Label   Selectivity Result
        -- ----------  -----   -----   ------------------------------
        A. Iupuls_U     int    spin	Spin index, only spin affected
        B. Iupuls_U    string  iso	Isotope label, type iso spins affected
        C. Iupuls_U     ---    ---	All spins in sys affected
        D. Iupuls_U     int*   flags	Spins with flags TRUE affected
        E. Iupuls_sp_U  ---             Spins with sys flags TRUE affected   */

// ******************* Pulse Propagators Along The X Axis *********************

gen_op Ixpuls_U(const spin_sys& sys, int spin,                 double beta)
  { return gen_op(Rx(sys, spin, beta)); }

gen_op Ixpuls_U(const spin_sys& sys, const std::string& iso,        double beta)
  { return gen_op(Rx(sys, iso, beta)); }

gen_op Ixpuls_U(const spin_sys& sys,                           double beta)
  { return gen_op(Rx(sys, beta)); }

gen_op Ixpuls_U(const spin_sys& sys, const flagvec& flags, double beta)
  { return gen_op(Rx(sys, flags, beta)); }

gen_op Ixpuls_sp_U(const spin_sys& sys,                        double beta)
  { return gen_op(Rx_sp(sys, beta)); }

// ******************* Pulse Propagators Along The Y Axis *********************

gen_op Iypuls_U(const spin_sys& sys, int spin,                  double beta)
  { return gen_op(Ry(sys, spin, beta)); }

gen_op Iypuls_U(const spin_sys& sys, const std::string& iso,         double beta)
  { return gen_op(Ry(sys, iso, beta)); }

gen_op Iypuls_U(const spin_sys& sys,                            double beta)
  { return gen_op(Ry(sys, beta)); }

gen_op Iypuls_U(const spin_sys& sys, const flagvec& flags, double beta)
  { return gen_op(Ry(sys, flags, beta)); }

gen_op Iypuls_sp_U(const spin_sys& sys,                         double beta)
  { return gen_op(Ry_sp(sys, beta)); }

// ******** Pulse Propagators In the XY Plane On Axis Phi Degrees From X ******

gen_op Ixypuls_U(const spin_sys& sys, int spin,        double phi, double beta)
  { return gen_op(Rxy(sys, spin, phi, beta)); }

gen_op Ixypuls_U(const spin_sys& sys, const std::string& I, double phi, double beta)
  { return gen_op(Rxy(sys, I, phi, beta)); }

gen_op Ixypuls_U(const spin_sys& sys,                  double phi, double beta)
  { return gen_op(Rxy(sys, phi, beta)); }

gen_op Ixypuls_U(const spin_sys& sys, const flagvec& flags,
                                                       double phi, double beta)
  { return gen_op(Rxy(sys, flags, phi, beta)); }

gen_op Ixypuls_U_sp(const spin_sys& sys,               double phi, double beta)
  { return gen_op(Rxy_sp(sys, phi, beta)); }
 
#endif 							//  PulseI.cc
