/* relaxExch.h  **************************************************
**								**
**                           G A M M A 				**
**                                				**
**      NMR Exchange                         Interface		**
**								**
**      Copyright (c) 1995					**
**      Scott A. Smith						**
**      Nikolai Skrynnikov					**
**                                				**
**      National High Magnetic Field Laboratory			**
**      1800 E. Paul Dirac Drive				**
**      Tallahassee, Florida, 32310				**
**                                				**
**      $Header: $
**                                				**
**                                				**
*****************************************************************/

/*****************************************************************
**								**
** Description							**
**								**
** The following functions herein provide easy access to many	**
** common exchange superoperators.  Typically, mutual exchange	**
** functions will demand a sys_dynamic variable as a argument.	**
** For non-mutual exchange functions a multi_sys variable is	**
** necessary.  							**
**								**
*****************************************************************/

//# include <super_op.h>

#ifndef   Relax_EXCH_h_			// Is this file already included?
#  define Relax_EXCH_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)

class gen_op;				// Know class gen_op
class super_op;				// Know class super_op
class sys_dynamic;			// Know class sys_dynamic

MSVCDLL super_op Rex(const sys_dynamic& sys);

	// Input		sys	: Dynamic spin system
	// Output		Rex     : Exchange superoperator
	// Note			        : The exchange information
	//				  matrix provided by sys
	//				  has an unique structure
	// Note				: The superoperator returned
	//				  from this function is for
	//				  intramolecular mutual exchange.

//           ^	                                  ^      ^
//           ^             ^     ^      ^   ^     ^      ^
//           X    = K  * [ K   x K    - E x E ] = K    - E
//            i,j    ex     i,j   i,j              i,j


MSVCDLL super_op Rex(const sys_dynamic& sys, gen_op Op);

	// Input		sys	: Dynamic spin system
	// 			Op	: Operator providing basis
	// Output		Rex     : Exchange superoperator
	// Note			        : The exchange information
	//				  matrix provided by sys
	//				  has an unique structure
	// Note				: The superoperator returned
	//				  from this function is for
	//				  intramolecular mutual exchange.

//           ^	                                  ^      ^
//           ^             ^     ^      ^   ^     ^      ^
//           X    = K  * [ K   x K    - E x E ] = K    - E
//            i,j    ex     i,j   i,j              i,j


#endif						// relaxExch.h

