/* nmr_tensor.h *********************************-*-c++-*-
**							**
**							**
** 	Project						**
**							**
**	GAMMA Library					**
**						 	**
**	NMR Library					**
**						 	**
**	Implementation   				**
**						 	**
**	Copyright (c) 1991			 	**
**	Scott Smith				 	**
**	Eidgenoessische Technische Hochschule	 	**
**	Labor fuer physikalische Chemie		 	**
**	8092 Zurich / Switzerland		 	**
**						 	**
**      $Header:
**						 	**
** Modifications:					**
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The NMR Library Provides Functions		**
**	for the Simulation of Magnetic Resonance	**
**	Experiments and Associated Mathematical		**
**	Capabilities.  This Particular Code 		**
**	Contains Functions for Tensor Operators.	**
**						 	**
*********************************************************/

#ifndef   Gnmr_tensor_h_		// Is file already included?
#  define Gnmr_tensor_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSys.h>
#include <Level1/SpinT.h>

// ____________________________________________________________________________
//                       DIPOLAR SPHERICAL TENSOR FUNCTIONS
// ____________________________________________________________________________


MSVCDLL spin_T T_D(const spin_sys &sys, int spin1, int spin2);

	// Input	sys   : spin system
	// 		spin1 : spin index
	//	 	spin2 : spin index
	// Output	SphT  : Irreducible spherical rank 2 spin ]  D
	//			tensor for the dipolar interaction] T (i,j)
	//			between the two spins specified   ]  2


MSVCDLL spin_T T_D(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
				spin_op &Im2, spin_op &Iz2, spin_op &Ip2);

	// Input	sys   : spin system
	// 		Im1   : spin operator I- for first spin
	// 		Iz1   : spin operator Iz for first spin
	// 		Ip1   : spin operator I+ for first spin
	// 		Im2   : spin operator I- for second spin
	// 		Iz2   : spin operator Iz for second spin
	// 		Ip2   : spin operator I+ for second spin
	// Output	SphT  : Irreducible spherical rank 2 spin ]  D
	//			tensor for the dipolar interaction] T (1,2)
	//			between the two spins specified   ]  2


MSVCDLL spin_op T_D(const spin_sys &sys, int spin1, int spin2, int m);

	// Input	sys   : spin system
	// 		spin1 : spin index
	//	 	spin2 : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of]  D
	//			rank 2 spin tensor for the dipolar] T  (i,j)
	//			interaction of specified 2 spins  ]  2m


// ______________________________________________________________________
//              CHEMICAL SHIFT SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


MSVCDLL spin_T T_CSA(const spin_sys &sys, int spin);

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : irreducible spherical rank 2 "spin"]  CSA
	//			tensor for the anisotropic part of ] T   (i)
	//			of thechemical shift interaction   ]  2
	//			for the spin specified
	// Note		      : Assumes the magnetic field is aligned


MSVCDLL spin_T T_CS2(const spin_sys &sys, int spin);

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : spherical rank 2 "spin" tensor     ]  CS
	//			for the chemical shift interaction ] T  (i)
	//			for the spin specified             ]  2
	// Note		      : Assumes the magnetic field is aligned
	//			with the frame z-axis


MSVCDLL spin_T T_CS2(const spin_sys &sys, int spin, coord &B);

	// Input	sys   : spin system
	// 		spin  : spin index
	//		B     : magnetic field vector (Cartesian)
	// Output	SphT  : spherical rank 2 "spin" tensor     ]  CS
	//			for the chemical shift interaction ] T  (i)
	//			for the spin specified             ]  2

MSVCDLL spin_op T_CS2(const spin_sys &sys, int spin, coord &B, int l, int m);

	// Input	sys   : spin system
	// 		spin  : spin index
	//		B     : magnetic field vector (Cartesian)
	//              l     : rank index
	//              m     : component index
	// Output	SphT  : irreducible spherical rank 2 "spin"    ]  CS
	//			tensor component for the chemical      ] T  (i)
	//			shift interaction for spin specified   ]  lm

MSVCDLL spin_T T_CS(const spin_sys &sys, int spin);

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : Irreducible spherical rank 1 spin ]  CS
	//			tensor for the chemical shift     ] T  (i)
	//			interaction for spin specified    ]  1

MSVCDLL spin_op T_CS(const spin_sys &sys, int spin, int m);

	// Input	sys   : spin system
	// 		spin  : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of  ]  CS
	//			rank 1 spin tensor for the chemical ] T  (i)
	//			shift interaction for specified spin]  1m

// ______________________________________________________________________
//               RANDOM FIELD SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


MSVCDLL spin_T T_RF(const spin_sys &sys, int spin);

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : Irreducible spherical rank 1 spin ]  RF
	//			tensor for the random field       ] T  (i)
	//			interaction for spin specified    ]  1


MSVCDLL spin_op T_RF(const spin_sys &sys, int spin, int l, int m);

	// Input	sys   : spin system
	// 		spin  : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of  ]  CS
	//			rank 1 spin tensor for the random   ] T  (i)
	//			field interaction for specified spin]  1m


// ______________________________________________________________________
//                QUADRUPOLAR SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


MSVCDLL spin_T T_Q(const spin_sys &sys, int spin);

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : Irreducible spherical rank 2 spin     ]  Q
	//			tensor for the quadrupolar interaction] T (i)
	//			for spin specified                    ]  2


MSVCDLL spin_op T_Q(const spin_sys &sys, int spin, int l, int m);

	// Input	sys   : spin system
	// 		spin  : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of    ]  Q
	//			rank 2 spin tensor for the quadrupolar] T  (i)
	//			interaction of specified spin         ]  2m

#endif						// nmr_tensor.h
