/* nmr_tensor.cc ********************************-*-c++-*-
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
**      $Header: $
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The NMR Library Provides Functions		**
**	for the Simulation of Magnetic Resonance	**
**	Experiments and Associated Mathematical		**
**	Capabilities.  These Functions Allow for	**
**	Use of the Tensors Associated with NMR		**
**	Hamiltonians.				 	**
**						 	**
*********************************************************/

#ifndef   Nmr_tensor_cc_		// Is file already included?
#  define Nmr_tensor_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Level1/nmr_tensor.h>
#include <HSLib/SpinOpCmp.h>

using std::cout;			// Using libstdc++ standard output

// ______________________________________________________________________
//                   DIPOLAR SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


spin_T T_D(const spin_sys &sys, int spin1, int spin2)

	// Input	sys   : spin system
	// 		spin1 : spin index
	//	 	spin2 : spin index
	// Output	SphT  : Irreducible spherical rank 2 spin ]  D
	//			tensor for the dipolar interaction] T (i,j)
	//			between the two spins specified   ]  2

  {
  if(spin1 == spin2)
    cout << "NMRLIB: Dipolar Tensor Component Between Same Spins Requested";
  spin_op Im1,Iz1,Ip1;			// Compute spin tensors for spin i
  Im1 = Im(sys,spin1);
  Iz1 = Iz(sys,spin1);
  Ip1 = Ip(sys,spin1);
  spin_op Im2,Iz2,Ip2;			// Compute spin tensors for spin j
  Im2 = Im(sys,spin2);
  Iz2 = Iz(sys,spin2);
  Ip2 = Ip(sys,spin2);
  return T22(sys,Im1,Iz1,Ip1,Im2,Iz2,Ip2);
  }


spin_T T_D(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
				spin_op &Im2, spin_op &Iz2, spin_op &Ip2)

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

  { return T22(sys,Im1,Iz1,Ip1,Im2,Iz2,Ip2); }


spin_op T_D(const spin_sys &sys, int spin1, int spin2, int m)

	// Input	sys   : spin system
	// 		spin1 : spin index
	//	 	spin2 : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of]  D
	//			rank 2 spin tensor for the dipolar] T  (i,j)
	//			interaction of specified 2 spins  ]  2m

  {
  if(spin1 == spin2)
    cout << "NMRLIB: Dipolar Tensor Component Between Same Spins Requested";
  return T22(sys, spin1, spin2, m);
  }


// ____________________________________________________________________________
//                  CHEMICAL SHIFT SPHERICAL TENSOR FUNCTIONS
// ____________________________________________________________________________


spin_T T_CSA(const spin_sys &sys, int spin) { return T22SSirr(sys,spin,UnitZ); }

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : irreducible spherical rank 2 "spin"]  CSA
	//			tensor for the anisotropic part of ] T   (i)
	//			of thechemical shift interaction   ]  2
	//			for the spin specified
	// Note		      : Assumes the magnetic field is aligned
//  {
//  coord B(0,0,1);
//  return T22SSirr(sys, spin, B);
//  }


spin_T T_CS2(const spin_sys &sys, int spin) { return T2(sys, spin, UnitZ); }

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : spherical rank 2 "spin" tensor     ]  CS
	//			for the chemical shift interaction ] T  (i)
	//			for the spin specified             ]  2
	// Note		      : Assumes the magnetic field is aligned
	//			with the frame z-axis

//{
//  coord B(0,0,1);
//  return T2(sys, spin, B);
//}


spin_T T_CS2(const spin_sys &sys, int spin, coord &B)

	// Input	sys   : spin system
	// 		spin  : spin index
	//		B     : magnetic field vector (Cartesian)
	// Output	SphT  : spherical rank 2 "spin" tensor     ]  CS
	//			for the chemical shift interaction ] T  (i)
	//			for the spin specified             ]  2

  { return T2(sys, spin, B); }


spin_op T_CS2(const spin_sys &sys, int spin, coord &B, int l, int m)

	// Input	sys   : spin system
	// 		spin  : spin index
	//		B     : magnetic field vector (Cartesian)
	//              l     : rank index
	//              m     : component index
	// Output	SphT  : irreducible spherical rank 2 "spin"    ]  CS
	//			tensor component for the chemical      ] T  (i)
	//			shift interaction for spin specified   ]  lm

{
// add a check that spins exist
  return T2(sys, spin, B, l, m);
}


spin_T T_CS(const spin_sys &sys, int spin)

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : Irreducible spherical rank 1 spin ]  CS
	//			tensor for the chemical shift     ] T  (i)
	//			interaction for spin specified    ]  1

{
// add a check that spins exist
return T11(sys, spin);
}


spin_op T_CS(const spin_sys &sys, int spin, int m)

	// Input	sys   : spin system
	// 		spin  : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of  ]  CS
	//			rank 1 spin tensor for the chemical ] T  (i)
	//			shift interaction for specified spin]  1m

{
// add a check that spins exist
return T11(sys, spin, m);
}


// ______________________________________________________________________
//               RANDOM FIELD SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


spin_T T_RF(const spin_sys &sys, int spin)

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : Irreducible spherical rank 1 spin ]  RF
	//			tensor for the random field       ] T  (i)
	//			interaction for spin specified    ]  1

{
// add a check that spins exist
return T11(sys, spin);
}


spin_op T_RF(const spin_sys &sys, int spin, int l, int m)

	// Input	sys   : spin system
	// 		spin  : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of  ]  CS
	//			rank 1 spin tensor for the random   ] T  (i)
	//			field interaction for specified spin]  1m

{
// add a check that spins exist
return T11(sys, spin, m);
l=0;				// Compiler likes this used
}


// ______________________________________________________________________
//                QUADRUPOLAR SPHERICAL TENSOR FUNCTIONS
// ______________________________________________________________________


spin_T T_Q(const spin_sys &sys, int spin)

	// Input	sys   : spin system
	// 		spin  : spin index
	// Output	SphT  : Irreducible spherical rank 2 spin     ]  Q
	//			tensor for the quadrupolar interaction] T (i)
	//			for spin specified                    ]  2

{
// add a check that spins exist
// add a check that the spin has I >= 1
  spin_op Im1,Iz1,Ip1;			// Compute spin tensors for spin j
  Im1 = Im(sys,spin);
  Iz1 = Iz(sys,spin);
  Ip1 = Ip(sys,spin);
  return T22(sys,Im1,Iz1,Ip1,Im1,Iz1,Ip1);
//return T22(sys, spin, spin);
}


spin_op T_Q(const spin_sys &sys, int spin, int l, int m)

	// Input	sys   : spin system
	// 		spin  : spin index
	//		m     : component index
	// Output	SOp   : Irreducible spherical component of    ]  Q
	//			rank 2 spin tensor for the quadrupolar] T  (i)
	//			interaction of specified spin         ]  2m

{
// add a check that spins exist
// add a check that the spin has I >= 1
return T22(sys, spin, spin, m);
l=0;				// Compiler likes this used
}

#endif						// nmr_tensor.cc

