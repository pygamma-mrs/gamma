/* HSauxil.h ************************************-*-c++-*-
**							**
**	              G A M M A 			**
**						 	**
**	NMR Library		Interface Definition	**
**						 	**
**	Copyright (c) 1991, 1992, 1993		 	**
**	Scott Smith				 	**
**	Eidgenoessische Technische Hochschule	 	**
**	Labor fuer physikalische Chemie		 	**
**	8092 Zurich / Switzerland		 	**
**						 	**
**      $Header: $
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
**	Capabilities. Most of the Library consists	**
**	in modules called nmr_"name", where "name"	**
**	indicates the type of functions contained	**
**	therein.  Include here are functions of a	**
**	more general nature not yet associated with	**
**	a group of routines.				**
**						 	**
*********************************************************/

#ifndef   Nmrlib_h_			// Is file already included?
#  define Nmrlib_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSystem.h>		// Know about spin systems
#include <HSLib/GenOp.h>		// Include operators
#include <Matrix/complex.h>		// Know about GAMMA complex #s

// ______________________________________________________________________
// A                     Denstity Matrix Functions
// ______________________________________________________________________

MSVCDLL gen_op sigma_eq(const spin_sys& sys);
MSVCDLL gen_op sigma_eq(const spin_sys& sys, const Isotope& I);

	///F_list sigma_eq	- high temperature equilibrium density matrix
	// Input		ss : spin system
	// Output		Op : high temperature equilibrium
	//			     density matrix


//gen_op sigma_eq(const spin_sys &sys, basis& bs);
	///F_list sigma_eq	- high temperature equilibrium density matrix
	// Input		bs : basis
	// Output		Op : high temperature equilibrium
	//			     density matrix in basis bs

// ______________________________________________________________________
// B                    Density Operator Coherence Selection
// ______________________________________________________________________


MSVCDLL void zero_mqc(const spin_sys &sys, gen_op &Op, int order, int type);

	///F_List zero_mqc	- Zeroes density matrix elements with specific
	///			  Coherences
	// Input	sys	: Spin system
	// Input	Op	: General operator (associated to sys)
	//		order	: Coherence order
	//			    <0 - populations (default)
	//			     0 - zero quantum coherence
	//		             1 - single quantum coherence
	//			     n - n quantum coherence		
	//		type	: Type of zeroing
	//			    <0 - zero coherences except "order" (default)
	//			     0 - zero coherences "order"
	//			     1 - zero coherences >= "order"
	//			    >1 - zero coherences <= "order"
	// Output	none	: Function is void.  Op has had
	//			  elements associated with one or
	//			  more coherence types zeroed


// ______________________________________________________________________
// ********************* Single Transition Operators ********************
// ______________________________________________________________________


MSVCDLL gen_op st_Op(gen_op &Ham, int lev1, int lev2, char axis);

	// Input		Ham  : General operator (Hamiltonian).
	// 			lev1 : Energy level 1
	// 			lev2 : Energy level 2
	// 			type : Transition operator type
	//				x = Ix 		p = I+
	//				y = Iy 		m = I-
	//				z = Iz
	// Return		Op   : General operator which is a
	//			       single transition operator
	//			       for the transtion between levels
	//			       lev1 and lev2.
	// Note			     : Result EXCLUSIVELY in EBR of Ham
	// Note			     : See EBW, page 34.


MSVCDLL void sqt_v(gen_op &Ham);

	// Input		Ham  : General operator (Hamiltonian).
	// Return		none : All transitions associated with
	//			       the input Hamiltonian are set to
	//			       standard output
	// Note			     : Sets Ham to its EBR


// ______________________________________________________________________
//                      VARIOUS USEFUL FUNCTIONS
// ______________________________________________________________________


// int read_units (std::string xstr, double &x);

	// Input		xstr	: string containing input value
	//			x       : actual value
	// Output		int	: TRUE if time interpreted O.K.
	//				  FALSE if not
	// Note				: time is set (if TRUE) as in timestr
	// Note				: error message if time not determined
	// Note				: default units are microseconds

// ----------------------------------------------------------------------
//                      Sorting the basis 
// ----------------------------------------------------------------------

MSVCDLL   int* sort_super_op_basis (const spin_sys& sys);
// ??? Shortened name, should remove this function

MSVCDLL   int* sort_LOp_basis (const spin_sys& sys);

	// Input		sys	: Spin system
	// Output		index	: An integer vector containing
	//				  the basis function order which
	//				  has coherence sorting.
	// Note				: Output order has populations 1st,
	//				  followed by ZQC, +SQC, -SQC, ...
	//				  up to the highest coherenece order
	//				  of the spin system
	// Note				: The default basis order is 
	//				  0,1,2,...,LS-1 and this is
	//				  usually not coherence ordered due
	//				  to the tensor products used in
	//				  formulating super operators
	///F_list sort_LOp_basis	- Returns coherence sorted LOp indices


MSVCDLL int* sort_Op_basis (const spin_sys& sys);

	// Input		sys	: Spin system
	// Output		index	: An integer vector containing
	//				  the basis function order from
	//				  from highest Fz to lowest
	// Note				: The default basis order is 
	//				  0,1,2,...,HS-1 and this is
	//				  usually not Fz ordered due to
	//				  single spin tensor products used
	//				  in formulating spin operators
	///F_list sort_Op_basis		- Returns Fz sorted Op indices


MSVCDLL void mqt_v(const spin_sys& sys, gen_op &Ham, int qn, int type, int ncols);

	// Input	sys	: A spin system
	// 		Ham	: General operator (associated to sys)
	//		qn	: Transition quantum number
	//			     0 - zero quantum transitions
	//		             1 - single quantum transitions (DEFAULT)
	//			     n - n quantum transitions		
	//		type	: Type of transition output
	//			    <0 - all transitions except qn quantum transitions
	//			     0 - only "qn" quantum transitions (DEFAULT)
	//			     1 - all transitions >= "qn" quantum transitions
	//			    >1 - all transitions <= "qn" quantum transitions
	// Return	none 	: All indicated transitions associated with
	//			  the input Hamiltonian are sent to standard output
	// Note		     	: Sets Ham to its EBR


MSVCDLL void wavefunction(const spin_sys& sys, gen_op &Op, int wf, int pbf);

	// Input	sys	: A spin system
	// 		Op	: General operator (associated to sys)
	//		wf	: Wavefunction number
	//		pbf	: Flag to specify how to print product basis function
	//			  -1 - designate product basis functions by total Fz
	//			   0 - designate product basis functions by number (DEFAULT)
	//			   1 - designate product basis functions by alpha(a) & beta(b)
// *** pbg=1 only works for systems with all I=1/2.
	// Return	none 	: Indicated wavefunction wf for the working basis
	//			  of operator Op is printed to standard output


MSVCDLL void wavefunctions(const spin_sys& sys, gen_op &Op, int pbf);

	// Input	sys	: A spin system
	// 		Op	: General operator (associated to sys)
	//			  -1 - designate product basis functions by total Fz
	//			   0 - designate product basis functions by number (DEFAULT)
	//			   1 - designate product basis functions by alpha(a) & beta(b)
// *** pbg=1 only works for systems with all I=1/2.
	// Return	none 	: All wavefunctions of the working basis
	//			  of operator Op are printed to standard output


MSVCDLL void eigensystem(std::ostream& ostr, gen_op Op);

	// Input 	Op	: General operator
	//			  -1 - designate product basis functions by total Fz
	//			   0 - designate product basis functions by number (DEFAULT)
	//			   1 - designate product basis functions by alpha(a) & beta(b)
	// Return	none 	: All wavefunctions of the working basis
	//			  of operator Op are printed to standard output

MSVCDLL double vecmax(row_vector &vx);

	// Input	vx       : Data vector
	//		i        : Empty integer
	//              max      : Empty double
	// Return		 : i & max are filled with
	//			   the vector values at it maximum


MSVCDLL complex integral(const row_vector& vx);

	// Input		vx  : Data vector
	// Return		z   : Integral, sum of all vector elements
	//			      between the limits specified.
	// Note			    : Range includes high & low [low, high]
	// Note			    : It is allowed that low = high


MSVCDLL void lwhh(row_vector &vx, int& i1, int& i2);

	// Input	vx       : Data vector
        //		ri	 : Flag for real versus imaginary
	//			   0=reals, non-zero=imaginaries
	// Return		 : Min and Max are returned


// ______________________________________________________________________
//                  PULSE SEQUENCE PARAMETER QUERY FUNCTIONS
// ______________________________________________________________________


MSVCDLL int query_isotope(const spin_sys& sys, std::string& Isotype);

	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type
	// Note       		Isotype : Isotope type

MSVCDLL int query_isotope(const spin_sys& sys, std::string& Isotype, const std::string& Query);

	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	//       		Query	: Question to be asked
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type
	// Note       		Isotype : Isotope type

MSVCDLL int query_isotope(int argc, char* argv[], int argn, const spin_sys& sys, std::string& Isotype);

	// Input		argc    : Number of arguments available
	//       		argv    : Vector of argument values
	// 			argn    : Argument index
	// 			sys     : Spin system
	//       		Isotype : Isotope type
	// Output		isoset  : Index of a spin having the
	//				  choosen isotope type


MSVCDLL double query_offset(spin_system& sys, int isoset, int askit=0);

	// Input		sys     : Spin system
	//       		isoset  : Index of a spin having the
	//				  choosen isotope type
	//			askit   : Flag to force query
	// Output		offset  : Carrier offset value (Hz) for
	//				  the choosen isotope type


MSVCDLL double query_offset(spin_system& sys, std::string& Isotype, int ask=0);


	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	//			ask     : Flag to force a query
	//				   0 - ask if non-zero center
	//				  >0 - always ask
	//				  <0 - never ask
	// Output		offset  : Carrier offset value (Hz) for
	//				  the choosen isotope type


MSVCDLL double query_Nyquist(spin_system& sys, int isoset, double lw=0, double fact=1.2);

	// Input		sys     : Spin system
	//       		isoset  : Index of a spin having the
	//				  choosen isotope type
	//			fact    : A scaling factor (1.2 = add 20%)
	// Output		Nyqf    : A choosen Nyquist frequency


MSVCDLL double query_Nyquist(const spin_system& sys, std::string& Isotype, double lw=0, double fact=1.2);

	// Input		sys     : Spin system
	//       		Isotype : Isotope type
	//			lw	: An expected linewidth
	//			fact    : A scaling factor (1.2 = add 20%)
	// Output		Nyqf    : An estimated Nyquist frequency


MSVCDLL void query_file1D(std::string& filename, int& type);

	// Input	filename: I/O filename
	//		type    : I/O file type
	// Output	none	: Function is void.  The string filename
	//		          is filled with an input name.  The type of
	//			  file is also set to one of the following:
	//			     1 - FrameMaker .mif format
	//			     2 - Felix .dat format
	//			     3 - NMRi format
	//			     4 - MatLab format


MSVCDLL void query_FelixFile1D(std::ostream& ostr, std::string& filename);

	// Input	filename: I/O filename
	// Output	none	: Function is void.  It asks for a filename
	//		          assuming the file will be a 1D Felix file


//void query_output1D(const spin_system& sys, std::string& filename, int& type, int FFT);

	// Input	filename: I/O filename
	// Input	filename: I/O filename
	//		type    : I/O file type
	//			     1 - FrameMaker .mif format
	//			     2 - Felix .dat format
	//			     3 - NMRi format
	//			     4 - MatLab format
	//		sw	: spectral width (in Hz)
	//		offset  : Spectral offset (in Hz)
	//		Omega	:
	//		FFT	: Flag for application of a Fourier Transform
	// Output	none	: Function is void.  The string filename
	//		          is filled with an input name.  The type of
	//			  file is also set to one of the following:


MSVCDLL void Felix1D_params(std::ostream& ostr, double Omega, double sw, int npts, double offset);

	// Input	ostr    : Output stream
	//			: Spectrometer frequency (MHz)
	//			: Spectral width (Hz)
	//			: Number of points
	//			: Offset frequency (Hz)
	// Output	none	: Function is void.  Felix parameters for 1D
	//		          spectral workup are sent into the output stream


MSVCDLL void Felix2D_params(std::ostream& ostr, double Omega, double sw, int npts, double offset, int PPM);

	// Input	ostr    : Output stream
	//			: Spectrometer frequency (MHz)
	//			: Spectral width (Hz)
	//			: Number of points
	//			: Offset frequency (Hz)
	// Output	none	: Function is void.  Felix parameters for 2D
	//		          2D homonuclear spectral workup are sent
        //                        into the output stream


MSVCDLL void Felix2D_params(std::ostream& ostr, double O2, double sw2, int npts2, double off2,
                                   double O1, double sw1, int npts1, double off1, int PPM);

	// Input	ostr    : Output stream
	//		O1	: Spectrometer frequency (MHz)
	//			: Spectral width (Hz)
	//			: Number of points
	//			: Offset frequency (Hz)
	// Output	none	: Function is void.  Felix parameters for 2D
	//		          2D homonuclear spectral workup are sent
    //                        into the output stream

 
#endif						// HSauxil.h 
