#include <gamma.h>
#include <vector>

void STJacobian3(const vector<double>& A, double delz, double eta,
                               const col_vector& X, col_vector& F, matrix& J)

	// Input		X	: Vector - {alpha,beta,gamma}
	//			A	: Vector - {Axx,Axy,Axz}
	//			F       : Vector - {f1,f2,f3}
	//			J	: Jacobian Matrix (3x3)
	// Output		void	: Given the values in X & A the
	//				  vector F and Jacobian matrix J
	//				  are generated
	// Note				: FOR FITTING {alpha,beta,gamma}
	//				  WHEN WE KNOW {delzz,eta}

  {
//			        Explicitly Get Our Input Variables 

  double Axx   = A[0];       double Axy  = A[1];       double Axz   = A[2]; 
  double alpha = X.getRe(0); double beta = X.getRe(1); double gamma = X.getRe(2);

//                                 Variables Involving Only Alpha

  double S2a    = sin(2.0*alpha);			// sin(2a)
  double C2a    = cos(2.0*alpha);			// cos(2a)

//                                 Variables Involving Only Beta

  double Cb        = cos(beta);				// sin(b)
  double Sb        = sin(beta);				// cos(b)
  double C2b       = cos(2.0*beta);			// sin(2b)
  double S2b       = sin(2.0*beta);			// cos(2b)
  double SbCb      = Sb*Cb;				// sin(b)*cos(b)
  double Ssqb      = Sb*Sb;				// sin(b)*sin(b)
  double Csqb      = Cb*Cb;				// cos(b)*cos(b)
  double CbSb      = Cb*Sb;				// cos(b*sin(b)
  double Csqbp1    = Csqb + 1.0;			// cos(b)*cos(b)+1
  double CsqbmSsqb = Cb*Cb-Sb*Sb;			// cos(b)*cos(b)-sin(b)*sin(b) 

//                                 Variables Involving Only Gamma

  double S2g    = sin(2.0*gamma);			// sin(2g)
  double C2g    = cos(2.0*gamma);			// cos(2g)

//                                 Variables Involving Only delzz
 
  double  delzzo4 =  0.25*delz;
  double mdelzzo4 = -0.25*delz;


  double term1, term2, term3;

//                       The 1st Function And Elements Of The 1st Jacobian Row
 
  term1 =  C2a*(3.0*Ssqb - eta*Csqbp1*C2g);		//  cos(2a)*{3.0*sin(b)*sin(b) + eta*[cos(b)*cos(b)+1.0]cos(2g)}
  term2 = -2.0*eta*S2a*Cb*S2g;				// -2.0*eta*sin(2a)*cos(b)*sin(2g)
  term3 = -3.0*Csqb + 1.0 - eta*Ssqb*C2g;		// -3.0*cos(b)*cos(b) + 1.0 - eta*sin(b)*sin(b)*cos(2g)
  double f1 =  delzzo4*(term1 + term2 + term3) - Axx;	// This is the function f1
  term1 = -6.0*S2a*Ssqb;				// -6.0*sin(2a)*sin(b)*sin(b)
  term2 = -2.0*eta*S2a*Csqbp1*C2g;			// -2.0*eta*sin(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 = -4.0*eta*C2a*Cb*S2g;				// -4.0*eta*cos(2a)*cos(b)*sin(2g)
  double df1da   =  delzzo4*(term1 + term2 + term3);	// This is @f1/@alpha
  term1 =  C2a*(6.0*SbCb - 2.0*eta*SbCb*C2g);		//  cos(2a)*{6.0*sin(b)*cos(b) - 2.0*eta*sin(b)*cos(b)*cos(2g)}
  term2 =  2.0*eta*S2a*Sb*S2g;				//  2.0*eta*sin(2a)*sin(b)*sin(2g)
  term3 =  6.0*CbSb - 2.0*eta*CbSb*C2g;			//  6.0*sin(b)*cos(b) - 2.0*eta*sin(b)*cos(b)*cos(2g)
  double df1db   =  delzzo4*(term1 + term2 + term3);	// This is @f1/@beta
  term1 = -2.0*eta*C2a*Csqbp1*S2g;			// -2.0*eta*cos(2a)*[cos(b)*cos(b)+1.0]*sin(2g)
  term2 = -4.0*eta*S2a*Cb*C2g;				// -4.0*eta*sin(2a)*cos(b)*cos(2g)
  term3 =  2.0*eta*Ssqb*S2g;				//  2.0*eta*sin(b)*sin(b)*sin(2g)
  double df1dg   =  delzzo4*(term1 + term2 + term3);	// This is @f1/@gamma

//                       The 2nd Function And Elements Of The 2nd Jacobian Row
 
  term1 =  3.0*S2a*Ssqb;;				//  3.0*sin(2a)*sin(b)*sin(b)
  term2 =  eta*S2a*Csqbp1*C2g;				//  eta*sin(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 = -2.0*eta*C2a*Cb*S2g;				// -2.0*eta*cos(2a)*cos(b)*sin(2g)
  double f2 =     delzzo4*(term1 + term2 + term3) - Axy;// This is the function f2
  term1 =  6.0*C2a*Ssqb;				//  6.0*cos(2a)*sin(b)*sin(b)
  term2 =  2.0*eta*C2a*Csqbp1*C2g;			//  2.0*eta*cos(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 =  4.0*eta*S2a*Cb*S2g;				//  4.0*eta*sin(2a)*cos(b)*sin(2g)
  double df2da   =  delzzo4*(term1 + term2 + term3);	// This is @f2/@alpha
  term1 =  6.0*S2a*SbCb;				//  6.0*sin(2a)*sin(b)*cos(b)		
  term2 = -2.0*eta*S2a*SbCb*C2g;			// -2.0*eta*sin(2a)*sin(b)*cos(b)*cos(2g)
  term3 =  2.0*eta*C2a*Sb*S2g;				//  2.0*eta*cos(2a)*sin(b)*sin(2g)
  double df2db   =  delzzo4*(term1 + term2 + term3);	// This is @f2/@beta
  term1 =  0.0;			
  term2 = -2.0*eta*S2a*Csqbp1*S2g;			// -2.0*eta*sin(2a)*[cos(b)*cos(b)+1]*sin(2g)
  term3 = -4.0*eta*C2a*Cb*C2g;				// -4.0*eta*cos(2a)*cos(b)*cos(2g)
  double df2dg   =  delzzo4*(term1 + term2 + term3);	// This is @f2/@gamma

//                       The 3rd Function And Elements Of The 3rd Jacobian Row
  term1 = -3.0*C2a*S2b;;				// -3.0*cos(2a)*sin(b)
  term2 =  2.0*eta*C2a*SbCb*C2g;			//  2.0*eta*cos(2a)*sin(b)*cos(b)*cos(2g)
  term3 = -1.0*eta*S2a*S2g;				// -1.0*eta*sin(2a)*sin(2g)
  double f3 =    mdelzzo4*(term1 + term2 + term3) - Axz;// This is function f3
  term1 = -6.0*S2a*S2b;					// -6.0*sin(2a)*sin(2b)
  term2 = -4.0*eta*S2a*SbCb*C2g;			// -4.0*eta*sin(2a)*sin(b)*cos(b)*cos(2g)
  term3 = -2.0*eta*C2a*S2g;				// -2.0*eta*cos(2a)*sin(2g)
  double df3da   = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@alpha
  term1 =  0.0;			
  term2 = -6.0*C2a*C2b;					// -6.0*cos(2a)*cos(2b)
  term3 =  2.0*eta*CsqbmSsqb*C2g;			//  2.0*eta*[sin(b)*sin(b)-cos(b)*cos(b)]*cos(2g)
  double df3db   = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@beta
  term1 =  0.0;			
  term2 = -4.0*eta*C2a*SbCb*S2g;			// -4.0*eta*cos(2a)*sin(b)*cos(b)*sin(2g)
  term3 = -2.0*eta*S2a*C2g;				// -2.0*eta*sin(2a)*cos(2g)
  double df3dg   = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@gamma

//                                  Build The Jacobian Matrix

  J = matrix(3,3);
  J.put(df1da,0,0); J.put(df1db,0,1); J.put(df1dg,0,2);
  J.put(df2da,1,0); J.put(df2db,1,1); J.put(df2dg,1,2);
  J.put(df3da,2,0); J.put(df3db,2,1); J.put(df3dg,2,2);

//					Build The F Vector

  F = col_vector(3);
  F.put(f1,0);
  F.put(f2,1);
  F.put(f3,2);

  return;
  }


void STJacobian5(const vector<double>& A, const col_vector& X, col_vector& F, matrix& J)

	// Input		X	: Vector - {delzz,eta,alpha,beta,gamma}
	//			A	: Vector - {Axx,Axy,Axz,Ayy,Ayz}
	//			F       : Vector + {f1,f2,f3,f4,f5}
	//			J	: Jacobian Matrix (5x5)
	// Output		void	: Given the values in X & A the
	//				  vector F and Jacobian matrix J
	//				  are generated
	// Note				: FOR FITTING {delzz,eta,alpha,beta,gamma}
	// Note				: The angles in X must be in radians
	//				: & the delzz value units in X should match
	//				  the units used in A

  {
//			        Explicitly Get Our Input Variables 

  double Axx  = A[0]; double Axy = A[1]; double Axz   = A[2]; double Ayy  = A[3]; double Ayz   = A[4];
  double delz = X.getRe(0); double eta = X.getRe(1); double alpha = X.getRe(2); double beta = X.getRe(3); double gamma = X.getRe(4);

//                                 Variables Involving Only Alpha

  double Sa     = sin(alpha);				// sin(a)
  double Ca     = cos(alpha);				// cos(a)
  double S2a    = sin(2.0*alpha);			// sin(2a)
  double C2a    = cos(2.0*alpha);			// cos(2a)

//                                 Variables Involving Only Beta

  double Cb        = cos(beta);				// sin(b)
  double Sb        = sin(beta);				// cos(b)
  double C2b       = cos(2.0*beta);			// sin(2b)
  double S2b       = sin(2.0*beta);			// cos(2b)
  double SbCb      = Sb*Cb;				// sin(b)*cos(b)
  double Ssqb      = Sb*Sb;				// sin(b)*sin(b)
  double Csqb      = Cb*Cb;				// cos(b)*cos(b)
  double CbSb      = Cb*Sb;				// cos(b*sin(b)
  double Csqbp1    = Csqb + 1.0;			// cos(b)*cos(b)+1
  double CsqbmSsqb = Cb*Cb-Sb*Sb;			// cos(b)*cos(b)-sin(b)*sin(b) 

//                                 Variables Involving Only Gamma

  double S2g    = sin(2.0*gamma);			// sin(2g)
  double C2g    = cos(2.0*gamma);			// cos(2g)

//                                 Variables Involving Only delzz
 
  double  delzzo4 =  0.25*delz;
  double mdelzzo4 = -0.25*delz;


  double term1, term2, term3;

//                       The 1st Function And Elements Of The 1st Jacobian Row
 
  term1 =  C2a*(3.0*Ssqb - eta*Csqbp1*C2g);		//  cos(2a)*{3.0*sin(b)*sin(b) + eta*[cos(b)*cos(b)+1.0]cos(2g)}
  term2 = -2.0*eta*S2a*Cb*S2g;				// -2.0*eta*sin(2a)*cos(b)*sin(2g)
  term3 = -3.0*Csqb + 1.0 - eta*Ssqb*C2g;		// -3.0*cos(b)*cos(b) + 1.0 - eta*sin(b)*sin(b)*cos(2g)
  double df1ddzz =     0.25*(term1 + term2 + term3);	// This is @f1/@delzz
  term1 =  C2a*Csqbp1*C2g;				// cos(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term2 = -2.0*S2a*Cb*S2g;				// -2.0*sin(2a)*cos(b)*sin(2g)
  term3 = -1.0*Ssqb*C2g;				// -1.0*sin(b)*sin(b)*cos(2g)
  double df1deta =  delzzo4*(term1 + term2 + term3);	// This is @f1/@eta
  term1 = -6.0*S2a*Ssqb;				// -6.0*sin(2a)*sin(b)*sin(b)
  term2 = -2.0*eta*S2a*Csqbp1*C2g;			// -2.0*eta*sin(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 = -4.0*eta*C2a*Cb*S2g;				// -4.0*eta*cos(2a)*cos(b)*sin(2g)
  double df1da   =  delzzo4*(term1 + term2 + term3);	// This is @f1/@alpha
  term1 =  C2a*(6.0*SbCb - 2.0*eta*SbCb*C2g);		//  cos(2a)*{6.0*sin(b)*cos(b) - 2.0*eta*sin(b)*cos(b)*cos(2g)}
  term2 =  2.0*eta*S2a*Sb*S2g;				//  2.0*eta*sin(2a)*sin(b)*sin(2g)
  term3 =  6.0*CbSb - 2.0*eta*CbSb*C2g;			//  6.0*sin(b)*cos(b) - 2.0*eta*sin(b)*cos(b)*cos(2g)
  double df1db   =  delzzo4*(term1 + term2 + term3);	// This is @f1/@beta

  term1 = -2.0*eta*C2a*Csqbp1*S2g;			// -2.0*eta*cos(2a)*[cos(b)*cos(b)+1.0]*sin(2g)
  term2 = -4.0*eta*S2a*Cb*C2g;				// -4.0*eta*sin(2a)*cos(b)*cos(2g)
  term3 =  2.0*eta*Ssqb*S2g;				//  2.0*eta*sin(b)*sin(b)*sin(2g)
  double df1dg   =  delzzo4*(term1 + term2 + term3);	// This is @f1/@gamma
  double f1      = delz*df1ddzz - Axx;			// This is function f1

//                       The 2nd Function And Elements Of The 2nd Jacobian Row
 
  term1 =  3.0*S2a*Ssqb;;				//  3.0*sin(2a)*sin(b)*sin(b)
  term2 =  eta*S2a*Csqbp1*C2g;				//  eta*sin(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 = -2.0*eta*C2a*Cb*S2g;				// -2.0*eta*cos(2a)*cos(b)*sin(2g)
  double df2ddzz =     0.25*(term1 + term2 + term3);	// This is @f2/@delzz
  term1 =  0.0;			
  term2 =  S2a*Csqbp1*C2g;				// sin(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 = -2.0*C2a*Cb*S2g;				// -2.0*cos(2a)*cos(b)*sin(2g)
  double df2deta =  delzzo4*(term1 + term2 + term3);	// This is @f2/@eta
  term1 =  6.0*C2a*Ssqb;				//  6.0*cos(2a)*sin(b)*sin(b)
  term2 =  2.0*eta*C2a*Csqbp1*C2g;			//  2.0*eta*cos(2a)*[cos(b)*cos(b)+1.0]*cos(2g)
  term3 =  4.0*eta*S2a*Cb*S2g;				//  4.0*eta*sin(2a)*cos(b)*sin(2g)
  double df2da   =  delzzo4*(term1 + term2 + term3);	// This is @f2/@alpha
  term1 =  6.0*S2a*SbCb;				//  6.0*sin(2a)*sin(b)*cos(b)		
  term2 = -2.0*eta*S2a*SbCb*C2g;			// -2.0*eta*sin(2a)*sin(b)*cos(b)*cos(2g)
  term3 =  2.0*eta*C2a*Sb*S2g;				//  2.0*eta*cos(2a)*sin(b)*sin(2g)
  double df2db   =  delzzo4*(term1 + term2 + term3);	// This is @f2/@beta
  term1 =  0.0;			
  term2 = -2.0*eta*S2a*Csqbp1*S2g;			// -2.0*eta*sin(2a)*[cos(b)*cos(b)+1]*sin(2g)
  term3 = -4.0*eta*C2a*Cb*C2g;				// -4.0*eta*cos(2a)*cos(b)*cos(2g)
  double df2dg   =  delzzo4*(term1 + term2 + term3);	// This is @f2/@gamma
  double f2      = delz*df2ddzz - Axy;			// This is function f2

//                       The 3rd Function And Elements Of The 3rd Jacobian Row
 
  term1 = -3.0*C2a*S2b;;				// -3.0*cos(2a)*sin(b)
  term2 =  2.0*eta*C2a*SbCb*C2g;			//  2.0*eta*cos(2a)*sin(b)*cos(b)*cos(2g)
  term3 = -1.0*eta*S2a*S2g;				// -1.0*eta*sin(2a)*sin(2g)
  double df3ddzz =    -0.25*(term1 + term2 + term3);	// This is @f3/@delzz
  term1 =  0.0;			
  term2 =  2.0*C2a*SbCb*C2g;				// 2.0*cos(2a)*sin(b)*cos(b)*cos(2g)
  term3 = -1.0*S2a*S2g;					// -1.0*sin(2a)*sin(2g)
  double df3deta = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@eta
  term1 = -6.0*S2a*S2b;					// -6.0*sin(2a)*sin(2b)
  term2 = -4.0*eta*S2a*SbCb*C2g;			// -4.0*eta*sin(2a)*sin(b)*cos(b)*cos(2g)
  term3 = -2.0*eta*C2a*S2g;				// -2.0*eta*cos(2a)*sin(2g)
  double df3da   = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@alpha
  term1 =  0.0;			
  term2 = -6.0*C2a*C2b;					// -6.0*cos(2a)*cos(2b)
  term3 =  2.0*eta*CsqbmSsqb*C2g;			//  2.0*eta*[sin(b)*sin(b)-cos(b)*cos(b)]*cos(2g)
  double df3db   = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@beta
  term1 =  0.0;			
  term2 = -4.0*eta*C2a*SbCb*S2g;			// -4.0*eta*cos(2a)*sin(b)*cos(b)*sin(2g)
  term3 = -2.0*eta*S2a*C2g;				// -2.0*eta*sin(2a)*cos(2g)
  double df3dg   = mdelzzo4*(term1 + term2 + term3);	// This is @f3/@gamma
  double f3      = delz*df3ddzz - Axz;			// This is function f3

//                       The 4th Function And Elements Of The 4th Jacobian Row
 
  term1 = C2a*(3.0*Ssqb + eta*C2g*Csqbp1);		// cos(2a)*{3.0*sin(b)*sin(b) + eta*cos(2g)*[1+cos(b)*cos(b)]}
  term2 = -2.0*eta*S2a*Cb*S2g;				// -2.0*eta*sin(2a)*cos(b)*sin(2g)
  term3 =  3.0*Csqb - 1.0 + eta*Ssqb*C2g;		// 3.0*cos(b)*cos(b) - 1.0 + eta*sin(b)*sin(b)*cos(2g)
  double df4ddzz =    -0.25*(term1 + term2 + term3);	// This is @f4/@delzz
  term1 = C2a*Csqbp1*C2g;				// cos(2a)*[cos(b)*cos(b) + 1.0]*cos(2g)
  term2 = -2.0*S2a*Cb*S2g;				// 2.0*sin(2a)*cos(b)*sin(2g)
  term3 = eta*Ssqb*C2g;					// eta*sin(b)*sin(b)*cos(2g)
  double df4deta = mdelzzo4*(term1 + term2 + term3);	// This is @f4/@eta
  term2 = -6.0*S2a*Ssqb;				// -6.0*sin(2a)*sin(b)*sin(b)
  term2 = -2.0*eta*S2a*Csqbp1*C2g;			// -2.0*eta*sin(2a)*[cos(b)*cos(b)+1]*cos(2g)
  term3 = -4.0*eta*C2a*Cb*S2g;				// -4.0*eta*cos(2a)*cos(b)*sin(2g)
  double df4da   = mdelzzo4*(term1 + term2 + term3);	// This is @f4/@alpha
  term2 =  6.0*C2a*SbCb - 2.0*eta*C2a*SbCb*C2g;		//  6.0*cos(2a)*sin(b)*cos(b)-2.0*eta*cos(2a)*S(b)*C(b)*cos(2g)
  term2 =  4.0*eta*S2a*Sb*S2g;				// -4.0*eta*sin(2a)*sin(b)*sin(2g)
  term3 = -6.0*SbCb + 2.0*eta*SbCb*C2g;			// -6.0*sin(b)*cos(b) + 2.0*eta*sin(b)*cos(b)*cos(2g)
  double df4db   = mdelzzo4*(term1 + term2 + term3);	// This is @f4/@beta
  term1 = -2.0*eta*C2a*Csqbp1*S2g;			// -2.0*eta*cos2a*[cos(b)*cos(b)-1]*sin(2g)
  term2 = -4.0*eta*S2a*Cb*S2g;				// -4.0*eta*sin2a*cos(b)*sin(2g)
  term3 = -2.0*eta*Ssqb*S2g;				// -2.0*eta*sin(b)*sin(b)*sin(2g)
  double df4dg   = mdelzzo4*(term1 + term2 + term3);	// This is @f4/@gamma
  double f4      = delz*df4ddzz - Ayy;			// This is function f4

//                       The 5th Function And Elements Of The 5th Jacobian Row
 
  term1 = -3.0*Sa*S2b;					// -3.0*sin(a)*sin(2b)
  term2 = 2.0*eta*Sa*SbCb*C2g;				//  2.0*eta*sin(a)*sin(b)*cos(b)*cos(2g)
  term3 = eta*Ca*Sb*S2g;				// eta*cos(a)*sin(b)*sin(2g)
  double df5ddzz =    -0.25*(term1 + term2 + term3);	// This is @f5/@delzz
  term1 = 0.0;
  term2 = 2.0*Sa*SbCb*C2g;				// 2.0*sin(a)*sin(b)*cos(b)*cos(2g)
  term3 = Ca*Sb*S2g;					// cos(a)*sin(b)*sin(2g)
  double df5deta = mdelzzo4*(term1 + term2 + term3);	// This is @f5/@eta
  term1 = -3.0*Ca*S2b + 2.0*eta*Sb*Cb*C2g;		// -3.0*cos(a)*sin(2b)
  term2 =  2.0*eta*Ca*SbCb*C2g;				//  2.0*eta*cos(a)*sin(b)*cos(b)*cos(2g)
  term3 = -1.0*eta*Sa*Sb*S2g;				// -1.0*eta*sin(a)*sin(b)*sin(2g)
  double df5da   = mdelzzo4*(term1 + term2 + term3);	// This is @f5/@alpha
  term1 = -6.0*Sa*C2b;					// -6.0*sin(a)*cos(2b)
  term2 = 2.0*eta*Sa*CsqbmSsqb*C2g;			//  2.0*eta*sin(a)*[cos(b)cos(b)-sin(b)sin(b)]*cos(2g)
  term3 = eta*Ca*Cb*S2g;				//  eta*cos(a)*cos(b)*sin(2g)
  double df5db   = mdelzzo4*(term1 + term2 + term3);	// This is @f5/@beta
  term1 =  0.0;
  term2 = -4.0*eta*Sa*SbCb*S2g;				// -4.0*eta*sin(a)*sin(b)*cos(b)*sin(2g)
  term3 =  2.0*eta*Ca*Sb*C2g;				//  2.0*eta*cos(a)*sin(b)*cos(2g)
  double df5dg   = mdelzzo4*(term1 + term2 + term3);	// This is @f5/@gamma
  double f5      = delz*df5ddzz - Ayz;			// This is function f5

//                                  Build The Jacobian Matrix

  J = matrix(5,5);
  J.put(df1ddzz,0,0); J.put(df1deta,0,1); J.put(df1da,0,2); J.put(df1db,0,3); J.put(df1dg,0,4);
  J.put(df2ddzz,1,0); J.put(df2deta,1,1); J.put(df2da,1,2); J.put(df2db,1,3); J.put(df2dg,1,4);
  J.put(df3ddzz,2,0); J.put(df3deta,2,1); J.put(df3da,2,2); J.put(df3db,2,2); J.put(df3dg,2,4);
  J.put(df4ddzz,3,0); J.put(df4deta,3,1); J.put(df4da,3,2); J.put(df4db,3,3); J.put(df4dg,3,4);
  J.put(df5ddzz,4,0); J.put(df5deta,4,1); J.put(df5da,4,2); J.put(df5db,4,3); J.put(df5dg,4,4);

//					Build The F Vector

  F = col_vector(5);
  F.put(f1,0);
  F.put(f2,1);
  F.put(f3,2);
  F.put(f4,3);
  F.put(f5,4);

  return;
  }


int main()
  {
//                  This Is A Typical Dipolar Spatial Cartesian Tensor (KHz)

/*
  matrix Dmx(3,3,h_matrix_type);
  Dmx.put  ( 0.79,0,0); Dmx.put_h( 0.00,0,1);
  Dmx.put_h( 0.00,0,2); Dmx.put  (-1.13,1,1);
  Dmx.put_h(-0.93,1,2); Dmx.put  ( 0.34,2,2);
*/
  matrix Dmx(3,3,d_matrix_type);
  Dmx.put  ( 0.79,0,0); 
  Dmx.put_h( 0.79,1,1);
  Dmx.put_h(-0.79-0.79,2,2);

gen_op Z(Dmx);
Z.set_EBR();
matrix DD = Z.get_mx();
double Dxx = DD.getRe(0,0);
double Dyy = DD.getRe(1,1);
double Dzz = DD.getRe(2,2);
double tmp;
if(fabs(Dxx) > fabs(Dzz)) { tmp=Dzz; Dzz=Dxx; Dxx=tmp; }
if(fabs(Dyy) > fabs(Dzz)) { tmp=Dzz; Dzz=Dyy; Dyy=tmp; }
if(fabs(Dyy) > fabs(Dxx)) { tmp=Dyy; Dyy=Dxx; Dxx=tmp; }
double Deta = (Dxx-Dyy)/Dzz;

  cout << "\n\n\tCartesian Tensor\n" << Dmx; 
  cout << "\n\n\tDiagonalized Tensor\n" << DD; 
  cout << "\n\n\tTrace Of Diagonalized Tensor: " << trace(DD);
  cout << "\n\n\tThe delxx value appears to be: " << Dxx;
  cout << "\n\n\tThe delyy value appears to be: " << Dyy;
  cout << "\n\n\tThe delzz value appears to be: " << Dzz;
  cout << "\n\n\tAsymmetry Appears To Be: " << Deta;


//                     Check The Spatial Tensor For Rank 2 Irreducibility
//                                   (Traceless and Symmetric)

  double cutoff = 1.e-15;
  bool TFT = (norm(Dmx.trace()) > cutoff);
  bool TFS = (!Dmx.is_hermitian(cutoff));
  if(TFT) cout << "\n\t\tBAD,Not A Traceless Tensor!";
  else   cout << "\n\t\tA OK, A Nice Traceless Tensor!";
  if(TFS)  cout << "\n\t\tBAD,Not A Symmetric Tensor!";
  else     cout << "\n\t\tA OK, A Nice Symmetric Tensor!";
  if(TFT || TFS)
    {
    cout << "\n\nCannot Operate On A Reducible Tensor....\n\n";
    exit(-1);
    }

//                     Set Spatial Tensor Components For Function Evaluations

  double Axx = Dmx.getRe(0,0);
  double Axy = Dmx.getRe(0,1);
  double Axz = Dmx.getRe(0,2);
  double Ayy = Dmx.getRe(1,1);
  double Ayz = Dmx.getRe(1,2);
  vector<double> A(5);
  A[0] = Axx; A[1] = Axy; A[2] = Axz; A[3] = Ayy; A[4] = Ayz;

//                Make An Initial Guess At The Spherical Components & Euler Angles

  double delzz = Dzz;				// Guess at delzz (KHz)
  double eta   = Deta;				// Guess at eta ([0,1], unitless)
  double alpha = 0.0;				// Guess at alpha (radians)
  double beta  = 0.0;				// Guess at beta  (radians)
  double gamma = 0.0;				// Guess at gamma (radians)

//  col_vector X(5);
//  X.put(delzz,0); X.put(eta,1); X.put(alpha,2); X.put(beta,3); X.put(gamma,4);
  col_vector X(3);
  X.put(alpha,0); X.put(beta,1); X.put(gamma,2);
  matrix J, Jinv;
  col_vector F, Y;
  for(int k=0; k<3; k++)
    {
    cout << "\n\n\n\n\t\t\tITERATION "        << k;
    cout << "\n\n\tStarting Values\n"         << X; 
    STJacobian3(A, delzz, eta, X, F, J);
//    STJacobian5(A, X, F, J);
    cout << "\n\n\tJacobian Matrix\n"         << J;
    cout << "\n\n\tFunction Evaluation\n"     << complex(1.e6)*F;
    Jinv = inv(J);
    Y = -1.0*(Jinv*F);
    cout << "\n\n\tJacobian Matrix Inverse\n" << Jinv;
    cout << "\n\n\tCorrection Values\n"       << Y;
    cout << "\n\n\tUpdated Values\n"          << X+Y; 
    X += Y;
    }
  }
