#include <gamma.h>
#include <vector>

void STJacobian(const col_vector& X, col_vector& F, matrix& J)

	// Input		X	: Vector - {alpha,beta,gamma}
	//			A	: Vector - {Axx,Axy,Axz}
	//			F       : Vector - {f1,f2,f3}
	//			J	: Jacobian Matrix (3x3)
	// Output		void	: Given the values in X & A the
	//				  vector F and Jacobian matrix J
	//				  are generated

  {
//			        Explicitly Get Our Input Variables 

  double x1 = X.getRe(0);
  double x2 = X.getRe(1);
  double x3 = X.getRe(2);

  double Jij, Fi;
  J = matrix(3,3);

  Jij = 3.0;                 J.put(Jij,0,0);
  Jij = x3*sin(x2*x3);       J.put(Jij,0,1); 
  Jij = x2*sin(x2*x3);       J.put(Jij,0,2);

  Jij = 2.0*x1;              J.put(Jij,1,0);
  Jij = -162.0*(x2 + 0.1);   J.put(Jij,1,1); 
  Jij = cos(x3);             J.put(Jij,1,2);

  Jij = -x2*exp(-1.0*x1*x2); J.put(Jij,2,0);
  Jij = -x1*exp(-1.0*x1*x2); J.put(Jij,2,1); 
  Jij = 20.0;  		     J.put(Jij,2,2);

//					Build The F Vector

  F = col_vector(3);
  Fi = 3.0*x1 - cos(x2*x3) - 0.5;                       F.put(Fi,0);
  Fi = x1*x1 - 81.0*(x2+0.1)*(x2+0.1) + sin(x3) + 1.06; F.put(Fi,1);
  Fi = exp(-1.0*x1*x2) + 20.0*x3 + (10.0*PI - 3.0)/3.0; F.put(Fi,2);

  return;
  }


int main()
  {
//                Make An Initial Guess At The Spherical Components & Euler Angles

  col_vector X(3);
  X.put(0.1,0); X.put(0.1,1); X.put(-0.1,2);

  matrix J, Jinv;
  col_vector F, Y;
  for(int k=0; k<6; k++)
    {
    cout << "\n\n\n\n\t\t\tITERATION "        << k;
    cout << "\n\n\tStarting Values\n"         << X; 
    STJacobian(X, F, J);
    cout << "\n\n\tJacobian Matrix\n"         << J;
    cout << "\n\n\tFunction Evaluation\n"     << F;
    Jinv = inv(J);
    Y = -complex1*(Jinv*F);
    cout << "\n\n\tJacobian Matrix Inverse\n" << Jinv;
    cout << "\n\n\tCorrection Values\n"       << Y;
    cout << "\n\n\tUpdated Values\n"          << X+Y; 
    X += Y;
    }
  }
