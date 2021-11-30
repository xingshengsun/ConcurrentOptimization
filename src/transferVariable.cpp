#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <vector>
#include "./JAMA/tnt_array2d.h"
#include "./JAMA/jama_eig.h"
#include "./JAMA/tnt_array2d_utils.h"

using namespace std;

void transferVariable(TNT::Array2D<double> &xc, TNT::Array2D<double> &W, 
											double r, double theta, TNT::Array2D<double> &x)
{
	// transfer variables from spherical coordinate system to Cartesian coordinate system
	TNT::Array2D<double> D(2,2);
  TNT::Array2D<double> V(2,2);
  JAMA::Eigenvalue<double> Eig(W);
  Eig.getD(D);
  Eig.getV(V);

	cout << D[0][0] << " " << D[0][1] << endl;
	cout << D[1][0] << " " << D[1][1] << endl;
	cout << endl;
	cout << V[0][0] << " " << V[0][1] << endl;
	cout << V[1][0] << " " << V[1][1] << endl;

	TNT::Array2D<double> q(2,1);
  TNT::Array2D<double> tmp(2,2);

	q[0][0] = r * cos(theta);
  q[1][0] = r * sin(theta);

	tmp = matpow(D,-0.5);
  tmp = matmult(V, tmp);
  x = xc + matmult(tmp, q);

	return;
}

