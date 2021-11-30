#include "../tnt_array1d.h"
#include "../tnt_array2d.h"
#include "../jama_eig.h"

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
	srand( (unsigned)time( NULL ) );

	{
		int i, j, k, n=3;

		TNT::Array2D<double> A(n,n);
		TNT::Array2D<double> B(n,n);

		for (i=0; i<3; i++)
			for (j=0; j<3; j++)
				B[i][j] = 2.0*((double)rand()/(double)RAND_MAX) - 1.0;

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				A[i][j] = 0.0;
				for (k=0; k<3; k++)
				{
					A[i][j] += B[i][k]*B[j][k];
				}
			}
		}

		JAMA::Eigenvalue<double> Eig(A);

//		A*V = V*D

		TNT::Array2D<double> D(n,n); Eig.getD(D);
/*
		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				cout << D[i][j] << ", ";
			}
			cout << endl;
		}
		cout << endl;
*/
		TNT::Array2D<double> V(n,n); Eig.getV(V);
		TNT::Array2D<double> AV(n,n);
		TNT::Array2D<double> VD(n,n);

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				AV[i][j] = 0.0; 
				VD[i][j] = 0.0;
				for (k=0; k<3; k++)
				{
					AV[i][j] += A[i][k]*V[k][j];
					VD[i][j] += V[i][k]*D[k][j];
				}
			}
		}

		double Error=0.0, Tolerance=1.0e-12;
		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				Error += pow(AV[i][j]-VD[i][j],2);
//				cout << AV[i][j] << ", " << VD[i][j] << endl;
			}
		}
		Error = sqrt(Error);
		try{if (Error > Tolerance) throw(1);}
		catch(int code){cout << "Jama::Eigenvalue Test " << code << endl;} 
	}

	{
		int i, j, k, n=3;

		TNT::Array1D<double> a(n);
		TNT::Array1D<double> b(n);
		TNT::Array1D<double> c(n);
		TNT::Array1D<double> d(n);
		TNT::Array2D<double> A(n,n);

		for (i=0; i<3; i++) 
			a[i] = 2.0*((double)rand()/(double)RAND_MAX) - 1.0;

		double anorm = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

		for (i=0; i<3; i++) a[i] /= anorm;

		for (i=0; i<3; i++) 
			d[i] = 2.0*((double)rand()/(double)RAND_MAX) - 1.0;

		b[0] = a[1]*d[2] - a[2]*d[1];
		b[1] = a[2]*d[0] - a[0]*d[2];
		b[2] = a[0]*d[1] - a[1]*d[0];

		double bnorm = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);

		for (i=0; i<3; i++) b[i] /= bnorm;

		c[0] = a[1]*b[2] - a[2]*b[1];
		c[1] = a[2]*b[0] - a[0]*b[2];
		c[2] = a[0]*b[1] - a[1]*b[0];

		for (i=0; i<3; i++)
		{
			A[i][0] = a[i];
			A[i][1] = b[i];
			A[i][2] = c[i];
		}

		JAMA::Eigenvalue<double> Eig(A);

//		A*V = V*D

		TNT::Array2D<double> D(n,n); Eig.getD(D);

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				cout << D[i][j] << ", ";
			}
			cout << endl;
		}
		cout << endl;

		TNT::Array2D<double> V(n,n); Eig.getV(V);
		TNT::Array2D<double> AV(n,n);
		TNT::Array2D<double> VD(n,n);

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				AV[i][j] = 0.0; 
				VD[i][j] = 0.0;
				for (k=0; k<3; k++)
				{
					AV[i][j] += A[i][k]*V[k][j];
					VD[i][j] += V[i][k]*D[k][j];
				}
			}
		}

		double Error=0.0, Tolerance=1.0e-12;
		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				Error += pow(AV[i][j]-VD[i][j],2);
				cout << AV[i][j] << ", " << VD[i][j] << endl;
			}
		}
		Error = sqrt(Error);
		try{if (Error > Tolerance) throw(2);}
		catch(int code){cout << "Jama::Eigenvalue Test " << code << endl;} 
	}

	return 0;
}
