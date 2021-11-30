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

void initializeVariables(double &m_tot, double &length, double &width, double &rho_str, 
												 double &rho_sft, double &t_mesh, TNT::Array2D<double> &xc1,
												 TNT::Array2D<double> &W1, TNT::Array2D<double> &xc2,
                         TNT::Array2D<double> &W2) 
{
	//m_tot = 88.5*1.5; // total mass [g]
  length = 10.0; // plate length [cm]
  width = 10.0; // plate width [cm]
  rho_str = 1.77; // mass density of strong material, AZ31B [g/cm^3]
  rho_sft = 1.1; // mass density of soft material, polyurea g/cm^3]
	t_mesh = 0.03; // single mesh thickness [cm]

	// center of first ellipsoid
	xc1[0][0] = 200.0;
  xc1[1][0] = 0.6;

	// shape of first ellipsoid
	W1[0][0] = 1.421849547750393e-04;
  W1[0][1] = 0.154898844120546;
  W1[1][0] = 0.154898844120546;
  W1[1][1] = 5.687485965163476e+02;

	// center of second ellipsoid
	xc2[0][0] = 250.0;
  xc2[1][0] = 2.5;

	// shape of second ellipsoid
	W2[0][0] = 5.687485965163474e-04;
  W2[0][1] = -0.007744942206027;
  W2[1][0] = -0.007744942206027;
  W2[1][1] = 0.355462386937598;

	return;
}
