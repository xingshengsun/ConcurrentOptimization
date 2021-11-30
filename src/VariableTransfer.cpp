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

void initializeVariables(double &m_tot, double &length, double &width, double &rho_str,
								         double &rho_sft, double &t_mesh, TNT::Array2D<double> &xc1,
												 TNT::Array2D<double> &W1, TNT::Array2D<double> &xc2,
												 TNT::Array2D<double> &W2);
using namespace std;


int main(int argc, char *argv[])
{
	int NumDesVar = 5;

	double t_mesh;
	int NumMesh1, NumMesh2;
	double t1, t2;

	double m_tot; // total mass [g]
	double length; // plate length [cm]
 	double width; // plate width [cm]
	double rho_str; // mass density of strong material, AZ31B [g/cm^3]
	double rho_sft; // mass density of soft material, polyurea g/cm^3]	

	TNT::Array2D<double> xc1(2,1);
	TNT::Array2D<double> W1(2,2);
	TNT::Array2D<double> xc2(2,1);
	TNT::Array2D<double> W2(2,2);

	initializeVariables(m_tot, length, width, rho_str, 
                      rho_sft, t_mesh, xc1,
                      W1, xc2, W2); 

	double area = length * width; // plate upper surface area [cm^2]

	TNT::Array2D<double> D1(2,2);
	TNT::Array2D<double> V1(2,2);
	JAMA::Eigenvalue<double> Eig1(W1);
	Eig1.getD(D1);
	Eig1.getV(V1);

/*
	cout << endl;
	cout << D1[0][0] << " " << D1[0][1] << endl;
	cout << D1[1][0] << " " << D1[1][1] << endl;
	cout << endl;
	cout << V1[0][0] << " " << V1[0][1] << endl;
  cout << V1[1][0] << " " << V1[1][1] << endl;
*/

	TNT::Array2D<double> D2(2,2);
  TNT::Array2D<double> V2(2,2);
  JAMA::Eigenvalue<double> Eig2(W2);
  Eig2.getD(D2);
  Eig2.getV(V2);

	/*
 	cout << endl;
  cout << D2[0][0] << " " << D2[0][1] << endl;
  cout << D2[1][0] << " " << D2[1][1] << endl;
	cout << endl;
  cout << V2[0][0] << " " << V2[0][1] << endl;
  cout << V2[1][0] << " " << V2[1][1] << endl;
	*/

 	vector<double> des_var;
 	des_var.assign(NumDesVar, 0.0);	
	
	/*
    Ensure that enough input files are specified
  */
  for(int i=1; i<=2; i++)
    {
      if( argv[i]==NULL )
        {
          std::cerr << argv[0] << ": ERROR: Must give the name of input file #" 
                    << i << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
  std::string _simulation_input  = argv[1];
  std::string _simulation_output = argv[2];
  std::string _tmpDir = _simulation_output + "Dir";
  
  /*
    Create run directory
  */
	const int dir_err = mkdir(_tmpDir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (-1 == dir_err)
  {
    std::cerr << "Error creating directory"<< std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*
    Read in simulation.in
  */
  fstream simulation_input(_simulation_input.c_str(), ios::in);
  if( !simulation_input.is_open() )
    {
      std::cerr << argv[0] << ": ERROR: Cannot open file "
                << _simulation_input << std::endl;
      std::exit(EXIT_FAILURE);
    }
  else
    {
      // Actually Read Values
			for(int i=0; i<NumDesVar; i++)
      	simulation_input >> des_var[i];

      // Sends error if any extra lines are found in the file
      //  (we'll keep it but I'm generally inclined to allow these for commenting)
      char tempchar[255];
      int inputLine=3;        
      while( true ){
          simulation_input >> tempchar;
          if(simulation_input.eof()) break;
          else
            {
              std::cerr << argv[0] << ": WARNING: Erroneous data at line " << inputLine 
                        << " of " << _simulation_input << ": " << tempchar << std::endl;
              inputLine++;
            }
        }
      simulation_input.close();
    }
	
	// Generate mesh file using Gmsh
	t1 = des_var[4];
	t2 = (m_tot / area - t1*rho_str) / rho_sft;

	//cout << t1 << " " << t2 << endl;

	NumMesh1 = round(t1/t_mesh);
	NumMesh2 = round(t2/t_mesh);

	NumMesh1 = max(NumMesh1, 1);
	NumMesh2 = max(NumMesh2, 1);


	// Creat file containing geometric parameters
	std::string _geofile= _tmpDir + "/GeoFile.geo";
	ofstream geofile(_geofile, ios::out);
	if( !geofile.is_open() ) {
		std::cerr << "geofile: ERROR: Cannot open file "
							<< _simulation_input << std::endl;
						     std::exit(EXIT_FAILURE);
	}

	geofile << "Plate_thick1 = " << t1 << " ;" << endl;
	geofile << "Plate_thick2 = " << t2 << " ;" << endl;
	geofile << "Num_mesh1 = " << NumMesh1 << " ;" << endl;
	geofile << "Num_mesh2 = " << NumMesh2 << " ;" << endl;
	geofile << "Include \"../GmshPlate.geo\";" << endl;
	geofile.close();

	// Run Gmsh with system call
	std::string RunGmsh = "cd " + _tmpDir+
		";module load gmsh/4.5.4;gmsh -3 -format key -o PlateElement.k GeoFile.geo";
	const char *commandGmsh = RunGmsh.c_str();
	system(commandGmsh);

  // For formatting exponents efficently:
  // Get exponent and coefficient and use "e" to express without wasting on decimal
  // Here we assume that exponent on A and B is less than 10 (and greater than zero)
  //
  // Notes: for JC in a standard formatting deck we get 80 colums for 8 paramters
  //  which  means we get 10 characters per paramter in fixed width.  I hate this format as I 
  //  find it roughly uneditable.  The fix is to use the comma delimited format, but
  //  the quirk of using it is that the fields still can't be wider than the original
  //  format column IF YOU INCLUDE THE COMMA IN CHARACTER COUNT.  So the formatting 
  //  for A and B are targeting a 9 character expression.

	TNT::Array2D<double> q1(2,1);
	TNT::Array2D<double> x1(2,1);
	TNT::Array2D<double> tmp1(2,2);

	q1[0][0] = des_var[0] * cos(des_var[1]);
	q1[1][0] = des_var[0] * sin(des_var[1]);	

	tmp1 = matpow(D1,-0.5);
	tmp1 = matmult(V1, tmp1);
	x1 = xc1 + matmult(tmp1, q1); 

	double A = x1[0][0];
	A *= 1.0e7;
	double DD2 = x1[1][0];

	/*
	cout << endl;
	cout << q1[0][0] << " " << q1[1][0] << endl;
	cout << endl;
	cout << A << " " << DD2 << endl;
	*/

	TNT::Array2D<double> q2(2,1);
  TNT::Array2D<double> x2(2,1);
  TNT::Array2D<double> tmp2(2,2);

  q2[0][0] = des_var[2] * cos(des_var[3]);
  q2[1][0] = des_var[2] * sin(des_var[3]);

  tmp2 = matpow(D2,-0.5);
	tmp2 = matmult(V2, tmp2);

  x2 = xc2 + matmult(tmp2, q2);

	double sig_max = x2[0][0];
  double C10 = x2[1][0];
	sig_max *= 1.0e7;
	C10 *= 1.0e7;

	//cout << C10 << endl;

	/*
  cout << endl;
  cout << q2[0][0] << " " << q2[1][0] << endl;
  cout << endl;
  cout << sig_max << " " << C10 << endl;
	*/

  int tenexp=6;
  double scale=pow(10,tenexp);

  int A_i1=floor(log10(A));
  double A_d=round(A/pow(10,(double)A_i1)*scale);
  int A_i2=(int)A_d;

	int sig_i1=floor(log10(sig_max));
  double sig_d=round(sig_max/pow(10,(double)sig_i1)*scale);
  int sig_i2=(int)sig_d;

	int C_i1=floor(log10(C10));
	double C_d=round(C10/pow(10,(double)C_i1)*scale);
	int C_i2=(int)C_d;

/*
	cout << endl;
	printf("RA,%ie%i\n", A_i2, A_i1-tenexp);
	cout << endl;
  printf("Rsig_max,%ie%i\n", sig_i2, sig_i1-tenexp);
*/

	std::string _matfile= _simulation_output.c_str();
	FILE * matfile;
	matfile = fopen(_matfile.c_str(),"w");
	fprintf(matfile, "*KEYWORD\n");
	fprintf(matfile, "*PARAMETER\n");
	fprintf(matfile, "RA,%ie%i\n", A_i2, A_i1-tenexp);
	fprintf(matfile, "RD2,%0.6f\n", DD2);
	fprintf(matfile, "Rsig_max,%ie%i\n", sig_i2, sig_i1-tenexp);
	fprintf(matfile, "RC10,%ie%i\n", C_i2, C_i1-tenexp);
	fprintf(matfile, "*INCLUDE\n");
	fprintf(matfile, "../../LsdynaMain.k\n");
	fprintf(matfile, "*end\n");

  fclose(matfile);
 /* 
//  Run Simulation with system call
	//  Need to change the directory of simulation solver on different cluster
  std::string RunDyna = "cd " + _tmpDir+
  ";module load ortiz/ls-dyna/7.0.0;ls-dyna i=MatFile.k memory=600000000";
  const char *commandDyna = RunDyna.c_str();
  system(commandDyna);

  // Read Velocities from nodeout file in subdirectory
  std::string _dyna_output = _tmpDir+"/nodout";
  double t,z,vz_r;
	int nCycle, maxCycle;
	maxCycle = 19900;
  // nodeout is an output file for displacements from every node from a specified node set
  // it records values in time associated text blocks like:
  // nodal print out for time step 1                              ( at time 0.0000000E+00 )

	// nodal point x-disp y-disp z-disp x-vel y-vel z-vel x-accl y-accl z-accl x-coor y-coor z-coor
	// 5 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 -1.0000E+05 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 2.3800E-01	 
  
	// Node 5 is the center of the ball, and the performance measure is the z-velocity at this point

	std::string lineDyna;

  // Open nodeout file
  fstream dyna_output(_dyna_output.c_str(), ios::in);
  if( !dyna_output.is_open() )
    {
      std::cerr << "nodeout: ERROR: Cannot open file "
                << _dyna_output << std::endl;
      std::exit(EXIT_FAILURE);
    }
  else
    {
			// XS: I have set the max number of cycles (30000) in the LSDYNA input file, in order to 
			// shorten LSDYNA calculations. Thus, the result with the max cycle number should be highlighted
			// in GA population.
			dyna_output.seekg(-370, dyna_output.end); // find the location of max cycle number
			getline(dyna_output,lineDyna);
			stringstream ssin1(lineDyna);
			ssin1 >> nCycle;

			if (nCycle < maxCycle) {
				dyna_output.seekg(-154, dyna_output.end); // go to the last line of nodout file
				getline(dyna_output,lineDyna); // read last line
				stringstream ssin(lineDyna);
				int i = 0;
				while (ssin.good() && i < 7){ // z-velicity is on the 7th column
					ssin >> z;
					++i;
				} 

				std::cout << "Residual velocity: " << z/100.0 << " m/s." << std::endl;
				if (z < 0.0) {
					std::cout << "Penetration happens." << std::endl; 
				}
				else {
					std::cout << "Penetration does not happen." << std::endl;
				}
				vz_r = -z; // Note: penetration happens if z<0
      	dyna_output.close();
			}
			else {
				std::cout << "Warning the cycle number is: " << nCycle << std::endl;
				vz_r = 1.0e8;
			}
    }
	
  //
  // Create designated simulation output file
  FILE * simulation_output;
  simulation_output=fopen(_simulation_output.c_str(),"w");
  if( simulation_output==NULL )
    {
      std::cerr << argv[0] << ": ERROR: Cannot open file "
                << _simulation_output << std::endl;
      std::exit(EXIT_FAILURE);
    }
  else
    {
     std::cout <<  "Write objective function " << vz_r;
     std::cout <<  " cm/s to "<< _simulation_output;
     std::cout <<  std::endl;
      fprintf(simulation_output,"%.20e", vz_r);
     }

  fclose(simulation_output);

  // Attempt to open KeepDir
  fstream clear_dir("../KeepDir", ios::in);
  if( !clear_dir.is_open() )
    {
      // No protector file, so clean up!
      std::cout<<"Working files removed" << endl;
      std::string ClearCmD = "rm -r " + _tmpDir;
      const char *command2 = ClearCmD.c_str();
      system(command2);
    }
  else
    {
     // File found, provide instructions if removal of directory desired
     std::cout<< "Files retained for posterity" << std::endl;
     std::cout<< "Remove file ../KeepDir for self cleaning" << std::endl;
    }*/
  return 0;
}
