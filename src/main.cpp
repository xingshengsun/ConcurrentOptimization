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
#include "initialize.h"
#include "readDynaOutput.h"
#include "runDyna.h"
#include "transferVariable.h"
#include "generateMesh.h"

using namespace std;

int main(int argc, char *argv[])
{
	int case_index; // 0: joint; 1: strong; 2: soft; 3: sequential
	double m_tot; // total mass [g]
	double vi; // inpact velocity [cm/s]
	double t_max; // total simulation time [s]

	std::string _case_file = "../case.in";

  fstream case_input(_case_file.c_str(), ios::in);
  if( !case_input.is_open() ) {
		std::cerr << "ERROR: Cannot open file "
              << _case_file << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else {
    // Actually Read Values
    case_input >> case_index;
		case_input >> m_tot;
		case_input >> vi;
		case_input >> t_max;	
    case_input.close();
	}

	int NumDesVar;
	if (case_index == 0)
			NumDesVar = 5;
	else if (case_index == 1)
			NumDesVar = 2;
	else if (case_index == 2)
			NumDesVar = 2;
	else if (case_index == 3)
			NumDesVar = 1;
	else {
			std::cerr << "ERROR: wrong index number in case.in" << std::endl;
			std::exit(EXIT_FAILURE);
	}

	double t_mesh;

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
	double area = length * width; // plate upper surface area [cm^2]
	double t1, t2;
	
	if (case_index == 0) {
		t1 = des_var[4] * m_tot / rho_str / area;		
		t2 = (m_tot / area - t1*rho_str) / rho_sft;
	}
	else if (case_index == 1) {
		t1 = 0.5 * m_tot / rho_str / area;
		t2 = t1;
	}
	else if (case_index == 2) {
		t1 = 0.5 * m_tot / rho_sft / area;
		t2 = t1;
	}
	else if (case_index == 3) {
		t1 = des_var[0] * m_tot / rho_str / area;  
		t2 = (m_tot / area - t1*rho_str) / rho_sft;
	}

	generateMesh(_tmpDir, t_mesh, t1, t2);

	// run ls-dyna calcualtions
	TNT::Array2D<double> x1(2,1);
	TNT::Array2D<double> x2(2,1);
	
	vector<double> param;

	if (case_index==0) {
		transferVariable(xc1, W1, des_var[0], des_var[1], x1);
		param.push_back(x1[0][0]);
		param.push_back(x1[1][0]);
		transferVariable(xc2, W2, des_var[2], des_var[3], x2);
    param.push_back(x2[0][0]);
		param.push_back(x2[1][0]);
  }
  else if (case_index==1) {
  	transferVariable(xc1, W1, des_var[0], des_var[1], x1);
    param.push_back(x1[0][0]);
    param.push_back(x1[1][0]);
	}
  else if (case_index == 2) {
 		transferVariable(xc2, W2, des_var[0], des_var[1], x2);
    param.push_back(x2[0][0]);
    param.push_back(x2[1][0]); 
	}

	runDyna(_tmpDir, case_index, param, vi, t_max);

	// Read information of energy from glstat file in subdirectory
	double maxInternal;
	double maxHourglass;
	double maxSliding;
	readGLSTAT(_tmpDir, maxInternal, maxHourglass, maxSliding);
	//cout << maxInternal << " " <<  maxHourglass << " " << maxSliding << endl;

  // Read Velocities from nodeout file in subdirectory
	double vz_r;
  readNODOUT(_tmpDir, vz_r);

	// Remove case with wrong hourglass and sliding energy
	if ((fabs(maxHourglass) > 0.1*fabs(maxInternal)) ||
			(fabs(maxSliding) > 0.1*fabs(maxInternal))) 
		vz_r = 1.0e8;

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

  return 0;
}
