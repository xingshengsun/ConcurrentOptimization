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

using namespace std;

void readGLSTAT(std::string _tmpDir, double &maxInternal, double &maxHourglass, double &maxSliding)
{
  // Read information of energy from glstat file in subdirectory
  maxInternal = 0.0;
	maxHourglass = 0.0;
	maxSliding = 0.0;
	double tmpInternal = 0.0;
	double tmpHourglass = 0.0;
	double tmpSliding = 0.0;

  std::string lineDyna;

	std::string _energy_output = _tmpDir+"/glstat";
  fstream energy_output(_energy_output.c_str(), ios::in);

  if( !energy_output.is_open() ) {
      std::cerr << "glstat: ERROR: Cannot open file "
                << _energy_output << std::endl;
      std::exit(EXIT_FAILURE);
  }
  else {
    while( getline(energy_output,lineDyna) ){
      if(energy_output.eof()) {
        break;
      }
      else if (strncmp(lineDyna.c_str(), " internal energy", 9) == 0) {
				//cout << lineDyna << endl;
				lineDyna.erase(1, 34);
      	tmpInternal = atof(lineDyna.c_str());
				//cout << tmpInternal << endl;
				maxInternal = max(maxInternal, fabs(tmpInternal));
			}
			else if (strncmp(lineDyna.c_str(), " hourglass energy", 9) == 0) {
        //cout << lineDyna << endl;
        lineDyna.erase(1, 34);
        tmpHourglass = atof(lineDyna.c_str());
        //cout << tmpHourglass << endl;
      	maxHourglass = max(maxHourglass, fabs(tmpHourglass));
			}
			else if (strncmp(lineDyna.c_str(), " sliding interface energy", 9) == 0) {
        //cout << lineDyna << endl;
        lineDyna.erase(1, 34);
        tmpSliding = atof(lineDyna.c_str());
        //cout << tmpSliding << endl;
				maxSliding = max(maxSliding, fabs(tmpSliding));
      }
    }
    energy_output.close();
  }
	
	return;
}


void readNODOUT(std::string _tmpDir, double &vz_r)
{
	// Read residual velocity from nodout file in subdirectory
	vz_r = 0.0;
	double t,z;
  int nCycle, maxCycle;
  maxCycle = 199000;
	// nodeout is an output file for displacements from every node from a specified node set
  // it records values in time associated text blocks like:
  // nodal print out for time step 1                              ( at time 0.0000000E+00 )

  // nodal point x-disp y-disp z-disp x-vel y-vel z-vel x-accl y-accl z-accl x-coor y-coor z-coor
  // 5 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 -1.0000E+05 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 0.0000E+00 2.3800E-01  

  // Node 5 is the center of the ball, and the performance measure is the z-velocity at this point

  std::string lineDyna;

  // Open nodeout file
	std::string _dyna_output = _tmpDir+"/nodout";
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

	return;
}
