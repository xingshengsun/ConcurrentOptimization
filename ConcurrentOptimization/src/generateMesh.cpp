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

void generateMesh(std::string _tmpDir, double t_mesh, double t1, double t2)
{
	// generate mesh file of plate using GMSH
	int NumMesh1 = (int) round(t1/t_mesh);
	int NumMesh2 = (int) round(t2/t_mesh);
	
	NumMesh1 = max(NumMesh1, 1);
	NumMesh2 = max(NumMesh2, 1);

	std::string _geofile= _tmpDir + "/GeoFile.geo";
	ofstream geofile(_geofile, ios::out);
  if( !geofile.is_open() ) {
    std::cerr << "geofile: ERROR: Cannot open file "
              << _tmpDir << std::endl;
                 std::exit(EXIT_FAILURE);
  }

	geofile << "Plate_thick1 = " << t1 << " ;" << endl;
  geofile << "Plate_thick2 = " << t2 << " ;" << endl;
  geofile << "Num_mesh1 = " << NumMesh1 << " ;" << endl;
  geofile << "Num_mesh2 = " << NumMesh2 << " ;" << endl;
  geofile << "Include \"../../GmshPlate.geo\";" << endl;
  geofile.close();

	// Run Gmsh with system call
  std::string RunGmsh = "cd " + _tmpDir+
		";module load gmsh/4.5.4;gmsh -3 -format key -o PlateElement.k GeoFile.geo";
	const char *commandGmsh = RunGmsh.c_str();
  system(commandGmsh);

	return;
}

