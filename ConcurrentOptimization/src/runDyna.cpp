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

void runDyna(std::string _tmpDir, int case_index, vector<double> param,
						 double vi, double t_max) 
{
	// run LS-DYNA calculations for joint design
	double A = 0.0;
	double sig_max = 0.0;
	double DD2 = 0.0;
	double C10 = 0.0;

	vi *= -1.0;

	if (case_index==0) {
		A = 1.0e7 * param[0];
		DD2 = param[1];
		sig_max = 1.0e7 * param[2];
    C10 = 1.0e7 * param[3];
	}
	else if (case_index==1) {
    A = 1.0e7 * param[0];
    DD2 = param[1];
  }
	else if (case_index == 2) {
    sig_max = 1.0e7 * param[0];
    C10 = 1.0e7 * param[1];
  }

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

	std::string _matfile= _tmpDir + "/MatFile.k";
  FILE * matfile;
  matfile = fopen(_matfile.c_str(),"w");
  fprintf(matfile, "*KEYWORD\n");
  fprintf(matfile, "*PARAMETER\n");
	if ((case_index==0) || (case_index==1)) {
  	fprintf(matfile, "RA,%ie%i\n", A_i2, A_i1-tenexp);
  	fprintf(matfile, "RD2,%0.6f\n", DD2);
  }
	if ((case_index==0) || (case_index==2)) {
		fprintf(matfile, "Rsig_max,%ie%i\n", sig_i2, sig_i1-tenexp);
  	fprintf(matfile, "RC10,%ie%i\n", C_i2, C_i1-tenexp);
	}
	fprintf(matfile, "Rvi,%0.0f\n", vi);
	fprintf(matfile, "Rt_max,%0.6f\n", t_max);
	fprintf(matfile, "*INCLUDE\n");
  fprintf(matfile, "../../LsdynaMain.k\n");
  fprintf(matfile, "*end\n");

  fclose(matfile);

	//  Run Simulation with system call
  //  Need to change the directory of simulation solver on different cluster
  std::string RunDyna = "cd " + _tmpDir+
  ";module load ortiz/ls-dyna/7.0.0;ls-dyna i=MatFile.k memory=600000000";
  const char *commandDyna = RunDyna.c_str();
  system(commandDyna);

	return;
}

