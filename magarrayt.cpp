//
// @(#)magarrayt.cpp  1.1  misha1-06oct108
//
// Copyright (c) 2008 of Mikhail Tchernychev
// All rights reserved.
//
// Created:  06oct108  by  misha
// Version:  06oct108  by  misha
//


#include "arrayinv.h"


#include <signal.h>

static int ibreak = 0;

void ctrl_c(int sig) {
  ibreak = 1;
  fprintf(stderr,"signal\n");
}

int progress_func(void * data, int percent, double discr) { // progress / cancel function

  fprintf(stderr,"iter=%d discrepancy=%lf\n", percent, discr); fflush(stderr);
  return ibreak;

}

void check_size();
void check_size2();


main(int argc, char **argv) {

  signal(SIGINT,ctrl_c);

  double positions[3], J[3], max_err, std_err;
  int i, j, jj;

  double depth   = 2.;
  double I       = 65.6;
  double D       = -12.27;
  double Az      = 90.;
  int n_channels = 2;
  double zmin    = 100;
  double zmax    = 100;
  int background = 1;
  double dm      = 1.e-3;
  int is_shift   = 0;
  double mag_ampl = 0.;
  int is_induced = 0;
  int is_fixed_z = 0;
  int max_samples = 10000;
  int n_samples   = 0;
  int n1, n2;
  int skip_tsynt = 0;

  double ** X, ** Y, ** Z, ** Alt, ** T, **Tsynt, *dTime;

  double rms_fit, max_fit, synt_max,  synt_min;
  int      max_iterations    = 2000;
  char  inv_message[1024];

	int dnls_size =  sizeof(DNLS1_STORAGE);
	int col_size   = sizeof(CMagDataCollection);

	check_size();
	check_size2();

  n1 = n2 = -1;

   X = Y = Z = Alt = T = Tsynt = NULL;
   dTime = NULL;

  if(argc == 1) {
    fprintf(stderr,"\nUsage: %s -I val -D val -Ch val -Bg val -Z val -dm val < input > output\n",argv[0]);
    fprintf(stderr,"-I    number  -  Inclination, degrees\n");
    fprintf(stderr,"-D    number  -  Declination, degrees\n");
    fprintf(stderr,"-Ch   number  - number of channels [2] \n");
    fprintf(stderr,"-Z    number  - initial depth below sea bottom, m [2 m]\n");
    fprintf(stderr,"-Bg   number  - background order (0, 1 or 2) [1]\n");
    fprintf(stderr,"-Zfix         - all data has fixed Z\n");
    fprintf(stderr,"-dm   number  - smooth depth order [1.e-3]\n");
    fprintf(stderr,"-S          - shift all profiles except first\n");
    fprintf(stderr,"-Az   number  - Local X axis azimuth [90]\n");
    fprintf(stderr,"-MaxSamples number  - Samples to reserve for the data[10000]\n");
    fprintf(stderr,"-Ind        - Induced only\n");
    fprintf(stderr,"-N1 number  - window start\n");
	fprintf(stderr,"-N2 number  - window end\n");
    fprintf(stderr,"-SkipTsynt  - work on file with X Y Z, ALT, T Tsynt\n"); 
    exit(0);
  }


  for(i=1; i<argc; i++) {
	  if(argv[i][0]=='-') {
		  if(strcmp(argv[i],"-I")==0) {
			  I = atof(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-D")==0) {
			  D = atof(argv[i+1]); i++;  continue;
		  }
		  if(strcmp(argv[i],"-Ch")==0) {
			  n_channels = atol(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-Z")==0) {
			  depth = atof(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-Bg")==0) {
			  background = atoi(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-dm")==0) {
			  dm = atof(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-S")==0) {
			  is_shift = 1; continue;
		  }
		  if(strcmp(argv[i],"-Zfix")==0) {
			  is_fixed_z = 1; continue;
		  }
		  if(strcmp(argv[i],"-Az")==0) {
			  Az = atof(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-MaxSamples")==0) {
			  max_samples = atoi(argv[i+1]); i++; continue;
		  }
		  if(strcmp(argv[i],"-N1")==0) {
			  n1 = atoi(argv[i+1])-1; i++; continue;
		  }
		  if(strcmp(argv[i],"-N2")==0) {
			  n2 = atoi(argv[i+1])-1; i++; continue;
		  }
		  if(strcmp(argv[i],"-Ind")==0) {
			  is_induced = 1; continue;
		  }
		  if(strcmp(argv[i],"-SkipTsynt")==0) {
			skip_tsynt = 1; continue;	  
		 } 
      }
  }

	// allocate memory

  X     = new double *[n_channels];
  Y     = new double *[n_channels];
  Z     = new double *[n_channels];
  T     = new double *[n_channels];
  Tsynt = new double *[n_channels];
  Alt   = new double *[n_channels];
  dTime = new double[max_samples];

  if(!X || !Y || !Z || !T || !Tsynt || !Alt) {
	  fprintf(stderr,"Cannot allocate memory for pointers. Stop\n");
	  exit(0);
  }

  for(i=0; i<n_channels; i++) {
	  X[i]       = new double[max_samples];
	  Y[i]       = new double[max_samples];
	  Z[i]       = new double[max_samples];
	  T[i]       = new double[max_samples];
	  Tsynt[i]   = new double[max_samples];
	  Alt[i]     = new double[max_samples];

	  if(!X[i] || !Y[i] || !Z[i] || !T[i]  || !Tsynt[i] || !Alt[i]) {
		  fprintf(stderr,"Cannot allocate memory for pointers. Stop\n");
		  exit(0);
	  }
  }


  // FILE * dat = fopen("field.dat", "rt");
  //#define stdin dat

  // read the data
  char buffer[1024], date_txt[20], time_txt[20], name[20];
  const char *IFS = " \t\n";

  int z_first = 1;
  double z_const;
  int is_z_const=1;

  j=0;

  while(!feof(stdin)) {
      if(!fgets(buffer, 1024, stdin)) break;

      char * p  = strtok (buffer, IFS);
      if(!p) continue;

	  if(buffer[0] == '#' || !isdigit(buffer[0])) continue;


      for(int i=0; i<n_channels; i++) {
          X[i][j]   =  atof(p); p = strtok(NULL, IFS);
          Y[i][j]   =  atof(p); p = strtok(NULL, IFS);
          Z[i][j]   =  atof(p); p = strtok(NULL, IFS);
          Alt[i][j] =  atof(p); p = strtok(NULL, IFS);
          T[i][j]   =  atof(p); p = strtok(NULL, IFS);

		  if(skip_tsynt) p = strtok(NULL, IFS); 

		  if(z_first) {
			  z_const = Z[i][j];
			  z_first = 0;
		  }
		  else {
			  if(z_const!=Z[i][j]) is_z_const = 0;
		  }

      }

      double dtime = 0.;
      strcpy(name,p);     p = strtok(NULL, IFS);
      if(p) {
            strcpy(date_txt,p); p = strtok(NULL, IFS);
            strcpy(time_txt,p);
            dTime[j] = date2sec(date_txt, time_txt, 1);
      }

      j++;
	  if(j>=max_samples) break;
   }

	//fclose(dat);

	n_samples = j;
	if(n1<0) n1 = 0;
	if(n2<0) n2 = n_samples-1;

	// invesion

	std::vector<double **> vX;
    std::vector<double **> vY;
	std::vector<double **> vZ;
	std::vector<double **> vA;
	std::vector<double **> vT;
	std::vector<double **> vTsynt;
	std::vector<double *> vTime;
	std::vector<int> vCh;
	std::vector<int> vN;
	std::vector<int> vN1;
	std::vector<int> vN2;

    std::vector<double> vX0;
    std::vector<double> vY0;
    std::vector<double> vZ0;
    std::vector<double> vJ;

	vX.push_back(X);
	vY.push_back(Y);
	vZ.push_back(Z);
	vT.push_back(T);
	vTsynt.push_back(Tsynt);
	vTime.push_back(dTime);

	vA.push_back(Alt);
	vCh.push_back(n_channels);
	vN.push_back(n_samples);
	vN1.push_back(n1);
	vN2.push_back(n2);

	vX0.push_back(NO_DATA);
	vY0.push_back(NO_DATA);
	vZ0.push_back(depth);
	vJ.push_back(NO_DATA);

	CMagCollection ObjectCollection;

	int ret = ArrayInversion<double>(ObjectCollection, vX, vY, vZ, vA, vT, vTsynt, vTime,
						  		     vCh, vN, vN1, vN2, vX0, vY0, vZ0, vJ,
							         &rms_fit, &max_fit, &synt_max, &synt_min,
									 I,D,Az, 0,
							         max_iterations,
							         progress_func,
							         inv_message);



	fprintf(stderr,"RMS fit %lf Max fit %lf\n", rms_fit, max_fit);
	fprintf(stderr,"Synt. max %lf Synt min %lf\n", synt_max, synt_min);

	fprintf(stderr,"X=%lf +/- %lf\n", vX0[0], vX0[1]);
	fprintf(stderr,"Y=%lf +/- %lf\n", vY0[0], vY0[1]);
	fprintf(stderr,"Z=%lf +/- %lf\n", vZ0[0], vZ0[1]);

	fprintf(stderr,"Jx=%lf +/- %lf\n", vJ[0], vJ[3]);
	fprintf(stderr,"Jy=%lf +/- %lf\n", vJ[1], vJ[4]);
	fprintf(stderr,"Jz=%lf +/- %lf\n", vJ[2], vJ[5]);



  // free memory

	 for(i=0; i<n_channels; i++) {
		 delete [] X[i];     X[i]   = NULL;
		 delete [] Y[i];     Y[i]   = NULL;
		 delete [] Z[i];     Z[i]   = NULL;
		 delete [] T[i];     T[i]   = NULL;
		 delete [] Alt[i];   Alt[i] = NULL;
		 delete [] Tsynt[i]; Tsynt[i] = NULL;
	 }

	 delete [] X;
	 delete [] Y;
	 delete [] Z;
	 delete [] T;
	 delete [] Tsynt;
	 delete [] Alt;
	 delete [] dTime;


}
