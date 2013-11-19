#include <iostream>
#include <stdlib.h>
#include <signal.h>


#include "magobject.h"
#include "magdipole.h"
#include "magcollection.h"
#include "maggriddata.h"
#include "magprofiles.h"
#include "maginv.h"
#include "magfunc.h"
#include "l2estimator.h"
#include "l1estimator.h"
#include "dnls1estimator.h"
#include "tn_estimator.h"
#include "magbackground.h"
#include "gradprofile.h"

using namespace std;

static int ibreak = 0;

void ctrl_c(int sig) {
  ibreak = 1;
  fprintf(stderr,"signal\n");
}

int progress_func(void * data, int percent, double discr) { // progress / cancel function

  fprintf(stderr,"iter=%d discrepancy=%lf\n", percent, discr); fflush(stderr);
  return ibreak;

}

int main(int argc, char *argv[])
{
 double I,D, T;    // magnetic field parameters
 int n_dipoles;    // number of dipoles
 double params[3];
 double x,y,z, Jt; // dipole parameters
 int n_lines;
 double x0,y0,z0, nx, ny, dl, length; // line parameters
 int n_data;
 char buffer[1024];
 char line_name[100];
 int is_induced = 1;
 int i,j;
 double max_fit, step, max_step;
 double dir_x, dir_y, dir_z;
 
 CMagCollection dipole_collection;
 CMagDataCollection lines_collection;


// ------------ get parameters ----------------------
 
 FILE * dat = fopen(argv[1],"rt");
 if(!dat) {
      fprintf(stderr,"Can not open parameter file\n");
      exit(0);
  }
  
  fgets(buffer, 1024, dat);   // total field parameters
  sscanf(buffer,"%lf%lf%lf", &I,&D,&T);
         
  fgets(buffer, 1024, dat);   // total dipoles
  sscanf(buffer,"%d", &n_dipoles);
  
  for(i=0; i<n_dipoles; i++) { // dipoles
    fgets(buffer, 1024, dat);   
    sscanf(buffer,"%lf%lf%lf%lf", &params[0],&params[1],&params[2], &Jt);     
    
    CMagDipole *pDipole = new CMagDipole(is_induced); 
    pDipole->SetNonLinearParams(params, 3);
    pDipole->SetLinearParams(&Jt, 1);
    dipole_collection.AddMagObject(pDipole);  
  }    
  	
  // total mag lines lines
  
  fgets(buffer, 1024, dat);  
  sscanf(buffer,"%d", &n_lines);   	
 
  // mag lines lines itself
  
  for(i=0; i<n_lines; i++) {
      fgets(buffer, 1024, dat);  
      sscanf(buffer,"%lf%lf%lf%lf%lf%lf%lf", 
            &x0, &y0, &z0, &dl, &nx, &ny, &length);  
      
      sprintf(line_name,"%d", i+1);
      CTestProfile * prof = new CTestProfile; 
      prof->InitProfile(line_name);
      
      int n_points = int(length/dl)+1;
      for(int j=0; j<n_points; j++)
          prof->Add(x0+nx*j*dl, y0+ny*j*dl, 0., 0., z0);
      
      lines_collection.AddData(prof);
  }    
      
      
   // total gradient lines lines
  
  fgets(buffer, 1024, dat);  
  sscanf(buffer,"%d", &n_lines);   	
 
  // grad lines lines itself
  FIELDTYPE ftype = totalmaggradient; 
  
  for(i=0; i<n_lines; i++) {
      fgets(buffer, 1024, dat);  
      sscanf(buffer,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &x0, &y0, &z0, &dl, &nx, &ny, &length, &dir_x, &dir_y, &dir_z);  
      
      sprintf(line_name,"grad%d", i+1);
      CGradProfile * prof = new CGradProfile; 

      prof->InitProfile(line_name);
      prof->SetFieldType(ftype);
      prof->SetVariableDirection();
           
      int n_points = int(length/dl)+1;
      for(int j=0; j<n_points; j++)
          prof->Add(x0+nx*j*dl, y0+ny*j*dl, 0., 0., z0, dir_x, dir_y, dir_z);
      
      lines_collection.AddData(prof);
  }    
       
      
   // max_fit, step, max_step
  
  fgets(buffer, 1024, dat);  
  sscanf(buffer,"%lf%lf%lf", &max_fit, &step, &max_step);
 
 fclose(dat);
  // ---------------------- end parameters
 
  CMagInv inversion;
  inversion.SetProgressFunc(progress_func);
  inversion.SetMagData(&lines_collection);
  inversion.SetMagCollection(&dipole_collection);
  inversion.SetEarthFieldDirections(I,D,90.);
  
  // compute field
  
  int ret = inversion.ComputeField(0,1);     
     
     
  // estimate range
  double fit;
  Z29RLS_PARAMS rls_params;
  rls_params.dm = 1.e+15;

  inversion.EstimateLinearParameters(Estimate_L2_rls, &fit, 1, &rls_params);
        
   DNLS1_STORAGE dstorage;
   dstorage.iopt =  2;
  
   char message[1024];
   inversion.EstimateAllParameters(Estimator_dnls1, &dstorage, message);
   inversion.EstimateNonLinearParametersRange2(Estimator_dnls1, max_fit, step, 
                                               max_step, &dstorage);  
     
   lines_collection.Dump(stdout);
   inversion.DumpObjects(stderr);  
        	
  return 0;
}
