//
// @(#)magarray.cpp  1.1  misha-11nov104
//
// Copyright (c) 2004 of Geometrics.
// All rights reserved.
//
// Created:  11nov104  by  misha@MISHA
// Version:  11nov104  by  misha@MISHA
//

#include "magarray.h"
#include "magfunc.h"

void CMagProfileArray::clear_all() {

    m_vX.clear();
    m_vY.clear();

    m_vXmiddle.clear();
    m_vYmiddle.clear();
    
    m_vZ.clear();
    m_vTime.clear();
    m_vField.clear();
    m_vAlt.clear();
    m_vComputedField.clear();
    m_nPos = 0;

	m_vPositionNames.clear();

    m_lfX1 = m_lfX2 = m_lfY1 = m_lfY2 = m_lfZ1 = m_lfZ2 = 0.;

    m_lfXaverage = 0.;
    m_lfYaverage = 0.;
    m_lfZaverage = 0.;

    m_nChannels = 0;
    m_iStart    = 0;
    m_iStop     = 0;
  
    m_iTotalStart = 0;  
    m_iTotalStop  = 0;

    m_lfMedianDepth = 0.;  
    m_llfMedianAlt  = 0.;

    m_lfNX = 1.;  // unit vector along line section 
    m_lfNY = 0.;

   m_lfMiddleX = 0.;
   m_lfMiddleY = 0.;
   m_lfLength  = 0.;

  }

void CMagProfileArray::InitProfile(char * name) {
    clear_all();
    SetObjectName(name);
}


void CMagProfileArray::Add(double * x, double * y, double time, double * field, double * z,
			     double * alt, double x_middle, double y_middle, char * name) {

  double z_in, alt_in;

  for(int i=0; i< m_nChannels; i++) {

  z_in = alt_in = 0.;
    
  if(z)   z_in   = z[i];  
  if(alt) alt_in = alt[i];  

  if(m_vField.size()) {
    m_lfX1 = MIN(m_lfX1, x[i]);
    m_lfX2 = MAX(m_lfX2, x[i]);
    m_lfY1 = MIN(m_lfY1, y[i]);
    m_lfY2 = MAX(m_lfY2, y[i]);
    m_lfZ1 = MIN(m_lfZ1, z_in);
    m_lfZ2 = MAX(m_lfZ2, z_in);
  }
  else {
    m_lfX1 = x[i];
    m_lfX2 = x[i];
    m_lfY1 = y[i];
    m_lfY2 = y[i];
    if(z) {
      m_lfZ1 = z[i];
      m_lfZ2 = z[i];
    }
  }

    m_vX.push_back(x[i]);
    m_vY.push_back(y[i]);
    m_vZ.push_back(z_in);
    m_vAlt.push_back(alt_in);

    m_vField.push_back(field[i]);
    m_vComputedField.push_back(0);
   
    m_lfXaverage = (m_lfX1 + m_lfX2) / 2.;
    m_lfYaverage = (m_lfY1 + m_lfY2) / 2.;
    m_lfZaverage = (m_lfZ1 + m_lfZ2) / 2.;
  }

  m_vTime.push_back(time-GetTime0());
  m_iStop++;
  m_iTotalStop  = m_nChannels*m_iStop;
  
  m_vXmiddle.push_back(x_middle);
  m_vYmiddle.push_back(y_middle);

  if(name) {
		std::string sName = name;
		m_vPositionNames.push_back(sName);
	}

 }



void CMagProfileArray::Correct() {
    int n =  m_vField.size();
    for(int i=0; i<n; i++) {
      m_vField[i] = m_vField[i] + m_vTime[i]* GetDrift() + GetDCShift(); 
    }
 }


int CMagProfileArray::SetStartStop(int start, int stop) {

  int i;
  if(start < 0) stop = 0;
  if(stop > int(m_vTime.size())) stop = m_vTime.size(); 

  if(start > stop || start >= int(m_vTime.size())) return -1;

  m_iStart = start;
  m_iStop  = stop;

  m_iTotalStart = m_nChannels*m_iStart;
  m_iTotalStop  = m_nChannels*m_iStop;

  m_nPos = 0;

  // compute median values for depth / alt

  int n = m_iTotalStop-m_iTotalStart;
  fprintf(stderr,"Total data points %d\n", n);
  double * d = new double[n];

  if(d) {
      for(i=0; i<n; i++)  d[i] = m_vZ[m_iTotalStart+i]; //GetZ(i);
      CalculateMode(d, n, &m_lfMedianDepth);
      for(i=0; i<n; i++)  d[i] = GetAlt(i);
      CalculateMode(d, n, &m_llfMedianAlt);
      delete [] d;
  }

  // calculate vector along section

  int n_pos = m_iStart;
    
  double x1, y1, x2, y2;
    
  x1 = y1 = x2 = y2 = 0.;
    
  for(i=0; i< m_nChannels; i++) {
        x1   += m_vX[m_nChannels*n_pos+i];    
        y1   += m_vY[m_nChannels*n_pos+i];  
    }
                   
  x1 = x1/m_nChannels; 
  y1 = y1/m_nChannels; 

 n_pos = m_iStop-1;    

 for(i=0; i< m_nChannels; i++) {
        x2   += m_vX[m_nChannels*n_pos+i];    
        y2   += m_vY[m_nChannels*n_pos+i];  
    }

  x2 = x2/m_nChannels; 
  y2 = y2/m_nChannels;

  double s = m_lfLength = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y1));
 
  m_lfNX = (x2-x1)/s;
  m_lfNY = (y2-y1)/s;

  // find middle point

  m_lfMiddleX = (x1+x2)/2.;
  m_lfMiddleY = (y1+y2)/2.;

  return 1;
}


double CMagProfileArray::GetDistanceFromMiddle(double x, double y) {

  return m_lfNX*(x-m_lfMiddleX) + m_lfNY*(y-m_lfMiddleY);

}

double CMagProfileArray::GetNearestTime(double x, double y) {

    double dist, dtime, dist0, dtime0;
    int k = 0;
    for(int i=m_iStart; i<m_iStop; i++) {
	dtime = m_vTime[i];
	for(int j=0; j<m_nChannels; j++) {
	    dist = sqrt(pow(m_vX[m_nChannels*i+j]-x,2) +
			pow(m_vY[m_nChannels*i+j]-y,2));
	    if(!k) { dist0 = dist; dtime0 = dtime;}
	    else {
		if(dist < dist0) {
		    dist0  = dist;
		    dtime0 = dtime;
		}
	    }
	    k++;
	}
    }

    return dtime0;
}


int CMagProfileArray::GetMiddlePointCoords(double *x, double *y, double *z, double *alt) {

    int n_pos = (m_iStart + m_iStop) /2;    
     
    double x0, y0, z0, alt0;
    
    x0 = y0 = z0 = alt0 = 0.;
    
    for(int i=0; i< m_nChannels; i++) {
        x0   += m_vX[m_nChannels*n_pos+i];    
        y0   += m_vY[m_nChannels*n_pos+i];  
        z0   += m_vZ[m_nChannels*n_pos+i];  
        alt0 += m_vAlt[m_nChannels*n_pos+i];          
    }
                   
     *x   = x0/m_nChannels;
     *y   = y0/m_nChannels;   
     *z   = z0/m_nChannels;   
     *alt = alt0/m_nChannels;
     
     return 1;            
}     


int CMagProfileArray::get_data_origin(double *x, double *y, double *z) { 

    *x = m_lfXaverage;
    *y = m_lfYaverage;
    *z = m_lfZaverage;

    return 1;

}


double CMagProfileArray::GetX(int n_pos) {

  n_pos += m_iTotalStart;
  if(n_pos < m_iTotalStart || n_pos >= m_iTotalStop) return 0.;
  return m_vX[n_pos]; 
}

double CMagProfileArray::GetY(int n_pos) {

  n_pos += m_iTotalStart;
  if(n_pos < m_iTotalStart || n_pos >= m_iTotalStop) return 0.;
  return m_vY[n_pos]; 
}

double CMagProfileArray::GetZ(int n_pos) {

  //return m_lfMedianDepth;
  n_pos += m_iTotalStart;
  if(n_pos < m_iTotalStart || n_pos >= m_iTotalStop) return 0.;
  return m_vZ[n_pos]; 
}

double CMagProfileArray::GetField(int n_pos) {

  n_pos += m_iTotalStart;
  if(n_pos < m_iTotalStart || n_pos >= m_iTotalStop) return 0.;
  return m_vField[n_pos]; 
}


double CMagProfileArray::GetAlt(int n_pos) {

  n_pos += m_iTotalStart;
  if(n_pos < m_iTotalStart || n_pos >= m_iTotalStop) return 0.;
  return m_vAlt[n_pos]; 
}


double CMagProfileArray::GetTime(int n_pos) {

  n_pos = int(n_pos/m_nChannels) + m_iStart; 
  if(n_pos < m_iStart || n_pos >= m_iStop) return 0.;

  return m_vTime[n_pos]; 
}

double CMagProfileArray::GetSyntField(int n_pos) {

  n_pos += m_iTotalStart;
  if(n_pos < m_iTotalStart || n_pos >= m_iTotalStop) return 0.;
  return m_vComputedField[n_pos];  
    
}

    
     

int CMagProfileArray::get_next_position(double *x, double *y, double *z, double * T, int *n_pos,
				     double *add_parameters, int size_in, int *size_out){ 

  if(m_nPos >=  m_iTotalStop) return 0;

  *x = GetX( m_nPos);
  *y = GetY( m_nPos);
  *z = GetZ( m_nPos);
  *T = GetField(m_nPos);
  *n_pos = m_nPos;
  GetAllPositionParams(add_parameters, size_in);
  
  m_nPos++;
 
  return 1;

}

int CMagProfileArray::get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
				  double *add_parameters, int size_in, int *size_out) {

  if(n_pos < 0 || n_pos >= m_iTotalStop-m_iTotalStart) return 0;

  *x = GetX(n_pos);
  *y = GetY(n_pos);
  *z = GetZ(n_pos);
  *T = GetField(n_pos);

  GetAllPositionParams(add_parameters, size_in);

  return 1;

}

int CMagProfileArray::set_synthetic_data(double *data, int replace) {

  int j=0;
  for(int i=m_iTotalStart; i<m_iTotalStop; i++, j++)  {
    m_vComputedField[i] = data[j];
    if(replace) {
        m_vField[i] =  data[j];
    }    
  }

  return 1;
} 

int CMagProfileArray::set_synthetic_data_channel(double *data, int channel,
                                                 int replace) {
  int j=0;
  for(int i=m_iStart; i<m_iStop; i++, j++)  {
    m_vComputedField[i*m_nChannels+channel] = data[j];
    if(replace) {
        m_vField[i*m_nChannels+channel] =  data[j];
    }    
  }

  return 1;
} 

int CMagProfileArray::get_position(int n, double *xd, double *yd, double *zd) { 

  double T;
  return get_at_position(n, xd, yd, zd, &T); 

}


int CMagProfileArray::GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			       double *z1, double *z2) { 

  *x1 = m_lfX1;
  *x2 = m_lfX2;
  *y1 = m_lfY1;
  *y2 = m_lfY2;
  *z1 = m_lfZ1;
  *z2 = m_lfZ2;
 
  return 1;

} 


void CMagProfileArray::Dump(FILE * dat) {

  int i, j, k=0;
  char date_txt[20], time_txt[20];

  int n_synt = GetSyntheticPoints();
 
  for(i=0; i<m_nChannels; i++) {
     j=i+1;  
     if(n_synt)
		 fprintf(dat,"#X%d   Y%d  Z%d  ALT%d  T%d Tsynt%d ", j,j,j,j,j,j);
     else
		 fprintf(dat,"#X%d   Y%d  Z%d  ALT%d  T%d ", j,j,j,j,j);         
  }

  if(GetTime(m_iStart)> 0.)
     fprintf(dat,"LINE DATE TIME\n");           
  else
	fprintf(dat,"LINE\n"); 


  for(i=m_iStart; i<m_iStop; i++) {

     date_txt[0] = 0;
	 time_txt[0] = 0;

     sec2date(GetTime(k), date_txt, time_txt, 1);

      for(j=0; j<m_nChannels; j++) {
          
       fprintf(dat,"%.3lf  %.3lf  %.3lf  %.3lf  %.3lf ",  GetX(k) + m_lfPos[X_OFFSET], 
	      GetY(k) + m_lfPos[Y_OFFSET],
	      GetZ(k) + m_lfPos[Z_OFFSET],
	      GetAlt(k) - m_lfPos[Z_OFFSET], GetField(k)); 
          if(n_synt) fprintf(dat," %.3lf ", GetSyntField(k));
          k++;   
	  }


	if(m_vPositionNames.size()) {
		std::string sName = m_vPositionNames.at(i);
		fprintf(dat,"%s %s %s\n", sName.c_str(), date_txt, time_txt);
	}
	else 
      fprintf(dat,"%s %s %s\n", GetName(), date_txt, time_txt);
           
   }

}


void CMagProfileArray::Dump(FILE * dat, char * head1, char *head2) {

  int i, j, k=0;
  char date_txt[20], time_txt[20];
  if(m_vPositionNames.size()==0) return;

  std::string sOldName ="";
  std::vector<int> v_data_points;
  v_data_points.clear();


  int nn = 0;

  for(i=m_iStart; i<m_iStop; i++) {
	std::string sName = m_vPositionNames.at(i);
	if(strcmp(sName.c_str(), sOldName.c_str())==0 || i==m_iStart) {
		nn++;
	}
	else {
		v_data_points.push_back(nn);
		nn=1;
	}
	sOldName = sName;
  }

  v_data_points.push_back(nn);

  fprintf(dat,"%s\n", head1);
  fprintf(dat,"#P %d %s %s\n", v_data_points.at(0), head2, m_vPositionNames.at(0).c_str());
  sOldName = m_vPositionNames.at(0);
  nn = 1;

  int n_synt = GetSyntheticPoints();
 
  for(i=0; i<m_nChannels; i++) {
     j=i+1;  
     if(n_synt)
		 fprintf(dat,"#X%d   Y%d  Z%d  ALT%d  T%d Tsynt%d ", j,j,j,j,j,j);
     else
		 fprintf(dat,"#X%d   Y%d  Z%d  ALT%d  T%d ", j,j,j,j,j);         
  }
 

 if(GetTime(m_iStart)> 0.)
     fprintf(dat,"LINE DATE TIME\n");           
  else
	fprintf(dat,"LINE\n"); 


         

  for(i=m_iStart; i<m_iStop; i++) {

	date_txt[0] = 0;
	time_txt[0] = 0;

     sec2date(GetTime(k), date_txt, time_txt, 1);

	 std::string sName = m_vPositionNames.at(i);

	 if(strcmp(sOldName.c_str(), sName.c_str()) && i!=m_iStart) {
		 fprintf(dat,"\n\n%s\n", head1);
		 fprintf(dat,"#P %d %s %s\n", v_data_points.at(nn), head2, sName.c_str());
		 for(int ik=0; ik<m_nChannels; ik++) {
			 j=ik+1;  
			 if(n_synt)
				 fprintf(dat,"X%d   Y%d  Z%d  ALT%d  T%d Tsynt%d ", j,j,j,j,j,j);
			 else
				 fprintf(dat,"X%d   Y%d  Z%d  ALT%d  T%d ", j,j,j,j,j);         
		 }
		 fprintf(dat,"LINE DATE TIME\n");
		 nn++;
	 }

      for(j=0; j<m_nChannels; j++) {
          
       fprintf(dat,"%.3lf  %.3lf  %.3lf  %.3lf  %.3lf ",  GetX(k) + m_lfPos[X_OFFSET], 
	      GetY(k) + m_lfPos[Y_OFFSET],
	      GetZ(k) + m_lfPos[Z_OFFSET],
	      GetAlt(k) - m_lfPos[Z_OFFSET], GetField(k)); 
          if(n_synt) fprintf(dat," %.3lf ", GetSyntField(k));
          k++;   
	  }


		fprintf(dat,"%s %s %s\n", sName.c_str(), date_txt, time_txt);
		sOldName = sName;

           
   }

}


void CMagProfileArray::DumpParameters(FILE * dat){

  assert(dat);  
  CMagObject::DumpParameters(dat);

  if(m_pPosfixed[X_OFFSET])
    fprintf(dat,"X-shift = %lg FIXED\n", m_lfPos[X_OFFSET]); 
  else
    fprintf(dat,"X-shift = %lg +/- %lg\n", m_lfPos[X_OFFSET], m_lfStdDevPos[X_OFFSET]); 

  if(m_pPosfixed[Y_OFFSET])
    fprintf(dat,"Y-shift = %lg FIXED\n", m_lfPos[Y_OFFSET]); 
  else
    fprintf(dat,"Y-shift = %lg +/- %lg\n", m_lfPos[Y_OFFSET], m_lfStdDevPos[Y_OFFSET]); 


  if(m_pPosfixed[Z_OFFSET])
    fprintf(dat,"Z-shift = %lg FIXED\n", m_lfPos[Z_OFFSET]);
  else
    fprintf(dat,"Z-shift = %lg +/- %lg\n", m_lfPos[Z_OFFSET], m_lfStdDevPos[Z_OFFSET]); 

}


void CMagProfileArray::DumpBinary(FILE * dat) {

  double a, dtime;
  int i, j, k=0;
  char line[20];

  int n = m_iStop-m_iStart;
  fwrite(&n, 1, sizeof(int), dat);
  fwrite(&m_nChannels, 1, sizeof(int), dat);
 
  strcpy(line, GetName());
 
  for(i=m_iStart; i<m_iStop; i++) {
      
      dtime = GetTime(k);
      
      for(j=0; j<m_nChannels; j++) {
       a = GetX(k);         fwrite(&a, 1, sizeof(double),dat); 
       a = GetY(k);         fwrite(&a, 1, sizeof(double),dat); 
       a = GetZ(k);         fwrite(&a, 1, sizeof(double),dat);
       a = GetAlt(k);       fwrite(&a, 1, sizeof(double),dat);
       a = GetField(k);     fwrite(&a, 1, sizeof(double),dat);
       a = GetSyntField(k); fwrite(&a, 1, sizeof(double),dat);
       k++;
      }        
      fwrite(&dtime, 1, sizeof(double),dat);
      fwrite(line, 1, 20, dat);           
   }

}



int CMagProfileArray::DepthSplineSmooth(double dm) {
    
   int ret = -1;
    
   int n =  m_vTime.size();
   int i,j;
   
   double * x = new double[n];
   double * u = new double[n];
   double * a = new double[n];
   double * b = new double[n];
   double * c = new double[n];
   double * d = new double[n];
   
   if(!x || !u || !a || !b || !c || !d) goto m_ret;
    
   for(i=0; i<n; i++) x[i] = double(i);
   
   for(i=0; i<m_nChannels; i++) {
       for(j=0; j<n; j++)  u[j] = m_vZ[j*m_nChannels+i];              
   
       int ret2 = splik(n, x, u, dm, a, b,c, d); 
    
       for(j=0; j<n; j++)  m_vZ[j*m_nChannels+i] = d[j];
   
   }
      ret = 0;       
       
 m_ret:
     if(x) delete [] x;
     if(u) delete [] u;      
     if(a) delete [] a;
     if(b) delete [] b;
     if(c) delete [] c;
     if(d) delete [] d;
 
    return ret;   
    
}    


int CMagProfileArray::CheckBackgroundL1(double * res) {
    
  int ret = -1;    
  int m,n, i,j, err =1;
  double *ct, *f, *a, *r;
  double prec, eps,z;
  int irank, iter, ind;
  double max_res = -1.;

  prec = 1.e-6;
  eps = 1.e-4;

  ct = f = a = r = NULL;

  n = m_iStop - m_iStart;
  m = 3;

  ct = new double[m*n];
  f  = new double[n];
  a  = new double[m];
  r  = new double[n];
  
  if(!ct || !r || !a || !r) goto m_ret; 
      
   for(i=0; i<m_nChannels; i++) {
       
       int i1, j1;
       for(i1=j1=0; i1<n; i1++,j1+=3){
            ct[j1]   = double(i1*i1);
            ct[j1+1] = (double)i1;
            ct[j1+2] = 1.;
       }    
       
       for(j=0; j<n; j++)  f[j] = m_vField[m_nChannels*(j+m_iStart)+i];               
       l1_driver_d(m, n, ct,  f, a, prec, eps, &irank, &iter, r, &z, &ind, NULL);   
   
       for(i1=0; i1<n; i1++) max_res = MAX(max_res, fabs(r[i1]));
   
   }   
      
   *res = max_res;      
   ret = 1;
    
m_ret:
    if(ct) delete [] ct;
    if(f)  delete [] f;
    if(a)  delete [] a;
    if(r)  delete [] r;
           

    return ret;   
    
}    


void CMagProfileArray::RecomputePositions90(double * shifts90){
    
	  int i,j;
      int n = m_vXmiddle.size();    
      double nx90[2], ny90[2], nx, ny, x, y, nx0, ny0;
      nx90[0] = ny90[0] = 0.;
      
      for(i=0; i<n; i++) {

      double s = 0.;
      for(j=i+1; j<n; j++) {
           nx = m_vXmiddle[j] - m_vXmiddle[i];
           ny = m_vYmiddle[j] - m_vYmiddle[i];          
           s = sqrt(nx*nx + ny*ny);
           if(s!=0.) break; 
       }       
       
       if(s!=0.) {    
           nx /= s; ny /= s;
           nx90[1] = -ny; ny90[1] = nx;
        }
        else {
            nx90[1] = nx90[0]; ny90[1] = ny90[0];                
        }    
    
        nx = (nx90[1] + nx90[0])/2.;
        ny = (ny90[1] + ny90[0])/2.;  
        
        s = sqrt(nx*nx + ny*ny);
        nx /= s; ny /= s;                    

        nx90[0] = nx90[1];
        ny90[0] = ny90[1];
        nx0 = ny;
        ny0 = -nx;
      
        for(int j=0; j<m_nChannels; j++) {
              x = m_vXmiddle[i] + nx*shifts90[j];
              y = m_vYmiddle[i] + ny*shifts90[j];
              if(i) {                      
               if((x-m_vX[(i-1)*m_nChannels+j])*nx0 + 
                  (y-m_vY[(i-1)*m_nChannels+j])*ny0 < 0.) {
                  x = (x + m_vX[(i-1)*m_nChannels+j])/2.;
                  y = (y + m_vY[(i-1)*m_nChannels+j])/2.;
                   m_vX[(i-1)*m_nChannels+j] = x;
                   m_vY[(i-1)*m_nChannels+j] = y;                                
                 }     
               }    
              m_vX[i*m_nChannels+j] = x;
              m_vY[i*m_nChannels+j] = y; 
        }
    }        
    
}    
