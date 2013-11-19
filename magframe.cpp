//
// @(#)magframe.cpp  1.1  misha-12jul2012
//
// Copyright (c) 2012 of Mikhail Tchernychev
// All rights reserved.
//

#include "magfunc.h"
#include "magframe.h"


CMagFrame::CMagFrame() {

  m_lfXShift =  m_lfYShift = m_lfZShift = 0;

  ClearFrame();

}

CMagFrame::~CMagFrame() {


}

int CMagFrame::ClearFrame() {
  m_mSensorMap.clear();
  m_lfnX =  m_lfnY =  m_lfnZ = 0.;
  m_iIsValid = 0;
  m_vPositions.clear();
  m_iCurrentPos = 0;
  return 0;
}


// insert or replace position

int CMagFrame::SetPosition(int id, double x, double y, double z) {

  std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.find(id);

  x -= m_lfXShift;
  y -= m_lfYShift;
  z -= m_lfZShift;

  if(pos == m_mSensorMap.end()) {
    RIGID_SENSOR sn;
    sn.x = x; sn.y = y; sn.z = z;
    m_mSensorMap.insert(std::make_pair<int,RIGID_SENSOR>(id, sn));
  }
    else {
      pos->second.x = x;
      pos->second.y = y;
      pos->second.z = z;
  }

  return 0;
}

int CMagFrame::SetPositions(std::string s_file) {

	FILE * dat = fopen(s_file.c_str(), "rt");
	if(!dat) return 0;

	 m_mSensorMap.clear();
	 reset();

	 char buffer[1024];
	 int id; 
	 double x, y, z;

	 // skip header
	 fgets(buffer,1022, dat);

	 while(!feof(dat)) {
		 if(!fgets(buffer,1022, dat)) break;
		 if(4==sscanf(buffer, "%d%lf%lf%lf", &id, &x, &y, &z)) {
			 SetPosition(id, x, y, z);
		 }
	 }

	 fclose(dat);

	 return 1;

}


// set field

int CMagFrame::SetField(int id, double field) {

  std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.find(id);

  if(pos == m_mSensorMap.end()) {
    return -1;
  }
    else {
      pos->second.F = field;
  }

  return 0;
}


int CMagFrame::SetFluxgates(double vx, double vy, double vz) {

  double s = sqrt(vx*vx + vy*vy + vz*vz);
  if(s<=0.) return -1;

  m_lfnX =  vx/s;
  m_lfnY =  vy/s;
  m_lfnZ =  vz/s;
  
  m_iIsValid =  1;

  return 0;
}

int CMagFrame::GetCosines(double *nx, double *ny, double *nz) {

  if(!m_iIsValid) return -1;

  *nx = m_lfnX;
  *ny = m_lfnY;
  *nz = m_lfnZ;

  return 0;
}

int CMagFrame::get_data_origin(double *x, double *y, double *z) {

  *x = *y = *z = 0.;
  return 0;

}

int CMagFrame::reset() {
  
  m_vPositions.clear();

  std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.begin();
  for(;pos!=m_mSensorMap.end(); ++pos) {
	  if(pos->second.F != NO_DATA) m_vPositions.push_back(pos);
  }

  m_iCurrentPos = 0;

  return 0;

}

int CMagFrame::get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				 double *add_parameters, int size_in, int *size_out) {

  if(m_iCurrentPos < m_vPositions.size() && m_iCurrentPos >=0) {
        std::map<int, RIGID_SENSOR>::iterator pos = m_vPositions[m_iCurrentPos];
	*xd    = pos->second.x;
	*yd    = pos->second.y;
	*zd    = pos->second.z;
	*Td    = pos->second.F;
	*n_pos = m_iCurrentPos;
	m_iCurrentPos++;
	GetAllPositionParams(add_parameters, size_in);
	return 1;
  }

 return 0;
 
}


int CMagFrame::get_at_position(int n_pos, double *x, double *y, double *z, double * T, 
			       double *add_parameters, int size_in, int *size_out) {

  if(n_pos <0 || n_pos >= m_vPositions.size()) return 0;

   std::map<int, RIGID_SENSOR>::iterator pos = m_vPositions[n_pos];

   *x   = pos->second.x;
   *y   = pos->second.y;
   *z   = pos->second.z;
   *T   = pos->second.F;
   
   GetAllPositionParams(add_parameters, size_in);

   return 1;

}

int CMagFrame::get_total_data_points() {

  std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.begin();
  int n = 0;
  for(;pos!=m_mSensorMap.end(); ++pos) {
	  if(pos->second.F != NO_DATA) n++;
  }

  return n;

}


int CMagFrame::get_total_positions() {

  return m_mSensorMap.size();

}


int CMagFrame::set_synthetic_data(double *data, int replace) {

	std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.begin();

	int i = 0;

	for(;pos!=m_mSensorMap.end(); ++pos) {

		if(pos->second.F == NO_DATA) continue;

		if(replace)
			pos->second.F = data[i];
		else 
			pos->second.Fsynt = data[i];

		i++;
	}

	return 1;

}

int CMagFrame::get_position(int n, double * xd, double * yd, double * zd) {

  double F;

  return get_at_position(n, xd, yd, zd, &F);


}

int CMagFrame::GetBoundingBox(double *x1, double *x2, double *y1, double *y2, 
			      double *z1, double *z2) {

  double xmin =  1.e+20;
  double xmax = -1.e+20;

  double ymin =  1.e+20;
  double ymax = -1.e+20;

  double zmin =  1.e+20;
  double zmax = -1.e+20;

  std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.begin();


  for(;pos!=m_mSensorMap.end(); ++pos) {
    if(pos->second.F == NO_DATA) continue;
    xmin = MIN(xmin, pos->second.x);
    xmax = MAX(xmax, pos->second.x);
    ymin = MIN(ymin, pos->second.y);
    ymax = MAX(ymax, pos->second.y);
    zmin = MIN(zmin, pos->second.z);
    zmax = MAX(zmax, pos->second.z);
  }

  *x1 = xmin; *x2 = xmax;
  *y1 = ymin; *y2 = ymax;
  *z1 = zmin; *z2 = zmax;

  return 1;

}


int CMagFrame::GetMinMax(int type, double * data_min, double * data_max) {

  double Tmin =  1.e+20;
  double Tmax = -1.e+20;


  std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.begin();


  for(;pos!=m_mSensorMap.end(); ++pos) {

   if(pos->second.F == NO_DATA) continue;

    if(type) {
      Tmin = MIN(Tmin, pos->second.Fsynt);
      Tmax = MAX(Tmax, pos->second.Fsynt);
    }
    else {
      Tmin = MIN(Tmin, pos->second.F);
      Tmax = MAX(Tmax, pos->second.F);
    }

  }

  *data_min = Tmin;
  *data_max = Tmax;

  return 1;

}

void CMagFrame::Dump(FILE * dat) {

	std::map<int, RIGID_SENSOR>::iterator pos =  m_mSensorMap.begin();

	if(!dat) return;
	fprintf(dat,"#ID X  Y Z T Tsynt DIFF\n");
	for(;pos!=m_mSensorMap.end(); ++pos) {
		if(pos->second.F != NO_DATA)
			fprintf(dat,"%d %.3lf %.3lf %.3lf  %.3lf %.3lf %lf\n",
			pos->first, pos->second.x, pos->second.y, pos->second.z,
			pos->second.F, pos->second.Fsynt, pos->second.F - pos->second.Fsynt);
	}

}

