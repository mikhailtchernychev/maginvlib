//
// @(#)magprofiles.cpp  1.1  misha-08sep103
//
// Copyright (c) 2003 Geometrics
//
//
// Created:  08sep103  by  misha@misha2.geometrics.com
// Version:  08sep103  by  misha@misha2.geometrics.com
//


#include "magprofiles.h"


CMagProfiles::CMagProfiles() :CMagData() {

  m_mapProfiles.clear();
  prof_pos = m_mapProfiles.begin();
  
  m_lfXaverage = 0.;
  m_lfYaverage = 0.;
  m_lfZaverage = 0.;
  m_nCurrentPoints = 0;
  m_iCurrentPos    = 0;
  m_iCurrentProfilePos = 0;
}


void CMagProfiles::destroy() {

  prof_pos = m_mapProfiles.begin(); 

  for(;prof_pos!=m_mapProfiles.end(); ++prof_pos) {
      if(prof_pos->second) delete prof_pos->second;     
    }

  m_mapProfiles.clear();
  prof_pos = m_mapProfiles.begin();

  m_lfXaverage = 0.;
  m_lfYaverage = 0.;
  m_lfZaverage = 0.;
  m_nCurrentPoints = 0;
  m_iCurrentPos    = 0;
  m_iCurrentProfilePos = 0;
}


CMagProfiles::~CMagProfiles() {

    destroy();

}


void CMagProfiles::Add(char * name, double x, double y, double time, double field, double z) {


  if(!m_nCurrentPoints) {
    m_lfXaverage = x;
    m_lfYaverage = y;
    m_lfZaverage = z;
    m_nCurrentPoints = 1;
  }
  else {
    int n = m_nCurrentPoints+1;
    m_lfXaverage =  (m_lfXaverage*m_nCurrentPoints + x) / n;
    m_lfYaverage =  (m_lfYaverage*m_nCurrentPoints + y) / n;
    m_lfZaverage =  (m_lfZaverage*m_nCurrentPoints + z) / n;
    m_nCurrentPoints++;
  }
    


   string key = name;

   map<string, CProfWrapper *>::iterator pos;

   pos = m_mapProfiles.find(key);

   if(pos != m_mapProfiles.end()) ((CTestProfile *)pos->second)->Add(x,y,time, field,z);
   

   CTestProfile * prof = new CTestProfile;
   if(!prof) return;

   prof->InitProfile(name);
   prof->Add(x,y,time,field,z);

   m_mapProfiles.insert(make_pair(key,prof));
  
 }



int CMagProfiles::get_data_origin(double *x, double *y, double *z) {
    
    *x = m_lfXaverage;
    *y = m_lfYaverage;
    *z = m_lfZaverage;

    return 1;
}

int CMagProfiles::reset(){

  prof_pos = m_mapProfiles.begin();
  m_iCurrentPos        = 0;
  m_iCurrentProfilePos = 0;
  return 1;

}


int CMagProfiles::is_valid(){

  if(get_total_data_points() > 0) return 1;
  return 0;
}


int CMagProfiles::get_total_data_points() {
 
  int n_total = 0;
  prof_pos = m_mapProfiles.begin(); 

  for(;prof_pos!=m_mapProfiles.end(); ++prof_pos) {
    if(prof_pos->second) n_total += prof_pos->second->GetTotalPoints();     
    }

  prof_pos = m_mapProfiles.begin(); 
  m_nCurrentPoints = n_total; 
  m_iCurrentPos        = 0;
  m_iCurrentProfilePos = 0;
  return n_total;

}




int CMagProfiles::get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				   double *add_parameters, int size_in, int *size_out) {


  if(prof_pos == m_mapProfiles.end()) return 0;

  CTestProfile *prof = (CTestProfile *)(prof_pos->second);

  if(m_iCurrentProfilePos >= prof->GetTotalPoints()) {
    prof_pos++;
    m_iCurrentProfilePos = 0;
    if(prof_pos == m_mapProfiles.end()) return 0;
    prof = (CTestProfile *)(prof_pos->second);
  }
  
  *xd = prof->GetX( m_iCurrentProfilePos);
  *yd = prof->GetY( m_iCurrentProfilePos);
  *zd = prof->GetZ( m_iCurrentProfilePos);
  *Td = prof->GetField(m_iCurrentProfilePos);
  *n_pos = m_iCurrentPos;

  m_iCurrentPos++;
  m_iCurrentProfilePos++;

  return 1;  

}




int CMagProfiles::get_at_position(int n_pos, double *x, double *y, double *z, double * T,
				   double *add_parameters, int size_in, int *size_out) {


  if(n_pos >=m_nCurrentPoints || n_pos < 0) return 0;

  map<string, CProfWrapper *>::iterator pos = m_mapProfiles.begin();
  if(pos == m_mapProfiles.end()) return 0;

  int n_start = 0;

  for(; pos != m_mapProfiles.end(); ++pos) {
    CTestProfile * prof = (CTestProfile *)(pos->second);
    if(!prof) continue;
    if(n_pos >= n_start && n_pos < n_start+prof->GetTotalPoints()) {
      int n = n_pos - n_start;
        *x = prof->GetX(n);
        *y = prof->GetY(n);
        *z = prof->GetZ(n);
        *T = prof->GetField(n);
	return 1;
    }
    n_start += prof->GetTotalPoints();
  }


  return 0;
}



int CMagProfiles::get_position(int n_pos, double *x, double *y, double *z) {

  double T;
  return get_at_position(n_pos, x, y, z, &T);

}



int CMagProfiles::set_synthetic_data(double *data) {

  map<string, CProfWrapper *>::iterator pos = m_mapProfiles.begin();
  if(pos == m_mapProfiles.end()) return 0;

  int n_start = 0;

  for(; pos != m_mapProfiles.end(); ++pos) {
    CTestProfile * prof = (CTestProfile *)(pos->second);
    if(!prof) continue;
    prof->set_synthetic_data(&data[n_start]);
    n_start += +prof->GetTotalPoints();
  }

  return 1;
}


void CMagProfiles::Dump(FILE * dat) {

  fprintf(dat,"X   Y    Z   T  Line Tsynt\n"); 

  map<string, CProfWrapper *>::iterator pos = m_mapProfiles.begin();
  if(pos == m_mapProfiles.end()) return;

  for(; pos != m_mapProfiles.end(); ++pos) {
    CTestProfile * prof = (CTestProfile *)(pos->second);
    if(!prof) continue;

    int n_synt = prof->GetSyntheticPoints();

    for(int i=0; i<prof->GetTotalPoints(); i++) { 
      fprintf(dat,"%lf  %lf  %lf  %lf %s",  prof->GetX(i), prof->GetY(i), prof->GetZ(i), prof->GetField(i), prof->GetName());  
      if(n_synt) fprintf(dat," %lf\n", prof->GetSyntField(i));
      else fprintf(dat,"\n");
    }

  }

}
