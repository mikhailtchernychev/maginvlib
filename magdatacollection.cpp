//
// @(#)magdatacollection.cpp  1.1  misha-10sep103
//
// Copyright (c) 2003 Geometrics
//
// Created:  10sep103  by  misha@misha2.geometrics.com
// Version:  10sep103  by  misha@misha2.geometrics.com
//

#include "magdatacollection.h"


CMagDataCollection::CMagDataCollection() {
  
  init();

}

void CMagDataCollection::init() {
  
  m_vData.clear();
  m_lfXaverage        = 0.;
  m_lfYaverage        = 0.;
  m_lfZaverage        = 0.;
  m_iCurrentPos       = 0;
  m_nTotalPoints      = 0;
  m_iCurrentTotalData = 0;
  m_lfX1 = m_lfX2 = m_lfY1 = m_lfY2 = m_lfZ1 = m_lfZ2 = 0.; 
  m_iCurrentData = 0;

}

void CMagDataCollection::destroy() {

  int n = m_vData.size();

  for(int i=0; i<n; i++) {
    CMagData * data = m_vData[i];

   	if(data) {
   	        data->DecrementRefCounter();
   	        if(data->GetRefCounter()<=0) delete data;
        }
    }     
 
  init();

}


CMagDataCollection::~CMagDataCollection() {

    destroy();

}


int  CMagDataCollection::AddData(CMagData * data) {

  double x1, y1, x2, y2, z1, z2;
  data->GetBoundingBox(&x1, &x2, &y1, &y2, &z1, &z2);

  if(m_vData.size()) {
    m_lfX1 = MIN(m_lfX1, x1);
    m_lfX2 = MAX(m_lfX2, x2);
    m_lfY1 = MIN(m_lfY1, y1);
    m_lfY2 = MAX(m_lfY2, y2);
    m_lfZ1 = MIN(m_lfZ1, z1);
    m_lfZ2 = MAX(m_lfZ2, z2);
  }
  else {
    m_lfX1 = x1;
    m_lfX2 = x2;
    m_lfY1 = y1;
    m_lfY2 = y2;
    m_lfZ1 = z1;
    m_lfZ2 = z2;
  }


  m_vData.push_back(data);
  data->IncrementRefCounter();
  
  m_nTotalPoints  += data->get_total_data_points();

  m_lfXaverage = (m_lfX1 + m_lfX2) / 2.;
  m_lfYaverage = (m_lfY1 + m_lfY2) / 2.;
  m_lfZaverage = (m_lfZ1 + m_lfZ2) / 2.;

  return 1;
}



int CMagDataCollection::ReplaceData(CMagData *p_old, CMagData *p_new) {
    
    int i, n;

    n = m_vData.size();
    for(i=0; i<n; i++) {
	if(m_vData[i] == p_old) {
        m_vData[i] = p_new;
		p_new->IncrementRefCounter();
		p_old->DecrementRefCounter();
		ReScan();
        return 0;
        }
	}
  
  return -1;
}



int CMagDataCollection::RemoveData(CMagData * pData) {

	int n = m_vData.size();
	m_vData.erase(std::remove(m_vData.begin(), m_vData.end(), pData), m_vData.end());

	if(n-1==m_vData.size()) { // was actually removed
		ReScan();
		pData->DecrementRefCounter();
		return 1;
	}

	return 0;
}


void  CMagDataCollection::ReScan() {

  m_lfXaverage        = 0.;
  m_lfYaverage        = 0.;
  m_lfZaverage        = 0.;
  m_iCurrentPos       = 0;
  m_nTotalPoints      = 0;
  m_iCurrentTotalData = 0;
  m_lfX1 = m_lfX2 = m_lfY1 = m_lfY2 = m_lfZ1 = m_lfZ2 = 0.; 
  m_iCurrentData = 0;

  double x1, y1, x2, y2, z1, z2;
 
  for(int i=0; i<int(m_vData.size()); i++) {
    CMagData * data = m_vData[i];
    data->GetBoundingBox(&x1, &x2, &y1, &y2, &z1, &z2);      
 
    if(i) {
        m_lfX1 = MIN(m_lfX1, x1);
        m_lfX2 = MAX(m_lfX2, x2);
        m_lfY1 = MIN(m_lfY1, y1);
        m_lfY2 = MAX(m_lfY2, y2);
        m_lfZ1 = MIN(m_lfZ1, z1);
        m_lfZ2 = MAX(m_lfZ2, z2);
    }
      else {
          m_lfX1 = x1;
          m_lfX2 = x2;
          m_lfY1 = y1;
          m_lfY2 = y2;
          m_lfZ1 = z1;
          m_lfZ2 = z2;
      }
      
    m_nTotalPoints  += data->get_total_data_points();
    m_lfXaverage = (m_lfX1 + m_lfX2) / 2.;
    m_lfYaverage = (m_lfY1 + m_lfY2) / 2.;
    m_lfZaverage = (m_lfZ1 + m_lfZ2) / 2.;              
      
  }    
}


int CMagDataCollection::get_data_origin(double *x, double *y, double *z) {
    
    *x = m_lfXaverage;
    *y = m_lfYaverage;
    *z = m_lfZaverage;

    //*x = m_lfX1;
    //*y = m_lfY1;
    //*z = m_lfZ1;

    return 1;
}

int CMagDataCollection::reset() {
  m_iCurrentPos       = 0;
  m_iCurrentTotalData = 0;
  m_iCurrentData      = 0;

  // now reset all data sets

  int n = m_vData.size();
  for(int i=0; i<n; i++) {
    if(m_vData[i]) m_vData[i]->reset();
  }

  return 1;

}


int CMagDataCollection::is_valid(){

  if(m_nTotalPoints > 0) return 1;
  return 0;
}



CMagData * CMagDataCollection::get_next_position(double *xd, double *yd, double *zd, double * Td, int *n_pos,
				                 int * data_pos, double *add_parameters, int size_in, int *size_out) {

  while(m_iCurrentPos < int(m_vData.size())) {
        CMagData * data =  m_vData[m_iCurrentPos];
	if(data) {
	  if(data->get_next_position(xd, yd, zd, Td, n_pos, add_parameters, size_in, size_out)) {
	    *data_pos = *n_pos;  
	    (*n_pos) += m_iCurrentTotalData;
	    return data;
	  }
	  else {
	    m_iCurrentTotalData += data->get_total_data_points();
	    m_iCurrentPos++;
	    continue;
	  }
	}
	m_iCurrentPos++;
    }

  return NULL;

}




CMagData * CMagDataCollection::get_at_position(int n_pos, double *x, double *y, double *z, double * T,
				   int * data_pos, double *add_parameters, int size_in, int *size_out) {


  if(n_pos >=  m_nTotalPoints || n_pos < 0) return NULL;

  int n       = m_vData.size();
  int n_start = 0;

  for(int i=0; i<n; i++) {
    CMagData * data = m_vData[i];
    if(!data) continue;
    if(n_pos >= n_start && n_pos < n_start+data->get_total_data_points()) {	
      int n = n_pos - n_start;
      data->get_at_position(n, x, y, z, T, add_parameters, size_in, size_out);
      *data_pos = n;
      return data;
    }
    n_start += data->get_total_data_points();
  }


  return NULL;
}



CMagData * CMagDataCollection::get_position(int n_pos, double *x, double *y, double *z) {

  double T;
  int data_pos;
  return get_at_position(n_pos, x, y, z, &T, &data_pos);

}



int CMagDataCollection::set_synthetic_data(double *data_in, int replace) {

  int n       = m_vData.size();
  int n_start = 0;

  int j = 0;
  for(int i=0; i<n; i++) {
    CMagData * data = m_vData[i];
    if(!data) continue;
    data->set_synthetic_data(&data_in[n_start], replace);
    n_start += data->get_total_data_points();  
  }

  return 1;
}

CMagData * CMagDataCollection::GetNextData() {

  if(m_iCurrentData >= int(m_vData.size())) return NULL;
  CMagData * data = m_vData[m_iCurrentData];
  m_iCurrentData++;
  return data;
}

void CMagDataCollection::Dump(FILE * dat) {

  int n       = m_vData.size();

  for(int i=0; i<n; i++) {
    if(m_vData[i]) m_vData[i]->Dump(dat);
    fprintf(dat,"\n\n");
  }

}


int CMagDataCollection::GetNearestDataPoint(double x0, double y0, double * d, double *x, double *y, double * dtime, 
	std::string & line_name) {

		std::string sLineName;
		double dist, x1, y1, dtime1;

		int n       = m_vData.size();

		int j = 0;

		for(int i=0; i<n; i++) {
			if(m_vData[i]) {
				if(-1!=m_vData[i]->GetNearestDataPoint(x0, y0, &dist, &x1, &y1, &dtime1,  sLineName)) {
					if(!j) {
						*d = dist; *x = x1; *y = y1; *dtime = dtime1; line_name = sLineName;
					}
					else {
						if(dist < *d) {
							*d = dist; *x = x1; *y = y1; *dtime = dtime1; line_name = sLineName;
						}
					}
					j++;
				}
			}
		}

		if(j) return 1;

		return -1;
}


void CMagDataCollection::DumpToXMLStream(std::ostream & str) {

  for(int i=0; i<int(m_vData.size()); i++) {
    if(m_vData[i]) m_vData[i]->DumpToXMLStream(str);
  }

}


int  CMagDataCollection::GetDataPoints(std::vector<int> & points) {

	 for(int i=0; i<int(m_vData.size()); i++) {
		 if(m_vData[i]) {
			 points.push_back(m_vData[i]->get_total_data_points());
		 }
	 }
   
	 return 0;
}


int CMagDataCollection::GetMinMax(int type, double * data_min, double * data_max) {

  double dmax, dmin, dmax_tmp, dmin_tmp;
  int n       = m_vData.size();
  int k = 0;

  for(int i=0; i<n; i++) {
    if(m_vData[i]) {
     if(m_vData[i]->GetMinMax(type, &dmin_tmp, &dmax_tmp)) {
          if(k) {
                dmin = MIN(dmin_tmp, dmin);
                dmax = MAX(dmax_tmp, dmax);
                }    
          else {
               dmin = dmin_tmp;                     
               dmax = dmax_tmp;
               }
           k++;
           }         
     }
   }   
   
   if(k) {
      *data_min = dmin;
      *data_max = dmax;
      return 1;
      }
      
     return 0;      
   }

















