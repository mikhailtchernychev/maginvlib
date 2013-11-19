//
// @(#)magobject.cpp  1.1  misha-09may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  09may103  by  misha@MISHA
// Version:  09may103  by  misha@MISHA
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "magobject.h"


CMagObject::CMagObject() {

  m_bIsValid         = 0;  // if object is valid
  m_bInduced         = 0;  // magnetization type
  m_nLinearParams    = 0;  // total number of linear (J) parameters 
  m_nNonLinearParams = 0;  // total number of non-linear (geometrics) parameters 

  m_lfJ   = NULL;          // linear parameters (magnetization)
  m_lfPos = NULL;          // non-linear parameters (positions)

  m_iLinearStart    = -1;  // valid only after calling CMagCollection::GetCurrentMagObject
  m_iNonLinearStart = -1;

  m_lfStdDevJ   = NULL;    // standard deviation for linear parameters estimation
  m_lfStdDevPos = NULL;    // standard deviation for non-linear parameters estimation

  m_pJfixed   = NULL;  
  m_pPosfixed = NULL;

  m_nLinearParamsReport    = 0;  // valid only after calling CMagCollection::GetCurrentMagObject
  m_nNonLinearParamsReport = 0;

  strcpy(m_sObjectType,"CMagObject");
  strcpy(m_sObjectName,"NoNameObject");
 
   m_iFieldType = anyfield;

  m_iObjectLinearAbilities    = 0;
  m_iObjectNonLinearAbilities = 0;

  m_lfJUp   = NULL;
  m_lfPosUp = NULL; 

  m_lfJLow   = NULL;
  m_lfPosLow = NULL;

  m_lfJMax   = NULL;
  m_lfPosMax = NULL;

  m_lfJMin   = NULL;
  m_lfPosMin = NULL;

  m_iRefCounter = 0;
  
  m_sID.clear();
}

CMagObject::CMagObject(int is_induced, int n_linear, int n_nonlinear, char *type, char *name) {

  int i;

  m_bIsValid         = 0;  
  m_bInduced         = is_induced;
  m_nLinearParams    = 0; 
  m_nNonLinearParams = 0;

  m_iLinearStart    = -1; // valid only after calling CMagCollection::GetCurrentMagObject
  m_iNonLinearStart = -1;

  m_lfJ     = new double[n_linear];
  m_pJfixed = new int[n_linear];
  if(m_lfJ && m_pJfixed) m_nLinearParams = n_linear;

  m_lfPos     = new double[n_nonlinear];
  m_pPosfixed = new int[n_nonlinear];
  if(m_lfPos && m_pPosfixed) m_nNonLinearParams = n_nonlinear;

  m_lfStdDevJ = new double[n_linear];
  if(m_lfStdDevJ) m_nLinearParams = n_linear;

  m_lfStdDevPos = new double[n_nonlinear];
  if(m_lfStdDevPos) m_nNonLinearParams = n_nonlinear;

  for(i=0; i<m_nLinearParams;    i++) { 
    m_pJfixed[i] = 0;  
    m_lfJ[i]     = m_lfStdDevJ[i] = 0.;
  }

  m_lfJUp    = new double[n_linear];
  m_lfPosUp  = new double[n_nonlinear];
  m_lfJLow   = new double[n_linear];
  m_lfPosLow = new double[n_nonlinear];

  m_lfJMax   = new double[n_linear];
  m_lfPosMax = new double[n_nonlinear];

  m_lfJMin   = new double[n_linear];   
  m_lfPosMin = new double[n_nonlinear];
  
  
  for(i=0; i<m_nNonLinearParams; i++) {
    m_pPosfixed[i] = 0;
    m_lfPos[i]     = m_lfStdDevPos[i] = 0.;   
    m_lfPosUp[i]   = 1.E38;
    m_lfPosLow[i]  = -1.E38;
    m_lfPosMax[i]  = -1.E38;
    m_lfPosMin[i]  = 1.E38;   
  }

  for(i=0; i<m_nLinearParams; i++) {
     m_lfJUp[i]  =  1.E38;
     m_lfJLow[i] = -1.E38;
     m_lfJMax[i] =  -1.E38;
     m_lfJMin[i] = 1.E38;   
  }

  m_nLinearParamsReport    = m_nLinearParams;  // valid only after calling CMagCollection::GetCurrentMagObject
  m_nNonLinearParamsReport = m_nNonLinearParams;

  if(type) strcpy(m_sObjectType, type);
  if(name) strcpy(m_sObjectName, name);

  m_bIsValid         = 1;  
  m_iFieldType = anyfield;

  m_iObjectLinearAbilities    = 0;
  m_iObjectNonLinearAbilities = 0;

  m_iRefCounter = 0;

  m_sID.clear();


}


void CMagObject::InitMinMaxParams() {

 int i;
    
 for(i=0; i<m_nNonLinearParams; i++) {
    m_lfPosMax[i]  = -1.E38;
    m_lfPosMin[i]  = 1.E38;   
  }

  for(i=0; i<m_nLinearParams; i++) {
     m_lfJMax[i] =  -1.E38;
     m_lfJMin[i] = 1.E38;   
  }    

}


void CMagObject::UpdateMinMaxParams() {

 int i; 
    
 for(i=0; i<m_nNonLinearParams; i++) {
    m_lfPosMax[i]  = MAX(m_lfPosMax[i], m_lfPos[i]);
    m_lfPosMin[i]  = MIN(m_lfPosMin[i], m_lfPos[i]);
  }

  for(i=0; i<m_nLinearParams; i++) {
     m_lfJMax[i] = MAX(m_lfJMax[i], m_lfJ[i]);
     m_lfJMin[i] = MAX(m_lfJMin[i], m_lfJ[i]); 
  }    

}
    
void CMagObject::DumpMinMaxParams(FILE * out) {

 int i;
 
 fprintf(out, "Non-linear min-max:\n");   
 for(i=0; i<m_nNonLinearParams; i++) {
    fprintf(out,"Param: %lg Min %lg Max %lg\n", 
        m_lfPos[i], m_lfPosMin[i], m_lfPosMax[i]);
  }
  
 fprintf(out, "Linear min-max:\n");
 for(i=0; i<m_nLinearParams; i++) {
 fprintf(out,"Param: %lg Min %lg Max %lg\n", 
        m_lfJ[i], m_lfJMin[i], m_lfJMax[i]);
  }    

}

void CMagObject::destroy() {

  m_nLinearParams    = 0; 
  m_nNonLinearParams = 0;
  m_iRefCounter = 0;
  m_sID.clear();

  if(m_lfJ)         { delete [] m_lfJ;         m_lfJ       = NULL;   }
  if(m_lfPos)       { delete [] m_lfPos;       m_lfPos     = NULL;   }
  if(m_pJfixed)     { delete [] m_pJfixed;     m_pJfixed   = NULL;   }
  if(m_pPosfixed)   { delete [] m_pPosfixed;   m_pPosfixed = NULL;   }
  if(m_lfJUp)       { delete [] m_lfJUp;       m_lfJUp     = NULL;   }
  if(m_lfPosUp)     { delete [] m_lfPosUp;     m_lfPosUp   = NULL;   }
  if(m_lfJLow)      { delete [] m_lfJLow;      m_lfJLow    = NULL;   }
  if(m_lfPosLow)    { delete [] m_lfPosLow;    m_lfPosLow  = NULL;   }
  if(m_lfJMax)      { delete [] m_lfJMax;      m_lfJMax    = NULL;   }
  if(m_lfPosMax)    { delete [] m_lfPosMax;    m_lfPosMax  = NULL;   }
  if(m_lfJMin)      { delete [] m_lfJMin;      m_lfJMin    = NULL;   }
  if(m_lfPosMin)    { delete [] m_lfPosMin;    m_lfPosMin  = NULL;   } 
  if(m_lfPosMin)    { delete [] m_lfPosMin;    m_lfPosMin  = NULL;   } 
  if(m_lfStdDevJ)   { delete [] m_lfStdDevJ;   m_lfStdDevJ  = NULL;   } 
  if(m_lfStdDevPos) { delete [] m_lfStdDevPos; m_lfStdDevPos  = NULL;   } 

}


void CMagObject::copy(CMagObject const &o) {

  destroy();

  m_sID = o.m_sID;

  m_bIsValid         = o.m_bIsValid;  
  m_bInduced         = o.m_bInduced;

  m_lfJ =  new double[o.m_nLinearParams];
  if(m_lfJ) {
    m_nLinearParams    = o.m_nLinearParams;
    memcpy(m_lfJ, o.m_lfJ, m_nLinearParams*sizeof(double));
  }

  m_lfPos =  new double[o.m_nNonLinearParams];
  if(m_lfPos) {
    m_nNonLinearParams    = o.m_nNonLinearParams;
    memcpy(m_lfPos, o.m_lfPos, m_nNonLinearParams*sizeof(double));
  }

  m_lfStdDevJ =  new double[o.m_nLinearParams];
  if(m_lfStdDevJ) {
    memcpy(m_lfStdDevJ, o.m_lfStdDevJ, m_nLinearParams*sizeof(double));
  }

  m_lfStdDevPos =  new double[o.m_nNonLinearParams];
  if(m_lfStdDevPos) {
    memcpy(m_lfStdDevPos, o.m_lfStdDevPos, m_nNonLinearParams*sizeof(double));
  }

  m_pJfixed = new int[o.m_nLinearParams];
  if(m_pJfixed) {
    memcpy(m_pJfixed, o.m_pJfixed, m_nLinearParams*sizeof(int));
  }

  m_pPosfixed =  new int[o.m_nNonLinearParams];
  if(m_pPosfixed) {
    memcpy(m_pPosfixed, o.m_pPosfixed, m_nNonLinearParams*sizeof(int));
  }

  m_lfJUp = new double[o.m_nLinearParams];
  if(m_lfJUp) {
    memcpy(m_lfJUp, o.m_lfJUp, m_nLinearParams*sizeof(double));
  }

  m_lfPosUp = new double[o.m_nNonLinearParams];
  if(m_lfPosUp) {
    memcpy(m_lfPosUp, o.m_lfPosUp, m_nNonLinearParams*sizeof(double));
  }

  m_lfJLow = new double[o.m_nLinearParams];
  if(m_lfJLow) {
    memcpy(m_lfJLow, o.m_lfJLow, m_nLinearParams*sizeof(double));
  }

  m_lfPosLow = new double[o.m_nNonLinearParams];
  if(m_lfPosLow) {
    memcpy(m_lfPosLow, o.m_lfPosLow, m_nNonLinearParams*sizeof(double));
  }

 m_lfJMax = new double[o.m_nLinearParams];
  if(m_lfJMax) {
    memcpy(m_lfJMax, o.m_lfJMax, m_nLinearParams*sizeof(double));
  }

  m_lfPosMax = new double[o.m_nNonLinearParams];
  if(m_lfPosMax) {
    memcpy(m_lfPosMax, o.m_lfPosMax, m_nNonLinearParams*sizeof(double));
  }


  m_lfJMin = new double[o.m_nLinearParams];
  if(m_lfJMin) {
    memcpy(m_lfJMin, o.m_lfJMin, m_nLinearParams*sizeof(double));
  }

  m_lfPosMin = new double[o.m_nNonLinearParams];
  if(m_lfPosMin) {
    memcpy(m_lfPosMin, o.m_lfPosMin, m_nNonLinearParams*sizeof(double));
  }

  strcpy(m_sObjectType, o.m_sObjectType);
  strcpy(m_sObjectName, o.m_sObjectName);

  m_iLinearStart    = o.m_iLinearStart; // valid only after calling CMagCollection::GetCurrentMagObject
  m_iNonLinearStart = o.m_iNonLinearStart;

  m_nLinearParamsReport    = o.m_nLinearParamsReport;
  m_nNonLinearParamsReport = o.m_nNonLinearParamsReport;

  m_iFieldType = o.m_iFieldType;

  m_iObjectLinearAbilities    =  o.m_iObjectLinearAbilities;
  m_iObjectNonLinearAbilities =  o.m_iObjectNonLinearAbilities;


}


// overloaded assignment

CMagObject const & CMagObject::operator=(CMagObject const &o) {

        if (this != &o)
        {
            destroy();
            copy(o);
        }
        // return (reference to) current object for
        // chain-assignments
        return (*this);
}


int CMagObject::GetLinearParams(double * params,    int size, int type) {

  double * pData = NULL;
  
  switch(type) {
      case 0: pData = m_lfJ; break;
      case 1: pData = m_lfStdDevJ; break;
      case 2: pData = m_lfJUp; break;
      case 3: pData = m_lfJLow; break;
      case 4: pData = m_lfJMax; break;
      case 5: pData = m_lfJMin; break;     
  default: break;
  }

  if(!pData) return 0;

  int i, j;
  i = j = 0;
  for(i=0; i<m_nLinearParams; i++) {
    if(!m_pJfixed[i]) {
      params[j] = pData[i]; j++;
      if(j >= size) return 0;
    }
  }


  return 1;
}


int CMagObject::GetNonLinearParams(double * params,    int size, int type) {
 
  double * pData = NULL;
  
  switch(type) {
      case 0: pData = m_lfPos; break;
      case 1: pData = m_lfStdDevPos; break;
      case 2: pData = m_lfPosUp; break;
      case 3: pData = m_lfPosLow; break;
      case 4: pData = m_lfPosMax; break;
      case 5: pData = m_lfPosMin; break;    
  default: break;
  }

  if(!pData) return 0;
 
  int i, j;
  i = j = 0;
  for(i=0; i<m_nNonLinearParams; i++) {
    if(!m_pPosfixed[i] && j<size) {
      params[j] = pData[i]; j++;
     }
  }

  return 1;
}


double CMagObject::GetNonLinearParam(int n_param, int type) {
 
  double * pData = NULL;
  
  switch(type) {
      case 0: pData = m_lfPos; break;
      case 1: pData = m_lfStdDevPos; break;
      case 2: pData = m_lfPosUp; break;
      case 3: pData = m_lfPosLow; break;
      case 4: pData = m_lfPosMax; break;
      case 5: pData = m_lfPosMin; break;       
  default: break;
  }

  if(!pData) return 0.;
 
  return pData[n_param];
}


double CMagObject::GetLinearParam(int n_param, int type) {

  double * pData = NULL;
  
  switch(type) {
      case 0: pData = m_lfJ; break;
      case 1: pData = m_lfStdDevJ; break;
      case 2: pData = m_lfJUp; break;
      case 3: pData = m_lfJLow; break;
      case 4: pData = m_lfJMax; break;
      case 5: pData = m_lfJMin; break;       
  default: break;
  }

  if(!pData) return 0.;

 return pData[n_param];
}


int CMagObject::GetAllPositionParams(double * params,    int size) {

  if(!m_lfPos) return 0;

  for(int i=0; i<m_nNonLinearParams && i<size; i++) {
      params[i] = m_lfPos[i];
    }

  return 1;
}



int CMagObject::SetLinearParams(double * params,    int size, int type) {


  double * pData = NULL;
  
  switch(type) {
      case 0: pData = m_lfJ; break;
      case 1: pData = m_lfStdDevJ; break;
      case 2: pData = m_lfJUp; break;
      case 3: pData = m_lfJLow; break;
      case 4: pData = m_lfJMax; break;
      case 5: pData = m_lfJMin; break;        
  default: break;
  }

  if(!pData) return 0;

  int i, j;
  i = j = 0;
  for(i=0; i<m_nLinearParams; i++) {
    if(!m_pJfixed[i]) {
      pData[i] = params[j]; j++;
      if(j >= size) return 0;
    }
  }

  return 1;
}


int CMagObject::SetNonLinearParams(double * params,    int size, int type) {

  double * pData = NULL;
  
  switch(type) {
      case 0: pData = m_lfPos; break;
      case 1: pData = m_lfStdDevPos; break;
      case 2: pData = m_lfPosUp; break;
      case 3: pData = m_lfPosLow; break;
      case 4: pData = m_lfPosMax; break;
      case 5: pData = m_lfPosMin; break;  
  default: break;
  }

  if(!pData) return 0;

  int i, j;
  i = j = 0;
  for(i=0; i<m_nNonLinearParams; i++) {
    if(!m_pPosfixed[i]) {
      pData[i] = params[j]; j++;
      if(j >= size) return 0;
    }
  }

  return 1;
}


int CMagObject::SetObjectType(char * type) {
 
  if(strlen(type) > OBJECT_TYPE_LEN) return 0;
  strcpy(m_sObjectType, type);
  return 1;	
 	
}


int CMagObject::SetObjectName(char * name) {

  if(strlen(name) > OBJECT_NAME_LEN) return 0;
  strcpy(m_sObjectName, name);
  return 1;	

}
  
int CMagObject::GetObjectType(char * type, int size) {
    
    if(size < OBJECT_TYPE_LEN) return 0;
    strcpy(type, m_sObjectType);	
    return 1;	

}

int CMagObject::GetObjectName(char * name, int size) {

    if(size < OBJECT_NAME_LEN) return 0;
    strcpy(name, m_sObjectName);	
    return 1;	

}

void CMagObject::DumpParameters(FILE * dat){
  
  assert(dat);
  fprintf(dat,"%s Name: %s\n", m_sObjectType, m_sObjectName); 

}


void CMagObject::UpdateReportedParams() {

  int i;

  m_nLinearParamsReport = m_nLinearParams;
  for(i=0; i<m_nLinearParams; i++) m_nLinearParamsReport -=  m_pJfixed[i];
  m_nNonLinearParamsReport = m_nNonLinearParams;
  for(i=0; i<m_nNonLinearParams; i++) m_nNonLinearParamsReport -=  m_pPosfixed[i]; 

} 



int  CMagObject::FixLinearParam(int n_param, int is_fixed, double val) {
   
  if(!m_lfJ || !m_pJfixed || !m_bIsValid) return 0; 

  if(n_param <0 || n_param >= m_nLinearParams) return 0;

  m_lfJ[n_param]     = val;
  m_pJfixed[n_param] = is_fixed;

  return 1;

}

int  CMagObject::FixNonLinearParam(int n_param, int is_fixed, double val) {
	
  if(!m_lfPos || !m_pPosfixed || !m_bIsValid) return 0; 
  if(n_param <0 || n_param >= m_nNonLinearParams) return 0;

  m_lfPos[n_param]     = val;
  m_pPosfixed[n_param] = is_fixed;
  return 1;

}


int CMagObject::IsLinearParamFixed(int n_param) {
    
  if(!m_lfJ || !m_pJfixed || !m_bIsValid) return -1; 
  return  m_pJfixed[n_param];
 
}  
         
int CMagObject::IsNonLinearParamFixed(int n_param) {
    
   if(!m_lfPos || !m_pPosfixed || !m_bIsValid) return -1;  
    
    return m_pPosfixed[n_param];
}    
  


int CMagObject:: FixLinearParamInterval(int n_param, double up, double low) {

  if(!m_lfJUp || !m_lfJLow || !m_bIsValid) return 0; 

  m_lfJUp[n_param]      = up;
  m_lfJLow[n_param]     = low;

  return 1;

}

  int  CMagObject::FixNonLinearParamInterval(int n_param, double up, double low) {

  if(!m_lfPosUp || !m_lfPosLow || !m_bIsValid) return 0; 

  m_lfPosUp[n_param]      = up;
  m_lfPosLow[n_param]     = low;

  return 1;

  }


  void CMagObject::GetFieldTypeText(std::string & str) {
	  
	  switch(m_iFieldType) {
	  case anyfield:           str = "ANYFIELD";             break;
	  case totalmagfield:      str = "TOTALMAGFIELD";        break;
	  case totalmaggradient:   str = "TOTALMAGGRADIENT";     break;
      case totalfinitemaggrad: str = "TOTALFINITEMAGGRAD";   break;
      case absmaggradient:     str = "ABSMAGGARIENT";        break;
      case absfinitegradient:  str = "ABSFINITEGRADIENT";    break;
	  default:
		  str = "UNDEFINED";
	  }
  }






