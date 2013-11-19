//
// @(#)magcollection.cpp  1.1  misha-11may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  11may103  by  misha@misha.local
// Version:  11may103  by  misha@misha.local
//


#include "magcollection.h"


CMagCollection::CMagCollection() {

	m_vObjects.clear();
	m_vCurrentLinearParams.clear();
	m_vCurrentNonLinearParams.clear();	

	m_iTotalLinearParameters    = 0;
        m_iTotalNonLinearParameters = 0;

        m_itLinear        = m_vCurrentLinearParams.begin();
        m_itNonLinear     = m_vCurrentNonLinearParams.begin();
        m_itCurrentObject = m_vObjects.begin();   

}


void CMagCollection::destroy() {

    vector<CMagObject *>::iterator p = m_vObjects.begin();

    for(;p!=m_vObjects.end(); ++p) {
	CMagObject * obj = (CMagObject *) *p;

	if(obj) {
   	        obj->DecrementRefCounter();
   	        if(obj->GetRefCounter()<=0) delete obj;
        }
    }
    
    m_vObjects.clear();
    m_vCurrentLinearParams.clear();
    m_vCurrentNonLinearParams.clear();	
    
    m_iTotalLinearParameters    = 0;
    m_iTotalNonLinearParameters = 0;
}


CMagCollection::~CMagCollection() {
    destroy();
}



int CMagCollection::AddMagObject(CMagObject * o) {

    m_vObjects.push_back(o);
    o->IncrementRefCounter();

    o->UpdateReportedParams();
    m_iTotalLinearParameters    += o->GetTotalLinearParams();
    m_iTotalNonLinearParameters += o->GetTotalNonLinearParams();

    m_vCurrentLinearParams.push_back(m_iTotalLinearParameters);
    m_vCurrentNonLinearParams.push_back(m_iTotalNonLinearParameters);

    return 1;
}

int CMagCollection::ResetPosition() {

    m_itLinear        = m_vCurrentLinearParams.begin();
    m_itNonLinear     = m_vCurrentNonLinearParams.begin();
    m_itCurrentObject = m_vObjects.begin(); 

    return 1;
}

CMagObject * CMagCollection::GetCurrentMagObject(int *nolin_start, int *lin_start) {
    
    if(m_itCurrentObject == m_vObjects.end()) return NULL;

    // note that non-linear parameters go first

    CMagObject * obj = (CMagObject *) *m_itCurrentObject;
    if(obj) {
          *lin_start = (int) *m_itLinear   + m_iTotalNonLinearParameters - obj->GetTotalLinearParams();
	  *nolin_start = (int)*m_itNonLinear - obj->GetTotalNonLinearParams();
	  obj->SetLinNonLinStarts(*lin_start, *nolin_start);
    }

     ++m_itCurrentObject;     
     ++m_itLinear;
     ++m_itNonLinear;
  
     return obj;

}

CMagObject * CMagCollection::GetObjectAt(int n) {
    
    if(n >= int(m_vObjects.size())) return NULL;
    
    return (CMagObject *)m_vObjects[n];
}
    

CMagObject * CMagCollection::GetObjectAt(char * name) {
    
    int i, n;
    char obj_name[80];

    n = m_vObjects.size();
    CMagObject * p;
    for(i=0; i<n; i++) {
	p = GetObjectAt(i);
	if(p) {
	    p->GetObjectName(obj_name, 80);	    
	    if(strcmp(name,  obj_name) == 0) return p;
	}
    }

    return NULL;
}


int CMagCollection::ReplaceObject(CMagObject *p_old, CMagObject *p_new) {
    
    int i, n;

    n = m_vObjects.size();
    for(i=0; i<n; i++) {
	if(GetObjectAt(i) == p_old) {
        m_vObjects[i] = p_new;
		p_new->IncrementRefCounter();
		p_old->DecrementRefCounter();
        return 0;
        }
	}
  
  return -1;
}

int CMagCollection::RemoveObject(CMagObject * pObject) {

	int n = m_vObjects.size();
	m_vObjects.erase(std::remove(m_vObjects.begin(), m_vObjects.end(), pObject), m_vObjects.end());

	if(n-1==m_vObjects.size()) { // was actually removed
		pObject->DecrementRefCounter();
		return 1;
	}

	return 0;
}



int CMagCollection::GetTotalLinearParams() {

   if(m_vCurrentLinearParams.size() == 0) return 0;
    vector<int>::iterator i = m_vCurrentLinearParams.end();
    --i;

    return *i;
}


int CMagCollection::GetTotalNonLinearParams() {

   if(m_vCurrentNonLinearParams.size() == 0) return 0;
    vector<int>::iterator i = m_vCurrentNonLinearParams.end();
    --i;

    return *i;

}


// --------------- get/set variables as signle vector ---------------------
//
// is_set      - 0 to get, 1 to set
// vars        - where to get variables
// size        - size of vars
// type        - 0 for linear varibales, 1 nor non-linear
// linear_only - shift linear variables in the fisrt part of vars
// is_errors   - get standard deviations instead of variables


int CMagCollection::GetSetVariables(int is_set, double * vars, int size, 
	int type, int linear_only, int is_errors) {


		typedef int (CMagObject::*ParametersFunc)(double * params,    int size, int type);
		typedef int (CMagObject::*ParameterNumberFunc)();

		ParameterNumberFunc pGetTotal;
		ParametersFunc      pGetParameters;

		if(!is_set) {
			switch(type) { // getting parameters
			case 0: pGetTotal      = &CMagObject::GetTotalLinearParams;
				pGetParameters = &CMagObject::GetLinearParams;
				break;
			case 1: pGetTotal      = &CMagObject::GetTotalNonLinearParams;
				pGetParameters = &CMagObject::GetNonLinearParams;
				break;
			default: return 0;	  
			} 
		}
		else {       // setting parameters
			switch(type) {
			case 0: pGetTotal      = &CMagObject::GetTotalLinearParams;
				pGetParameters = &CMagObject::SetLinearParams;
				break;
			case 1: pGetTotal      = &CMagObject::GetTotalNonLinearParams;
				pGetParameters = &CMagObject::SetNonLinearParams;
				break;
			default: return 0;	  
			} 
		}

		int nonlin_start, lin_start;
		int m1  = 0;
		int ret = 0;	
		int n_nonlin = 0;
		int k_start = 0;
		if(linear_only && type==0) n_nonlin = GetTotalNonLinearParams();

		ResetPosition();

		CMagObject * p = GetCurrentMagObject(&nonlin_start, &lin_start);

		while(p) {
			lin_start -= n_nonlin;
			k_start = (type == 0) ? lin_start : nonlin_start;
			m1 = (p->*pGetTotal)();
			(p->*pGetParameters)(&vars[k_start], size-k_start, is_errors);
			p = GetCurrentMagObject(&nonlin_start, &lin_start);
			ret += m1;
		}

		return ret;

}



void CMagCollection::DumpObjects(FILE * dat) {

  assert(dat);

  vector<CMagObject *>::iterator p = m_vObjects.begin();
  int i =0;

    for(;p!=m_vObjects.end(); ++p) {
      fprintf(dat,"\n\nObject #%d\n", i+1);
	CMagObject * obj = (CMagObject *) *p;
	if(obj) {
	  obj->DumpParameters(dat);
	  i++;
	}
    }
}

 
void CMagCollection::InitMinMaxForObjects() {

  vector<CMagObject *>::iterator p = m_vObjects.begin();
   for(;p!=m_vObjects.end(); ++p) {
       	CMagObject * obj = (CMagObject *) *p;
	    if(obj) obj->InitMinMaxParams();
    }
}

void CMagCollection::UpdateMinMaxForObjects() {

  vector<CMagObject *>::iterator p = m_vObjects.begin();
   for(;p!=m_vObjects.end(); ++p) {
       	CMagObject * obj = (CMagObject *) *p;
	    if(obj) obj->UpdateMinMaxParams();
    }
}

void CMagCollection::DumpMinMaxForObjects(FILE * dat) {
  assert(dat);

  vector<CMagObject *>::iterator p = m_vObjects.begin();
  int i =0;

    for(;p!=m_vObjects.end(); ++p) {
      fprintf(dat,"\n\nMin Max dump: Object #%d\n", i+1);
      	CMagObject * obj = (CMagObject *) *p;
      		if(obj) {
      		obj->DumpMinMaxParams(dat);
      		i++;
	     }
    }
    
}    
    




int CMagCollection::ReScanAllObjects(int * total_linear, int * total_non_linear) {

  vector<CMagObject *>::iterator p = m_vObjects.begin();

  m_iTotalLinearParameters    = 0;
  m_iTotalNonLinearParameters = 0;
  m_vCurrentLinearParams.clear();
  m_vCurrentNonLinearParams.clear();	

  for(;p!=m_vObjects.end(); ++p) {
	CMagObject * o = (CMagObject *) *p;
	if(o) {
          o->UpdateReportedParams();
	  m_iTotalLinearParameters    += o->GetTotalLinearParams();
	  m_iTotalNonLinearParameters += o->GetTotalNonLinearParams();
	  m_vCurrentLinearParams.push_back(m_iTotalLinearParameters);
	  m_vCurrentNonLinearParams.push_back(m_iTotalNonLinearParameters);
	}
    }
   
  ResetPosition();

  *total_linear     = m_iTotalLinearParameters;
  *total_non_linear = m_iTotalNonLinearParameters;

  return 1;

}


int CMagCollection:: GetAbilities(int *linear_ab, int *non_linear_ab) {

  vector<CMagObject *>::iterator p = m_vObjects.begin();
  int i =0;

  *linear_ab = 0;
  *non_linear_ab = 0;

  for(;p!=m_vObjects.end(); ++p) {
	CMagObject * o = (CMagObject *) *p;
	if(o) {
	  if(!o->CanComputeField()) continue;
	  if(!i) {
	    *linear_ab     = o->GetObjectLinearAbilities(); 
	    *non_linear_ab = o->GetObjectNonLinearAbilities(); 
	    i++;
            continue; 
	  }
           *linear_ab     = MIN(*linear_ab,     o->GetObjectLinearAbilities()); 
           *non_linear_ab = MIN(*non_linear_ab, o->GetObjectNonLinearAbilities());
	   i++;
	}
    }
   
   return 1;

}




void CMagCollection::DumpToStream(std::ostream & str, std::string comment,  char * close_str,
								  char * delimit, std::string offset) {

   str << comment.c_str() << "Objects" << std::endl;
   str << comment.c_str()  << "TOTAL_OBJECTS=" << m_vObjects.size() << std::endl;

   std::string comment1 = comment;
   comment1.insert(0, offset.c_str());

   vector<CMagObject *>::iterator p = m_vObjects.begin();
    int i =0;

    for(;p!=m_vObjects.end(); ++p) {
	   CMagObject * obj = (CMagObject *) *p;
	   if(obj) {
	    obj->DumpToStream(str, comment1, close_str, delimit, offset);
	  i++;
	   }
    }

	str << comment.c_str()  << close_str << "Objects" << std::endl;


}


void CMagCollection::DumpToXMLStream(std::ostream & str) {

  char quot = '\"';
  str << "<Objects>" << std::endl;
  str << "<params " << "TOTAL_OBJECTS=" << quot << m_vObjects.size() << quot << "/>" << std::endl;

  vector<CMagObject *>::iterator p = m_vObjects.begin();
    int i =0;

    for(;p!=m_vObjects.end(); ++p) {
	   CMagObject * obj = (CMagObject *) *p;
	   if(obj) {
	    obj->DumpToXMLStream(str);
	  i++;
	   }
    }

	str << "</Objects>" << std::endl;
}



