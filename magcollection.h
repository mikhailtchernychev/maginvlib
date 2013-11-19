//
// @(#)magcollection.h  1.1  misha-11may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  11may103  by  misha@misha.local
// Version:  11may103  by  misha@misha.local
//

#ifndef MAGCOLLECTION_H
#define MAGCOLLECTION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <vector>
#include <algorithm>

#include "magobject.h"

using namespace std;

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) (((x) < (y)) ? (y) : (x))
#endif

class CMagCollection {

protected:

    int m_iTotalLinearParameters;
    int m_iTotalNonLinearParameters;
    	
    vector<int> m_vCurrentLinearParams;
    vector<int> m_vCurrentNonLinearParams;	
    vector<CMagObject *> m_vObjects;

    vector<int>::iterator          m_itLinear;
    vector<int>::iterator          m_itNonLinear;
    vector<CMagObject *>::iterator m_itCurrentObject;  

    void destroy();

public:
    CMagCollection();
    ~CMagCollection();
    int AddMagObject(CMagObject * o);
    int ResetPosition();
    CMagObject * GetCurrentMagObject(int *nolin_start, int *lin_start);
    CMagObject * GetObjectAt(int n);
    CMagObject * GetObjectAt(char *name);
    int ReplaceObject(CMagObject *p_old, CMagObject *p_new);
    int GetTotalLinearParams();
    int GetTotalNonLinearParams();
    int GetTotalObjects() { return m_vObjects.size(); }
    void DumpObjects(FILE * dat);
    int GetSetVariables(int is_set, double * vars, int size, 
			int type, int linear_only, int is_errors);
    int ReScanAllObjects(int *total_linear, int * total_non_linear);

    // 0 - can do nothing
    // 1 - can compute field
    // 2 - can compute second derivatives    
    int GetAbilities(int *linear_ab, int *non_linear_ab);

    // min-max
    void InitMinMaxForObjects();
    void UpdateMinMaxForObjects();
    void DumpMinMaxForObjects(FILE * dat);
	int  RemoveObject(CMagObject * pObject);

    virtual void DumpToStream(std::ostream & str, std::string comment="#", 
							  char * close_str="/",char * delimit=",",
							  std::string offset = "  ");

	virtual void DumpToXMLStream(std::ostream & str);
};


#endif
