//
// @(#)inversion.h  1.1  misha-08aug103
//
// Copyright (c) 2003 Geometrics
//
// Created:  08aug103  by  misha@misha2.geometrics.com
// Version:  08aug103  by  misha@misha2.geometrics.com
//

#ifndef INVERSION_H
#define INVERSION_H

#include "magdatacollection.h"
#include "magfunc.h"

#include <iostream>

class CInversion {

public:
  virtual CMagDataCollection * GetMagData()     { return 0; }
  virtual int GetDataVector(double *) { return 0; }
  virtual int GetWholeInversionMatrix(double * a, int packing = 1) { return 0; }
  virtual int GetWholeInversionMatrixRow(double * a, int n_pos, int size) { return 0; }
  virtual int GetAllVariables(double * vars, int size, int is_errors = 0) { return 0; }
  virtual int SetAllVariables(double * vars, int size, int is_errors = 0) { return 0; }
  virtual int ComputeField(double *T_synt, int size, int type=0) {return 0; }
  virtual double ComputeFieldAtPosition(int n_pos, int type=0) { return 0.; }

  virtual PROGRESS_FUNC GetProgressFunc() { return NULL; }
  virtual void SetProgressFunc(PROGRESS_FUNC pfunc) { }
  virtual int GetTotalAbilities() { return 0; }
  virtual ~CInversion() {};
  virtual void DumpToStream(std::ostream & str, std::string comment="#", char * close_str="/", 
 	           char * delimit=",", std::string offset = "  " ) {};
  virtual void DumpToXMLStream(std::ostream & str) {};
  virtual void SetInvID(std::string id) {};
  virtual void AddDoubleProperty(char * name, double prop) { return;}
  virtual double  GetDoubleProperty(char * name) { return 0.;}
  virtual void SetTheoryStdDevFlag(int flag) {};
  virtual int  GetTheoryStdDevFlag()         { return 0;};
};

#endif
