//
// @(#)magobject.h  1.1  misha-09may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  09may103  by  misha@MISHA
// Version:  09may103  by  misha@MISHA
//

#include <stdio.h>
#include <stdlib.h>

#ifndef MAGOBJECT_H
#define MAGOBJECT_H

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x, y) (((x) < (y)) ? (y) : (x))
#endif

#ifndef M_PI
# define M_PI     3.14159265358979324
#endif

#include "fieldtype.h"
#include <string>
#include <map>
#include <iostream>

#define OBJECT_TYPE_LEN 20
#define OBJECT_NAME_LEN 80


// different values in add_params fiels

// additonal offsets for observation point

#define X_OFFSET 0
#define Y_OFFSET 1
#define Z_OFFSET 2

// gradient direction for gradient computation

#define DIRECTION_X  3
#define DIRECTION_Y  4
#define DIRECTION_Z  5

// separation for total gradient computaion

#define FINITE_SEPARATION 6


#define MAX_FIELD_TYPES  4
#define MAX_FIELD_PARAMS 30

// Important note about  Magnetic Moment (J) units:
//
// Package is using internal J scaling, which is 1.
//
// Distance units are **meters**
//
// For forward modeling:
// 1. If J is in A*m^2,  multiply J by 100 before calling
//    the program (or multiply result by 100).
// 2. If J is in cgs, divide it by 10 before calling the program
//
// For inversion:
// 1. To obtain J in A*m^2, divide results by 100
// 2. To obtain J in cgs, multiply by 10
//
// To get approx. mass in kg, divide by 5*10e-2
//

#include <string>

class CMagObject;

typedef int (CMagObject::*GetDerivativesFunc)(FIELDTYPE type, double *add_params, double *, int, double *, double,
					      double, double, double, double, double, double, double, double,
					      int, int *, int *);

typedef double (CMagObject::*GetFieldFunc)(FIELDTYPE type, double * add_params,
					   double x_obs, double y_obs, double z_obs,
					   double x0,    double y0,    double z0,
					   double A,     double B,     double C,
					   int data_pos);


class CMagObject {

protected:


  int m_bIsValid;          // if object is valid
                           // 0 - no both Lin/nonlin 1 - lin only 2 -both
  int m_bInduced;          // magnetization type
  int m_nLinearParams;     // total number of linear (J) parameters
  int m_nNonLinearParams;  // total number of non-linear (geometrics) parameters

  double * m_lfJ;          // linear parameters (magnetization)
  double * m_lfPos;        // non-linear parameters (positions)

  double * m_lfStdDevJ;    // standard deviation for linear parameters estimation
  double * m_lfStdDevPos;  // standard deviation for non-linear parameters estimation

  double * m_lfJUp;        // upper limits for J;
  double * m_lfPosUp;      // upper limits for positions;

  double * m_lfJLow;        // upper limits for J;
  double * m_lfPosLow;      // upper limits for positions;

  double * m_lfJMax;        // max. possible value for J;
  double * m_lfPosMax;      // max. possible value for position;

  double * m_lfJMin;        // min. possible value for J;
  double * m_lfPosMin;      // min. possible value for position;

  char m_sObjectType[OBJECT_TYPE_LEN];  // type and name of the object
  char m_sObjectName[OBJECT_NAME_LEN];

  int * m_pJfixed;    // holds 0 for variable value and 1 for fixed - linear variables
  int * m_pPosfixed;  // holds 0 for variable value and 1 for fixed - non-linear variables

  int m_iLinearStart;      // where linear or non-linear parts start in common matrix.
  int m_iNonLinearStart;   // these are valid ONLY right after calling CMagCollection::GetCurrentMagObject

  int m_nLinearParamsReport;
  int m_nNonLinearParamsReport;

  int m_iObjectLinearAbilities;    // 0 - can do nothing
  int m_iObjectNonLinearAbilities; // 1 - can compute field
                                   // 2 - can compute second derivatives

  FIELDTYPE m_iFieldType;   // field type to compute only for this object;
                            // if "anyfield" than compute for all fields possible


  int m_iRefCounter;        // reference counter; tremove object only when it is 0

  std::string m_sID;

  // protected methods

  virtual void destroy();
  virtual void copy(CMagObject const &o);

public:
  CMagObject();
  CMagObject(int is_induced, int n_linear, int n_nonlinear, char *type=NULL, char *name=NULL);
  virtual ~CMagObject() { destroy(); }

  CMagObject (CMagObject const &o);                  // copy constructor;
  CMagObject const &operator=(CMagObject const &o); // overloaded assignment

  int IsValid() { return m_bIsValid==2; }

  FIELDTYPE GetFieldType() { return m_iFieldType; }
  void SetFieldType(FIELDTYPE field) { m_iFieldType = field; }
  void GetFieldTypeText(std::string & str);


  virtual void SetInduced(int is_induced) { m_bInduced = is_induced; }

  virtual int GetTotalLinearParams()    { return m_nLinearParamsReport;   }
  virtual int GetTotalNonLinearParams() { return m_nNonLinearParamsReport;}

  virtual int GetAllLinearParams()    { return m_nLinearParams;   }
  virtual int GetAllNonLinearParams() { return m_nNonLinearParams;}

  // params - input / output array
  // size   - size of "params"
  // type: 0 for values, 1 for standard errors

  virtual int GetLinearParams(double * params,    int size, int type = 0);
  virtual int GetNonLinearParams(double * params, int size, int type = 0);
  virtual int GetAllPositionParams(double * params, int size);

  virtual int SetLinearParams(double * params,    int size, int type = 0);
  virtual int SetNonLinearParams(double * params, int size, int type = 0);

  double GetNonLinearParam(int n_param, int type = 0);
  double GetLinearParam(int n_param, int type = 0);

  // x_obs, y_obs, z_obs - position of the observation point
  // x0, y0, z0          - position of the origin
  // A, B, C             - Earth's magnetic field direction coefficients

  virtual double ComputeField(FIELDTYPE type, double * add_params,  double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
			              int data_pos) {return 0.; }

  virtual int GetLinearDerivatives(FIELDTYPE type, double * add_params,  double * deriv, int size,
				   double * pos_deriv,
				   double x_obs, double y_obs, double z_obs,
				   double x0,    double y0,    double z0,
				   double A,     double B,     double C,
				   int data_pos,
				   int * n_params, int * start = NULL) { *n_params = 0; return 0;}

  virtual int GetNonLinearDerivatives(FIELDTYPE type, double * add_params, double * deriv, int size,
				      double * pos_deriv,
				      double x_obs, double y_obs, double z_obs,
				      double x0,    double y0,    double z0,
				      double A,     double B,     double C,
				      int data_pos,
				      int * n_params, int * start = NULL) { *n_params = 0; return 0;}

  void SetLinNonLinStarts(int lin_start, int nonlin_start) {
     m_iLinearStart = lin_start;
     m_iNonLinearStart = nonlin_start;
  }

  void SetID(std::string id) { m_sID = id; };
  const char * GetID() { return (m_sID.empty() ? NULL : m_sID.c_str()); }

  void UpdateReportedParams();

  // service: get / save name and type

  virtual int SetObjectType(char * type);
  virtual int SetObjectName(char * type);

  virtual int GetObjectType(char * type, int size);
  virtual int GetObjectName(char * type, int size);

  virtual void DumpParameters(FILE * dat);

  // 0 - can do nothing
  // 1 - can compute field
  // 2 - can compute  derivatives
  int GetObjectLinearAbilities() { return m_iObjectLinearAbilities; }
  int GetObjectNonLinearAbilities() { return m_iObjectNonLinearAbilities; }

  virtual void SetObjectLinearAbilities( int ab) { m_iObjectLinearAbilities = ab; }
  virtual void SetObjectNonLinearAbilities(int ab) { m_iObjectNonLinearAbilities = ab; }

  virtual int FixLinearParam(int n_param, int is_fixed = 0, double val = 0.);
  virtual int FixNonLinearParam(int n_param, int is_fixed = 0, double val= 0.);

  virtual int FixLinearParamInterval(int n_param, double up=1.e+38, double low=-1.e+38);
  virtual int FixNonLinearParamInterval(int n_param, double up=1.e+38, double low=-1.e+38);

  virtual int IsLinearParamFixed(int n_param);
  virtual int IsNonLinearParamFixed(int n_param);

  // for correction object

  virtual int CanComputeField() { return 0; }

 // for error forward estimation

 virtual void InitMinMaxParams();
 virtual void UpdateMinMaxParams();
 virtual void DumpMinMaxParams(FILE * out);

 // reference counter

 void IncrementRefCounter() { m_iRefCounter++; }
 void DecrementRefCounter() { m_iRefCounter--; }

 int GetRefCounter() { return m_iRefCounter; }

 // get header information
 virtual int GetInfoHeader(char * str, int size) { return 0; }
 virtual int GetDataHeader(char * str, int size) { return 0; }

 virtual void DumpToStream(std::ostream & str, std::string comment="#",
						   char * close_str="/", char * delimit=",",
						   std::string offset = "  ") {};

 virtual void DumpToXMLStream(std::ostream & str) {};
};





#endif






