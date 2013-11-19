//
// @(#)maginv.h  1.1  misha-14may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  14may103  by  misha@misha.local
// Version:  14may103  by  misha@misha.local
//

#ifndef MAGINV_H
#define MAGINV_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#include "magdatacollection.h"
#include "magcollection.h"
#include "magfunc.h"
#include "inversion.h"

#ifndef M_PI
# define M_PI     3.14159265358979324
#endif



class CMagInv : public CInversion {

protected:

    // Earth's fieled part
    double I,  D;                 //  inclination and dec. of main field
    double Az;                    //  azimuth of X axes
    double A,B,C;                 //  unit vector comp. for main mag. field
    double I_rad, D_rad;          //  inclination and declination in radian

    // progress and general part
    PROGRESS_FUNC m_ProgressFunc;
    void *data;                   //  any data

    // synthetic / observed difference

    double m_lfLinearDiscrepancy;    // RMS observed / synthetic field
    double m_lfLinearMaxDiff;        // max. difference  observed / synthetic after linear phase
    double m_lfNonLinearDiscrepancy; // RMS observed / synthetic field after non-linear phase
    double m_lfNonLinearMaxDiff;     // max. difference  observed / synthetic after linear phase


    // Linear inversion part
    //double eps;                   // parameter to discard small singular numbers for z29rls
    //double epsm;                  // parameter to check matrix V (if less them epsm, V[i,j]=0) (z29rls);
    //double dm;                    // parameter to find regularization (z29rls)
    //double am;                    // estimated reg. parameter after z29rls
    //double linear_fit;            // fit after linear inversion


   // non-linear inversion part
    //int n_iter;                   // non-linear iterations made
    //int iopt;                     // how to compute jacobian (1 - num, 2, 3 - analytical)
    //int maxfev;                   // max. non-linear iterations
    //double fit;                   // RMS between observed and modeled fields
    //double * std_err;             // standatd errors computed based of fit and cov. matrix
    //double max_diff;              // maximum difference between observed and computed fields
    //int info;                     // return from dnls1


    // Note! data, collection and background are not destroyed aftre program terminates.
    // if you need to replace either object make sure previous one is destroyed!

    // data

    CMagDataCollection * m_pData;          // data (measurements)

    // models (objects)

    CMagCollection * m_pCollection; // objects to be used

	std::string m_sInvID;

	std::vector<std::string> m_vPropNames;
	std::vector<double>      m_sPropVals;	

	// only estimate std. dev theoretically

	int m_iTheoryStdDevFlag;

public:

    // constructor / destructor

    CMagInv();
    virtual ~CMagInv();

    // data / collection / background

    CMagDataCollection * GetMagData() { return m_pData;}
    CMagCollection * GetMagCollection() { return m_pCollection;}

    void SetMagData(CMagDataCollection * pData) {  m_pData = pData;}
    void SetMagCollection(CMagCollection * pCollection) { m_pCollection = pCollection; }

    void SetEarthFieldDirections(double Incl, double Decl, double Azimuth);
    void SetEarthCoeff();
    void GetEarthCoeff(double *Aout, double *Bout, double *Cout, double *Azout);
	void SetEarthCosines(double nx, double ny, double nz);

    // obtain matrixes for inversion

    int GetTotalParameters(int *n_data, int *n_linear_vars, int *n_nonlin_vars);
    int GetLinearInversionMatrix(double * a, double *y, int packing = 1);
    int GetNonLinearInversionMatrix(double * a, double *y, int packing = 1);
    int GetWholeInversionMatrix(double * a, int packing = 1);
    int GetWholeInversionMatrixRow(double * a, int n_pos, int size);

    int GetGeneralInversionMatrix(double * a, GetDerivativesFunc df, int m_vars,
				  int lin_offset, int packing);

    int GetGeneralInversionMatrixRow(int n_pos, double * a, GetDerivativesFunc df, int lin_offset);

    int GetDataVector(double *y);

    int GetLinearVariables(double * vars, int size, int linear_only = 1, int is_errors = 0);
    int SetLinearVariables(double * vars, int size, int linear_only = 1, int is_errors = 0);

    int GetAllVariables(double * vars, int size, int is_errors = 0);
    int SetAllVariables(double * vars, int size, int is_errors = 0);

    int GetNonLinearVariables(double * vars, int size, int is_errors = 0);
    int SetNonLinearVariables(double * vars, int size, int is_errors = 0);

    // set / get progress function

    PROGRESS_FUNC GetProgressFunc() { return m_ProgressFunc; }
    void SetProgressFunc(PROGRESS_FUNC pfunc) { m_ProgressFunc = pfunc; }


    int ComputeField(double *T_synt, int size, int type=0);
    int ComputeField(int type=0, int replace=0);

    double ComputeFieldAtPosition(int n_pos, int type=0);

    int EstimateLinearParameters(LINEAR_ESTIMATOR_FUNC func, double * fit, int packing=1, void * params = NULL,
				 char * message = NULL, int store_results = 1);
    int EstimateAllParameters(NONLINEAR_ESTIMATOR_FUNC func, void * func_param = NULL,
                              char * message=NULL, int store_results = 1);

    int EstimateNonLinearParametersRange(LINEAR_ESTIMATOR_FUNC func, double fit, double step, double max_step,
                                         int packing = 1, void * params = NULL);

    int EstimateNonLinearParametersRange2(NONLINEAR_ESTIMATOR_FUNC func, double fit, double step,
                                         double max_step, void * params = NULL);

	int EstimateStandardDeviations();

    double GetNonLinearDiscrepancy() { return  m_lfNonLinearDiscrepancy;}
    double GetNonLinearMaxDiff()     { return  m_lfNonLinearMaxDiff;}
    double GetLinearDiscrepancy()    { return  m_lfLinearDiscrepancy;}
    double GetLinearMaxDiff()        { return  m_lfLinearMaxDiff;}

    void DumpObjects(FILE * dat);

    virtual int GetTotalAbilities();

	virtual void DumpToStream(std::ostream & str,   std::string comment="#", 
		                      char * close_str="/", char * delimit=",",
							  std::string offset = "  ");

	virtual void DumpToXMLStream(std::ostream & str);

	virtual void SetInvID(std::string id) { m_sInvID = id; }

    virtual void    AddDoubleProperty(char * name, double prop);
    virtual double  GetDoubleProperty(char * name);

	 int GetNearestDataPoint(double x0, double y0, double * d, double *x, double *y, double *dtime,
		 	              std::string & line_name) {
						  if(m_pData) return m_pData->GetNearestDataPoint(x0, y0, d, x, y, dtime, line_name);
						  return -1;
	 }

	virtual void SetTheoryStdDevFlag(int flag) { m_iTheoryStdDevFlag = flag; }
    virtual int  GetTheoryStdDevFlag()         { return m_iTheoryStdDevFlag; }

};


#endif







