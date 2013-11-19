//
// @(#)profwrapper.h  1.1  Owner-18jul103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  18jul103  by  misha@MISHA
// Version:  18jul103  by  misha@MISHA
//

#ifndef PROFWRAPPER_H
#define PROFWRAPPER_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>

class CProfWrapper {

private:
  double m_lfDCShift;  // DC shift to be estimated
  double m_lfDrift;    // drift to be estimated
  double m_lfTime0;    // time of the profile start

public:
  CProfWrapper() { m_lfDCShift = 0.; m_lfDrift=0.; m_lfTime0=0.; }
  virtual ~CProfWrapper() {};

  // dealing with parameters

  void SetParameters(double shift, double drift) {
    m_lfDCShift = shift; m_lfDrift=drift;
  }

  double GetDCShift() { return m_lfDCShift; }
  double GetDrift()   { return m_lfDrift; }
  double GetTime0()   { return m_lfTime0; }
  void   SetTime0(double time0) { m_lfTime0 = time0; }

  // replace these functions with something meaningful

  virtual int GetTotalPoints()   { return 0;}
  virtual double GetX(int n)     { return 0.;}
  virtual double GetY(int n)     { return 0.;}
  virtual double GetTime(int n)  { return 0.;}
  virtual double GetAbsoluteTime(int n)  { return 0.;}
  virtual double GetField(int n) { return 0.;}
  virtual void   Correct()       { return; }

  virtual int  GetFineIndex(int n) { return 0;}
  virtual void SetFineMode(int mode) { };
  virtual void AddConstant(double b) { };
  virtual double GetDepth(int n) { return 0.; }
  virtual double  GetAltitude(int n) { return 0.;}
  virtual int GetLineName(std::string & sName) { sName = "no_line"; return 0;}
  virtual int GetChannel() { return 0; }
};


#endif

