//
// @(#)fieldtype.h  1.1  misha-12sep103
//
// Copyright (c) 2003 Geometrics, Inc
//
// Created:  12sep103  by  misha@misha2.geometrics.com
// Version:  12sep103  by  misha@misha2.geometrics.com
//



#ifndef FIELDTYPE_H
#define FIELDTYPE_H

typedef enum _fieldtype {
  anyfield,
  totalmagfield,
  totalmaggradient,
  totalfinitemaggrad,
  absmaggradient,
  absfinitegradient
} FIELDTYPE;


#endif
