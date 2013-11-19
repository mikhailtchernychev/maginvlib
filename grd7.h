//
// @(#)grd7.h  1.1  misha-21apr103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  21apr103  by  misha@MISHA
// Version:  21apr103  by  misha@MISHA
//

#ifndef SURFER_7_GRID_H
#define SURFER_7_GRID_H

#define SURFER7_HEAD    0x42525344  // "DSRB"
#define SURFER7_GRD     0x44495247  // "GRID"
#define SURFER7_DAT     0x41544244  // "DBTA"
#define SURFER7_FAULT   0x49544c46  // "FLTI"

//#pragma pack(push,1)
//#pragma pack(1)

#pragma pack(1) // gcc


typedef struct _surfer7_tag {
  long id;   //  The type of data in the following section.  
             //  See the next table for a list of valid values.
  long Size; //  The number of bytes in the section (not including this tag).
             //  Skipping this many bytes after reading the tag will align 
             // the file pointer on the next tag.
} SURFER7_TAG;


typedef struct _surfer7_header {
  long version;  // Version number of the file format.  
                 // Currently must be set to 1. 
} SURFER7_HEADER;

typedef struct _surfer_7_grid {

  long nRow;         // number of rows in the grid
  long nCol;         // number of columns in the grid
  double xLL;        // X coordinate of the lower left corner of the grid
  double yLL;        // Y coordinate of the lower left corner of the grid
  double xSize;      // spacing between adjacent nodes in the X direction (between columns)
  double ySize;	     // spacing between adjacent nodes in the Y direction (between rows)
  double zMin;       // minimum Z value within the grid
  double zMax;       // maximum Z value within the grid
  double Rotation;   // not currently used
  double BlankValue; // nodes are blanked if  greater or equal to this value

} SURFER7_GRID;


typedef struct _surfer7_fault_info {
  long nTraces;    // 	number of fault traces (polylines)
  long nVertices;  //	total number of vertices in all the traces
  // data section: variable-sized data block consisting of an array of
  // Trace structures immediately followed by the array of vertices
} SURFER7_FAULT_INFO;

// A Data section containing an array of Trace structures and an 
// array of Vertex structures must immediately follow a Fault Info section.
// The number of Trace structures in the array is nTraces, and the number
//  of Vertex structures is nVertices.


typedef struct _surfer7_trace {
  long iFirst; // 0-based index into the vertex array for the first vertex of this trace
  long nPts;   // number of vertices in this trace
} SURFER7_TRACE;


typedef struct _surfer7_vertex {
  double x; //  X coordinate of the vertex
  double y; //	Y coordinate of the vertex
} SURFER7_VERTEX;


//#pragma pack(pop)

#pragma pack() // gcc


void * LoadSurfer7Grid(char * file_name, long * nx, long *ny, 
			 double * x1, double *y1, 
			 double *dx, double *dy,
		         double *zmin, double *zmax, int is_float=0);



#endif









