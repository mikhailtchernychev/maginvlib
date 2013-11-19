//
// @(#)magfunc.cpp  1.1  misha-16may103
//
// Copyright (c) 2003 of Geometrics.
// All rights reserved.
//
// Created:  16may103  by  misha@misha.local
// Version:  16may103  by  misha@misha.local
//



// functons used in mag processing

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "grd7.h"
#include "magfunc.h"


// ------------------------ functions -------------------------------

// n - number of rows
// m - number of colums
// a - matrix m x n, columns packed
// b - vector m
// c = vector n

void MatMultVect(double *a, double *b, int n, int m, double *c) {

  for(int i=0; i<n; i++) {
    c[i] = 0.;
    for(int j=0; j<m; j++) {
      c[i] += (a[j*n+i]*b[j]);
    }
  }
}


// matrix is packed in rows

void MatMultVectRow(double *a, double *b, int n, int m, double *c) {

  for(int i=0; i<n; i++) {
    c[i] = 0.;
    for(int j=0; j<m; j++) {
      c[i] += (a[j+ m*i]*b[j]);
    }
  }
}



// compute discrepancy(fit) between two vectors

void compute_fit(int n, double * y, double * synt, double * fit, double * max_fit) {

   double d     = 0.;
   double d_max = 0.;
   double dtmp;

   for(int i=0; i<n; i++) {
	  dtmp =  y[i] - synt[i];
	  d +=  (dtmp*dtmp);
	  d_max = MAX(d_max, fabs(dtmp));
	}

	*fit     = sqrt(d/n);
	*max_fit = d_max;

}


void * LoadSurfer7Grid(char * file_name, long * nx, long *ny,
			 double * x1, double *y1,
			 double *dx, double *dy,
			 double *zmin, double *zmax, int is_float) {

  double * z = NULL;
  float * z1 = NULL;
  FILE * dat = fopen(file_name, "rb");
  if(!dat) return NULL;

  long lver;
  SURFER7_TAG tag;
  //SURFER7_HEADER header;
  SURFER7_GRID grid;

  fread(&tag, sizeof(SURFER7_TAG), 1, dat);

  if(tag.id != SURFER7_HEAD) {
    fclose(dat);
    return NULL;
  }

  fread(&lver, sizeof(long), 1, dat);

  while(1) {
    if(1!=fread(&tag, sizeof(SURFER7_TAG), 1, dat)) break;

    switch(tag.id) {

    case SURFER7_GRD: {
      fread(&grid, sizeof(SURFER7_GRID), 1, dat); // read the header
      *ny   = grid.nRow;
      *nx   = grid.nCol;
      *x1   = grid.xLL;
      *y1   = grid.yLL;
      *dx   = grid.xSize;
      *dy   = grid.ySize;
      *zmin = grid.zMin;
      *zmax = grid.zMax;
      fread(&tag, sizeof(SURFER7_TAG), 1, dat); // read data tag;
      if(!is_float) {
	z = (double *) malloc(tag.Size);
	if(!z) goto end_func;
	fread(z, tag.Size, 1, dat);
      }
      else {
	int n =  *nx * *ny;
	z1 = (float *) malloc(sizeof(float)*n);
	if(!z1) goto end_func;
	double dtmp;
	int index;
	for(int i=0; i<*ny; i++) {
	  for(int j=0; j<*nx; j++) {
	  fread(&dtmp, sizeof(double), 1, dat);
	  index = ((*ny)-i-1)*(*nx)+j;
	  z1[index] = (float)dtmp;
	  }
	}
	z = (double *)z1;
      }

      break;
    }

    default: // skip all other tags
      fseek(dat, tag.Size, SEEK_CUR);
      break;
    }
  }



end_func:
  fclose(dat);
  return (void *)z;

}



/* ------------------------------- time functions --------------------------- */

void parse_time(char *time_txt, double * lfTime)
{
 char *ptr1 =  time_txt;
 int i,j, max_size;

 max_size = strlen(time_txt);
 lfTime[0] = lfTime[1] = lfTime[2] = 0.;

 char * buf = new char[strlen(time_txt)+1];
 i=0;

 for(j=0; j<max_size; j++)
    {
     ptr1 = strpbrk(ptr1, "0123456789");
     if ( ptr1 == NULL) break;

     strcpy(buf, ptr1);

     if(i<2) {
       strcpy(buf, ptr1);
       buf[strspn(buf, "0123456789")] = 0;
       lfTime[i++] = atof(buf);
     }
     else { // sec. field can have comma
       char * p=buf;
       while(*p) {
	 if(*p==',') *p='.';
	 p++;
       }
       lfTime[i++] = atof(buf);
       break;
     }

     ptr1 = ptr1+strspn(ptr1, "0123456789");
    }

 delete [] buf;
}


double date2sec(char * date_txt, char * time_txt, int format){

  double lfTime[3];
  long year, month, day, hh, mm,  ss;
  unsigned long sec;
  double dsec;

  parse_time(time_txt, lfTime);

  hh   = (int)lfTime[0];
  mm   = (int)lfTime[1];
  dsec = lfTime[2];

  parse_time(date_txt, lfTime);

  year = int(lfTime[2]);

  // check if year is first
  if(lfTime[0] > 1000) 	format = 3;

  if(format <=2) {
    if(lfTime[0] > 12.) format = 2;
    else if(lfTime[1] > 12.) format = 1;
  }

  switch(format) {
      case 2:
        month = (int)lfTime[1];
        day   = (int)lfTime[0];
       break;
       case 3:
        month = (int)lfTime[1];
        day   = (int)lfTime[2];
        year  = (int)lfTime[0];
        break;
       default:
        month = (int)lfTime[0];
        day   = (int)lfTime[1];
  }


  //if(format==1) sscanf(date_txt, "%ld/%ld/%ld", &month, &day, &year);
  //       else   sscanf(date_txt, "%2.2ld/%2.2ld/%2.2ld", &day, &month, &year);

  //sscanf(time_txt, "%ld:%ld:%lf", &hh, &mm, &dsec);

  // check

  if( !((month>0 && month < 13) &&
      (day  >0 && day <   32) &&
      (hh  >=0 && hh < 25)    &&
      (mm  >=0 && mm <60)     /*&&
	  (dsec < 60.)*/ )) return 0.;



  ss = (int) dsec;

  if(year<1000) {
    if(year>80) year +=1900;
    else        year +=2000;
  }

  year = year-1980;

  month = month-3;

  if ( month<0 ) {month+=12; year-=1;}
  day=day+(month*306+4)/10+((year*1461)>>2);
  sec=ss+60*(mm+60*(hh+ 24*(day)));
  return (double) sec+ (dsec - (double)ss);;
}

int sec2date (double dsec, char * date_txt, char * time_txt, int format ){

  if(!format || dsec <=0.) {
    if(date_txt) date_txt[0] = 0;
    if(time_txt) time_txt[0] = 0;
    return 0;
  }

  unsigned long sec;

  sec=(unsigned long)dsec;

  // prevent 1000 ms
  double dms= (dsec-(double)sec);
  if(int(1000.*(dsec-sec)+0.5) >= 1000) {
    dsec = sec+1;
    sec++;
  }

  long year, month, day, hh, mm, ss;
  unsigned long  Temp,rest,quot;

  ss  =sec%60;
  mm  =(sec/60)%60;
  hh=(sec/3600)%24;
  Temp      =(sec/3600)/24*4-1;
  rest=Temp%1461;
  quot=Temp/1461;
  if (rest<0){ rest+=1461 ; quot-=1; }
  year =quot+1980;
  Temp      =(rest/4+1)*4+(rest/4+1)-3;
  month=Temp/153+3;
  day  =Temp%153/5+1;
  if(month>12){ year+=1; month-=12;}

  dsec = ss + (dsec-(double)sec);

  if(year > 1900 && year < 2000) year-=1900;
  if(year >=2000) year -= 2000;

  if(date_txt) {
    if(format == 1) sprintf(date_txt,"%2.2ld/%2.2ld/%2.2ld", month, day, year);
            else  sprintf(date_txt,"%2.2ld/%2.2ld/%2.2ld", day, month, year);
  }

  if(time_txt)
    sprintf(time_txt,"%2.2ld:%2.2ld:%2.2ld.%3.3ld", hh, mm, ss, int(1000.*(dsec-ss)+0.5));

  return 0;
}



/* -------------------- compute median value ---------------------- */

void SortSumData(double *dataset,int numrow)
{
   int j;
   int k;
   double TempX;
   char abort;

   if ( numrow <= 1 )
      return;
   for ( j = 0; j <= numrow - 1; ++j ) {
      abort = 0;
      TempX = dataset[j];
      k = j - 1;
      while ( ! (abort) && (k >= 0) ) {
         if ( TempX < dataset[k] ) {
            dataset[k + 1] = dataset[k];
            k = k - 1;
         }
         else {
            abort = 1;
         }
      }
      dataset[k + 1] = TempX;
   }

}

void CalculateMode(double *dataset, int numrow, double *mode)
{
   char Even;
   int Num;

   if(numrow==1) { *mode = *dataset; return; }

   Even = 1;
   if ( numrow % 2 != 0 ) {
      Even = 0;
   }
  SortSumData(dataset,numrow);
   if ( Even ) {
      Num = numrow / 2 - 1;
      (*mode) = (dataset[Num] + dataset[Num + 1]) / 2;
   }
   else {
      Num = numrow / 2 - 1;
      (*mode) = dataset[Num];
   }
}


