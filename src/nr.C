/* Program file nr.C - 
* Vector and matrix allocation and deallocation routines 
* from Numerical Recipes
* Last change: KL 01.08.02
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 3.0, August 2015
* 
============================================
* allocates vectors/matrices with subscript range v[nl..nh]
*/

/* #include"gcc.h" */
//#include "/usr/include/sys/types.h"
#include <sys/types.h>
#include<complex>
#include<stdio.h>
#include<stdlib.h>
#include<cnr.h>
// #include<malloc.h>
#include<RLexact.h>

void nrerror(const char*);


void nrerror(const char *error_text){
  
  printf("Numerical Recipes run-time error ...\n");
  printf("%s\n",error_text);
  printf("... now exiting to system.\n");
  exit(1);
}

long long *ivector(long long nl, long long nh){
  long long *di;

  di=(long long*) malloc((unsigned) (nh-nl+1)*sizeof(long long));
  if(di==NULL){
    fprintf(stderr,"Allocation error in ivector()\n");
    exit(-1);
  }
  di-=nl;
  return(di);
}

void freeivector(long long *di, long long nl, long long nh){
 free(di+nl);
}

long long **imatrix(long long ml, long long mh, long long nl, long long nh){
  long long **im;
  long long i;
  
  im=(long long**) malloc((unsigned) (mh-ml+1)*sizeof(long long*));
  if(im==NULL) nrerror("Allocation error 1 in imatrix()\n");
  im-=ml;

  for(i=ml;i<=mh;i++){
    im[i]=(long long*) malloc((unsigned) (nh-nl+1)*sizeof(long long));
    if(im[i]==NULL) nrerror("Allocation error 2 in imatrix()\n");
    im[i] -= nl;
  }
  
  return(im);
}

void freeimatrix(long long **im, long long ml, long long mh, long long nl, long long nh){
  long long i;

  for(i=ml;i<=mh;i++)
    free(im[i]+nl);

  free(im+ml);
}

unsigned long long *lvector(long long nl, long long nh){
  unsigned long long *di;

  di=(unsigned long long*) malloc((unsigned) (nh-nl+1)*sizeof(unsigned long long));
  if(di==NULL){
    fprintf(stderr,"Allocation error in lvector()\n");
    exit(-1);
  }
  di-=nl;
  return(di);
}

void freelvector(unsigned long long *di, long long nl, long long nh){
 free(di+nl);

}

unsigned long long **lmatrix(long long ml, long long mh, long long nl, long long nh)
{
  unsigned long long **dm;
  long long i;
  
  dm=(unsigned long long**) malloc((unsigned) (mh-ml+1)*sizeof(unsigned long long*));
  if(dm==NULL) nrerror("Allocation error 1 in imatrix()\n");
  dm-=ml;

  for(i=ml;i<=mh;i++){
    dm[i]=(unsigned long long*) malloc((unsigned) (nh-nl+1)*sizeof(unsigned long long));
    if(dm[i]==NULL) nrerror("Allocation error 2 in imatrix()\n");
    dm[i] -= nl;
  }
  
  return(dm);
}

void freelmatrix(unsigned long long **dm, long long ml, long long mh, long long nl, long long nh){
  long long i;

  for(i=ml;i<=mh;i++)
    free(dm[i]+nl);

  free(dm+ml);
}

float *vector(long long nl, long long nh){
  float *dv;

  dv=(float*) malloc((unsigned) (nh-nl+1)*sizeof(float));
  if(dv==NULL) {
    fprintf(stderr,"nl: %lld  nh: %lld\n",nl,nh);
    fprintf(stderr,"Allocation error in vector()\n");
    exit(-1);
  }
  dv-=nl;
  return(dv);
}

void freevector(float *dv, long long nl, long long nh){
  free(dv+nl);
}

float **matrix(long long ml, long long mh, long long nl, long long nh)
{
  float **dm;
  long long i;

  dm=(float**) malloc((unsigned) (mh-ml+1)*sizeof(float*));
  if(!dm) nrerror("Allocation error 1 in matrix()\n");
  dm -= ml;

  for(i=ml;i<=mh;i++){
    dm[i]=(float*) malloc((unsigned) (nh-nl+1)*sizeof(float));
    if(!dm[i])  nrerror("Allocation error 2 in matrix\n");
    dm[i] -= nl;
  }
  
  return(dm);
}

void freematrix(float **dm, long long ml, long long mh, long long nl, long long nh){
  long long i;

  for(i=ml;i<=mh;i++)
    free(dm[i]+nl);

  free(dm+ml);
}

double *dvector(long long nl, long long nh){
  double *dv;

  dv=(double*) malloc((unsigned) (nh-nl+1)*sizeof(double));
  if(dv==NULL) {
    fprintf(stderr,"nl: %lld  nh: %lld\n",nl,nh);
    fprintf(stderr,"Allocation error in dvector()\n");
    exit(-1);
  }
  dv-=nl; //move pointer, so that i.e. dv[nl] accessing the first element of dvec.
  return(dv);
}

void freedvector(double *dv, long long nl, long long nh){
  free(dv+nl); //pointers must be freed at the start of their actual memory location
}

double **dmatrix(long long ml, long long mh, long long nl, long long nh)
{
  double **dm;
  long long i;
  
  dm=(double**) malloc((unsigned) (mh-ml+1)*sizeof(double*));
  if(dm==NULL) nrerror("Allocation error 1 in dmatrix()\n");
  dm -= ml;

  for(i=ml;i<=mh;i++){
    dm[i]=(double*) malloc((unsigned) (nh-nl+1)*sizeof(double));
    if(dm[i]==NULL) nrerror("Allocation error 2 in dmatrix()");
    dm[i] -= nl;
  }
  
  return(dm);
}

void freedmatrix(double **dm, long long ml, long long mh, long long nl, long long nh){
  long long i;

  for(i=ml;i<=mh;i++)
    free(dm[i]+nl);

  free(dm+ml);
}

komplex *cvector(long long nl, long long nh){ //(kvector = cvector as per RLexact.h)
  komplex *cv;

  cv=(komplex*) malloc((unsigned) (nh-nl+1)*sizeof(komplex));
  if(cv==NULL) nrerror("Allocation error in cvector()\n");
  cv-=nl;
  return(cv);
}

void freecvector(komplex *cv, long long nl, long long nh){
  free(cv+nl);
}

komplex **cmatrix(long long ml, long long mh, long long nl, long long nh) // cmatrix = kmatrix
{
  komplex **cm;
  long long i;
  
  cm=(komplex**) malloc((unsigned) (mh-ml+1)*sizeof(komplex*));
  if(cm==NULL) nrerror("Allocation error 1 in cmatrix()\n");
  cm -= ml;

  for(i=ml;i<=mh;i++){
    cm[i]=(komplex*) malloc((unsigned) (nh-nl+1)*sizeof(komplex));
    if(cm[i]==NULL) nrerror("Allocation error 2 in cmatrix()");
    cm[i] -= nl;
  }
  
  return(cm);
}

void freecmatrix(komplex **cm, long long ml, long long mh, long long nl, long long nh){
  long long i;

  for(i=ml;i<=mh;i++)
    free(cm[i]+nl);

  free(cm+ml);
}











