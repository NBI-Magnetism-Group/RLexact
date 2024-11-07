/* Program file RLhamil.C - 
* Applying the Hamilton operator to a given state
* Last change: AR 25.09.12
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 2.4, August 2002
* 
============================================
*/

//#include "/usr/include/sys/types.h"
#include <sys/types.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <nr.h>
#include <cnr.h>
#include <RLexact.h>

/* Functions defined in this file */
void Eigenvector_test(int *, komplex *, komplex *);
void FillHamilton(int *, komplex **);
void Hamilton(komplex *, komplex *, int *);
double HamDiag();
void Hamil2(int *, komplex, komplex *);
void matrixelement(komplex, int *, komplex, komplex *);

/* Functions defined elsewhere */
extern int LookUpU(unsigned long);
extern int Count(unsigned long);
void fatalerror(char *, int);
unsigned long FindUnique(unsigned long, int *);   
extern void WriteGSEnergy(komplex);
extern void WriteState(char*, komplex *);
extern void time_stamp(time_t*, int, char*);

/* Global variables defined in RLexact.c */
extern unsigned long unique[];
extern int    Nunique, Nuniq_k, Nsym, Nspins;
extern int    Nocc[], Nsymvalue[], uniq_k[];
extern int    hamil_coup[NCOUP][2];
extern double Jzz[],Jxy[],Janis[];
extern double Jdip[],geom_13[],r_vector[NCOUP][3];
#ifdef M_SYM
extern int    m;
#else
extern double h,field[3];
#endif /* M_SYM */
extern double sine[],cosine[],sqroot[];
extern int    Nunique,Ncoup;

/* Regional variables in this file */
unsigned long bitmap, new_state;
int           n_2, u_occ;
unsigned long index1, index2;
komplex this_;


#ifdef TEST_EIG
void Eigenvector_test(int k[NSYM], komplex* evec, komplex* tmp)
{
  int i;
  komplex product;

  printf("Entering Eigenvector_test \n");
  Hamilton(evec,tmp,k);
  product = zero; 
  printf("Middle of Eigenvector_test \n");
  for(i=0;i<Nunique;i++) 
    product += evec[i]*conj(tmp[i]);
  WriteGSEnergy(product);

#ifdef WRITE_STATES
  WriteState("Eigenvector:",evec);
  WriteState("H|e>:",tmp);
#endif  /* WRITE_STATES */

  return;
}
#endif   /* TEST_EIG */


void FillHamilton(int k[], komplex **hamilton)
/* Fills the Hamiltonian matrix, used only with MATRIX */
{
  double diag;
  komplex *next;
  int i,j;
  time_t time_fill;

  time_stamp(&time_fill,START,"filling hamiltonian");
  next = kvector(0,Nunique-1);

  /* Find smallest larger power of 2 */
  for (n_2=1; n_2<Nuniq_k; n_2=n_2<<1);

  /* Fill matrix vector by vector */
  for(i=0;i<Nuniq_k;i++) 
  {    
    index1 = uniq_k[i];
    /* Reset result vector */
    for(index2=0;index2<Nunique;index2++)
      next[index2]=zero;
    bitmap = unique[index1];
    diag = 0;
    if(u_occ=Nocc[index1]) 
    {
      diag=HamDiag();     
#ifdef TEST_HAMILTON
      printf(" Diagonal-element: %g\n",diag);
#endif /* TEST_HAMILTON */
      Hamil2(k, one, next);
    }
    else
    {
      fatalerror("Non-allowed unique encountered in FillHamilton",index1);
    }  /* if Nocc[] */ 
    /* Transfer result vector to matrix */ 
    for(j=0;j<Nuniq_k;j++)
    {
      index2 = uniq_k[j];
      hamilton[i+1][j+1]=next[index2];
    }
    hamilton[i+1][i+1] += diag;

  }  /* for (i=0.. */

  freekvector(next,0,Nunique-1);
  time_stamp(&time_fill,STOP," ");

#ifdef TEST_FILLHAM
  printf(" FillHamilton() returning matrix: \n");
  for(i=0;i<Nuniq_k;i++) 
  {
   for(j=0;j<Nuniq_k;j++) 
    printf(" ( %lg + %lg i)", real(hamilton[i+1][j+1]),imag(hamilton[i+1][j+1]));
   printf("\n");
  }
#endif /* TEST_FILLHAM */

  return;
}

/* Hamilton applies the Hamilton operator to one state vector. */
/* Mostly used with LANCZOS algorithm */
void Hamilton(komplex *this_v, komplex *next_v, int k[])
{
  double diag;

#ifdef TEST_HAMILTON
  printf(" Hamilton() called with vector: \n");
  for(index1=0;index1<Nunique;index1++) 
    printf(" re: %g, im: %g\n",real(this_v[index1]),imag(this_v[index1]));
#endif
  /* Find smallest larger power of 2 */
  for (n_2=1; n_2<Nunique; n_2=n_2<<1);
  /* Initialize for all unique */
  for(index1=0;index1<Nunique;index1++) 
    next_v[index1] = zero;  // *** index added h***
  printf("test_ham 1 \n");
  for(index1=0;index1<Nunique;index1++) 
  {    
    this_ = this_v[index1];
    bitmap = unique[index1];
    diag = 0;
    printf("test_ham 2 \n");
    if(u_occ=Nocc[index1] && (this_ != zero)) 
    {
      diag=HamDiag();   
      printf("test_ham 3: index1= %i \n",index1);
      next_v[index1] += diag*this_; /* for off-diagonals this 
                        update is done in matrix_element, called by Hamil2() */
#ifdef TEST_HAMILTON
      printf(" Diagonal-element: %g\n",diag);
#endif /* TEST_HAMILTON */
      Hamil2(k, this_, next_v);
    }  /* if Nocc[] ... */
  }  /* for (index1=0.. */
#ifdef TEST_HAMILTON
  printf(" HAMILTON() returning vector: \n");
  for(index1=0;index1<Nunique;index1++) 
    printf(" ( %g + %g i)\n",real(next_v[index1]),imag(next_v[index1]));
#endif /* TEST_HAMILTON */

  return;
}

double HamDiag() 
/* Calculates the diagonal value of the Hamiltonian */
{
  double sz=0, diagonal=0;
  int j;
  unsigned long s0,s1,mask0, mask1;

#ifdef TEST_HAMDIAG
  printf(" bitmap: %ld ",bitmap);
#endif
  /* Field Sz term */
#ifndef M_SYM
  sz = Count(bitmap)-Nspins/2;
#ifdef TEST_HAMDIAG
  printf(" Sz = %g ",sz);
#endif /* TEST_HAMDIAG */
#endif /* M_SYM */

  /* Sz1 Sz2 term: run through spin pairs */
  for (j=0;j<Ncoup;j++) 
  {  
    mask0 = ((unsigned long) 1)<<hamil_coup[j][0];
    mask1 = ((unsigned long) 1)<<hamil_coup[j][1];
    s0 = (bitmap & mask0)!=0;
    s1 = (bitmap & mask1)!=0;
#ifdef TEST_HAMDIAG
    printf(" coup: (%d,%d). s0, s1: %d, %d ",
	   hamil_coup[j][0],hamil_coup[j][1],s0,s1);
#endif

    if( (s0+s1)==1 )       /* Spins are of opposite sign */
#ifdef DIPOLE
      diagonal -= Jzz[j]-Jdip[j]*geom_13[j];
#else
      diagonal -= Jzz[j];
#endif /* DIPOLE */
    else
#ifdef DIPOLE
      diagonal += Jzz[j]+Jdip[j]*geom_13[j];
#else
      diagonal += Jzz[j];
#endif /* DIPOLE */
  }

#ifdef M_SYM
  return diagonal/4.0;
#else
  return diagonal/4.0-h*sz;  /* The field is always in the Z direction */
#endif /* M_SYM */
}

/* Hamil2() deals with two-spin interactions */
void Hamil2(int k[], komplex coof, komplex *next)
{
 int j;
 unsigned long mask0,mask1,s0,s1;
 double sz;

 /* Run through all couplings */
 for(j=0;j<Ncoup;j++) 
 {  
  mask0 = ((unsigned long) 1)<<hamil_coup[j][0];
  mask1 = ((unsigned long) 1)<<hamil_coup[j][1];
  s0 = bitmap & mask0;
  s1 = bitmap & mask1;
#ifdef TEST_HAM2
  printf(" Coupling between %d and %d, (s0,s1)= %d %d \n"
          ,hamil_coup[j][0],hamil_coup[j][1],s0,s1);
#endif 
  if (s0==0)  /* case s0 down; S+.. terms (and Sz.. for DIPOLE) */ 
  {
#ifdef DIPOLE  /* S+Sz */
   if (s1==0) sz=-0.5; else sz=0.5;
   new_state = (bitmap | mask0);
   matrixelement(((-1.5*sz*Jdip[j]*r_vector[j][Z]*r_vector[j][X]) + 
    I*(1.5*sz*Jdip[j]*r_vector[j][Z]*r_vector[j][Y])),k,coof,next);
#endif
   if (s1==0) /* down down: S+S+ term */ 
   {
#ifndef M_SYM
    new_state = (bitmap|mask0|mask1);
    matrixelement(komplex (Janis[j]/2.0,0.0),k,coof,next);
#ifdef DIPOLE   /* S+S+ and SzS+ */
    matrixelement((-0.75*Jdip[j]*(SQR(r_vector[j][X])-SQR(r_vector[j][Y]))) + 
     I*(1.5*Jdip[j]*r_vector[j][X]*r_vector[j][Y]),k,coof,next);
    new_state = (bitmap | mask1);
    matrixelement((0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][X]) +
     I*(-0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][Y]),k,coof,next);
#endif
#endif   /* not M_SYM */
   }
   else /* down up: S+S- terms */
   {
    new_state = ((bitmap | mask0) & ~mask1);
    matrixelement(komplex (Jxy[j]/2.0,0.0),k,coof,next);
#ifdef DIPOLE   /* S+S- and SzS- */
    matrixelement(komplex (-0.25*Jdip[j]*geom_13[j],0.0),k,coof,next);
    new_state = (bitmap & ~mask1);
    matrixelement((0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][X]) +
     I*(0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][Y]),k,coof,next);
#endif
   }  /* if s1==0. */
  } 
  else  /* case s0 up; S-.. terms (and Sz.. for DIPOLE) */
  {
#ifdef DIPOLE   /* S-Sz term */
   if (s1==0) sz=-0.5; else sz=0.5;
   new_state = (bitmap & ~mask0);
   matrixelement((-1.5*sz*Jdip[j]*r_vector[j][Z]*r_vector[j][X]) + 
    I*(-1.5*sz*Jdip[j]*r_vector[j][Z]*r_vector[j][Y]),k,coof,next);
#endif
   if (s1==0) /* up down: S-S+ term (and SzS+) */ 
   {
    new_state = ((bitmap | mask1) & ~mask0);
    matrixelement(komplex (Jxy[j]/2.0,0.0),k,coof,next);
#ifdef DIPOLE  /* S-S+ and SzS+ */
    matrixelement(komplex(-0.25*Jdip[j]*geom_13[j],0.0),k,coof,next);
    new_state = (bitmap | mask1);
    matrixelement((-0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][X])+
     I*(0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][Y]),k,coof,next);
#endif
   }
   else /* up up: S-S- terms (and SzS-) */
   {
#ifndef M_SYM
    new_state = (bitmap & ~mask0) & ~mask1;
    matrixelement(komplex (Janis[j]/2.0,0.0),k,coof,next);
#ifdef DIPOLE   /*  S-S- and SzS-  */
    matrixelement((-0.75*Jdip[j]*(SQR(r_vector[j][X])-SQR(r_vector[j][Y])))+
     I*(-1.5*Jdip[j]*r_vector[j][X]*r_vector[j][Y]),k,coof,next);
    new_state = (bitmap &  ~mask1);
    matrixelement((-0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][X])+
     I*(-0.75*Jdip[j]*r_vector[j][Z]*r_vector[j][Y]),k,coof,next);
#endif /* DIPOLE */
#endif /* M_SYM */
   }  /* if s1==0.. */
  } /* if s0==0.. */
 }  /* for(j=0.. */
  return;
}


void matrixelement(komplex Jval, int k[], 
                   komplex coof, komplex *next)
{
  long l;
  int i,new_occ,j,T[NSYM];
  unsigned long uniq;
  double norm;
  komplex C;

#ifdef TEST_MAT_ELEM
  printf("J: (%g + %g i)", real(Jval),imag(Jval));
#endif
  uniq=FindUnique(new_state,T); /* Unique after spin-flip */
  l=LookUpU(uniq);  /* Find position in table */
  /* Check for existence of new state with this k[] */
  if(new_occ=Nocc[l]) 
  {
    norm = sqroot[new_occ]/sqroot[u_occ];
    for (i=0, j=0; i<(Nsym);i++)
    {
#ifdef TEST_MAT_ELEM
      printf("i=%i, T(i)=%i, k(i)=%i, phase=%g \n",
	     i,T[i],k[i],T[i]*k[i]*(Nspins/(double)Nsymvalue[i]));
#endif
      j += T[i]*k[i]*(Nspins/Nsymvalue[i]);
    }
    j = j%Nspins;

#ifdef TEST_HAMILTON
    printf(" Coupling from state %d to %d -> unique %d.\n"
           ,bitmap,new_state,uniq);
    printf(" J: (%lg +i %lg) j: %d, cos: %g\n",real(Jval),imag(Jval),j,cosine[j]);
#endif
    C = norm*Jval*coof*(cosine[j]+I*sine[j]); 
    next[l] += C;
#ifdef TEST_MAT_ELEM
    printf("coof: ( %g + %g i), contribution: (%g + %g i)",
	    real(coof), imag(coof), real(C), imag(C));
    printf("l: %d next: (%g + %g i)\n", 
	    l, real(next[l]), imag(next[l]));
#endif
  } 
  return;
}


