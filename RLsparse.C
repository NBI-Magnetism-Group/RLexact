/* Program file RLsparse.C - 
* Generating and applying sparse Hamiltonian file
* Last change: SJ 20.12.16
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 4.0, September 2017
* 
============================================
*/

//#include "/usr/include/sys/types.h"
#include <sys/types.h>

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<errno.h>
#include<complex>
#include"RLexact.h"
#include <cnr.h>

// Global variables defined in RLexact.c
extern unsigned long long *unique;
extern long long Ncoup;
extern long long Nunique;
extern long long twom;
extern long long Nelem;
extern long long Nspins, Nsym, Ncoup;
extern long long hamil_coup[][2];
extern long long Nsymvalue[NSYM];
extern long long *Nocc;
extern double Jzz[], Jxy[], Janis[];
#ifdef RING_EXCHANGE
  extern double Jr[];
  extern long long ring_coup[][4];
  extern long long Nring, Ncoup;
#endif /* RING_EXCHANGE */
extern double sine[],cosine[];
extern double h;
extern double field[3];

// Functions defined elsewhere
extern long long LookUpU(unsigned long long);
void fatalerror(const char *, long long);
unsigned long long FindUnique(unsigned long long, int *);
long long Count(unsigned long long);
extern void LogMessageChar(const char *);
extern void LogMessageInt(long long);
extern void LogMessageCharDouble(const char *, double);
extern void LogMessageCharInt(const char *, long long);
extern void LogMessageChar3Vector(const char *, double, double, double);
#ifndef M_SYM
extern void Hamil_Zeeman(unsigned long long, unsigned long long *, long long, int *, long long *, int *, komplex *, double *, FILE*, FILE*, FILE*);
#endif
extern void Hamil2_sparse(unsigned long long, unsigned long long *, long long, long long, int *, long long *, int *, komplex *, double *, FILE*, FILE*, FILE*);
#ifdef RING_EXCHANGE
extern void Hamil4_sparse(unsigned long long, unsigned long long *, long long, long long, int *, long long *, int *, komplex *, double *, FILE*, FILE*, FILE*);
#endif /*RING_EXCHANGE*/

// Functions defined in this file 
void MakeSparse();
void ApplySparse(komplex *vectin, komplex *vectout, long long *k);
void FillHamilSparse(komplex **hamil, long long *k);
void WriteCouplingFiles(unsigned long long, unsigned long long, int *, komplex, long long, long long, int *, long long *, FILE*, FILE*, FILE*);
//extern void WriteCouplingFiles(unsigned long long, unsigned long long, long long *, double, long long, long long, long long *, long long *, FILE*, FILE*, FILE*);
// File buffers
char Indexbuf[BUFFERSIZE];
char Jbuf[BUFFERSIZE];
char Tbuf[BUFFERSIZE];
char Nbuf[BUFFERSIZE];
char Diagbuf[BUFFERSIZE];
char Ntotbuf[BUFFERSIZE];

/* MakeSparse writes a series of 6 files containing all information on the q-independent matrix elements */
/* TODO: isolate file opening/closing in separate functions in RLio.C  */

void MakeSparse()
{
  long long i,strength,l,index,j,u_cycle,totcount=0;
  int nelem;
  unsigned long long mask0,mask1,s0,s1,bitmap,downup,updown,new_state,u;
  double diag,sz;
  komplex J(0,0);
  char filename[150];
  int T[NSYM];

  FILE *indexfile,*Jfile,*Tfile,*nelemfile,*diagfile, *ntotfile;
  
/* Open files, begin */
#ifdef M_SYM
  sprintf(filename,"%sIndex%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sIndex%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  indexfile=fopen(filename,"w");
  if (indexfile== NULL)
   {
    fatalerror("Cannot open indexfile, sorry!",errno);
   }
  setvbuf(indexfile,Indexbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sJ%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sJ%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  Jfile=fopen(filename,"w");
  if (Jfile== NULL)
   {
    fatalerror("Cannot open Jfile, sorry!",errno);
   }
  setvbuf(Jfile,Jbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sT%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sT%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  Tfile=fopen(filename,"w");
  if (Tfile== NULL)
   {
    fatalerror("Cannot open Tfile, sorry!",errno);
   }
  setvbuf(Tfile,Tbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sN%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sN%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  nelemfile=fopen(filename,"w");
  if (nelemfile== NULL)
   {
    fatalerror("Cannot open nelemfile, sorry!",errno);
   }
  setvbuf(nelemfile,Nbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sDiag%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sDiag%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  diagfile=fopen(filename,"w");
  if (diagfile== NULL)
   {
    fatalerror("Cannot open diagfile, sorry!",errno);
   }
  setvbuf(diagfile,Diagbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sNtot%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sNtot%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  ntotfile=fopen(filename,"w");
  if (ntotfile== NULL)
   {
    fatalerror("Cannot open ntotfile, sorry!",errno);
   }
  setvbuf(diagfile,Ntotbuf,_IOFBF,BUFFERSIZE);


/* Open files, end */
  #ifdef TEST_MAKESPARSE
    LogMessageChar("all files opened in makesparse\n");
  #endif

/* Loop over uniques, begin */
  for(i=0;i<Nunique;i++) 
  {    
    #ifdef TEST_MAKESPARSE
      LogMessageCharInt("\n Handling unique ",i);
      LogMessageCharInt(", bitmap = ",unique[i]);
      LogMessageChar("\n");
    #endif
    bitmap = unique[i];
    diag = 0;
    nelem = 0;

    #ifndef M_SYM /*The Zeeman term*/
    new_state = bitmap;
    J=0;
    Hamil_Zeeman(bitmap,&new_state,i,&nelem,&totcount,T,&J,&diag,indexfile,Tfile,Jfile);
    #endif

    /*Calculate on- and off-diagonal elements for normal Heisenberg */
    /* Loop over couplings, begin */
    for(j=0;j<Ncoup;j++)
    {  
      /* we must set new_state and J (=J_anis/2 or J_xy/2) to some
   default values to avoid bugs. If appropriate, they will be
   reset to true values later. - treue 20080204 */
      new_state = bitmap;
      J=0;
      Hamil2_sparse(bitmap,&new_state,i,j,&nelem,&totcount,T,&J,&diag,indexfile,Tfile,Jfile);  
    } /* loop over couplings, end */

    #ifdef RING_EXCHANGE
     /*loop over ring couplings*/
    for(j=0;j<Nring;j++) 
     {  
          new_state = bitmap;
          J=0;
          Hamil4_sparse(bitmap,&new_state,i,j,&nelem,&totcount,T,&J,&diag,indexfile,Tfile,Jfile);
     }  /*loop over ring couplings, end*/
        
    #endif /* RING_EXCHANGE */

    /* Write diagonal elements and non-diagonal counter to files */
    errno=0;
    if(fwrite(&diag,sizeof(double),1,diagfile)<=0) {
      LogMessageCharInt("Error in writing diagfile for unique ",i);
      fatalerror("Error was: ",errno);
    }
    errno=0;
    if(fwrite(&nelem,sizeof(int),1,nelemfile)<=0) {
      LogMessageCharInt("Error in writing nelemfile for unique ",i);
      fatalerror("Error was: ",errno);
    }
    #ifdef TEST_MAKESPARSE
      LogMessageCharInt(" i = ",i);
      LogMessageCharInt(", nelem = ",nelem);
      LogMessageCharDouble(", diag = ",diag);
      LogMessageChar("\n");
    #endif
  } /* Loop over uniques, end */
  
  Nelem = totcount;

  //write total count to file
  if(fwrite(&Nelem,sizeof(long long),1,ntotfile)<=0) {
    LogMessageChar("Error in writing Ntotfile");
  fatalerror("Error was: ",errno);
    }    

  
#ifdef TEST_MAKESPARSE
  LogMessageCharInt(" MakeSparse() recorded a total of ",Nelem);
  LogMessageChar(" off-diagonal elements \n");
#endif // TEST_MAKESPARSE

  /* Close files */
  fclose(nelemfile);
  fclose(diagfile);
  fclose(indexfile);
  fclose(Tfile);
  fclose(Jfile);
  fclose(ntotfile);
}


void WriteCouplingFiles(unsigned long long bitmap, unsigned long long new_state, int *T, komplex J, long long i, 
  long long j, int *nelem, long long *totcount, FILE* indexfile, FILE* Tfile, FILE* Jfile)
{
   long long l;
   unsigned long long u;
   u=FindUnique(new_state,T); /* Unique */
    l=LookUpU(u);  /* Find position in table */
  
    /* For elements in lower triangle of H, begin write couplings to file */
  #ifdef TEST_MAKESPARSE
  LogMessageCharInt("l = ",l);
  LogMessageCharInt("and i = ",i);
  LogMessageCharDouble(", J = ",real(J));
  LogMessageCharDouble("+ i",imag(J));
  LogMessageChar("\n");
  #endif

   if(i<=l && (real(J)!=0 || imag(J)!=0)) {
#ifdef TEST_MAKESPARSE
   LogMessageCharInt(" Chg. spins", hamil_coup[j][0]);
   LogMessageCharInt("and", hamil_coup[j][1]);
   LogMessageCharInt(", from bitmap", bitmap);
   LogMessageCharInt("to bitmap",new_state);
   LogMessageCharInt(", meaning unique: ",u);
   LogMessageCharDouble(", J = ",real(J));
   LogMessageCharDouble("+ i",imag(J));
   LogMessageChar("\n");

#ifdef TEST_HAM4
  LogMessageCharInt("   new_state = ",new_state);
  LogMessageCharInt(" unique = ",u);
  LogMessageCharDouble(", J = ",J);
  LogMessageChar("\n");
#endif

#endif
    errno=0;
        if(fwrite(&l,sizeof(long long),1,indexfile)<=0) {
            LogMessageCharInt("Error in writing indexfile for unique ",i);
      fatalerror("Error was: ",errno);
    }     
    errno=0;
          //if(fwrite(&J,sizeof(double),1,Jfile)<=0) {
          if(fwrite(&J,sizeof(komplex),1,Jfile)<=0) {
            LogMessageCharInt("Error in writing Jfile for unique ",i);
      fatalerror("Error was: ",errno);   
    }
    errno=0;
          if(fwrite(T,sizeof(int),Nsym,Tfile)<=0) {
            LogMessageCharInt("Error in writing Tfile for unique ",i);
      fatalerror("Error was: ",errno);
    }
    *totcount+=1;
    
    *nelem+=1;
    
#ifdef VERBOSE_MAKESPARSE //this flag doesnt exist anymore - cant see the relevance either!
          if(totcount%100000==00) {
            LogMessageCharInt("", totcount);
      LogMessageCharInt(" recorded, i = ",i);
      LogMessageChar("\n");
          }
#endif
      } 
      /* write coupling to file, end */

#ifdef TEST_MAKESPARSE
      LogMessageCharInt("nelem = ", *nelem);
      LogMessageChar("\n");
#endif 
}


// ApplySparse applies a sparse Hamilton operator read from series of 5 files to a state vector. CR 170300
// This function in RLexact uses the most computational time!
// TODO: Consider if the reading from files could be done more efficiently, e.g. should they be defined globally and kept open?
// TODO: Test if BUFFERSIZE is optimal

void ApplySparse(komplex *vectin, komplex *vectout, long long *k)
{
  int T[NSYM];
  long long l,s,index,j,phi,u_occ,elemcount=0;
  int nelem;
  long long i,r;
  double norm,diag;
  komplex elemJ;
  komplex coup;
  char filename[150];

  FILE *indexfile,*Jfile,*Tfile,*nelemfile,*diagfile, *ntotfile;

  /* Open files, begin */
#ifdef M_SYM
  sprintf(filename,"%sIndex%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sIndex%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  indexfile=fopen(filename,"r");
  if (indexfile== NULL)
   {
    fatalerror("Cannot open indexfile, sorry! Error was: ",errno);
   }
  setvbuf(indexfile,Indexbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sJ%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sJ%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  Jfile=fopen(filename,"r");
  if (indexfile== NULL)
   {
    fatalerror("Cannot open Jfile, sorry! Error was: ",errno);
   }
  setvbuf(Jfile,Jbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sT%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sT%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  Tfile=fopen(filename,"r");
  if (indexfile== NULL)
   {
    fatalerror("Cannot open Tfile, sorry! Error was: ",errno);
   }
#ifdef TEST_APPLYSPARSE
  for(i=0;i<Nelem;i++)  {
    fread(T,sizeof(long long),Nsym,Tfile);
    LogMessageCharInt("Coupling ",i);
    LogMessageCharInt(", T = ",T[0]);
    for (long long j=1;j<Nsym;j++) {
      LogMessageCharInt(",",T[j]);
    }
    LogMessageChar("\n");
  }
  fclose(Tfile);
  Tfile=fopen(filename,"r");
#endif
  setvbuf(Tfile,Tbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sN%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sN%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  nelemfile=fopen(filename,"r");
  if (indexfile== NULL)
   {
    fatalerror("Cannot open nelemfile, sorry! Error was: %s",errno);
   }
  setvbuf(nelemfile,Nbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sDiag%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sDiag%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  diagfile=fopen(filename,"r");
  if (indexfile== NULL)
   {
    fatalerror("Cannot open indexfile, sorry! Error was: ",errno);
   }
  setvbuf(diagfile,Diagbuf,_IOFBF,BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename,"%sNtot%llu_%llu.bin",MATRIXFILENAME,Nspins,twom);
#else
  sprintf(filename,"%sNtot%llu.bin",MATRIXFILENAME,Nspins);
#endif /* M_SYM */
  errno=0;
  ntotfile=fopen(filename,"r");
  if (ntotfile== NULL)
   {
    fatalerror("Cannot open ntotfile, sorry! Error was: ",errno);
   }
  setvbuf(ntotfile,Ntotbuf,_IOFBF,BUFFERSIZE);
  
   /*Read total number of elements */
    
    errno=0;
    if(fread(&Nelem,sizeof(long long),1,ntotfile)<=0) {
      LogMessageChar("Error in reading ntotfile for unique ");
      fatalerror("Error was: ",errno);
    }

/* Open files, end */

#ifdef TEST_APPLYSPARSE_READ
      LogMessageCharInt("\n Total number of elements, Nelem =",Nelem);
#endif
  
#ifdef TEST_APPLYSPARSE
  LogMessageChar(" ApplySparse() called with vector\n");
  for(i=0;i<Nunique;i++) 
    LogMessageCharDouble(" ",real(vectin[i]));
    LogMessageCharDouble(" + i ",imag(vectin[i]));
    LogMessageChar("\n");
#endif

  /* Apply sparse H to input vector (vectin) */
  
  /* Initialise output vector (vectout) */
  for(i=0;i<Nunique;i++) 
    vectout[i] = zero;

  for(i=0;i<Nunique;i++) /* Run through unique */
  {    
    errno=0;
    if(fread(&nelem,sizeof(int),1,nelemfile)<=0) {
      LogMessageCharInt("Error in reading nelemfile for unique ",i);
      fatalerror("Error was: ",errno);
    }
    errno=0;
    if(fread(&diag,sizeof(double),1,diagfile)<=0) {
      LogMessageCharInt("Error in reading diagfile for unique ",i);
      fatalerror("Error was: ",errno);
    }

#ifdef TEST_APPLYSPARSE_READ
      LogMessageCharInt("\n unique no =",i);
      LogMessageCharInt("correspondig to unique ",unique[i]);
      LogMessageCharInt(" number of couplings nelem =",nelem);
      LogMessageCharDouble(" and Diag =", diag);
      LogMessageChar("\n");
      //LogMessageCharInt("and T=",T);
#endif

    u_occ = Nocc[i];
#ifdef TEST_APPLYSPARSE
  LogMessageCharInt(" i = ",i);
  LogMessageCharInt(", nelem = ",nelem);
  LogMessageCharDouble(", diag = ",diag);
  LogMessageCharInt(", Nocc = ",u_occ);
  LogMessageChar("\n");
#endif

    /* Diagonal term */
    vectout[i] += diag*vectin[i];
    
#ifdef TEST_APPLYSPARSE_DEEP
    for (r=0;r<Nunique;r++) {
      LogMessageCharInt("vectout[",r);
      LogMessageCharDouble("] = ",real(vectout[r]));
      LogMessageCharDouble(" + i ",imag(vectout[r]));
      LogMessageChar("\n");
    }
#endif // TEST_APPLYSPARSE_DEEP


    /* Loop over interactions concerning the particular unique (i) */
    for(j=0;j<nelem;j++) {
      errno=0;
      if(fread(&l,sizeof(long long),1,indexfile)<=0) {
        LogMessageCharInt("Error in reading indexfile for unique ",i);
  fatalerror("Error was: ",errno);
      }
      errno=0;
      if(fread(&elemJ,sizeof(komplex),1,Jfile)<=0) {
        LogMessageCharInt("Error in reading Jfile for unique ",i);
  fatalerror("Error was: ",errno);
      }
      errno=0;
      if(fread(T,sizeof(int),Nsym,Tfile)<=0) {
        LogMessageCharInt("Error in reading Tfile for unique ",i);
  fatalerror("Error was: %s",errno);
      }
      elemcount++;

#ifdef TEST_APPLYSPARSE_READ
      LogMessageCharDouble("\n nelem =",j);
      LogMessageCharInt(" Index =",l);
      LogMessageCharInt(" equals Nocc[l] =",Nocc[l]);
      LogMessageCharDouble("and elemJ =", real(elemJ));
      LogMessageCharDouble("+i", imag(elemJ));
      LogMessageChar("\n");
      //LogMessageCharInt("and T=",T);
#endif

      
      // If state is allowed for present k[], begin 
      if((Nocc[l]>0) && (u_occ>0)) {
        if(elemcount>Nelem) 
    fatalerror("Too few elements in H, i=",i);
        /* Calculate phase from T[] and K[] */
  phi=0;
        for(s=0;s<Nsym;s++) {
          phi += T[s]*k[s]*Nspins/Nsymvalue[s];
#ifdef DEBUG_APPLYSPARSE
    if (s==0) {
      for (long long testsym=0;testsym<Nsym;testsym++) {
        LogMessageCharInt("sym ",testsym);
        LogMessageCharInt(" has T = ",T[testsym]);
        LogMessageCharInt(" q = ",k[testsym]);
        LogMessageChar("\n");
      }
    }
    LogMessageCharInt("From state ",i);
    LogMessageCharInt(" to ", l);
    LogMessageCharInt(" s = ",s);
    LogMessageCharInt(" T = ",T[s]);
    //LogMessageCharInt(" q = ",q);
    LogMessageCharInt(" k = ",k[s]);
    LogMessageCharInt(" got phase ",phi);
    LogMessageChar("\n");
#endif /* DEBUG_APPLYSPARSE */
        }
        phi = phi%Nspins;
        coup = elemJ*(cosine[phi]+I*sine[phi])*sqrt((double)Nocc[l]/u_occ);
#ifdef TEST_APPLYSPARSE
    LogMessageCharInt(" k1 = ", k[3]);
  LogMessageCharInt(" k2 = ", k[2]);
  LogMessageCharInt(" Coupling to state ",l);
  LogMessageCharInt(" phase = ",phi);
  LogMessageCharDouble(" coup = ",real(coup));
  LogMessageCharDouble(" + i ",imag(coup));
  LogMessageChar("\n");
#endif // TEST_APPLYSPARSE

        vectout[l] += coup*vectin[i];
        /* Add the other side of the diagonal, utilize that H is Hermitean */
        if(i!=l) {
    vectout[i] += real(coup)*vectin[l]-I*imag(coup)*vectin[l];  // TODO: replace with conjucate function ??
#ifdef TEST_APPLYSPARSE
//    printf("conjugate gives (from %llu to %llu): %g+i*%g (invector=%g+i*%g)\n",l,i,real(real(coup)*vectin[l]-I*imag(coup)*vectin[l]),imag(real(coup)*vectin[l]-I*imag(coup)*vectin[l]),real(vectin[l]),imag(vectin[l]));
#endif // TEST_APPLYSPARSE
  }
      }  /* If state is allowed for present k[], end */
#ifdef TEST_APPLYSPARSE_DEEP
    for (r=0;r<Nunique;r++) {
      LogMessageCharInt("vectout[",r);
      LogMessageCharDouble("] = ",real(vectout[r]));
      LogMessageCharDouble(" + i ",imag(vectout[r]));
      LogMessageChar("\n");
    }
#endif // TEST_APPLYSPARSE_DEEP
    }  /* Loop over interactions, end */

  } /* Loop over uniques, end */
    
  if(elemcount != Nelem) {
    LogMessageCharInt(" Used: ",elemcount);
    LogMessageCharInt(" elements of H, avaliable: ",Nelem);
    fatalerror("Exiting",0);
  }

#ifdef TEST_APPLYSPARSE
  LogMessageChar(" ApplySparse() returning vector: \n");
  for(i=0;i<Nunique;i++) {
    LogMessageCharDouble(" ",real(vectout[i]));
    LogMessageCharDouble(" + i ",imag(vectout[i]));
    LogMessageChar("\n");
  }
#endif

/* Close files */
  fclose(nelemfile);
  fclose(diagfile);
  fclose(indexfile);
  fclose(Tfile);
  fclose(Jfile);
  fclose(ntotfile);
}

// FillHamilSparse fills Hamiltonian for exact matrix diagonalization with help of 
// sparse matrix file. SJ 12-01-17
void FillHamilSparse(komplex **hamilI, long long *k)
{
  int T[NSYM];
  long long l, s, index, j, phi, u_occ, elemcount = 0;
  int nelem;
  
  long long i, r;
  long long is=10000;
  long long countZ = 0;
  long long SecondC;
  double norm, diag;
  long long *Nocc_c;
  Nocc_c = (long long*)malloc(sizeof(long long)*(Nunique+1));
  komplex coup, elemJ;
  char filename[150];
  FILE *indexfile, *Jfile, *Tfile, *nelemfile, *diagfile;

  /* Open files */
#ifdef M_SYM
  sprintf(filename, "%sIndex%llu_%llu.bin", MATRIXFILENAME, Nspins, twom);
#else
  sprintf(filename, "%sIndex%llu.bin", MATRIXFILENAME, Nspins);
#endif /* M_SYM */
  errno = 0;
  indexfile = fopen(filename, "r");
  if (indexfile == NULL)
  {
    fatalerror("Cannot open indexfile, sorry! Error was: ", errno);
  }
  setvbuf(indexfile, Indexbuf, _IOFBF, BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename, "%sJ%llu_%llu.bin", MATRIXFILENAME, Nspins, twom);
#else
  sprintf(filename, "%sJ%llu.bin", MATRIXFILENAME, Nspins);
#endif /* M_SYM */
  errno = 0;
  Jfile = fopen(filename, "r");
  if (indexfile == NULL)
  {
    fatalerror("Cannot open Jfile, sorry! Error was: ", errno);
  }
  setvbuf(Jfile, Jbuf, _IOFBF, BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename, "%sT%llu_%llu.bin", MATRIXFILENAME, Nspins, twom);
#else
  sprintf(filename, "%sT%llu.bin", MATRIXFILENAME, Nspins);
#endif /* M_SYM */
  errno = 0;
  Tfile = fopen(filename, "r");
  if (indexfile == NULL)
  {
    fatalerror("Cannot open Tfile, sorry! Error was: ", errno);
  }
  setvbuf(Tfile, Tbuf, _IOFBF, BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename, "%sN%llu_%llu.bin", MATRIXFILENAME, Nspins, twom);
#else
  sprintf(filename, "%sN%llu.bin", MATRIXFILENAME, Nspins);
#endif /* M_SYM */
  errno = 0;
  nelemfile = fopen(filename, "r");
  if (indexfile == NULL)
  {
    fatalerror("Cannot open nelemfile, sorry! Error was: %s", errno);
  }
  setvbuf(nelemfile, Nbuf, _IOFBF, BUFFERSIZE);
#ifdef M_SYM
  sprintf(filename, "%sDiag%llu_%llu.bin", MATRIXFILENAME, Nspins, twom);
#else
  sprintf(filename, "%sDiag%llu.bin", MATRIXFILENAME, Nspins);
#endif /* M_SYM */
  errno = 0;
  diagfile = fopen(filename, "r");
  if (indexfile == NULL)
  {
    fatalerror("Cannot open indexfile, sorry! Error was: ", errno);
  }
  setvbuf(diagfile, Diagbuf, _IOFBF, BUFFERSIZE);
  SecondC = 0;

  for (i = 0; i < Nunique; i++) /* Run through unique */
  {
    if (Nocc[i] > 0) {
      Nocc_c[i] = countZ;
    }
    else
    {
      Nocc_c[i] = countZ;
      countZ++;
    }
  }

  for (i = 0; i<Nunique; i++) /* Run through unique */
  {
    errno = 0;

    if (fread(&nelem, sizeof(int), 1, nelemfile) <= 0) {
      LogMessageCharInt("Error in reading nelemfile for unique ", i);
      fatalerror("Error was: ", errno);
    }
    errno = 0;
    if (fread(&diag, sizeof(double), 1, diagfile) <= 0) {
      LogMessageCharInt("Error in reading diagfile for unique ", i);
      fatalerror("Error was: ", errno);
    }

    u_occ = Nocc[i];

    /* Diagonal terms */
    if (Nocc[i] > 0) {
      hamilI[i - Nocc_c[i] + 1][i - Nocc_c[i] + 1] += diag;
    }


    /* Loop over interactions */
    for (j = 0; j<nelem; j++) {
      errno = 0;
      if (fread(&l, sizeof(long long), 1, indexfile) <= 0) {
        LogMessageCharInt("Error in reading indexfile for unique ", i);
        fatalerror("Error was: ", errno);
      }
      errno = 0;
      if (fread(&elemJ, sizeof(komplex), 1, Jfile) <= 0) {
        LogMessageCharInt("Error in reading Jfile for unique ", i);
        fatalerror("Error was: ", errno);
      }
      errno = 0;
      if (fread(T, sizeof(int), Nsym, Tfile) <= 0) {
        LogMessageCharInt("Error in reading Tfile for unique ", i);
        fatalerror("Error was: %s", errno);
      }

      elemcount++;

      // If state is allowed for present k 
      if ((Nocc[l]>0) && (u_occ>0)) {
        if (elemcount > Nelem){
          LogMessageCharInt("\nl = ", l);
          LogMessageCharInt("\nNocc[l] = ",Nocc[l] );
          LogMessageCharInt("\nNelem= ", Nelem);
          LogMessageCharInt("\nelemcount ", elemcount);
          fatalerror("Too few elements in H, i=", i);
        }
        /* Calculate phase */
        phi = 0;
        for (s = 0; s<Nsym; s++) {
          phi += T[s] * k[s] * Nspins / Nsymvalue[s];
        }
        phi = phi%Nspins;
        coup = one*elemJ*(cosine[phi] + I*sine[phi])*sqrt((double)Nocc[l] / u_occ);
        /* Off diagonal */
        if (i != l) {
          hamilI[i - Nocc_c[i] + 1][l - Nocc_c[l] + 1] += coup;
            /* Other side of diagonal */
          hamilI[l - Nocc_c[l] + 1][i - Nocc_c[i] + 1] += conj(coup);
        }
        else{
          hamilI[i - Nocc_c[i] + 1][i - Nocc_c[i] + 1] += coup;
        }
      }
      else {
        SecondC++;

      }
    }

  }

  if (elemcount != Nelem) {
    LogMessageCharInt(" Used: ", elemcount);
    LogMessageCharInt(" elements of H, avaliable: ", Nelem);
    fatalerror("Exiting", 0);
  }
  fclose(nelemfile);
  fclose(diagfile);
  fclose(indexfile);
  fclose(Tfile);
  fclose(Jfile);
  free(Nocc_c); //SJ 010616
}
