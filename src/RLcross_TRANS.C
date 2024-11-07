/* Program file RLcross.C - 
* Calculating cross sections 
* Last change: SJ 23.06.16
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 3.1, 2016
*  NEW
============================================
*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <complex>
#include <RLexact.h>
#include <nr.h>

#ifdef FIND_CROSS

/* Functions declared elsewhere */
extern unsigned long long SymOp(long long, unsigned long long);
extern double LowestLanczos(long long*, komplex*, long long*, long long);
extern void WriteCross(long long, long long*, long long);
extern unsigned long long FindUnique(unsigned long long, long long *);
extern long long LookUpU(unsigned long long);
extern void LogMessageChar(const char *);
extern void LogMessageInt(long long);
extern void LogMessageImag(const double, const double);
extern void LogMessageCharDouble(const char *, double);
extern void LogMessageCharInt(const char *, long long);
extern void LogMessageChar3Vector(const char *, double, double, double);
extern void fatalerror(const char*, long long);

/* Functions declared in this file */
void ApplySzq(long long*);
#ifdef M_SYM
void ApplySmp(long long *, long long);
#endif
double lengthofvector(komplex*);
#ifdef LANCZOS
void CrossLanczos(long long *); 
#endif /* LANCZOS */

#ifdef NEVER //doesnt work, SJ 270616
#ifdef MATRIX
void CrossMatrix();
int sym;
#endif /* MATRIX */
bool NonZero(unsigned long long, long long*);
#endif /*NEVER*

/* Global variables defined in RLexact.c */
#ifdef M_SYM
extern long long twom;
#endif /* M_SYM */
extern long long Nspins,Nunique,hamil_coup[NCOUP][2];
extern long long Nuniq_k, Nsym;
extern long long uniq_k[];
extern long long *Nocc;
extern long long Nsymvalue[NSYM];
extern double **position;
extern double cosine[],sine[],sqroot[];
extern unsigned long long *unique;
extern komplex *gs;
extern double gs_energy;
extern komplex *szmpgs,*tmp;
extern double szqlength;
extern long long symlist[NSYM];
extern long long q_gs[NSYM];
// double *energies;
//FILE *outfilezz, *outfilepm, *outfilemp;

/* Regional variables defined here */
long long k[NSYM];


#ifdef LANCZOS
//void CrossLanczos(int symvalue[NSYM])
void CrossLanczos(long long *symvalue) //(Note: symvalue =qvector)
{
 long long j,r,Nener;
 long long flag[3];
 flag[0] = SZZ;
 flag[1] = SPM;
 flag[2] = SMP;
 #ifdef TEST_CROSS
     LogMessageCharInt("\nCalling CrossLanczos for m =",(twom)/2);
     LogMessageChar("\n");
#endif // TEST_CROSS

  /* Sort out which symmetry elements are in real-space */
  /* For a start: Consider only periodic boundary conditions */
  /* Apply S+-(k),S-+(k),Sz(k) to groundstate */

for (int i = 0; i < 3; i++)
{ 
  Nener =0;
  
  switch (flag[i]) {
    case SZZ: 
      ApplySzq(symvalue);
      break;
    case SPM: 
      #ifdef M_SYM
        ApplySmp(symvalue,0);
      #endif
      break;
    case SMP:
      #ifdef M_SYM 
        ApplySmp(symvalue,1);
      #endif
      break;
  } /*switch*/

  szqlength=lengthofvector(szmpgs); //for output in cross files

#ifdef TEST_LANCCROSS
     LogMessageChar("Calling LowestLanczos\n");
#endif // TEST_LANCCROSS


  LowestLanczos(symvalue,NULL,&Nener,CROSS); // This will do the actual work! When done, global variables *energies and *Cross contains the relevant values for cross-section
  
  WriteCross(Nener,symvalue,flag[i]);


} /*for flag*/
  
 return;
}

#endif /* LANCZOS */


/* ApplySzq applies the Sz(q) operator to a state vector */
/* ToDo: Must ensure that the unique is valid under the new q-symmetry */
// ReWritten by Kim, 14.07.00
void ApplySzq(long long *q)
{
long long i,n,phase,cycle;
unsigned long long state,gsstate,new_state;
komplex res,factor;
long long sym, T[NSYM],diffQ[NSYM];

LogMessageCharInt("Nsym =",Nsym);
for (sym=0;sym<Nsym;sym++) 
{
  diffQ[sym]=q[sym]-q_gs[sym];
  if (diffQ[sym]<0) 
  {
    diffQ[sym]+=Nsymvalue[sym];
  }
   if (symlist[sym]==SPIN_FLIP && diffQ[sym]==0)
   {
     for(i=0;i<Nunique;i++)
     {
     szmpgs[i]=zero;
     return;
     }
  }

  #ifdef TEST_APPLYSZQ
    LogMessageCharInt("ApplySzq: sym=",sym);
    LogMessageCharInt(", diffQ=",diffQ[sym]);
    LogMessageCharInt(", q=",q[sym]);
    LogMessageCharInt(", q_gs=",q_gs[sym]);
    LogMessageCharInt(", Nsym=",Nsym);
    LogMessageCharInt(", symlist=",symlist[sym]);
    LogMessageChar("\n");
  #endif // TEST_APPLYSZQ
} 

for (i=0; i<Nunique; i++)
{
  res=zero;
  gsstate = unique[i];
  #ifdef TEST_APPLYSZQ
  LogMessageCharInt("Nunique loop, i =",i);
  LogMessageCharInt("Groundstate coefficients for unique",gsstate);
  LogMessageCharDouble(" are",real(gs[i]));
  LogMessageCharDouble("+ i",imag(gs[i]));
  LogMessageChar("\n");
  #endif // TEST_APPLYSZQ
// *** TODO: This works only with spin flip, or with an identity symmetry as 1st symmetry. Too specific solution. FIX THIS !!! ***
  LogMessageChar("\n");

  if (Nocc[i] != 0) //
  { 
      state = 1; //bitmap for transloop
      LogMessageChar("\nIn Transloop");
      TRANSLOOP_BEGIN
        phase=0;
        for (long long j=1; j<Nsym;j++) //sum over all symmetries except spin-flip/identity
        {  // Dirty non-general solution !!!
          phase += diffQ[j]*T[j]*Nspins/Nsymvalue[j]; //T[j] in fact means r_j?!! And diffQ[j]*Nspins/Nsumvalue[j] is q normalized to length of system!
          //one uses diffQ because of theorem - to be calculated!
          #ifdef TEST_APPLYSZQ
          LogMessageCharInt("\n j =",j);
          LogMessageCharInt(", Nsymvalue = ",Nsymvalue[j]);
          LogMessageCharInt("phase =",phase);
          #endif
        }

        #ifdef TEST_APPLYSZQ
        LogMessageCharInt("\nNunique i=",i);
        LogMessageCharInt(" T[0]=",T[0]);
        LogMessageCharInt(", T[1]=",T[1]);
        LogMessageCharInt("\nNew_state =",new_state);
        LogMessageCharInt(", gsstate =",gsstate);
        LogMessageCharDouble(" with phase =",phase);
        LogMessageCharInt(" Also, new_state & gsstate =",(new_state & gsstate));
        #endif

          phase = phase%Nspins;
             
        //something with spin up or down?
        if (new_state & gsstate) 
          res += cosine[phase] + I*sine[phase]; //Note: cosine[i] = cos(2*PI*i/Nspins); 
        else 
          res -= cosine[phase] + I*sine[phase];
        #ifdef TEST_APPLYSZQ
        LogMessageCharInt("\nAfter if/else: phase =",phase);
        LogMessageCharDouble(". Res is now: ",real(res));
        LogMessageCharDouble("+ i",imag(res));
        LogMessageChar("\n");
        #endif // TEST_APPLYSZQ
    TRANSLOOP_END
    LogMessageChar("Exit Transloop\n");

    res /= (2*sqroot[Nspins]);
    #ifdef TEST_APPLYSZQ
      LogMessageCharDouble("After Transloop res =",real(res));
      LogMessageCharDouble("+ i",imag(res));
      LogMessageChar("\n");
    #endif // TEST_APPLYSZQ

    szmpgs[i] = gs[i]*res;
  }
  else //Nocc == 0
  {
    LogMessageChar("State is forbidden!\n");
    szmpgs[i] = zero;
  }

  #ifdef TEST_APPLYSZQ
  LogMessageCharInt ("ApplySzq output: m =",twom/2);
  LogMessageCharInt(", unique ",i);
  LogMessageCharInt(" corresponding to state ", unique[i]);
  LogMessageCharDouble("with weight: ",real(szmpgs[i]));
  LogMessageCharDouble(" + i",imag(szmpgs[i]));
  LogMessageChar("\n");
  #endif
} 

return;
}


#ifdef M_SYM
/* ApplySmp applies the S-S+(q) operator (to the ground state vector).
Last Change:  SJ 27.07.16
The which_q parameter determines if it is the S-S+(q) operator 
(which_q=0) or the S-S+(q) operator (which_q=1). */
void ApplySmp(long long *q, long long which_q)
{
  unsigned long long i,l,j,k;
  long long T[Nsym];
  unsigned long long gsstate,downup,updown,new_state,u,state,smp_state;
  long long phase,n_flip,u_cycle,new_cycle;
  double norm;
  long long u_occ;
  long long sym,diffQ[NSYM];
  komplex res;

#ifdef TEST_APPLYSMP
  LogMessageCharInt("ApplySMP q =(",q[0]);
  LogMessageCharInt(",",q[1]);
  LogMessageChar(")\n");
#endif


for(i=0;i<Nunique;i++) /* Run through unique */
  {    
    gsstate = unique[i];
    res = zero;
    u_occ = Nocc[i];

   /* Check if unique is component of ground state. */
    if(Nocc[i] != 0) 
    {
      n_flip=0;
      /* Run through possible spin-flips */
      for(k=0;k<Nspins;k++) {
        for(j=0;j<Nspins;j++) {
        
          downup = ((unsigned long long) 1)<<k; //spin to raise
          updown = ((unsigned long long) 1)<<j; //spin to lower
#ifdef TEST_APPLYSMP
  LogMessageCharInt("\ngsstate:",gsstate);
  LogMessageCharInt(", downup:", downup);
  LogMessageCharInt(", updown:",updown);
  LogMessageCharInt(" - !(gsstate&downup)) ->",!(gsstate&downup));
  LogMessageCharInt(" gsstate&updown ->",gsstate&updown);
#endif 
/* Check if up-down operation is possible, including operation on the same spin */
        if( ( !(gsstate & downup)) && (( gsstate & updown) || (k==j)) ) //if=true as long as condition != 0
        { 
          n_flip++;
          if (k==j)
            smp_state = gsstate;
          else
            smp_state = (gsstate | downup)&(~updown);

#ifdef TEST_APPLYSMP
  LogMessageCharInt("\nk=",k);
  LogMessageCharInt(" j=",j);
  LogMessageCharInt("smp_state=",smp_state);
#endif

          u = FindUnique(smp_state,T); // Unique after spin-flip
          l = LookUpU(u);  // Find position in table 
          LogMessageCharInt("SMP: Nocc[l]=",Nocc[l]);
          LogMessageChar("\n");
          if (Nocc[l] != 0)
          {
            LogMessageCharInt("Nsym =",Nsym);
            for (sym=0;sym<Nsym;sym++) 
            {
              diffQ[sym]=q[sym]-q_gs[sym];
              if (diffQ[sym]<0) 
              {
                diffQ[sym]+=Nsymvalue[sym];
              }
               if (symlist[sym]==SPIN_FLIP && diffQ[sym]==0)
               {
                 for(i=0;i<Nunique;i++)
                 {
                 szmpgs[i]=zero;
                 return;
                 }
              }

              #ifdef TEST_APPLYSMP
                LogMessageCharInt("ApplySMP: sym=",sym);
                LogMessageCharInt(", diffQ=",diffQ[sym]);
                LogMessageCharInt(", q=",q[sym]);
                LogMessageCharInt(", q_gs=",q_gs[sym]);
                LogMessageCharInt(", Nsym=",Nsym);
                LogMessageCharInt(", symlist=",symlist[sym]);
                LogMessageChar("\n");
              #endif // TEST_APPLYSMP
            } 

            state = 1; //bitmap for transloop
            LogMessageChar("\nIn Transloop");
            TRANSLOOP_BEGIN
              phase=0;
              for (long long s=1; s<Nsym;s++) //sum over all symmetries except spin-flip/identity
              {  // Dirty non-general solution !!!
                phase += diffQ[s]*T[s]*Nspins/Nsymvalue[s]; //T[j] in fact means r_j?!! And diffQ[j]*Nspins/Nsumvalue[j] is q normalized to length of system!
                //one uses diffQ because of theorem - to be calculated!
                #ifdef TEST_APPLYSMP
                LogMessageCharInt("\n j =",s);
                LogMessageCharInt(", Nsymvalue = ",Nsymvalue[s]);
                LogMessageCharInt("phase =",phase);
                #endif
              }

              #ifdef TEST_APPLYSMP
              LogMessageCharInt("\nNunique i=",i);
              LogMessageCharInt(" T[0]=",T[0]);
              LogMessageCharInt(", T[1]=",T[1]);
              LogMessageCharInt("\nNew_state =",smp_state);
              LogMessageCharInt(", gsstate =",gsstate);
              LogMessageCharDouble(" with phase =",phase);
              LogMessageCharInt(" Also, new_state & gsstate =",(new_state & gsstate));
              #endif

              phase = phase%Nspins;
            
              res += cosine[phase] + I*sine[phase];
              #ifdef TEST_APPLYSMP
              LogMessageCharInt("\nAfter if/else: phase =",phase);
              LogMessageCharDouble(". Res is now: ",real(res));
              LogMessageCharDouble("+ i",imag(res));
              LogMessageChar("\n");
              #endif // TEST_APPLYSMP
            TRANSLOOP_END
            LogMessageChar("Exit Transloop\n");

           #ifdef TEST_APPLYSMP
              LogMessageCharDouble("After Transloop res =",real(res));
              LogMessageCharDouble("+ i",imag(res));
              LogMessageChar("\n");
            #endif // TEST_APPLYSMP

            norm = sqroot[Nocc[l]]/sqroot[u_occ]; //should these reversed? meh

            if (which_q==0) //Meaning S+S-(q)
              norm /= sqroot[2*(twom/2+1)*Nspins];  // Splus operator effect 
            else if (which_q==1) //S-S+(q)
              norm /= sqroot[2*(twom/2+1)*Nspins];  // Sminus operator effect 

            #ifdef TEST_APPLYSMP
              LogMessageCharInt("u_occ=",u_occ);
              LogMessageCharInt(", Nocc=",Nocc[l]);
              LogMessageCharDouble(" -> norm =",(norm));
            #endif

            szmpgs[l] = gs[i]*res*norm;                        
            /*
            if (which_q==1)
            {
              phase = (2*Nspins*Nspins - k*q + (*s)*(q+kGNonZero*Nspins/2) )%Nspins;
                    if (phase > 0 && kGNonZero == 1) 
                       phase = Nspins - phase;  // PATCH for complex conjugation 
            }
            else if (which_q==0)
            {
               phase = (2*Nspins*Nspins - j*q + (*s)*(q-kGNonZero*Nspins/2) )%Nspins;
                    if (phase > 0 && kGNonZero == 1)
                       phase = Nspins - phase;  // PATCH for complex conjugation 
            }
                    
           */

          } // Nocc 
          else {szmpgs[l]=zero;}
        }   // Interchangable 
      }     // Double k,j loop 
    }
    }       // Nocc 
  }         // Unique loop 


  return;
}
#endif //MSYM

#ifdef NEVER
#ifdef MATRIX
// CrossMatrix calculates now only S^zz(q) on the whole set of eigenstates
// WARNING: all symmetries are considered to be spatial periodic translations !!
// Written by Kim, 14.07.00
void CrossMatrix(long long symvalue[NSYM])
 {
  long long i,j,q[NSYM];
  komplex res,basis_vector[NUNIQUE], szq_vector[NUNIQUE];

  QLOOP_BEGIN  //* These are the Q's in S(Q), perhaps do this differently? 
    for (i=0; i<Nuniq_k; basis_vector[i++]=1);
//    ApplySzq(basis_vector,symvalue,q); // basis_vector is the trace of a diagonal matrix 
    for (i=0; i<Nuniq_k; i++)
     {
      res=0;
      for (j=0; j<Nuniq_k; j++)
//        res+=basis_vector[j]*eigenstates[i][j]; // Multiplying with a diagonal matrix 
//      sqz_vector=res*conj(res)
;
     }
  QLOOP_END
  
  return;
 }
#endif /* MATRIX */

bool NonZero(unsigned long long state, long long *q) {
  // simple implementation of routine to determine if state is value in q-space q
  unsigned long long new_state;
  long long eksponent;
  long long sym, T[NSYM];

  LogMessageCharInt("\nIn NonZero:  state =",state);
  LogMessageCharInt(" q =(",q[0]);
  LogMessageCharInt(",",q[1]);
  LogMessageChar(")\n");

TLOOP_BEGIN // goes through all symmetries defined, setting new_state and T[]
  if (new_state==state) { // if we got back to the unique
    eksponent =0;
    for (long long i=0;i<Nsym;i++) {
      eksponent+=q[i]*T[i]*Nspins/Nsymvalue[i];  // q_i*d_i*N/n_i

      LogMessageCharInt("i =",i);
      LogMessageCharInt("T[i] =",T[i]);
      LogMessageCharInt(", eksponent%Nspins =",eksponent%Nspins);
      LogMessageChar("\n");

    }
    if (eksponent%Nspins != 0) { // equiv. "if coefficient not 1"
      LogMessageChar(" - Non Zero false!\n");
      return false; // coefficient will add up to 0
    }
  }
TLOOP_END
  LogMessageChar(" - Non Zero true!\n");
  return true; // all coefficients were equal to one.
}

#endif /*NEVER*/

 double lengthofvector(komplex *v) {
  double length=0.0;
  for (long long i=0;i<Nunique;i++) {
    length+=((real(v[i])*real(v[i]))+(imag(v[i])*imag(v[i])));
  }
  return length;
}


#endif /* FIND_CROSS */