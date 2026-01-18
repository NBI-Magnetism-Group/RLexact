/* Program file RLhamil.C -
* Applying the Hamilton operator to a given state
* Last change: KL 25.09.12
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94
* Version 4.0, September 2017
*
============================================
*/

// #include "/usr/include/sys/types.h"
#include <sys/types.h>

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// #include <nr.h>
#include <cnr.h>
#include <RLexact.h>
#include <complex>

/* Functions defined in this file */
void Hamil2_sparse(unsigned long long, unsigned long long *, long long, long long, int *, long long *, int *, komplex *, double *, FILE *, FILE *, FILE *);
void Hamil_Zeeman(unsigned long long, unsigned long long *, long long, int *, long long *, int *, komplex *, double *, FILE *, FILE *, FILE *);
void Hamil4_sparse(unsigned long long, unsigned long long *, long long, long long, int *, long long *, int *, komplex *, double *, FILE *, FILE *, FILE *);

void Eigenvector_test(long long *, komplex *, komplex *);

/* Functions defined elsewhere */
extern long long LookUpU(unsigned long long);
extern long long Count(unsigned long long);
void fatalerror(const char *, long long);
unsigned long long FindUnique(unsigned long long, int *);
extern void WriteGSEnergy(komplex);
extern void WriteState(const char *, komplex *);
extern void time_stamp(time_t *, long long, const char *);

extern void LogMessageChar(const char *);
extern void LogMessageInt(long long);
extern void LogMessageImag(long long);
extern void LogMessageCharDouble(const char *, double);
extern void LogMessageCharInt(const char *, long long);
extern void LogMessageChar3Vector(const char *, long long, long long, long long);
// extern void WriteCouplingFiles(unsigned long long, unsigned long long, long long *, double, long long, long long, long long *, long long *, FILE*, FILE*, FILE*);
void WriteCouplingFiles(unsigned long long, unsigned long long, int *, komplex, long long, long long, int *, long long *, FILE *, FILE *, FILE *);

/* Global variables defined in RLexact.c */
extern unsigned long long *unique;
extern long long Nsym, Nspins;
extern long long Nunique, Nuniq_k;
extern long long *uniq_k;
// extern long long *Nocc;
extern long long Nsymvalue[];
extern long long hamil_coup[NCOUP][2];
extern double Jzz[], Jxy[], Janis[];
extern double Jdip[], geom_13[], r_vector[NCOUP][3];
extern long long m;
extern double h, field[3];
extern double sine[], cosine[], sqroot[];
extern long long Ncoup;
#ifdef RING_EXCHANGE
extern double Jr[NRING];
extern long long ring_coup[NRING][4];
#endif /* RING_EXCHANGE */

/* Regional variables in this file */
unsigned long long bitmap, new_state;
long long n_2, u_occ;
unsigned long long index1, index2;
komplex this_;

void Hamil_Zeeman(unsigned long long bitmap, unsigned long long *new_state, long long i, int *nelem, long long *totcount, int *T, komplex *J, double *diag, FILE *indexfile, FILE *Tfile, FILE *Jfile)
{ // only supports fields along one coordinate axis
#ifdef TEST_HAMZEE
  LogMessageCharDouble("\nIn Hamil_Zeeman. h =", h);
  LogMessageCharInt(" bitmap =", bitmap);
#endif
  if (h != 0) // if there is no field, this doesnt matter
  {

#ifdef TEST_HAMZEE
    LogMessageCharDouble("\nfield[0] =", field[0]);
    LogMessageCharDouble(", field[1] =", field[1]);
    LogMessageCharDouble(", field[2] =", field[2]);
#endif
    unsigned long long mask0, s0;
    if (field[2] > 0) /* Mag field along z*/
    {
      double sz;
      sz = Count(bitmap) - Nspins / 2;
#ifdef TEST_HAMZEE
      LogMessageCharInt("\nHamzee Sz: Count(bitmap) =", Count(bitmap));
      LogMessageCharInt(" , and sz =", sz);
      LogMessageChar("\n");
#endif
      *diag -= h * sz * field[2];
    }

    if (field[0] > 0 || field[1] > 0) // Magnetic field "transverse", i.e. along x or y
    {                                 // Here J is simply the coupling strength
      for (int k = 0; k < Nspins; k++)
      {
        mask0 = ((unsigned long long)1) << hamil_coup[k][0];
        s0 = (bitmap & mask0) != 0;
#ifdef TEST_HAMZEE
        LogMessageCharInt("k =", k);
        LogMessageCharInt(" and s0 =", s0);
        LogMessageChar("\n");
#endif

        if (s0 == 0) // S+, spin can be raised
        {
          *new_state = (bitmap | mask0);
          *J = -h / 2 * field[0] - (h / 2 * I) * field[1];
#ifdef TEST_HAMZEE
          LogMessageCharInt("Hamzee: S+: From state ", bitmap);
          LogMessageCharInt("to state ", *new_state);
          LogMessageCharDouble("with J=", real(*J));
          LogMessageCharDouble("+ i", imag(*J));
          LogMessageChar("\n");
#endif
        }
        else // then s0==1, S-, spin can be lowered
        {
          *new_state = (bitmap & ~(mask0));
          *J = -h / 2 * field[0] + (h / 2 * I) * field[1];

#ifdef TEST_HAMZEE
          LogMessageCharInt("Hamzee: S-: From state ", bitmap);
          LogMessageCharInt("to state ", *new_state);
          LogMessageCharDouble("with J=", real(*J));
          LogMessageCharDouble("+ i", imag(*J));
          LogMessageChar("\n");
#endif
        }
        WriteCouplingFiles(bitmap, *new_state, T, *J, i, i, nelem, totcount, indexfile, Tfile, Jfile);
      }
    }
  }
}

void Hamil2_sparse(unsigned long long bitmap, unsigned long long *new_state,
                   long long i, long long j, int *nelem,
                   long long *totcount, int *T, komplex *J,
                   double *diag, FILE *indexfile, FILE *Tfile, FILE *Jfile,
                   struct FLAGS *input_flags)
{
#ifdef TEST_HAM2
  LogMessageChar("Now entering Hamil2_sparse function \n");
#endif
  unsigned long long mask0, mask1, s0, s1;

  mask0 = ((unsigned long long)1) << hamil_coup[j][0];
  mask1 = ((unsigned long long)1) << hamil_coup[j][1];
  s0 = (bitmap & mask0) != 0;
  s1 = (bitmap & mask1) != 0;
#ifdef TEST_HAM2
  LogMessageCharInt(" Coupling between ", hamil_coup[j][0]);
  LogMessageCharInt(" and ", hamil_coup[j][1]);
  LogMessageCharInt(", (s0,s1)= (", s0);
  LogMessageCharInt(", ", s1);
  LogMessageChar(")\n");
  LogMessageCharInt("Also, mask0 = ", mask0);
  LogMessageCharInt(" and mask1 = ", mask1);
#endif
  if ((s0 + s1) == 1) /* Spins are of opposite sign */
    *diag -= Jzz[j] * 0.25;
  else
    *diag += Jzz[j] * 0.25;

  if (s0 == 0) /* case s0 down; S+.. terms (and Sz.. for DIPOLE) */
  {

    if (s1 == 0) /* down down: S+S+ term */
    {
      if (input_flags->m_sym)
      {
        *new_state = (bitmap | mask0 | mask1);
        *J = Janis[j] / 2;
      }
    }
    else /* down up: S+S- terms */
    {
      *new_state = ((bitmap | mask0) & ~mask1);
      *J = Jxy[j] / 2.0;
    } /* if s1==0. */
  }
  else /* case s0 up; S-.. terms (and Sz.. for DIPOLE) */
  {

    if (s1 == 0) /* up down: S-S+ term (and SzS+) */
    {
      *new_state = ((bitmap | mask1) & ~mask0);
      *J = Jxy[j] / 2.0;
    }
    else /* up up: S-S- terms (and SzS-) */
    {
      if (!input_flags->m_sym)
      {

        *new_state = (bitmap & ~mask0) & ~mask1;
        *J = Janis[j] / 2;
      }
    } /* if s1==0.. */
  } /* if s0==0.. */

#ifdef TEST_HAM2
  LogMessageCharDouble("with J=", real(*J));
  LogMessageCharDouble("+ i", imag(*J));
  LogMessageChar("\n");
#endif

  WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
}

#ifdef RING_EXCHANGE
void Hamil4_sparse(unsigned long long bitmap, unsigned long long *new_state, long long i, long long j, int *nelem, long long *totcount, int *T, komplex *J, double *diag, FILE *indexfile, FILE *Tfile, FILE *Jfile)
{
  unsigned long long mask0, mask1, mask2, mask3, s0, s1, s2, s3;

  mask0 = ((unsigned long long)1) << ring_coup[j][0];
  mask1 = ((unsigned long long)1) << ring_coup[j][1];
  mask2 = ((unsigned long long)1) << ring_coup[j][2];
  mask3 = ((unsigned long long)1) << ring_coup[j][3];

  s0 = (bitmap & mask0) != 0;
  s1 = (bitmap & mask1) != 0;
  s2 = (bitmap & mask2) != 0;
  s3 = (bitmap & mask3) != 0;

#ifdef TEST_HAM4
  LogMessageChar(" Ring coupling: ");
  LogMessageInt(s0);
  LogMessageInt(s1);
  LogMessageInt(s2);
  LogMessageInt(s3);
  LogMessageChar("\n");
#endif

  if (s0 + s1 + s2 + s3 == 0 || s0 + s1 + s2 + s3 == 4 || s0 + s1 + s2 + s3 == 2) /* Hzzzz */
  {
    *diag += Jr[j] / 16.0;
#ifdef TEST_HAM4
    LogMessageCharDouble("  Hzzz term, J = ", Jr[j] / 16.0);
    LogMessageChar("\n");
#endif
  }
  else
  {
    *diag -= Jr[j] / 16.0;
#ifdef TEST_HAM4
    LogMessageCharDouble("  Hzzz term, J = ", -Jr[j] / 16.0);
    LogMessageChar("\n");
#endif
  }

  if (s0 + s1 + s2 + s3 == 2) /* m=0 rings*/
  {
    if (s0 == 0 && s1 == 1 && s2 == 0 && s3 == 1) /*Hpmpm, only alternating rings*/
    {
#ifdef TEST_HAM4
      LogMessageChar("  Hpmpm term");
      LogMessageChar("\n");
#endif
      *new_state = ((((bitmap | mask0) & ~mask1) | mask2) & ~mask3);
      *J = Jr[j] / 2.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 1 && s1 == 0 && s2 == 1 && s3 == 0)
    {
#ifdef TEST_HAM4
      LogMessageChar("  Hpmpm term");
      LogMessageChar("\n");
#endif
      *new_state = ((((bitmap | mask1) & ~mask0) | mask3) & ~mask2);
      *J = Jr[j] / 2.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 0 && s1 == 0 && s2 == 1 && s3 == 1) /*Hzpm down down up up*/
    {
#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask1) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 1 && s1 == 1 && s2 == 0 && s3 == 0) /*Hzpm up up down down*/
    {
#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask3) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 0 && s1 == 1 && s2 == 1 && s3 == 0) /*Hzpm down up up down*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif
      *new_state = ((bitmap | mask3) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask2);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 1 && s1 == 0 && s2 == 0 && s3 == 1) /*Hzpm up down down  up*/
    {
#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif
      *new_state = ((bitmap | mask1) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 1 && s1 == 0 && s2 == 1 && s3 == 0) /*Hzpm up down up down */
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif
      *new_state = ((bitmap | mask3) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s0 == 0 && s1 == 1 && s2 == 0 && s3 == 1) /*Hzpm down up down up */
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask0) & ~mask1);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
  } /* if m=0 ring */

  if (s0 + s1 + s2 + s3 == 1) /* if m=-1 ring */
  {
    if (s0 == 1) /*Hzpm up down down down*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask1) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask0);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s1 == 1) /*Hzpm down up down down*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask0) & ~mask1);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask1);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask1);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s2 == 1) /*Hzpm down down up down*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif
      *new_state = ((bitmap | mask1) & ~mask2);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s3 == 1) /*Hzpm down down down up*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask0) & ~mask3);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
  } /* if m=-1*/

  if (s0 + s1 + s2 + s3 == 3) /* if m=1 ring */
  {
    if (s0 == 0) /*Hzpm down up up up*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask0) & ~mask1);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask0) & ~mask2);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s1 == 0) /*Hzpm up down up up*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask1) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask1) & ~mask3);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s2 == 0) /*Hzpm up up down up*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask2) & ~mask1);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask3);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask2) & ~mask0);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
    if (s3 == 0) /*Hzpm up up up down*/
    {

#ifdef TEST_HAM4
      LogMessageChar("  Hzzpm term");
      LogMessageChar("\n");
#endif

      *new_state = ((bitmap | mask3) & ~mask0);
      *J = Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask2);
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);

      *new_state = ((bitmap | mask3) & ~mask1);
      *J = -Jr[j] / 8.0;
      WriteCouplingFiles(bitmap, *new_state, T, *J, i, j, nelem, totcount, indexfile, Tfile, Jfile);
    }
  } /*for m=1 rings*/
}
#endif /*RING_EXCHANGE*/

#ifdef TEST_EIG
void Eigenvector_test(long long k[NSYM], komplex *evec, komplex *tmp)
{
  long long i;
  komplex product;

  LogMessageChar("Entering Eigenvector_test.\n");
  Hamilton(evec, tmp, k);
  product = zero;
  LogMessageChar("Middle of Eigenvector_test.\n");
  for (i = 0; i < Nunique; i++)
    product += evec[i] * conj(tmp[i]);
  WriteGSEnergy(product);

#ifdef WRITE_STATES
  WriteState("Eigenvector:", evec);
  WriteState("H|e>:", tmp);
#endif /* WRITE_STATES */

  return;
}
#endif /* TEST_EIG */
