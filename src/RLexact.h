/* Include file RLexact.h - 
* Main include file for diagonalization program RLexact
* Last change: SJ 20.11.2017
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94  
* Version 4.0, September 2017
* 
============================================
*/

// Phyxsics of problem 
#define M_SYM
//#define DIPOLE     
//#define RING_EXCHANGE

// Calculation of problem and output requests 
#define LANCZOS 
//#define MATRIX   
#define FIND_EIGENSTATE
#define WRITE_ENERGIES
#define FIND_CROSS //Find cross in terms of S^zz (q,w), S^xx (q,w) and S^yy(q,w) 
//#define FIND_CROSS_PM //Find cross in S^+-(q,w) and S^-+(q,w) instead of S^xx (q,w) and S^yy(q,w). Requires MSYM and FIND_CROSS.
#define WRITE_STATES
//#define FIND_MAG //Debugging required! Should only be used WITHOUT MSYM SJ 20/11/17
//#define WRITE_MAGNETISATION //Should only be used WITHOUT MSYM , Works only for MATRIX-mode, SJ 31/5/16
//#define STRUCTURE  //doesnt currently work, SJ, 2/3/16 //Rename to Sqw_Q_VALUES AP 4/11/24 - Kim er næsten sikker
#define MOTIVE //spin positions

// Dimensions of problem 
#define NCOUP 400
#define NCOUPSTR 10
#define NSPINS 100
#define NSYM 20
#define NSYMADD 10
#define NUNIQUE 4500
#define MAXARRAYSIZE 200 // max number of entries on a given line in the input file
#ifdef RING_EXCHANGE
  #define NRING 100
  #define NRINGSTR 10
#endif

// Highest possible state
// #ifdef LANCZOS
//  #define MAX_STATE ((((unsigned long) 1)<<(Nspins-1)) -1) 
// #else
  #define MAX_STATE ((unsigned long long)1<<(Nspins)) 
// #endif

// Set up size of buffer for sparse matrix read/write
#ifdef MATRIX
  #define BUFFERSIZE BUFSIZ
#else
  #define BUFFERSIZE 4194304
#endif
 
#define USE_COMPLEX

#ifdef LANCZOS
// Generated vectors with squared lengths below the following 
// cutoff are considered to be zero in Lanczos calculation because of numerical precision.
#define ZERO_VEC_LENGTH 1e-9
#define SHORT_VEC_LENGTH (double) 1/NUNIQUE
#define MAX_LANCZOS 600
#endif  // LANCZOS 
#define right(i) (((i&1)<<(Nspins-1)) + (i>>1))
#define inverse(i) (i^(MAX_STATE-1))

// Symmetry related definitions 
#define MASK_L1 0x77777777
#define MASK_L2 0x88888888
#define MASK_L3 0xF000F000
#define MASK_L4 0x0FFF0FFF
#define MASK_P1 0x0000FFFF
#define MASK_P2 0xFFFF0000

// Hard coded symmetry flags
#define SPIN_FLIP    0
#define FCC32Tx      1
#define FCC32Ty      2
#define FCC32Tz      3
#define FCC32MIRROR  4
#define FCC32R2      5
#define FCC32R4      6
#define FCC32R3      7
#define IDENTITY     9
#define LINEAR_T    10

#define SYMNUM(i,j) ( (long long) (log( (double)SymOp( (long long)i,((unsigned long long)1)<<j)+1E-6) / log((double)2) )) 

// Deal with complex numbers 
#ifdef USE_COMPLEX
//typedef std::complex<double> komplex;
#define komplex std::complex<double>
//  #define komplex _Complex double
//  #define komplex std::complex<double>;
  #define kvector cvector
  #define freekvector freecvector
  #define kmatrix cmatrix
  #define freekmatrix freecmatrix
  #define I komplex (0.0, 1.0) //_Complex_I
  #define zero komplex (0.0, 0.0)
  #define one komplex (1.0, 0.0)
  #define Arg(a) (a==zero ? 0.0 : arg(a))
  #define arg(a) atan2(imag(a),real(a))
  #define abs(a) sqrt(real(a)*real(a) + imag(a)*imag(a))
  #define sqrabs(a) (real(a)*real(a) + imag(a)*imag(a))
  #define eksp(a) exp(real(a))*(cos(imag(a))+I*sin(imag(a))) 
  #define conj(a)  real(a)-I*imag(a)
  #define skrt(a) (sqrt(abs(a))*(cos(Arg(a)/2.0)+I*sin(Arg(a)/2.0) ) )
//#define skrt(a) (sqrt(real(a))*(cos(Arg(a)/2.0)+I*sin(Arg(a)/2.0) ) )

#else
  #define komplex double
  #define skrt(a) sqrt(a)
  #define real(a) (a)
  #define imag(a) 0.0
  #define conj(a) (a)
  #define abs(a) fabs(a)
  #define eksp(a) exp(a)
  #define Arg(a) 0.0
  #define I 0.0
  #define zero 0.0
  #define sqrabs(a) a*a
  #define kvector dvector
  #define freekvector freedvector
  #define kmatrix dmatrix
  #define freekmatrix freedmatrix
#endif

// Various definitions 
#define SPIN_0_UP ((unsigned long long) 1)
#define SQR(a) ((a)*(a))
#define RANDOM ( (double) 2*rand()/RAND_MAX - 1 )
#define PI 3.14159265358979323846
#define LARGE_NUMBER 9.999E+99
#define SMALL_NUMBER 1E-6
#define float double

/* Output definitions */
#define FILEEND ".dat"
#define LOGFILEEND ".log"
#define SZZEND ".szz" //only used of FIND_CROSS
#define SXXEND ".sxx" //only used of FIND_CROSS
#define SYYEND ".syy" //only used of FIND_CROSS
#define SPMEND ".spm" //only used of FIND_CROSS and FIND_CROSS_PM
#define SMPEND ".smp" //only used of FIND_CROSS and FIND_CROSS_PM

#define MATRIXFILENAME "SqMat"
#define COEND ".gs"
#define GSCOEND ".gs"
#define UNIEND ".un"
#define UNIOBSEND ".uno"

// Flags 
#define START 1
#define STOP 0
#define X 0
#define Y 1
#define Z 2
#define REAL 0
#define IMAG 1
#define NORMAL 0
#define RECONSTRUCT 1
#define CROSS 2
#define SZZ 0
#define SPM 1
#define SMP 2
//parallelization flags
#define MODEN 0 //Normal mode
#define MODEW 1 //Write Hamiltonian to file
#define MODEGS 2 //Groundstate mode
#define MODERC 3//Reconstruct mode
#define MODEQ 4 //Cross section mode
#define UNIMODEN 0 //Normal mode
#define UNIMODEW 1 //Write unique
#define UNIMODER 2 //Read unique

// Program pieces 
#define TLOOP_BEGIN for(sym=0;sym<Nsym;sym++) T[sym]=0;       \
                    sym=0, new_state=state;                         \
                    while(T[Nsym-1]<Nsymvalue[Nsym-1]) {
#define TLOOP_END   new_state=SymOp(0,new_state);                         \
		    for (sym=0; ++T[sym]==Nsymvalue[sym]; )   \
	              if (sym<Nsym-1) { T[sym++]=0;           \
                                new_state=SymOp(sym,new_state);}}      /* while */
 /*TRANSLOOP: for 1D this will be a loop over Nspins. statrs at sym=1 to avoid spin-flip/identity symmetry*/
#define TRANSLOOP_BEGIN T[0]=1;     \
                    for(sym=1;sym<Nsym;sym++) T[sym]=0;       \
                    sym=1, new_state=state;                         \
                    while(T[Nsym-1]<Nsymvalue[Nsym-1]) {
#define TRANSLOOP_END   new_state=SymOp(1,new_state);                         \
		    for (sym=1; ++T[sym]==Nsymvalue[sym]; )   \
	              if (sym<Nsym-1) { T[sym++]=0;           \
                                new_state=SymOp(sym,new_state);}}      /* while */
#define QLOOP_BEGIN for(sym=0;sym<Nsym;sym++) q[sym]=0;       \
                    sym=0;                                    \
                    while(q[Nsym-1]<Nsymvalue[Nsym-1]) {
#define QLOOP_END   for (sym=0; ++q[sym]==Nsymvalue[sym]; ) { \
	              if (sym<Nsym-1) q[sym++]=0; } }  /* while */

// General debugging requests from RLexact.C

#define VERBOSE
#define VERBOSE_TIME_LV1
#define VERBOSE_TIME_LV2


// Debugging output requests from RLexact.c 
//#define TEST_GS_SEARCH 
//#define TEST_ALLOCATE
//#define MAIN_LOOP_MESSAGES

// Debugging output requests from RLhamil.c 
//#define TEST_HAMILTON  /*this tests an old version of Hamilton, no longer in use*/
//#define TEST_MAT_ELEM 
//#define TEST_HAM2  
//#define TEST_HAMZEE
//#define TEST_FILLHAM 

// Debugging output requests from RLsparse.C
//#define TEST_MAKESPARSE
//#define TEST_APPLYSPARSE
//#define TEST_APPLYSPARSE_READ
//#define TEST_APPLYSPARSE_DEEP
//#define DEBUG_APPLYSPARSE
//#define TEST_HAM4


// Debugging output requests from RLmatrix.c 
//#define MATRIX_MESSAGES 
//#define TEST_EIG 
//#define TEST_EVEC
//#define TEST_MATRIXMAG  
//#define TEST_CALC_M 
//#define TEST_CALC_SZZQ 

// Debugging output requests from RLlanczos.c 
//#define TEST_SEED 
//#define TEST_EIGENVECTORS
//#define TEST_FINDGROUND
//#define TEST_STATES 
//#define TEST_ENERGIES
//#define TEST_LANCZOS_VECTOR
//#define LANCZOS_MESSAGES
//#define VERBOSE_LANCZOS
//#define TEST_FINDMAG
//#define DEBUG_SEED // WARNING: May cause premature stopping of lanczos algorithm
//#define DEBUG_LANCSTEP
#define TEST_LANCCROSS

// Debugging output requests from RLtables.c 
//#define TEST_TABLES    
//#define TEST_INVERTMATRIX
//#define TEST_OCC
//#define TEST_ISUNIQUE
//#define TEST_FINDUNIQUE  
//#define TEST_FINDUNIQUE_DETAIL
#define TEST_FILLUNIQUE 
//#define TEST_READUNIQUE 
//#define TEST_FILLUNIQUE_LIST
//#define TEST_FILLUNIQUEOBSERVABLES
//#define TEST_COUNT 
//#define TEST_LOOKUP   
//#define TEST_CHECKPHASE 

// Debugging output requests from RLsymm.c 
//#define TEST_SYM
//#define TEST_SYMCOUPLING   

// Debugging output requests from RLio.c and RLutil.c
#define TEST_WRITECROSS
//#define TEST_SPINFLIP
#define TEST_INPUT
//#define TEST_ROTATION

// Debugging output requests from RLcross.c 
//#define PRINT_STATES 
//#define PRINT_ENERGIES  
//#define TEST_APPLYSMP
//#define TEST_APPLYSZQ 
#define TEST_CROSS


// Debugging output requests from Diagonal.C 
//#define TEST_DIAGONALIZE 

// Debugging output requests from regcc.cpp
//#define TEST_MULTIMATCH
//#define TEST_FILEREAD
//

// Found in RLexact.h Previously named MEMORY_DEBUG, but not called anywhere.
//#define TEST_MEMORY
