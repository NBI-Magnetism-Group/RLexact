/* Program file RLio.c - 
* Processes all user input/output
* and some geometric conversions
*
* Last change: SJ 19.12.16
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

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex>
#include <cnr.h>
#include <RLexact.h>

/* Functions defined in this file */
long long intro();
long long ReadCoupPattern(char *);
#ifndef M_SYM
void TransformCoup(long long);
#endif /* M_SYM */
#ifdef FIND_MAG
void WriteMaggs(long long *);
#endif //FIND_MAG
void time_stamp(time_t *, long long, const char *);
void outro();
void fatalerror(const char *, long long);
void Warning(const char *, long long);
void LogMessageChar(const char *);
void OutMessageChar(const char *);
void LogMessageInt(long long);
void LogMessageImag(long long);
void LogMessageCharDouble(const char *, double);
void LogMessageCharInt(const char *, long long);
void OutMessageCharInt(const char *, long long);
void LogMessageChar3Vector(const char *, double, double, double);
void WriteState(const char *, komplex *);
void WriteStates(komplex **);
void WritehmQ(long long *);
void WriteResults(long long);
void WriteCross(long long, long long*, long long);
void WriteGSEnergy(komplex);
void WriteEnergy(double);
void WriteQvalue(long long*);
void WriteGSdata(double,long long*);
void WriteGSstate(komplex*);
void ReadGSdata(double*,long long*,komplex*);
void ReadGSenergy(double*, long long*);

/* Functions in RLutils.C */
extern void NormalizeVector(double *);
extern void RotateVector(double *);
extern void Bubblesort(double *, double *, long long);
extern void FillRotationMatrix(double *);

/* Functions defined elsewhere */
extern void MakeSymCoup();
extern void itoa(long long,char[]);
extern void filereader(char*,char*,long long);
extern long long filesizer(char*);
extern long long multimatch(char*,long long,const char*,double**,long long*,long long);
extern long long multimatch(char*,long long,const char*,long long**,long long*,long long);
extern long long matchlines(char*,const char*,double*,bool);
extern long long matchlines(char*,const char*,long long*,bool);
extern double* dvector(long long,long long);

/* Global variables defined in RLexact.c */
extern long long    Nspins,Nsym,Nsymadd,Nunique,Nuniq_k;
extern long long    **symadd;
extern long long    Ndimensions;
extern long long    *TransIds;
extern long long    Trans_Qmax[3];
extern long long    hamil_coup[NCOUP][2], Ncoup; 
#ifdef RING_EXCHANGE
extern long long    ring_coup[NCOUP][4], Nring; 
#endif /* RING_EXCHANGE */
extern komplex*    gs;
#ifdef M_SYM
extern long long    twom;
extern double  mstart, mend;
#else
extern double h, hstart, hend, hstep; 
extern double field[3];
#endif /* M_SYM */
extern double Jzz[NCOUP],Jxy[NCOUP],Janis[NCOUP];
#ifdef RING_EXCHANGE
extern double       Jr[NCOUP];
#endif /* RING_EXCHANGE */
extern double Jdip[NCOUP],geom_13[NCOUP],r_vector[NCOUP][3];
extern long long symlist[NSYM];   
extern double *energies;
#ifdef FIND_MAG
extern double *magnetisation;
extern double maggs;
#endif /* FIND_MAG */
extern long long Nq_choice;
extern long long **q_choice;
#ifdef FIND_CROSS
extern FILE *outfilezz;
#ifndef FIND_CROSS_PM
extern FILE *outfilexx, *outfileyy;
#endif
#ifdef FIND_CROSS_PM
extern FILE  *outfilepm, *outfilemp;
#endif
extern double *cross;
#endif /* FIND_CROSS */
#ifdef MOTIVE
extern long long Nspins_in_uc;
extern float **spin_positions;
#endif //MOTIVE
extern FILE *gscoinfile;
extern FILE *outfile;
extern FILE *logfile;
extern long long q_gs[NSYM];
extern double gs_energy;
extern char *infile_name;
extern bool name_on_commandline;
extern double szqlength;
extern long long mode;
extern long long unimode;
extern int spinflip_number;
extern int spinflip_present;
/* Lanczos variables defines in RLlanz.C */
extern double Ritz_conv;

/* Regional variables defined here */
FILE *infile;
FILE *gscofile;
char* filedata;
long long filesize;

/************************************************/
long long intro()
 {
  /* Introduce, read input and perform file handling */
  char outfile_name[30];
  char logfile_name[30];
  char gscofile_name[30];
  char gscoinfile_name[30];
  char qstr[2];

  OutMessageChar("Welcome to the Exact Diagonalization Program, RLexact \n");
#ifdef MATRIX
  OutMessageChar(" Matrix"); 
#endif  /* MATRIX */
#ifdef LANCZOS
  OutMessageChar(" Lanzcos "); 
#endif  /* LANCZOS */ 
  OutMessageCharInt(" diagonalization. BUFFERSIZE =",BUFFERSIZE);
#ifdef M_SYM
  OutMessageChar(" M-symmetry present. \n");
#else
 OutMessageChar(" M-symmetry absent. \n");
#endif /* M_SYM */
  OutMessageChar(" Observables:");
#ifdef FIND_MAG
  OutMessageChar(" Magnetization,");
#endif /* FIND_MAG */
#ifdef CROSS
  OutMessageChar(" S^zz(q,w),");
  #ifndef M_SYM
  #ifndef FIND_CROSS_PM
  OutMessageChar(" S^xx(q,w), S^yy(q,w),");
  #endif
  #ifdef FIND_CROSS_PM
  OutMessageChar(" S^+-(q,w), S^-+(q,w),");
  #endif
  #endif
  #endif  /* CROSS */
  OutMessageChar(" Energy.\n");
  OutMessageChar(" For more information, see the manual.\n");
  if (!name_on_commandline) {
    OutMessageChar(" Please type the name of the input file : ");
    scanf("%s",infile_name);
    OutMessageChar(" ... \n");
  }

  strcpy(outfile_name,infile_name);
  strcat(outfile_name,FILEEND);
  errno=0;
  outfile=fopen(outfile_name,"w");
  if (outfile == NULL)
    {
      fatalerror("Cannot open general output file, sorry!",errno);
      return -1;
    }
  fflush(outfile);
  strcpy(logfile_name,infile_name);
  strcat(logfile_name,LOGFILEEND);
  errno=0;
  logfile=fopen(logfile_name,"w");
  if (logfile == NULL)
    {
      fatalerror("Cannot open logfile, sorry!",errno);
      return -1;
    }

  ReadCoupPattern(infile_name); // this function does the actual file-input handling

#ifdef FIND_CROSS

  if (mode==MODEQ || mode==MODERC) {
    strcpy(gscoinfile_name,infile_name);
    strcat(gscoinfile_name,GSCOEND);
    errno=0;
    gscoinfile=fopen(gscoinfile_name,"r");
    if (gscoinfile == NULL)
      {
	fatalerror("Cannot open gs file, sorry!",errno);
	return -1;
      }
    fflush(gscoinfile);
  }

  if (mode==MODEGS) {
    strcpy(gscofile_name,infile_name);
    strcat(gscofile_name,COEND);
    errno=0;
    gscofile=fopen(gscofile_name,"w");
    if (gscofile == NULL)
      {
	fatalerror("Cannot open carry over file, sorry!",errno);
	return -1;
      }
    fflush(gscofile);
  }
  if (mode==MODERC) {
    strcpy(gscofile_name,infile_name);
    strcat(gscofile_name,COEND);
    errno=0;
    gscofile=fopen(gscofile_name,"a");
    if (gscofile == NULL)
      {
	fatalerror("Cannot open carry over file, sorry!",errno);
	return -1;
      }
    fflush(gscofile);
  }

  strcpy(outfile_name,infile_name);
  strcat(outfile_name,SZZEND);
  errno=0;
  if ((outfilezz=fopen(outfile_name,"w")) == NULL)
   {
    fatalerror("Cannot open output file for Szz, sorry!",errno);
    return -1;
   }
  fflush(outfilezz);

#ifndef FIND_CROSS_PM
  strcpy(outfile_name,infile_name);
  strcat(outfile_name,SXXEND);
  errno=0;
  if ((outfilexx=fopen(outfile_name,"w")) == NULL)
   {
    fatalerror("Cannot open output file for S+-, sorry!",errno);
    return -1;
   }
  fflush(outfilexx);

  strcpy(outfile_name,infile_name);
  strcat(outfile_name,SYYEND);
  errno=0;
  if ((outfileyy=fopen(outfile_name,"w")) == NULL)
   {
    fatalerror("Cannot open output file for S-+, sorry!",errno);
    return -1;
   }
  fflush(outfileyy);
#endif

#ifdef FIND_CROSS_PM
  strcpy(outfile_name,infile_name);
  strcat(outfile_name,SPMEND);
  errno=0;
  if ((outfilepm=fopen(outfile_name,"w")) == NULL)
   {
    fatalerror("Cannot open output file for S+-, sorry!",errno);
    return -1;
   }
  fflush(outfilepm);

  strcpy(outfile_name,infile_name);
  strcat(outfile_name,SMPEND);
  errno=0;
  if ((outfilemp=fopen(outfile_name,"w")) == NULL)
   {
    fatalerror("Cannot open output file for S-+, sorry!",errno);
    return -1;
   }
  fflush(outfilemp);
#endif

#endif /* FIND_CROSS */
#ifdef TEST_INPUT
  LogMessageChar("All filenames are OK. \n");
#endif /* TEST_INPUT */


 
 return 0;  /* All is OK */
 }

/* ----------------------------------------------------------------------- */

long long ReadCoupPattern(char *filename)
 {
 /* Read all couplings (pairs, types and strengths) 
    and symmetry numbers from file */
  long long Ncoupstr;        /* Number of different coupling strengths */
#ifdef RING_EXCHANGE
  long long Nringstr;
  double hamring[NRINGSTR];
#endif /* RING_EXCHANGE */
  long long n1,n2,str,i,j, symconstruct;
  double hamzz[NCOUPSTR],hamxy[NCOUPSTR],hamanis[NCOUPSTR];
#ifdef DIPOLE
  double r,rx,ry,rz;
  double Dip[NCOUPSTR];
#endif /* DIPOLE */
#ifndef M_SYM
  double B2;
#endif /* M_SYM */
  char test[80];
  long long filesize; // of inputfile

  filesize = filesizer(filename);
  char *filedata = (char*)malloc(filesize*sizeof(char));
  filereader(filename,filedata,filesize); // the entire file is now in filedata

#ifdef TEST_INPUT
  LogMessageChar("Input file opened...\n");
#endif /* TEST_INPUT */

  matchlines(filedata, "Number of spins", &Nspins, true);
#ifdef TEST_INPUT
  LogMessageCharInt(" Nspins:",Nspins);
  LogMessageChar("\n");
#endif /* TEST_INPUT */

// *********** Mandatory input: Mode *************
  matchlines(filedata, "Mode", &mode, true);
#ifdef TEST_INPUT
  switch(mode) {
  case MODEN:
    LogMessageChar("Mode is: NORMAL MODE. \n");
    break;
  case MODEGS:
    LogMessageChar("Mode is: GS MODE. \n");
    break;
  case MODEQ:
    LogMessageChar("Mode is: Q MODE. \n");
    break;
  default:
    fatalerror("Unknown mode",mode);
  }
#endif /* TEST_INPUT */

// *********** Mandatory input: Unique mode *************
  matchlines(filedata, "Unimode", &unimode, true);
#ifdef TEST_INPUT
  switch(unimode) {
  case UNIMODEN:
    LogMessageChar("Unique mode is: NORMAL MODE. \n");
    break;
  case UNIMODEW:
    LogMessageChar("Unique mode is: WRITE MODE. \n");
    break;
  case UNIMODER:
    LogMessageChar("Unique mode is: READ MODE. \n");
    break;
  default:
    fatalerror("Unknown unique mode",unimode);
  }
#endif /* TEST_INPUT */


// ******** Mandatory input: Read number of lines to input **********
  matchlines(filedata, "Number of couplings", &Ncoup, true);
  matchlines(filedata, "Number of coupling strength", &Ncoupstr, true);
#ifdef RING_EXCHANGE
  matchlines(filedata, "Number of rings", &Nring, true);
  matchlines(filedata, "Number of ringstrength", &Nringstr, true);
#endif /* RING_EXCHANGE */
#ifdef TEST_INPUT
  LogMessageCharInt("Ncoup:", Ncoup);
  LogMessageCharInt("Ncoupstr:", Ncoupstr);
  LogMessageChar("\n");
#ifdef RING_EXCHANGE
  LogMessageCharInt("Nring:", Nring);
  LogMessageCharInt("Nringstr:", Nringstr);
  LogMessageChar("\n");
#endif /* RING_EXCHANGE */
#endif /* TEST_INPUT */

  matchlines(filedata, "Number of hardcoded symmetries", &Nsym, true);
  matchlines(filedata, "Number of custom symmetries", &Nsymadd, true);
  matchlines(filedata, "Construct symmetries", &symconstruct, true);
  if ((Nsym+Nsymadd) <=0)
   {
    fatalerror("At least one symmetry must be defined; use for instance IDENTITY",IDENTITY);
   }
#ifdef TEST_INPUT
  LogMessageCharInt("Number of hardcoded symmetries: Nsym=", Nsym);
  LogMessageCharInt("\nNumber of custom symmetries: Nsymadd=", Nsymadd);
  LogMessageCharInt("\nsymconstruct:",symconstruct);
  LogMessageChar("\n");
#endif /* TEST_INPUT */ 
  
// ********* Input symmetry info ********************
  matchlines(filedata, "Hardcoded symmetries", symlist, false);
#ifdef TEST_INPUT
  LogMessageCharInt("Number of hardcoded symmetries scanned:",Nsym); 
  LogMessageChar("\n");
#endif /* TEST_INPUT */   
  long long* dummy = (long long*) malloc(Nsymadd*sizeof(long long));
  symadd = (long long**) malloc(Nsymadd*sizeof(long long*));
 multimatch(filedata,filesize,"Custom symmetry", symadd, dummy, Nsymadd);
 
 matchlines(filedata, "Number of dimensions", &Ndimensions, true);
 TransIds = (long long*) malloc(3*sizeof(long long));
 for (int i = 0;i<3;i++) TransIds[i] = 0; //Ugly,but needed for loop in RLcross.
 matchlines(filedata, "Translation indices", TransIds, 1);


#ifdef TEST_INPUT
 LogMessageCharInt("Number of dimensions:", Ndimensions);
 LogMessageChar("\n");
 LogMessageChar("Translations at indices:");
 for (int i=0; i<Ndimensions;i++) LogMessageInt(TransIds[i]);
 LogMessageChar("\n");
#endif
 

  #ifdef TEST_SPINFLIP
    LogMessageCharInt("Testing for spin flip sym. Matching symlist with SPIN_FLIP=",SPIN_FLIP);
    LogMessageChar("\n symlist=");
  #endif

  for (i=0; i<Nsymadd; i++) // Give numbers to custom symmetries
  {  
    symlist[Nsym+i]=NSYM+i;

    #ifdef TEST_SPINFLIP
      LogMessageCharInt("\nsymlist[Nsym+i]=",symlist[Nsym+i]);
    #endif
  }  

  #ifdef TEST_SPINFLIP
    LogMessageChar("\n");
  #endif


  spinflip_number = -1;
  spinflip_present = 0;
  for (i=0; i<Nsym; i++) // Look for spin flip
  {  
    #ifdef TEST_SPINFLIP
      LogMessageCharInt("\nMatch SPIN_FLIP with symlist[Nsym+i]=",symlist[Nsym+i]);
    #endif
    if (symlist[i] == SPIN_FLIP){
      spinflip_present = 1;
      spinflip_number = i;
        #ifdef TEST_SPINFLIP
          LogMessageCharInt("Spinflip found! at i =", i);
          LogMessageChar("\n");
        #endif
    }

  }
    #ifdef TEST_SPINFLIP
      LogMessageCharInt("\nRLio: spinflip_number=",spinflip_number);
    #endif  
  free(dummy);
// ************ Input GS q-values ******************

  matchlines(filedata, "Number of chosen GS q-values", &Nq_choice, true);
#ifdef TEST_INPUT
  LogMessageCharInt("Nq_choice:",Nq_choice); 
#endif /* TEST_INPUT */   
  dummy = (long long*) malloc(Nq_choice*sizeof(long long));
  q_choice = (long long**) malloc(Nq_choice*sizeof(long long*));

 multimatch(filedata,filesize,"Chosen GS q-value", q_choice, dummy, Nq_choice);
#ifdef TEST_INPUT
  for (i=0; i<Nq_choice; i++)
    {
      LogMessageChar("q_choice: ( ");
      for (j=0; j<Nsym+Nsymadd; j++){
	LogMessageCharInt(" ",q_choice[i][j]); 
      }
    LogMessageChar(") \n");
    }
#endif /* TEST_INPUT */   
  free(dummy);
  
#ifdef MOTIVE
  
  matchlines(filedata, "Number of spins in unit cell", &Nspins_in_uc, true);

#ifdef TEST_INPUT
  LogMessageCharInt("Number of spins in unit cell: ", Nspins_in_uc);
  LogMessageChar("\n");
#endif //TEST_INPUT
  
  dummy = (long long*) malloc(Nspins_in_uc*sizeof(long long));
  spin_positions = (double**)malloc(Nspins_in_uc*sizeof(double*));
  multimatch(filedata, filesize, "Relative position", spin_positions, 
      dummy, Nspins_in_uc);
  free(dummy);

#ifdef TEST_INPUT
  for (i=0; i<Nspins_in_uc;i++){
      LogMessageChar3Vector("Spin pos.",spin_positions[i][X],
                                        spin_positions[i][Y],
                                        spin_positions[i][Z]);
  }
#endif //TEST_INPUT


  matchlines(filedata, "Qmax translation", Trans_Qmax, true);
#ifdef TEST_INPUT
 LogMessageChar3Vector("Qmax translation", Trans_Qmax[X],
                                           Trans_Qmax[Y],
                                           Trans_Qmax[Z]);

#endif



#endif //MOTIVE


#ifdef M_SYM
  matchlines(filedata, "M start", &mstart, true);
  matchlines(filedata, "M end", &mend, true);
#ifdef TEST_INPUT
  LogMessageCharDouble("mstart:", mstart);
  LogMessageCharDouble("mend:", mend);
  LogMessageChar("\n");
#endif /* TEST_INPUT */        
#else
  matchlines(filedata, "H start", &hstart, true);
  matchlines(filedata, "H end", &hend, true);
  matchlines(filedata, "H step", &hstep, true);
#ifdef TEST_INPUT
  LogMessageCharDouble("hstart:", hstart);
  LogMessageCharDouble("hend:", hend);
  LogMessageCharDouble("hstep:", hstep);
  LogMessageChar("\n");
#endif /* TEST_INPUT */    
  matchlines(filedata, "Hx", &field[X], true);
  matchlines(filedata, "Hy", &field[Y], true);
  matchlines(filedata, "Hz", &field[Z], true);
#ifdef TEST_INPUT
  LogMessageChar3Vector("Field direction:", field[X], field[Y], field[Z]);
  LogMessageChar("\n");
#endif /* TEST_INPUT */     
  for (B2=0, j=X; j<=Z; j++)
    B2+=SQR(field[j]);
  if (B2 > 0)
    for (j=X; j<=Z; field[j++]/=sqrt(B2));
  else
    fatalerror(" Field direction not well defined",0);
#ifdef TEST_INPUT
  LogMessageChar3Vector("Normalized field direction:",  field[X],
                                                        field[Y],
                                                        field[Z]);
  LogMessageChar("\n");
#endif /* TEST_INPUT */    
  //FillRotationMatrix(field); //doesnt work and unnecessary, SJ 20.02.17
#endif /* M_SYM */

#ifdef LANCZOS
// ********* Lanzcos-Mandatory input : Read Lanzcos numbers *********
  matchlines(filedata, "Ritz_conv", &Ritz_conv, true);
  #ifdef TEST_INPUT
  LogMessageCharDouble("Ritz_conv:", Ritz_conv);
  LogMessageChar("\n");
#endif /* TEST_INPUT */
#endif /* LANCZOS */

  dummy= (long long*) malloc(Ncoupstr*sizeof(long long));
  double** dummyresdouble = (double**) malloc(Ncoupstr*sizeof(double*));
  multimatch(filedata,filesize,"Coupling strength vector", dummyresdouble,
                                                           dummy, Ncoupstr);

  for (long long k=0;k<Ncoupstr;k++) {
    hamzz[k]=dummyresdouble[k][0];
    hamxy[k]=dummyresdouble[k][1];
    hamanis[k]=dummyresdouble[k][2];
#ifdef DIPOLE
    Dip[k]=dummyresdouble[k][3];
#endif /* DIPOLE */
#ifdef TEST_INPUT
    LogMessageChar3Vector(" Jzz, Jxy, Janis:",hamzz[k],hamxy[k],hamanis[k]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }
  free(dummy);
  free(dummyresdouble);

#ifdef RING_EXCHANGE
  dummy= (long long*) malloc(Nringstr*sizeof(long long));
  double** dummyresdouble1 = (double**) malloc(Nringstr*sizeof(double*));
  multimatch(filedata,filesize,"Ring strength", dummyresdouble1, dummy,
                                                                 Nringstr);

  for (long long k=0;k<Nringstr;k++) {
    hamring[k]=dummyresdouble1[k][0];
#ifdef TEST_INPUT
    LogMessageCharDouble(" Jr:", hamring[k]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }
  free(dummy);
  free(dummyresdouble1);
#endif /* RING_EXCHANGE */

  dummy = (long long*) malloc(Ncoup*sizeof(long long));
  long long **dummyres = (long long**) malloc(Ncoup*sizeof(long long*));
  multimatch(filedata,filesize,"Coupling vector",dummyres , dummy, Ncoup);

  for (long long k=0;k<Ncoup;k++) {
    hamil_coup[k][0]=dummyres[k][0];
    hamil_coup[k][1]=dummyres[k][1];
    Jzz[k]=hamzz[dummyres[k][2]];
    Jxy[k]=hamxy[dummyres[k][2]];
    Janis[k]=hamanis[dummyres[k][2]];    
#ifdef TEST_INPUT
    LogMessageCharInt(" Coupling from ",dummyres[k][0]);
    LogMessageCharInt(", to ",dummyres[k][1]);
    LogMessageCharInt(", strength ",dummyres[k][2]);
    LogMessageChar("\n");
    LogMessageCharDouble(" Jzz ",Jzz[k] );
    LogMessageCharDouble(", Jxy ", Jxy[k]);
    LogMessageCharDouble(", Janis ",Janis[k]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
    /*#ifndef M_SYM //NEVER: rotating J is the wrong way!
    TransformCoup(k);
#endif /* M_SYM */

#ifdef DIPOLE
    Jdip[k]=Dip[dummyres[k][2]];
    r_vector[k][X]=dummyres[k][3];
    r_vector[k][Y]=dummyres[k][4];
    r_vector[k][Z]=dummyres[k][5];
    NormalizeVector(r_vector[k]);
    RotateVector(r_vector[k]);
    geom_13[k] = 1.0-3.0*SQR(r_vector[k][Z]);
#ifdef TEST_ROTATION
  LogMessageChar3Vector("   transformed direction: ", r_vector[k][X],
                                                      r_vector[k][Y],
                                                      r_vector[k][Z]);
  LogMessageChar("\n");
#endif /* TEST_ROTATION */
#endif /* DIPOLE */
  }

#ifdef RING_EXCHANGE
  dummy = (long long*) malloc(Nring*sizeof(long long));
  long long **dummyresring = (long long**) malloc(Nring*sizeof(long long*));
  multimatch(filedata,filesize,"Coupling ring vector", dummyresring, dummy,
                                                                     Nring);

  for (long long k=0;k<Nring;k++) {
    ring_coup[k][0]=dummyresring[k][0];
    ring_coup[k][1]=dummyresring[k][1];
    ring_coup[k][2]=dummyresring[k][2];
    ring_coup[k][3]=dummyresring[k][3];
    Jr[k]=hamring[dummyresring[k][4]];
#ifdef TEST_INPUT
    LogMessageCharInt(" Coupling from:", ring_coup[k][0]);
    LogMessageCharInt(", over:", ring_coup[k][1]);
    LogMessageCharInt(", and:", ring_coup[k][2]);
    LogMessageCharInt(", to:", ring_coup[k][3]);
    LogMessageCharInt(", and back to:", ring_coup[k][0]);
    LogMessageCharDouble(" Jr: ",Jr[k]);
    LogMessageChar("\n");
#endif /* TEST_INPUT */
  }
  free(dummy);
  free(dummyresring);
#endif /* RING_EXCHANGE */

  if (symconstruct==1)
    MakeSymCoup();

  LogMessageChar("RLio.C successful. \n");
  free(filedata);
  return 0;
 }

#ifndef NEVER //it doesnt work currently: rotating J is the wrong way!
void TransformCoup(long long j)
{
  double vec[3];

  vec[Z]=Jzz[j];
  vec[X]=Jxy[j]+Janis[j];
  vec[Y]=Jxy[j]-Janis[j];
  RotateVector(vec);
  Jzz[j]=vec[Z];
  Jxy[j]=(vec[X]+vec[Y])/2.0;
  Janis[j]=(vec[X]-vec[Y])/2.0;
  LogMessageCharDouble("Now Jz =",Jzz[j]);
  LogMessageCharDouble(", Jxy =",Jxy[j]);
  LogMessageCharDouble(", Janis =",Janis[j]);
  LogMessageChar("\n");
  return;
}
#endif /* M_SYM */

#ifdef FIND_CROSS
#ifdef LANCZOS
long long SortCross (long long Nener) {
  /* Finds and sums crosssections for same energies, and sorts energies and
   * cross by energy.
   * Could be optimized greatly, but MAX_LANCZOS is small: no point */
  long long XS = 0; // number of unique crosssections
  bool newener;

  double *tmpener = dvector(0,MAX_LANCZOS-1);
  double *tmpcross = dvector(0,MAX_LANCZOS-1);

  for (long long i=0;i<Nener;i++) {
    newener = true;
    for (long long j=0;j<XS;j++) {
      if (energies[i]==tmpener[j]) {
      	newener = false;
      	tmpcross[j]+=cross[i];
      }
    }
    if (newener) {
      tmpener[XS]=energies[i];
      tmpcross[XS++]=cross[i];
    }
  }

  Bubblesort(tmpener,tmpcross,XS); // Maybe not necessary, but nice

  for (long long i=0;i<XS;i++) {
    cross[i]=tmpcross[i];
    energies[i]=tmpener[i];
  }
  freedvector(tmpener,0,MAX_LANCZOS-1);
  freedvector(tmpcross,0,MAX_LANCZOS-1);
  return XS;
}
#endif /* LANZCOS */
#endif /* FIND_CROSS */



void time_stamp(time_t *tim, long long flag, const char *string)
 { 
  if (flag == START)
  {
   *tim = -clock();
   fprintf(logfile,"Starting %s ... ", string);
  }
  if (flag == STOP)
  {
   *tim += clock();
   fprintf(logfile," %s done, %lg seconds used\n",
            string, *tim/(double)CLOCKS_PER_SEC);
  } 

  return;
 }

void outro()
 {
  fclose(outfile);
#ifdef FIND_CROSS
  fclose(outfilezz);

  #ifndef FIND_CROSS_PM
  fclose(outfilexx);
  fclose(outfileyy);
  #endif
  #ifdef FIND_CROSS_PM
  fclose(outfilepm);
  fclose(outfilemp);
  #endif
#endif /* FIND_CROSS */
  LogMessageChar("\n End of diagonalization program RLexact.\n");
  printf("\n End of diagonalization program RLexact.\n");

  return;
 }

void fatalerror(const char *str, long long i)
 {
  fprintf(logfile, "*** FATAL ERROR *** %s %lld %s \n",str,i,strerror(i));
  printf("*** FATAL ERROR *** %s %lld %s \n",str,i, strerror(i));
  exit(-1);
 }

void Warning(const char *str, long long i)
 {
  fprintf(logfile, "*** WARNING *** %s %lld \n",str,i);
  fflush(logfile);
  return;
 }

void LogMessageInt(const long long i)
 {
  fprintf(logfile, " %lld ",i);
  fflush(logfile);
  return;
 }

void LogMessageImag(const double a, const double b)
 {
  fprintf(logfile, " (  %lg + %lg i )",a,b);
  fflush(logfile);
  return;
 }

void LogMessageKomplex(const komplex z)
 {
  fprintf(logfile, " (  %lg + %lg i )",real(z),imag(z));
  fflush(logfile);
  return;
 }

void LogMessageChar(const char *str)
 {
  fprintf(logfile, " %s ",str);
  fflush(logfile);
  return;
 }
void OutMessageChar(const char *str)
 {
  fprintf(stdout, " %s ",str);
  return;
 }

void LogMessageCharDouble(const char *str, double d)
 {
  fprintf(logfile, " %s %g ",str,d);
  fflush(logfile);
  return;
 }

void LogMessageCharKomplex(const char *str, komplex z)
{
  fprintf(logfile, " %s (  %lg + %lg i )",str, real(z),imag(z));
  fflush(logfile);
}

void LogMessageCharInt(const char *str, long long i)
 {
  fprintf(logfile, " %s %lld ",str,i);
  fflush(logfile);
  return;
 }
void OutMessageCharInt(const char *str, long long i)
 {
  fprintf(stdout, " %s %lld ",str,i);
  return;
 }

void LogMessageChar3Vector(const char *str, double d1, double d2, double d3)
 {
  fprintf(logfile, " %s (%g, %g, %g ) \n",str,d1,d2,d3);
  fflush(logfile);
  return;
 }

/* Overload for longs */
void LogMessageChar3Vector(const char *str, 
                          long long l1,long long l2, long long l3)
{
  fprintf(logfile, " %s (%lld, %lld, %lld ) \n",str,l1,l2,l3);
  fflush(logfile);
  return;
 }

void WriteState(const char *msg, komplex *state)
 {
   /* Output one state vector -- with a text message */
  long long i;

  fprintf(outfile,"%s \n",msg);
  for(i=0;i<Nunique;i++)
    fprintf(outfile,"  %g + i * %g \n",real(state[i]),imag(state[i]));

  return;
 }

void WriteStates(komplex **hamil)
 {
   long long i,j;

   fprintf(outfile,"{\n");
   for (i=0;i<Nuniq_k;i++)
   {
    fprintf(outfile,"( ");
    for (j=0; j<Nuniq_k; j++)
      fprintf(outfile,"( %lg +i %lg ), ", real(hamil[j+1][i+1]),
                                          imag(hamil[j+1][i+1]));
    fprintf(outfile,") \n");
   }
  
   fprintf(outfile,"}\n");
  return;
 }

void WriteGSdata(double energy, long long *symmetry)
{
  fwrite(&energy,sizeof(double),1,gscofile);
  fwrite(symmetry,sizeof(long long),Nsym,gscofile);
  return;
}
void WriteGSstate(komplex *state)
{
  fwrite(state,sizeof(komplex),Nunique,gscofile);
  return;
}

void ReadGSdata(double* energy, long long *symmetry, komplex *state)
{
  //  gscoinfile=fopen("./12-0.gs","r");
  fread(energy,sizeof(double),1,gscoinfile);
  fread(symmetry,sizeof(long long),Nsym,gscoinfile);
  fread(state,sizeof(komplex),Nunique,gscoinfile);
  return;
}

void ReadGSenergy(double* energy, long long *symmetry)
{
  //  gscoinfile=fopen("./12-0.gs","r");
  fread(energy,sizeof(double),1,gscoinfile);
  fread(symmetry,sizeof(long long),Nsym,gscoinfile);
  return;
}

void WritehmQ(long long *q)
 {
  long long i;

  fprintf(outfile,"!");
#ifdef M_SYM
  fprintf(outfile," m= %g ; ",double(twom)/2);
#else
  fprintf(outfile," h= %g ; ",h);
#endif
#ifdef TEST_WRITEHMQ
  printf("\n no of q values in dat files, Nsym=%d\n",Nsym);
#endif 
  fprintf(outfile,"q= ( ");
  for (i=0; i<Nsym; i++)
    fprintf(outfile,"%lld ",q[i]);
  fprintf(outfile," ) \n");

  return;
 }

void WriteResults(long long N)
 {
/* Output the energies and other observables of the eigenstates */
  long long i,q;
  #ifndef WRITE_MAGNETISATION 
  //otherwise the energy and magnetisation pairs get mixed up
    Bubblesort(energies,NULL,N);
  #endif

  fprintf(outfile, "[ \n");
  for (i=0; i<N; i++)
    {

      // throw away rounding errors
      if (fabs(energies[i])<SMALL_NUMBER) {
	energies[i]=0;
      }

      fprintf(outfile,"energy= %9.6g ",energies[i]);
#ifdef FIND_MAG
      // throw away rounding errors
      if (fabs(magnetisation[i])<SMALL_NUMBER) {
	     magnetisation[i]=0;
      }
      #ifdef WRITE_MAGNETISATION
      #ifdef MATRIX //only meaningful for matrix-calculations
      fprintf(outfile,", mag_z= %9.6g ",magnetisation[i]); 
      #endif //MATRIX
      #endif //WRITE_MAGNETISATION
#endif  /* FIND_MAG */

      fprintf(outfile, "\n");
    }
  fprintf(outfile, "] \n");
  fflush(outfile);

  return;
 }

#ifdef FIND_CROSS
#ifdef LANCZOS
void WriteCross(long long Nener, long long *symvalue, long long flag)
{
/* Output the cross sections of the ground state */

  FILE * crossfile;
  if (flag == 0 )  { crossfile = outfilezz; } //SZZ
  #ifndef FIND_CROSS_PM
  else if (flag == 1 )  { crossfile = outfilexx; } //SXX
  else if (flag == 2 )  { crossfile = outfileyy; } //SYY
  #endif
  #ifdef FIND_CROSS_PM
  else if (flag == 1 )  { crossfile = outfilemp; } //SMP
  else if (flag == 2 )  { crossfile = outfilepm; } //SPM
  #endif

  else {LogMessageChar("\nError in WriteCross flag!\n");}

#ifdef TEST_WRITECROSS
  LogMessageChar("In WriteCross, ");
  #ifdef M_SYM
  LogMessageCharInt(", m =",twom/2);
  #endif
  LogMessageCharDouble(", gs_energy=",gs_energy);

  LogMessageCharInt(", q_gs = (",q_gs[0]);
  for (long long i=1;i<Nsym;i++) LogMessageCharInt(",",q_gs[i]);
  LogMessageChar(")");
  for (int i = 0; i < Nsym; i++)
  {
    LogMessageCharInt("\n i =",i);
    LogMessageCharInt(", symvalue[i] =",symvalue[i]);
  }
  LogMessageChar("\n");
#endif

#ifdef CSVOUT
// Print header
#ifdef M_SYM
  fprintf(crossfile, "m,");
#else
  fprintf(crossfile, "h,");
#endif

  for (int i=0;i<Nsym;i++) fprintf(crossfile, "q%d,", i);
  fprintf(crossfile, "E,S\n");

//This doesn't sort the array as below, but this shouldn't be needed as
//CSV-out is mainly for loading into other software and plotting/looking
//at data there.
  for (long long i = 0;i<Nener;i++) 
  {
#ifdef M_SYM
  fprintf(crossfile,"%g,",double(twom)/2);
#else
  fprintf(crossfile,"%g,",double(h));
#endif
  for (int i=0;i<Nsym;i++) fprintf(crossfile,"%lld,", symvalue[i]-q_gs[i]);

  fprintf (crossfile,"%g,%g\n",energies[i]-gs_energy,cross[i]);
  }
#else

#ifdef M_SYM
  fprintf(crossfile, "[\n m= %g, q= (%lld",double(twom)/2,symvalue[0]-q_gs[0]);
  //fprintf(crossfile, "[ \n m= %g, q= (%lld",double(twom)/2,symvalue[0]);
#else
  fprintf(crossfile, "[ \n h= %g, q= (%lld",double(h),symvalue[0]-q_gs[0]);
  //fprintf(crossfile, "[ \n h= %g, q= (%lld",double(h),symvalue[0]-q_gs[0]);
#endif /* M_SYM */

  for (long long i=1;i<Nsym;i++) //print spin-wave q: q'=symvalue[i]-q_gs[i]
  {
    fprintf(crossfile, " %lld",symvalue[i]-q_gs[i]); //
    //fprintf(crossfile, ",%2lld",symvalue[i]); //
  }
  
  fprintf(crossfile, ")\n length of Szpm(q)|gs> = %g",szqlength);
  fprintf(crossfile, " S(q) = %g\n\n",2*PI*szqlength);

  #ifdef TEST_WRITECROSS
  for (int i = 0; i < Nener+1; i++)
  {
    LogMessageCharDouble("\nenergies[i] =",energies[i]);
    LogMessageCharDouble(", cross[i] =",cross[i]);
  }
  #endif

  long long Xsections = SortCross(Nener); 
  //destroys *energies and *Cross. Maybe not so smart?
  
  for (long long i = 0;i<Nener;i++) 
  {
    // for (long long i = 0;i<Nener;i++) {
    // fprintf (outfilezz," omega = %9.6g Szz = %9.6g\n",energies[i]-gs_energy,cross[i]);
    fprintf (crossfile," %9.6g %9.6g\n",energies[i]-gs_energy,cross[i]);
    //fprintf (crossfile," %9.6g %9.6g\n",energies[i],cross[i]);

  }
  fprintf (crossfile,"]\n\n");
#endif //CSVOUT
  return;
}
#endif /* LANCZOS */
#endif /* FIND_CROSS */

void WriteGSEnergy(komplex E)
 {
   printf("<e|H|e> = %g +i %g\n",real(E),imag(E));

  return;
 }

void WriteEnergy(double re)
 {
  printf(" E = %g \n",re);

  return;
 }

#ifdef FIND_MAG
void WriteMaggs(long long *q){
  int i; 
  WriteQvalue(q);
  fprintf(outfile, "\nMagnetization of groundstate %g\n",maggs);
  fprintf(outfile, "\n");
}
#endif  //FIND_MAG

void WriteQvalue(long long *qvec){
  int i;
  fprintf(outfile, "q-vector: (");
  fprintf(outfile," %lld",qvec[0]);
  for(i=0; i<Nsym; i++){
    fprintf(outfile,", %lld",qvec[i]);
  }
fprintf(outfile,")");
}
