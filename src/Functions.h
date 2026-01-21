
#include "RLexact.h"



#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Functions in RLio
void WritehmQ(long long *, struct FLAGS *);
void WriteResults(long long, struct FLAGS *);
void WriteCross(long long, long long *, long long, struct FLAGS*);



// Functions in RLcross
void CrossLanczos(long long *, struct FLAGS*);
void ApplySzq(long long *);
double lengthofvector(komplex *);
void ApplySmp(long long *, long long, komplex *);



// Functions in RLsparse
void MakeSparse(struct FLAGS*);
void ApplySparse(komplex *vectin, komplex *vectout, long long *k, struct FLAGS*);
void FillHamilSparse(komplex **hamil, long long *k, struct FLAGS*);
void WriteCouplingFiles(unsigned long long, unsigned long long, int *, 
    komplex, long long, long long, int *, long long *, 
    FILE *, FILE *, FILE *, struct FLAGS*);


// Functions in RLlancz
double NextLanczos(komplex *, komplex *, komplex *,
                   unsigned long long, long long *, struct FLAGS*);
long long LanczosLoop(long long, long long *, komplex *, struct FLAGS*);
double LowestLanczos(long long *, komplex *, long long *, long long, struct FLAGS*);



// Functions in RLmatrix
double Matrix_gs(komplex **, long long *, long long *, komplex *, struct FLAGS*);


// Functions in RLhamil
void Hamil2_sparse(unsigned long long, unsigned long long *, long long, long long, 
    int *, long long *, int *, komplex *, 
    double *, FILE *, FILE *, FILE *, struct FLAGS*);
void Hamil_Zeeman(unsigned long long, unsigned long long *, long long, 
    int *, long long *, int *, komplex *, double *, 
    FILE *, FILE *, FILE *, struct FLAGS*);



#endif