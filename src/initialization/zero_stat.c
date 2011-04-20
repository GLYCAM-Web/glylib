#include <stats.h>
#include <mylib.h>
//#include "../inc/stats.h"
//#include "../inc/mylib.h"

//meanvar zero_meanvar();
meanvar zero_meanvar(){
meanvar M;
M.t='\0'; // type population (p) or sample (s)
M.n=0; // number of units in sample
M.m=0; // mean
M.v=0; // variance
M.s=0; // standard deviation
return M;
}

statsarray zero_statsarray(){ // if only one allocation
statsarray S;
S.t='\0'; // type, population (p) or sample (s)
S.n=0; // number in sample/population
return S;
} 

statsarray init_statsarray(){ // for dynamic allocations
statsarray S;
S.t='\0'; // type, population (p) or sample (s)
S.n=0; // number in sample/population
S.d=(double*)calloc(1,sizeof(double));
return S;
} 

autocorr zero_autocorr(){ // if only one allocation
autocorr A;
A.k=0; // number in autocorrelation function
return A;
} 
autocorr init_autocorr(){ // for dynamic allocations
autocorr A;
A.k=0; // number in autocorrelation function
A.a=(double*)calloc(1,sizeof(double)); //the autocorr function
return A;
} 


