#include <stats.h>
#include <mylib.h>
//#include "../inc/stats.h"
//#include "../inc/mylib.h"

/********* get_autocorr_est_array *********/
//autocorr get_autocorr_est_array(statsarray S,meanvar M);
autocorr get_autocorr_est_array(statsarray S,meanvar M){
autocorr A;
int k=0,j=0;
double *Dmu;

// initialize
A.k=M.n;
A.a=(double*)calloc((A.k),sizeof(double));
Dmu=(double*)calloc(S.n,sizeof(double));

for(j=0;j<S.n;j++){
	Dmu[j]=S.d[j]-M.m;
	}
for(k=0;k<(A.k);k++){
	for (j=0;j<(A.k-k);j++){
		A.a[k]+=Dmu[j+k]*Dmu[j];
		}
	A.a[k]/=((A.k-k)*M.v); // assumes whole population or biased sample
	}
free(Dmu);
return A;
}
