#include <stats.h>
#include <mylib.h>
//#include "../inc/stats.h"
//#include "../inc/mylib.h"

/********* get_meanvar_array *********/
// written by BLFoley May 2008
//meanvar get_meanvar_array(statsarray S);
meanvar get_meanvar_array(statsarray S){
meanvar M;
int i=0;

// initialize
M=zero_meanvar();
M.t=S.t;
M.n=S.n;

for(i=0;i<S.n;i++){
	M.m+=S.d[i];
	M.v+=S.d[i]*S.d[i];
	}
M.v-=(M.m*M.m/(double)M.n);
M.m/=(double)M.n;
switch(M.t){ // calculate standard deviation
	case 's': 
		  M.v/=(double)(M.n-1); // for sample / unbiased pop estimate
		  break;
	case 'p': 
		  M.v/=(double)(M.n); // for whole population or biased sample
		  break;
	default: mywhine("unknown type for statistical set in get_meanvar_array"); 
	} 
M.s=sqrt(M.v);
return M;
}
