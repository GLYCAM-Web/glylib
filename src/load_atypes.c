// File load_atypes.c begun on 20071128 by BLFoley
// Purpose: load database of atom types into atype
//	substructure within a types superstructure
// NOTE!  T needs to be initialized -before- calling
//	this function!
#include <mylib.h>
#include <general.h>
#include <molecules.h>
#include <declarations.h>
//#include "../inc/mylib.h"
//#include "../inc/general.h"
//#include "../inc/molecules.h"
//#include "../inc/declarations.h"

void load_atypes(fileset F, types *T){
char line[10001],*entry;
int test=0,la=0,numbnd=0;
F.F=myfopen(F.N,"r"); // 
while(fgets(line,10000,F.F)!=NULL){// while not end of file
	if(line[0]!='#'){//check that the line does not start with a hash 
		if(strchr(line,EOF)!=NULL){read_eek("EOF",F.N);}
		T[0].na++;// increment and allocate
		T[0].a=(atype*)realloc(T[0].a,T[0].na*sizeof(atype));
		// read data into structure by splitting on the : token
		//strcpy(entry,strtok(line,":")); // get atomic number
		test=sscanf(strtok(line,":"),"%d",&T[0].a[T[0].na-1].n);
		if(test!=1){mywhine("load_atypes: error reading atomic number");}
		T[0].a[T[0].na-1].N=strdup(strtok(NULL,":")); // get atom name
		T[0].a[T[0].na-1].NT=strdup(strtok(NULL,":")); // get element name
		//strcpy(entry,strtok(NULL,":")); // get mass
		test=sscanf(strtok(NULL,":"),"%lf",&T[0].a[T[0].na-1].m);
		if(test!=1){mywhine("load_atypes: error reading atomic mass");}
		T[0].a[T[0].na-1].desc=strdup(strtok(NULL,":")); // get description 
		entry=strdup(strtok(NULL,":")); // get bonding info 
		test=sscanf(strtok(NULL,":"),"%d",&T[0].a[T[0].na-1].nlp);// get number lone pairs
		if(test!=1){mywhine("load_atypes: error reading number of lone pairs");}
		// split bonding info
		if(strchr(entry,'(')==NULL){mywhine("load_atypes: bonding field does not contain bond orders");}
		test=sscanf(strtok(entry,"("),"%d",&T[0].a[T[0].na-1].nb);// get number typical bonds
		if(test!=1){mywhine("load_atypes: error reading number of typical bonds");}
		numbnd=T[0].a[T[0].na-1].nb;
		if(numbnd==0){
			numbnd=1;
			T[0].a[T[0].na-1].bo=(double*)calloc(1,sizeof(double));
			test=sscanf(strtok(NULL,",)"),"%lf",&T[0].a[T[0].na-1].bo[0]);// get bond order
			if(test!=1){mywhine("load_atypes: error reading bond order");}
			if(T[0].a[T[0].na-1].bo[0]!=0) T[0].a[T[0].na-1].nb=-1;
			}
		else{
			T[0].a[T[0].na-1].bo=(double*)calloc(numbnd,sizeof(double));
			for(la=0;la<numbnd;la++){
				test=sscanf(strtok(NULL,",)"),"%lf",&T[0].a[T[0].na-1].bo[la]);// get bond order
				if(test!=1){mywhine("load_atypes: error reading bond order");}
				} 
			}
		} // close if line is not a comment
	}// close while not EOF
// return; 
return;
}
