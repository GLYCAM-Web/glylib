// Function written by B. Lachele Foley, 2007
#include <mylib.h>
//#include "../inc/mylib.h"

/************** myfopen(const char,const char) ****************/ 
FILE* myfopen(const char *myfilename,const char *myopentype){
FILE *MYFP;
MYFP=fopen(myfilename,myopentype);
if(MYFP==NULL){
	printf("Error opening file.\n");
	printf("Expected location is:\n\n");
	printf("\t%s\n\n",myfilename);
	printf("Attempted file open for %s\n",myopentype);
	printf("Exiting.\n");
	exit(1);
	}
return MYFP;
}


/************** myfreopen(const char,const char,FILE *) ****************/ 
/* This function is useful when the open status of a file is not know.
The freopen function will attempt to "close" the file, but will ignore
any close errors (such as the file not being open) and then (re)open it. */
FILE* myfreopen(const char *myfilename,const char *myopentype, FILE *alreadyopen){
FILE *MYFP;
MYFP=freopen(myfilename,myopentype,alreadyopen);
if(MYFP==NULL){
	printf("Error opening file.\n");
	printf("Expected location is:\n\n");
	printf("\t%s\n\n",myfilename);
	printf("Attempted file open for %s\n",myopentype);
	printf("Exiting.\n");
	exit(1);
	}
return MYFP;
}
