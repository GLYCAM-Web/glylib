// Function written by B. Lachele Foley, 2007

// prints and error message and a usage statement then exits
#include <mylib.h>

void mywhineusage(const char *msg,const char *usg){ 
	printf("\nThis program encountered the following error and will exit:\n\t%s\n\n",msg);
	printf("Usage:\n%s\n",usg); 
        exit(1);
}
