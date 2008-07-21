// Function written by B. Lachele Foley, 2007

// prints an error message and exits
#include <mylib.h>
//#include "../inc/mylib.h"

void mywhine(const char *msg){ 
	printf("\nThis program encountered the following error and will exit:\n\t%s\n\n",msg);
        exit(1);
}
