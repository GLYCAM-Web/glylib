// Function written by B. Lachele Foley, 2007
#include <mylib.h>
//#include "../inc/mylib.h"
/**************************  read_eek() *******************************/
/* This exits if there is a read problem */
void read_eek(const char *EEK, const char *fdesc){
printf("Unexpected %s during read of file: %s\n",EEK,fdesc);
printf("Exiting.\n");
exit(1);
return;
}

