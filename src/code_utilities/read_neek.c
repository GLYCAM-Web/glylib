/**************************  read_neek() *******************************/
/* This exits if there is a read problem 
   Author: Michael Tessier (slight mods by BLFoley)*/
//#include <load_pdb.h>
//#include "../inc/load_pdb.h"
#include "../inc/mylib.h"
void read_neek(const char *NEEK, int x, int y){
printf("Unexpected %s when reading line %d at field f[%d]\n",
        NEEK,x,y);
printf("Exiting.\n");
exit(1);
return;
}
/**************************  read_fneek() *******************************/
/* This exits if there is a read problem and gives more inf
   Author: Michael Tessier plus Lachele Foley */
//#include "../inc/mylib.h"
void read_fneek(const char *NEEK, int x, int y, const char *FILENAME){
printf("Unexpected: %s \n\treading line %d at field %d\n\tfile:%s\n",
        NEEK,x,y,FILENAME);
printf("Exiting.\n");
exit(1);
return;
}
