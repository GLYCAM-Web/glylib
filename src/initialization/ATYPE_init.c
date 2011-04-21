/* ATYPE_init.c -- begun 20071015 by BLFoley
 * This function initializes an atom type list with several
 * commonly-used atom types */
#include <mylib.h>
#include <molecules.h>
#include <PDB.h>
//#include "../inc/mylib.h" 
//#include "../inc/molecules.h" 
//#define NUMAT 15 // change this if you add more types
//atype *ATYPE_init();
/* For convenience:
	typedef struct {
		int n; // atomic number
		char N[21]; // element name
		char NT[21]; // name for element type
		double m; // mass of element type
		char desc[301]; // brief free-form description field
		int nb; // number of typical bonds for this atom type
		double *bo; // typical bond orders (nb of these)
	} atype; // atom types, with other info
*/
atype *ATYPE_init(){
atype *A;

NUMAT=15;

A=(atype*)calloc(NUMAT,sizeof(atype));

A[0].n=6;
strcpy(A[0].N,"CG");
strcpy(A[0].NT,"C");
A[0].m=12.011;
strcpy(A[0].desc,"sp3");
A[0].nb=4;
//
A[1].n=6;
strcpy(A[1].N,"CY");
strcpy(A[1].NT,"C");
A[1].m=12.011;
strcpy(A[1].desc,"sp3");
A[1].nb=4;
//
A[2].n=1;
strcpy(A[2].N,"H");
strcpy(A[2].NT,"H");
A[2].m=1.00794;
strcpy(A[2].desc,"sp3");
A[2].nb=1;
//
A[3].n=1;
strcpy(A[3].N,"H1");
strcpy(A[3].NT,"H");
A[3].m=1.00794;
strcpy(A[3].desc,"sp3");
A[3].nb=1;
//
A[4].n=1;
strcpy(A[4].N,"H2");
strcpy(A[4].NT,"H");
A[4].m=1.00794;
strcpy(A[4].desc,"sp3");
A[4].nb=1;
//
A[5].n=1;
strcpy(A[5].N,"HC");
strcpy(A[5].NT,"H");
A[5].m=1.00794;
strcpy(A[5].desc,"sp3");
A[5].nb=1;
//
A[6].n=1;
strcpy(A[6].N,"HO");
strcpy(A[6].NT,"H");
A[6].m=1.00794;
strcpy(A[6].desc,"sp3");
A[6].nb=1;
//
A[7].n=1;
strcpy(A[7].N,"HW");
strcpy(A[7].NT,"H");
A[7].m=1.00794;
strcpy(A[7].desc,"sp3");
A[7].nb=1;
//
A[8].n=7;
strcpy(A[8].N,"N"); 
strcpy(A[8].NT,"N");
A[8].m=14.01;
strcpy(A[8].desc,"sp2");
A[8].nb=3;
//
A[9].n=8;
strcpy(A[9].N,"OH");
strcpy(A[9].NT,"O");
A[9].m=15.994;
strcpy(A[9].desc,"sp3");
A[9].nb=2;
//
A[10].n=8;
strcpy(A[10].N,"OS");
strcpy(A[10].NT,"O");
A[10].m=15.994;
strcpy(A[10].desc,"sp3");
A[10].nb=2;
//
A[11].n=8;
strcpy(A[11].N,"O"); 
strcpy(A[11].NT,"O");
A[11].m=15.994;
strcpy(A[11].desc,"sp2");
A[11].nb=1;
//
A[12].n=8;
strcpy(A[12].N,"O2");
strcpy(A[12].NT,"O");
A[12].m=15.994;
strcpy(A[12].desc,"sp2");
A[12].nb=1;
//
A[13].n=8;
strcpy(A[13].N,"OW");
strcpy(A[13].NT,"O");
A[13].m=15.994;
strcpy(A[13].desc,"sp3");
A[13].nb=2;
//
A[14].n=11;
strcpy(A[14].N,"IP");
strcpy(A[14].NT,"Na+");
A[14].m=15.994;
strcpy(A[14].desc,"sp3");
A[14].nb=0;

//theoretically, something like this should work...
//A[]={
//{6,"C",  "C", 12.011, "sp2", 3, {1,1,2}},
//{6,"CG", "C", 12.011, "sp3", 4, {1,1,1,1}},
//{6,"CY", "C", 12.011, "sp3", 4}
//{1,"H",  "H", 1.00794, "sp3", 1}
//{1,"H1", "H", 1.00794, "sp3", 1}
//{1,"H2", "H", 1.00794, "sp3", 1}
//{1,"HC", "H", 1.00794, "sp3", 1}
//{1,"HO", "H", 1.00794, "sp3", 1}
//{1,"HW", "H", 1.00794, "sp3", 1}
//{7,"N",  "N", 14.01,  "sp2", 3}
//{8,"OH", "O", 15.994, "sp3", 2}
//{8,"OS", "O", 15.994, "sp3", 2}
//{8,"O",  "O", 15.994, "sp2", 1}
//{8,"O2", "O", 15.994, "sp2", 1}
//{8,"OW", "O", 15.994, "sp3", 2}
//{11,"IP", "Na+", 15.994, "sp3", 0}
//};

return A;
}
