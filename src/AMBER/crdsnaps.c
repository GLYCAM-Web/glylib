// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <AMBER/amber_coords.h>
/*********************** crdsnaps(const char cstop, const char cscrd) *******************/
int crdsnaps(const char *cstop,const char *cscrd,int csdot){
int csa=0,csnat=0,csnsnp=0,csstat=0,cssig=0;
char cstst[501];
double csEa=0,csEb=0,csEc=0;
FILE *CST,*CSC;

// Probe top file for number of atoms
CST=myfopen(cstop,"r");
csstat=fscanf(CST,"%s",cstst);
if(csstat!=1) mywhine("Error reading topology file");
//printf("cstst is >>>%s<<<\n",cstst);
while(csstat==1){
	//printf("\t in loop cstst is >>>%s<<<\n",cstst);
	if(strcmp(cstst,"POINTERS")==0){
		csstat=fscanf(CST,"%s",cstst);
//printf("cstst is >>>%s<<<\n",cstst);
		if(csstat!=1) mywhine("Error reading topology file"); 
		csstat=fscanf(CST,"%d",&csnat);
//printf("csnat is >>>%d<<<\n",csnat);
		if(csstat!=1) mywhine("Error reading number of atoms in topology file."); 
		//break;
		}
	if(strcmp(cstst,"BOX_DIMENSIONS")==0) csnat++;
	//if(csnat!=0) break; 
	csstat=fscanf(CST,"%s",cstst);
	}
if(csnat==0) mywhine("Got zero for number of atoms.");
//printf("crdsnaps after CST loop\n");
fclose(CST);
//printf("crdsnaps after close CST\n");
// Probe crd file for total number of snapshots 
//printf("crdsnaps opening CSC\n");
CSC=myfopen(cscrd,"r");
//printf("crdsnaps fgets CSC\n");
fgets(cstst,201,CSC);
//printf("cstst is >>>%s<<<\n",cstst);
if(cstst==NULL) mywhine("Error reading coordinate file."); 
//if(strlen(cstst)==0) mywhine("Error reading coordinate file."); 
//printf("cstst about to fscanf \n");
csstat=fscanf(CSC,"%lf %lf %lf",&csEa,&csEb,&csEc);;
//printf("cstst after fscanf  csstat is >>>%d<<<\n",csstat);
// START HERE -- fix this for loop structure
printf(" ");
while(csstat==3){
	csa=1;
	while(csa<csnat){
		csstat=fscanf(CSC,"%lf %lf %lf",&csEa,&csEb,&csEc);
//printf("csa is %d; csnat is %d, csstat is %d\n",csa,csnat,csstat);
		if(csstat!=3) break;
		csa++;
		//switch(cssig){	
			//case 0:	printf("\b-");
				//cssig=1;
				//break;
			//case 1:	printf("\b\\");
				//cssig=2;
				//break;
			//case 2:	printf("\b|");
				//cssig=3;
				//break;
			//case 3:	printf("\b/");
				//cssig=0;
				//break;
			//default: mywhine("something wrong in cssig in crdsnaps.c");
			//}
		} 
	if(csa==csnat){csnsnp++;}
	else if (csa!=1){mywhine("Mismatch between NATOMs in top and crd files.\n\
\tPerhaps there is a mismatch between BOX info in crd and top."); }
//printf("csa is %d; csnat is %d; csnsnp is %d\n",csa,csnat,csnsnp);
	csstat=fscanf(CSC,"%lf %lf %lf",&csEa,&csEb,&csEc);
	if((csnsnp%csdot)==0){printf("%d  ",csnsnp);}
	switch(cssig){	
		case 0:	printf("\b-");
			cssig=1;
			break;
		case 1:	printf("\b\\");
			cssig=2;
			break;
		case 2:	printf("\b|");
			cssig=3;
			break;
		case 3:	printf("\b/");
			cssig=0;
			break;
		default: mywhine("something wrong in cssig in crdsnaps.c");
		}
//printf("csstat is %d\n",csstat);
	}
printf("\n");
return csnsnp;
}
