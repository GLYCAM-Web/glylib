/* File read_prep.c begun on 20071120 by BLFoley

Purpose:  Read prep file contents into a molecule structure

NOTE:	At present, this function only reads in the prep file.
	It will not set bonding, etc., in the main residue.

NOTE:	Much information is read into the void pointer space.

Since prep files are generated on a per-residue basis, each 
entry will be considered part of the molecule structure.

This function will rely on read_prepres, which will read in
a single residue.  That function will rely on read_prepatom
as well.  Both read_prepres and read_prepatom will return
enture residues/atoms.

In order to allow for the addition of many prep files to an
expanding database of residues, this function will add residues
to an existing molecule structure (rather than returning a
new structure).  

NOTE:   THE PARENT MOLECULE STRUCTURE MUST BE PREINITIALIZED
	THIS INCLUDES AN INITIAL ALLOCATION OF THE RESIDUE
	SUB STRUCTURES (M.nr should be set to zero if there are
	no existing residues in the molecule).
*/ 
#include <mylib.h>
#include <molecules.h>
#include <general.h>
#include <AMBER/read_prep.h>
#include <declarations.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
//#include "../inc/general.h"
//#include "../inc/read_prep.h"
//#include "../inc/declarations.h"



void read_prep(molecule *M, fileset F, types *TYP){
int rpa=0,stopflag=1,rpuse=0,rpsz=0,rpbins=0;
char line[501];
ambermprepinfo *amp;
fpos_t f_linestart;

M[0].nVP=1;
amp=(ambermprepinfo*)calloc(1,sizeof(ambermprepinfo));
F.F=myfreopen(F.N,"r",F.F); // don't depend on calling function to have opened
amp[0].FN=strdup(F.N);
fgets(line,500,F.F);
sscanf(line,"%d %d %d",&amp[0].IDBGEN,&amp[0].IREST,&amp[0].ITYPF);
fgets(line,500,F.F);
amp[0].NAMDBF=strdup(line);
M[0].VP=amp;

rpuse=malloc_usable_size(M[0].r);
rpsz=sizeof(residue);
rpbins=rpuse/rpsz;
if(rpbins<1){mywhine("read_prep: no space in residue -- is molecule not initialized?");}

//printf("about to read through the file. Current line is\n\t>>%s<<\n",line);
fgetpos(F.F,&f_linestart);// get current position
while(fgets(line,500,F.F)!=NULL){// while not end of file
//printf("\treading a line ...  line is\n>>%s<<\n",line);
	stopflag=1;
	for(rpa=0;rpa<strlen(line);rpa++){// if the first word is not "STOP" (and only that??)
		if((line[rpa]!='\t')&&(line[rpa]!=' ')){
//printf("rpa is %d and line[rpa] is >>%c<< and strlen(line) is %d\n",rpa,line[rpa],strlen(line));
			if(rpa+4>strlen(line)){break;}
			else{ 
				if((line[rpa]=='S')&&(line[rpa+1]=='T')&&(line[rpa+2]=='O')&&(line[rpa+3]=='P')){ 
					stopflag=0;
					}
				break;
				}
			}
		}
//printf("stopflag is %d\n",stopflag);
	if(stopflag==0){break;} // break loop
//printf("about to allocate space for the residue...\n");
	M[0].nr++;// allocate new memory for a residue
	M[0].r=(residue*)realloc(M[0].r,M[0].nr*sizeof(residue));
	fsetpos(F.F,&f_linestart); // rewind to start of line
	M[0].r[M[0].nr-1]=read_prepres(F,TYP);// call read_prepres to read in the residue information
	fgetpos(F.F,&f_linestart);// get current position 
	} // end while not end of file
return;
}

/******************** read_prepres ****************/
// START HERE -- add in the void pointer read stuff
/* this function reads in residue information from a prep file
	-- currently, it only reads in number, name, type and charge
	-- connectivity, etc., shall be defined later, though it will
		read in the information */
residue read_prepres(fileset F, types *TYP){ // will the position be what is expected?
amberrprepinfo *arp;
residue R;
char line[501],desc[1001],ctmp[51],ctmp2[51]; // hope no prep file line is longer than this...
char dumat[21]; // the symbol used for dummy atoms 
char *linestat; // did the line read work?
int itmp=0,rra=0,rline=0;
char whinetext[1001];

// initialize R
R.n=0; // residue number given in input file
R.N=strdup(""); // reside name -- can live with residue names of 200 chars or less???  :-)
R.na=0; // number of atoms in residue
R.m=0; // molecular weight
R.COM=zero_coord(); // center of mass for residue
R.a=(atom*)calloc(1,sizeof(atom)); // atom structures (na of these)
R.nbs=0; // number of bond sets
R.bs=(molbondset*)calloc(1,sizeof(molbondset)); // (consecutive bonds, use these for plotting, etc.)
//int nr; // number of simple rings (no cage structures, etc.)
//int nrbs; // number of ring bondsets defined
//bondset *rbs; // bondsets for rings
//int nrc; // number of ring/reference coordinate sets defined
//coord_3D *rc; // coordinates for ring/reference centers
//int nrp; // number of ring planes defined
//plane *rp; // equations for average/approximate/exact/etc. ring planes (where useful)
R.ni=1; // number of other indices
R.i=(int*)calloc(R.ni,sizeof(int)); // other indices, as needed (ni of these)
R.nd=1; // number of double-precision parameters
R.d=(double*)calloc(R.nd,sizeof(double)); // other parameters, as needed (nd of these)

// allocate space for prepinfo
arp=(amberrprepinfo*)calloc(1,sizeof(amberrprepinfo));
arp[0].IOPR=(llcharset*)calloc(1,sizeof(llcharset));

// Read residue header information into residue structure
/* A sample header is:
0-[ALPHA-D-ALLO-] terminal, RESP 0.010 HF/6-31Gstar//HF/6-31Gstar
0na.dat
0NA   INT 0
CORRECT OMIT DU BEG
 0.1940

By line, they are:
1-2] R.D
3] R.N, error if not INT (for now), error if not zero (for now)
4] error if not "CORRECT", ignore, read in dummy atom name, ignore
5] R.d[0]
*/
// first two lines (description)
linestat=fgets(desc,500,F.F); // first line
//printf("line 1 looks like: >>%s<<\n",desc);
if(linestat==NULL){read_eek("EOF during residue header read, line 1",F.N);}
arp[0].TITLE=strdup(desc); // add description to arp (destined for void pointer)
R.D=strdup(desc); // copy to residue structure
linestat=fgets(line,500,F.F); // second line
//printf("line 2 looks like: >>%s<<\n",line);
if(linestat==NULL){read_eek("EOF during residue header read, line 2",F.N);}
arp[0].NAMF=strdup(line);
// next line (name and other)
linestat=fgets(line,500,F.F); // third line
if(linestat==NULL){read_eek("EOF during residue header read, line 3",F.N);}
sscanf(line,"%s %s %d",desc,ctmp,&itmp);
//printf("line 3 looks like: >>%s<<\n",line);
R.N=strdup(desc); // set the residue name
arp[0].NAMRES=strdup(desc);
if(strcmp(ctmp,"INT")!=0){read_eek("entry not INT on line 3 of residue header",F.N);}
if(itmp!=0){read_eek("KFORM not zero on line 3 of residue header",F.N);}
arp[0].INTX=strdup(ctmp);
arp[0].KFORM=itmp;
// fourth line (dummy atom symbol & check)
linestat=fgets(line,500,F.F); // fourth line (dummy atom symbol & check)
if(linestat==NULL)read_eek("EOF during residue header read, line 4",F.N);
sscanf(line,"%s %s %s %s",desc,ctmp,dumat,ctmp2);
if(strcmp(desc,"CORRECT")!=0){read_eek("entry not CORRECT on line 4 of residue header",F.N); }
arp[0].IFIXC=strdup(desc);
arp[0].IOMIT=strdup(ctmp);
arp[0].ISYMDU=strdup(dumat);
arp[0].IPOS=strdup(ctmp2);
// fifth line (charge on residue)
linestat=fgets(line,500,F.F); // fourth line (dummy atom symbol & check)
if(linestat==NULL)read_eek("EOF during residue header read, line 5 (residue charge)",F.N);
sscanf(line,"%lf",&arp[0].CUT);
// scan first three dummy atom lines, make sure that the third field
// 	of each is dumat and ignore otherwise (for now, at least)
for(rra=0;rra<3;rra++){
	linestat=fgets(line,500,F.F); // fourth line (dummy atom symbol & check)
	if(linestat==NULL){
		sprintf(whinetext,"EOF during residue header read, dummy atom line %d",(rra+1));
		read_eek(whinetext,F.N); }
	sscanf(line,"%s %s %s ",desc,ctmp,ctmp2);
	if(strcmp(ctmp2,dumat)!=0){
//printf("\t dumat is >>%s<<, and line looks like: >>%s<<\n",dumat,line);
		sprintf(whinetext,"dummy atom symbol, line %d, does not match declared value",(rra+1));
		read_eek(whinetext,F.N); }
	}

while(fgets(line,500,F.F)!=NULL){ // while not end of file
	rline=1;
//printf("line is: %s\n",line);
	for(rra=0;rra<strlen(line);rra++){ // check for an empty line
		if((line[rra]!='\t')&&(line[rra]!=' ')&&(line[rra]!='\n')){ // if this character is NOT a space or tab
			rline=0;
			break; } }
	if(rline==1) break;// if an empty line, break loop
//printf("calling read atom for line: %s\n",line);
	R.na++;//	allocate space
	R.a=(atom*)realloc(R.a,R.na*sizeof(atom));
	R.a[R.na-1]=read_prepatom(line,TYP);//read_prepatom
	// update molecular weight for residue
	R.m+=TYP[0].a[R.a[R.na-1].t].m;
	} // end read of the residue's atoms

// while not "DONE", read any after-atoms residue information
while(fgets(line,500,F.F)!=NULL){ // while not end of file
	rline=1;
	for(rra=0;rra<strlen(line);rra++){ // check for an empty line
		if((line[rra]!='\t')&&(line[rra]!=' ')){ // if this character is NOT a space or tab
			rline=0;
			break; } }
	if(rline==1) continue;// if an empty line, check next line 
	if((line[0]=='D')&&(line[1]=='O')&&(line[2]=='N')&&(line[3]=='E')){break;}
	arp[0].nIOPR++;
	arp[0].IOPR=(llcharset*)realloc(arp[0].IOPR,arp[0].nIOPR*sizeof(llcharset));
	arp[0].IOPR[arp[0].nIOPR-1].T=strdup(line);
	} // end scan of post-atom residue info

//printf("Residue %s contains %d atoms\n",R.N,R.na);
// for future:
//   scan atom info to set bonding
R.nVP=1;
R.VP=arp;

return R;	
}


/******************* read_prepatom *******************/

atom read_prepatom(const char *line, types *TYP){
amberaprepinfo *aap;
atom A;
char N[21],T[21],C[21]; // if any entry is longer than 21 chars, something is wrong...
char whinetext[201]; // error message text
int raa=0;

// Most of this will go into an amberaprepinfo structure.
// Some bits also go into an atom structure
// Here is a typical prep-file atom line:
//	>> 4 C1   CG  M  3  2  1  1.400   113.9      60.0     0.4780<<
// 	   1  2   3   4  5  6  7    8      9          10       11
//	The atom structure bits:
//	1: A.n
//	2: A.N=strdup(2)
//	3: A.T=strdup(3)  --and-- A.t=number for type 3

aap=(amberaprepinfo*)calloc(1,sizeof(amberaprepinfo));

sscanf(line,"%d %s %s %s %d %d %d %lf %lf %lf %lf",\
	&aap[0].I,N,T,C,&aap[0].NA,&aap[0].NB,&aap[0].NC,&aap[0].R,&aap[0].THETA,&aap[0].PHI,&aap[0].CHG);
A.N=strdup(N);
aap[0].IGRAPH=strdup(N);
A.T=strdup(T);
aap[0].ISYMBL=strdup(T);
aap[0].ITREE=strdup(C);
A.n=aap[0].I-3; 
// set the atom's type according to the data in TYP
A.t=-1;
for(raa=0;raa<TYP[0].na;raa++){
//printf("TYP[0].a[raa].N is >>%s<< ; A.T is >>%s<<\n",TYP[0].a[raa].N,A.T);
	if(strcmp(TYP[0].a[raa].N,A.T)==0){ 
		A.t=raa; 
		break;} }
if(A.t==-1){
	sprintf(whinetext,"read_prepatom:\n\tMatch not found for atom type %s.  The prep file line follows:\n\n%s\n\n",A.T,line);
	mywhine(whinetext);
	}
//printf("The entire line is:\n");
/*printf("%d %s %s %s %d %d %d %f %f %f %f\n",aap[0].I,aap[0].IGRAPH,aap[0].ISYMBL,\
        aap[0].ITREE,aap[0].NA,aap[0].NB,aap[0].NC,aap[0].R,aap[0].THETA,aap[0].PHI,aap[0].CHG);*/
//
A.nVP=1;
A.VP=aap;

return A;
}

