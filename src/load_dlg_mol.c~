/* File load_dlg_mol.c begun on 20071015 by BLFoley
 * Function load_dlg_mol: loads AutoDock dlg output into a
 * 	molecule structure, placing the initial structure in
 * 	the main coordinate set and the docked structures into
 * 	alternate coordinate sets.  Also calculates COM for
 * 	the main molecule and each residue (initial and docked).
*/
#include "../inc/mylib.h"
#include "../inc/general.h"
#include "../inc/molecules.h"
#include "../inc/declarations.h"

//dockinfo *load_dlg_mol(fileset F,atype *AT);

dockinfo *load_dlg_mol(fileset F,types *T){
int a=0,b=-1000,anum=0,rnum=0,ndock=0,ai=0,di=0,ri=0,ti=0;
int scntst=0,*AI,nat=0,localdebug=1;
int fRes = -1;
fpos_t* start = (fpos_t*)calloc(1,sizeof(fpos_t));
char line[501],dum1[50],dum2[50],at[10],res[10]; // maybe no lines longer than 500...
char whinetext[501];
char* cid = (char*)calloc(5,sizeof(char));
double x=0,y=0,z=0,xl=0,xh=0,yl=0,yh=0,zl=0,zh=0;
atype* AT = T[0].a;
dockinfo *D;

if(localdebug>0){printf("load_dlg_mol: At top.\n");}
// initialize some stuff
D=(dockinfo*)calloc(1,sizeof(dockinfo));
D[0].M.nr=1;
D[0].M.na=0;
D[0].M.r=(residue*)calloc(1,sizeof(residue));
D[0].M.r[0].na=D[0].M.r[0].m=0;
D[0].M.r[0].a=(atom*)calloc(1,sizeof(atom));
D[0].M.COM.i=D[0].M.COM.j=D[0].M.COM.k=D[0].M.m=0;
D[0].M.r[0].COM.i=D[0].M.r[0].COM.j=D[0].M.r[0].COM.k=D[0].M.r[0].m=0;
//if(localdebug>1){dprint_molecule(&D[0].M,10);}
if(localdebug>0){printf("load_dlg_mol: After init.\n");}
rewind(F.F);
// find the start of the input PDBQ file
while(fgets(line,500,F.F)!=NULL){
if(localdebug>3){printf("line is >>>%s<<<\n",line);}
	if(strstr(line,"INPUT PDBQ FILE")!=NULL){
		fgets(line,500,F.F);
		fgets(line,500,F.F);
		break;
		}}
if(localdebug>0){printf("load_dlg_mol: After found PDBQ.\n");}
//Find the smallest residue number
fgetpos(F.F,start);
while(strstr(line,"__________________________________") == NULL){
	if(strstr(line,"INPUT-PDBQ: ATOM")!=NULL ||
           strstr(line,"INPUT-PDBQ: HETATM")!=NULL){
		scntst=sscanf(line,"%s %s %d %s %s %d",dum1,dum2,&anum,at,res,&rnum);
		if(rnum < fRes || fRes < 0)
			fRes = rnum;
	}
	fgets(line,500,F.F);
}
fRes -= 1;
fsetpos(F.F,start);
free(start);
// read in the PDBQ info to the structure
while(fgets(line,500,F.F)!=NULL){
if(localdebug>3){printf("load_dlg_mol: Reading file.\n");}
if(localdebug>3){printf("line is >>>%s<<<\n",line);}
	if(strstr(line,"INPUT-PDBQ")!=NULL){
if(localdebug>3){printf("load_dlg_mol: Reading PDBQ.\n");}
		if(strstr(line,"INPUT-PDBQ: ATOM")!=NULL ||
		   strstr(line,"INPUT-PDBQ: HETATM")!=NULL){
			scntst=sscanf(line,"%s %s %d %s %s %d %lf %lf %lf",\
				dum1,dum2,&anum,at,res,&rnum,&x,&y,&z);
if(localdebug>3){printf("line is >>>%s<<<\n",line);}
			rnum-=fRes;
			//***START EDIT MNT Aug,14 2008
			if(scntst!=9){
				scntst=sscanf(line,"%s %s %d %s %s %c %d %lf %lf %lf",\
					dum1,dum2,&anum,at,res,cid,&rnum,&x,&y,&z);
				if(scntst!=10){read_eek("end of line while reading INPUT-PDBQ: ATOM",F.N);}
			}//***END EDIT
			if(b==-1000){
				xl=xh=x; 
				yl=yh=y; 
				zl=zh=z; 
				b=0;
				}
			else{
				if(x<xl){xl=x;}
				if(x>xh){xh=x;}
				if(y<yl){yl=y;}
				if(y>yh){yh=y;}
				if(z<zl){zl=z;}
				if(z>zh){zh=z;}
				}
			if(rnum>D[0].M.nr){
if(localdebug>3){printf("\tload_dlg_mol: need more space for rnum -- assigning\n");}
				ri=D[0].M.nr; // number already allocated
				D[0].M.nr=rnum;
				D[0].M.r=(residue*)realloc(D[0].M.r,rnum*sizeof(residue));
				for(b=ri;b<rnum;b++){
					D[0].M.r[b].m=0;
					D[0].M.r[b].na=0;
					D[0].M.r[b].N = (char*)calloc(11,sizeof(char));
					D[0].M.r[b].D = NULL;
					D[0].M.r[b].bs = NULL;
					D[0].M.r[b].rbs = NULL;
					D[0].M.r[b].rc = NULL;
					D[0].M.r[b].rp = NULL;
					D[0].M.r[b].i = NULL;
					D[0].M.r[b].d = NULL;
					D[0].M.r[b].VP = NULL;
					D[0].M.r[b].N[0]='\0';
					D[0].M.r[b].nd = 0;
					D[0].M.r[b].ni = 0;
					D[0].M.r[b].nr = 0;
					D[0].M.r[b].nrp = 0;
					D[0].M.r[b].nrc = 0;
					D[0].M.r[b].nrbs = 0;
					D[0].M.r[b].nbs = 0;
					D[0].M.r[b].nVP = 0;
					D[0].M.r[b].a=(atom*)calloc(1,sizeof(atom)); 
					D[0].M.r[b].a[0].ni=0;
					D[0].M.r[b].a[0].nb=0;
					D[0].M.r[b].a[0].nmb=0;
					D[0].M.r[b].a[0].nalt=0;
					D[0].M.r[b].a[0].nvec=0;
			       		}
//if(localdebug>1){dprint_molecule(&D[0].M,10);}
				}
			ri=rnum-1;
			D[0].M.r[ri].n=rnum+fRes;//Set the residue number
if(localdebug>3){printf("\tload_dlg_mol: ri is %d ; D[0].M.r[ri].N is >>%s<<\n",ri,D[0].M.r[ri].N);}
			//if(D[0].M.r[ri].na==0){ strcpy(D[0].M.r[ri].N,res); }
			if(D[0].M.r[ri].na==0){ D[0].M.r[ri].N = strdup(res); }
			else{if(strcmp(D[0].M.r[ri].N,res)!=0){
				sprintf(whinetext,"residue name mismatch reading initial PDBQ\n\
between read residue >>%s<< and comparison residue (D[0].M.r[%d].N >>%s<<)",res,ri,D[0].M.r[ri].N);
				mywhine(whinetext);}}
			ai=D[0].M.r[ri].na; // atom index, for convenience
			D[0].M.r[ri].na++;
			D[0].M.na++;
			if(anum>nat){nat=anum;}; // overall atom count
			D[0].M.r[ri].a=(atom*)realloc(D[0].M.r[ri].a,(ai+1)*sizeof(atom));
			D[0].M.r[ri].a[ai].ni=0;
			D[0].M.r[ri].a[ai].nb=0;
			D[0].M.r[ri].a[ai].nmb=0;
			D[0].M.r[ri].a[ai].nalt=0;
			D[0].M.r[ri].a[ai].nvec=0;
			D[0].M.r[ri].a[ai].n=anum;
			//strcpy(D[0].M.r[ri].a[ai].N,at);
			D[0].M.r[ri].a[ai].N = strdup(at);
if(localdebug>3){printf("\tload_dlg_mol: ri is %d ; ai is %d ; D[0].M.r[ri].a[ai].N is >>%s<<\n",ri,ai,D[0].M.r[ri].a[ai].N);}
			D[0].M.r[ri].a[ai].x.i=x;
			D[0].M.r[ri].a[ai].x.j=y;
			D[0].M.r[ri].a[ai].x.k=z;
			//***START EDIT MNT Aug,14 2008 
			if(scntst == 10)
				D[0].M.r[ri].a[ai].cID = cid[0];
			if(strstr(line,"INPUT-PDBQ: HETATM")!=NULL)
				D[0].M.r[ri].a[ai].D = strdup("HETATM");
			else
				D[0].M.r[ri].a[ai].D = strdup("ATOM");//***END EDIT
			for(b=0;b<T[0].na;b++){ // might be a problem, the NUMAT(fixed)
				if(at[0]==AT[b].NT[0]){
					D[0].M.r[ri].a[ai].t=b;
					ti=b; // for convenience
					break;
					}
				}
			// add to the COM determinations
			D[0].M.m+=AT[ti].m; // add to molecule MW
			D[0].M.COM.i+=D[0].M.r[ri].a[ai].x.i*AT[ti].m; // molecule COM
			D[0].M.COM.j+=D[0].M.r[ri].a[ai].x.j*AT[ti].m;
			D[0].M.COM.k+=D[0].M.r[ri].a[ai].x.k*AT[ti].m;
			D[0].M.r[ri].m+=AT[ti].m; // add to residue MW
			D[0].M.r[ri].COM.i+=D[0].M.r[ri].a[ai].x.i*AT[ti].m; // residue COM
			D[0].M.r[ri].COM.j+=D[0].M.r[ri].a[ai].x.j*AT[ti].m;
			D[0].M.r[ri].COM.k+=D[0].M.r[ri].a[ai].x.k*AT[ti].m;
if(localdebug>1){printf("\tload_dlg_mol: ri = %d ; ai = %d (%s) type %s mass %f  (running= %f )\n",ri,ai,D[0].M.r[ri].a[ai].N,AT[ti].N,AT[ti].m,D[0].M.r[ri].m);}
if(localdebug>1){printf("\tload_dlg_mol: res main x y z is %20.15e  %20.15e  %20.15e \n",D[0].M.r[ri].a[ai].x.i,D[0].M.r[ri].a[ai].x.j,D[0].M.r[ri].a[ai].x.k);}
if(localdebug>1){printf("\tload_dlg_mol: running COM is %20.15e  %20.15e  %20.15e \n",D[0].M.r[ri].COM.i,D[0].M.r[ri].COM.j,D[0].M.r[ri].COM.k);}
			}
//if(localdebug>1){dprint_molecule(&D[0].M,10);}
		}
	if(strstr(line,"Number of atoms found in molecule")!=NULL){break;} 
	}
// ONE DAY:  add in any missing atoms
// -- open the prep database
// -- scan for the residue name
// -- record number/name/type of atoms
// -- mark as present or not, etc...
//
if(localdebug>0){printf("load_dlg_mol: About to start COM calc.\n");}
// finish COM calculation for molecule
D[0].M.COM.i/=D[0].M.m;
D[0].M.COM.j/=D[0].M.m;
D[0].M.COM.k/=D[0].M.m;
//printf("center of mass (MW=%20.12e) is at: %20.12e %20.12e %20.12e \n",D[0].M.COM.i,D[0].M.COM.j,D[0].M.COM.k);
// also for the residues -- plus record atom number correlation
if(localdebug>0){printf("load_dlg_mol: continuing residue-level COM calc.\n");}
AI=(int*)calloc(nat,sizeof(int));
for(a=0;a<D[0].M.nr;a++){ 
	if(D[0].M.r[a].na==0){
		sprintf(whinetext,"lack of atoms in residue >>%s<<, number %d",D[0].M.r[a].N,(a+1));
		read_eek(whinetext,F.N);
		}
	D[0].M.r[a].COM.i/=D[0].M.r[a].m;
	D[0].M.r[a].COM.j/=D[0].M.r[a].m;
	D[0].M.r[a].COM.k/=D[0].M.r[a].m;
if(localdebug>0){printf("\tload_dlg_mol: about to set ai indices.\n");}
	for(b=0;b<D[0].M.r[a].na;b++){
		if(D[0].M.r[a].a[b].n>nat){mywhine("problem counting atoms (n) in load_dlg");}
		AI[D[0].M.r[a].a[b].n-1]=b;
		}
	}

if(localdebug>0){printf("load_dlg_mol: Scanning to dockings.\n");}
// read in the alternate coordinate sets
while(fgets(line,500,F.F)!=NULL){
	if(strstr(line,"Number of requested LGA dockings")!=NULL){
if(localdebug>3){printf("load_dlg_mol: line is >>%s<<\n",line);} 
		sscanf(line,"Number of requested LGA dockings = %d ",&ndock);
if(localdebug>0){printf("load_dlg_mol: ndock is >>%d<<\n",ndock);} 
		break; }}
// allocate the memory -- use the internal ring coordinates to hold
// the extra centers of mass
if(localdebug>0){printf("load_dlg_mol: About to allocate; ndock is %d.\n",ndock);}
D[0].i=0;
D[0].n=ndock;
D[0].TR=(coord_3D*)calloc(ndock,sizeof(coord_3D));//translation for docked structure
D[0].Q=(vectormag_3D*)calloc(ndock,sizeof(vectormag_3D));// Quaternion x,y,z,w
D[0].QN=(vectormag_3D*)calloc(ndock,sizeof(vectormag_3D));// Quaternion nx,ny,nz,angle
D[0].QN0=(vectormag_3D*)calloc(ndock,sizeof(vectormag_3D));// quat0 nx,ny,nz,angle
D[0].eFEB=(double*)calloc(ndock,sizeof(double));// Est. Free Energy of Binding, kcal/mol [=(1)+(3)]
D[0].eKi=(double*)calloc(ndock,sizeof(double));// Est. Inhibition Const. Ki 
D[0].Tmp=(double*)calloc(ndock,sizeof(double));// temperature, Kelvin
D[0].fDE=(double*)calloc(ndock,sizeof(double));// Final Docked Energy, kcal/mol [=(1)+(2)]
D[0].fIE=(double*)calloc(ndock,sizeof(double));// (1) Final Intermolecular Energy
D[0].fIEL=(double*)calloc(ndock,sizeof(double));// (2) Final Internal Energy of Ligand
D[0].TFE=(double*)calloc(ndock,sizeof(double));// (3) Torsional Free Energy 
D[0].M.nrc=ndock;
D[0].M.rc=(coord_3D*)calloc(ndock,sizeof(coord_3D)); // for alt struct mol COM's
for(a=0;a<D[0].M.nr;a++){
	D[0].M.r[a].nrc=ndock;
	D[0].M.r[a].rc=(coord_3D*)calloc(ndock,sizeof(coord_3D)); // for alt struct res COM's
	for(b=0;b<D[0].M.r[a].na;b++){
		D[0].M.r[a].a[b].nalt=ndock;
		D[0].M.r[a].a[b].xa=(coord_3D*)calloc(ndock,sizeof(coord_3D));
		}
	}


if(localdebug>0){printf("load_dlg_mol: About to read alternate sets.\n");}
// read in the docked structures
di=-1; // index to check that all docked structures were found
while(fgets(line,500,F.F)!=NULL){
	if(strstr(line,"FINAL LAMARCKIAN GENETIC ALGORITHM DOCKED STATE")!=NULL){
if(localdebug>3){printf("\tload_dlg_mol: found a set; about to read energy info.\n");}
	di++;
	while(fgets(line,500,F.F)!=NULL){
		if(strstr(line,"DOCKED: USER                           _______")!=NULL){break;}
		else{
			if(strstr(line,"Quaternion nx,ny,nz,angle")!=NULL){
				sscanf(line,"%s %s %s %lf %lf %lf %lf",dum1,dum1,dum1,\
					&D[0].QN[di].i,&D[0].QN[di].j,&D[0].QN[di].k,&D[0].QN[di].d); 
				}
			if(strstr(line,"Quaternion x,y,z,w ")!=NULL){
				sscanf(line,"%s %s %s %lf %lf %lf %lf",dum1,dum1,dum1,\
					&D[0].Q[di].i,&D[0].Q[di].j,&D[0].Q[di].k,&D[0].Q[di].d); 
				}
			if(strstr(line,"Estimated Free Energy of Binding ")!=NULL){
				sscanf(line,"%s %s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,\
					dum1,dum1,dum1,dum1,&D[0].eFEB[di]);
				}
			if(strstr(line,"Estimated Inhibition Constant, Ki")!=NULL){
				sscanf(line,"%s %s %s %s %s %s %s %lf %s %s %lf",dum1,dum1,dum1,dum1,\
					dum1,dum1,dum1,&D[0].eKi[di],dum1,dum1,&D[0].Tmp[di]);
				}
			if(strstr(line,"Final Docked Energy")!=NULL){
				sscanf(line,"%s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,dum1,dum1,\
					&D[0].fDE[di]);
				}
			if(strstr(line,"Final Intermolecular Energy")!=NULL){
				sscanf(line,"%s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,dum1,\
					dum1,dum1,&D[0].fIE[di]);
				}
			if(strstr(line,"Final Internal Energy of Ligand")!=NULL){
				sscanf(line,"%s %s %s %s %s %s %s %s %s %lf",dum1,dum1,dum1,\
					dum1,dum1,dum1,dum1,dum1,dum1,&D[0].fIEL[di]);
				}
			if(strstr(line,"Torsional Free Energy")!=NULL){
				sscanf(line,"%s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,\
					dum1,dum1,dum1,&D[0].TFE[di]);
				}
			if((strstr(line,"NEWDPF about")!=NULL)&&(di==0)){
				sscanf(line,"%s %s %s %s %lf %lf %lf",dum1,dum1,dum1,dum1,\
					&D[0].RC.i,&D[0].RC.j,&D[0].RC.k); 
				}
			if(strstr(line,"NEWDPF tran0 ")!=NULL){
				sscanf(line,"%s %s %s %s %lf %lf %lf",dum1,dum1,dum1,dum1,\
					&D[0].TR[di].i,&D[0].TR[di].j,&D[0].TR[di].k);
				}
			if(strstr(line,"NEWDPF quat0")!=NULL){
				sscanf(line,"%s %s %s %s %lf %lf %lf %lf",dum1,dum1,dum1,dum1,\
					&D[0].QN0[di].i,&D[0].QN0[di].j,&D[0].QN0[di].k,&D[0].QN0[di].d); 
				}
			}
		}
	nat=0;
if(localdebug>3){printf("\tload_dlg_mol: found a set; about to read coordinates.\n");}
	while(fgets(line,500,F.F)!=NULL){
		if(strstr(line,"DOCKED: ATOM")!=NULL){
if(localdebug>3){printf("\tload_dlg_mol: found a set; in atom scan.\n");}
			nat++;
			scntst=sscanf(line,"%s %s %d %s %s %d %lf %lf %lf",\
				dum1,dum2,&anum,at,res,&rnum,&x,&y,&z);
			//if(rnum>3020){rnum%=3020;}//Added to handle LSTC_1918 (Temporary)
			rnum-=fRes;
			//***START EDIT MNT Aug,14 2008
			if(scntst!=9){
				scntst=sscanf(line,"%s %s %d %s %s %c %d %lf %lf %lf",\
					dum1,dum2,&anum,at,res,cid,&rnum,&x,&y,&z);
				if(scntst!=10){
					sprintf(whinetext,"end of line while reading DOCKED: ATOM (run %d)",(di+1));
					read_eek(whinetext,F.N);}
			}//***END EDIT
			if(x<xl){xl=x;}
			if(x>xh){xh=x;}
			if(y<yl){yl=y;}
			if(y>yh){yh=y;}
			if(z<zl){zl=z;}
			if(z>zh){zh=z;}
			ri=rnum-1;
			if(strcmp(D[0].M.r[ri].N,res)!=0){
				sprintf(whinetext,"residue name mismatch reading DOCKED: ATOM (run %d)",(di+1));
				mywhine(whinetext); }
			ai=AI[anum-1]; // atom index
if(localdebug>3){printf("\tload_dlg_mol: anum=%d ; ri=%d ; ai=%d ; di=%d.\n",anum,ri,ai,di);}
			D[0].M.r[ri].a[ai].xa[di].i=x;
			D[0].M.r[ri].a[ai].xa[di].j=y;
			D[0].M.r[ri].a[ai].xa[di].k=z;
			ti=D[0].M.r[ri].a[ai].t; // for convenience
			// add to the COM determinations
			D[0].M.rc[di].i+=D[0].M.r[ri].a[ai].xa[di].i*AT[ti].m; // molecule COM
			D[0].M.rc[di].j+=D[0].M.r[ri].a[ai].xa[di].j*AT[ti].m;
			D[0].M.rc[di].k+=D[0].M.r[ri].a[ai].xa[di].k*AT[ti].m;
			D[0].M.r[ri].rc[di].i+=D[0].M.r[ri].a[ai].xa[di].i*AT[ti].m; // residue COM
			D[0].M.r[ri].rc[di].j+=D[0].M.r[ri].a[ai].xa[di].j*AT[ti].m;
			D[0].M.r[ri].rc[di].k+=D[0].M.r[ri].a[ai].xa[di].k*AT[ti].m;
			}
		if(strstr(line,"DOCKED: ENDMDL")!=NULL){ 
			// finish COM calculation for molecule
if(localdebug>3){printf("\tload_dlg_mol: found a set; in COM calc.\n");}
			if(nat!=D[0].M.na){
				sprintf(whinetext,"didn't find enough atoms in run %d",(di+1));
				mywhine(whinetext);}
			D[0].M.rc[di].i/=D[0].M.m;
			D[0].M.rc[di].j/=D[0].M.m;
			D[0].M.rc[di].k/=D[0].M.m;
//printf("center of mass (MW=%20.12e) is at: %20.12e %20.12e %20.12e \n",D[0].M.COM.i,D[0].M.COM.j,D[0].M.COM.k);
			// also for the residues -- plus record atom number correlation
			for(a=0;a<D[0].M.nr;a++){ 
				D[0].M.r[a].rc[di].i/=D[0].M.r[a].m;
				D[0].M.r[a].rc[di].j/=D[0].M.r[a].m;
				D[0].M.r[a].rc[di].k/=D[0].M.r[a].m;
				}
			break; // leave this while read loop
			}
		} // close while read DOCKED info
		} // close if FINAL LAMARCKIAN etc.
	} // close overall read through file while-fgets line

/// Add BOX information to the molecule structure.
D[0].M.nBOX=1;
D[0].M.BOX=(boxinfo*)calloc(1,sizeof(boxinfo));
D[0].M.BOX[0].STYPE=strdup("non-periodic");
D[0].M.BOX[0].GTYPE=strdup("rectangular, defined by minimum and maximum x,y,z coordinate values");
D[0].M.BOX[0].nC=D[0].M.BOX[0].nCD=2; 
D[0].M.BOX[0].C=(coord_nD*)calloc(D[0].M.BOX[0].nC,sizeof(coord_nD));
D[0].M.BOX[0].CD=(char**)calloc(D[0].M.BOX[0].nCD,sizeof(char*));
D[0].M.BOX[0].C[0].nD=D[0].M.BOX[0].C[1].nD=3;
D[0].M.BOX[0].C[0].D=(double*)calloc(3,sizeof(double));
D[0].M.BOX[0].C[1].D=(double*)calloc(3,sizeof(double));
D[0].M.BOX[0].CD[0]=strdup("Lowest x, y and z values in the data set.");
D[0].M.BOX[0].C[0].D[0]=xl;
D[0].M.BOX[0].C[0].D[1]=yl;
D[0].M.BOX[0].C[0].D[2]=zl;
D[0].M.BOX[0].CD[1]=strdup("Highest x, y and z values in the data set.");
D[0].M.BOX[0].C[1].D[0]=xh;
D[0].M.BOX[0].C[1].D[1]=yh;
D[0].M.BOX[0].C[1].D[2]=zh;

//D[0].M.boxl.i=xl; // removed on 20080813 by BLFoley
//D[0].M.boxl.j=yl;
//D[0].M.boxl.k=zl;
//D[0].M.boxh.i=xh;
//D[0].M.boxh.j=yh;
//D[0].M.boxh.k=zh;

if(localdebug==2){
	dXprint_molecule(&D[0].M,5);
	//dXprint_coord_3D(&D[0].M.boxl); // removed on 20080813 by BLFoley
	//dXprint_coord_3D(&D[0].M.boxh);
	}
if(localdebug>2){dXprint_molecule(&D[0].M,150);}

if(localdebug>1){
printf("Molecule COM (main) is %20.15e %20.15e %20.15e\n",D[0].M.COM.i,D[0].M.COM.j,D[0].M.COM.k);
printf("\tThere are %d alternate COM's:\n",D[0].M.nrc);
for(b=0;b<D[0].M.nrc;b++){
printf("\t\tCOM (alt %d) is %20.15e %20.15e %20.15e\n",b,D[0].M.rc[b].i,D[0].M.rc[b].j,D[0].M.rc[b].k);
}
printf("\tThere are %d Residues:\n",D[0].M.nr); 
for(a=0;a<D[0].M.nr;a++){
printf("\t\tmain COM (RES %d) is %20.15e %20.15e %20.15e\n",a,D[0].M.r[a].COM.i,D[0].M.r[a].COM.j,D[0].M.r[a].COM.k);
}
printf("\t ...each with %d alternate COM's:\n",D[0].M.nrc); 
for(a=0;a<D[0].M.nr;a++){
for(b=0;b<D[0].M.r[a].nrc;b++){
printf("\t\tres %d alt %d : %20.15e %20.15e %20.15e\n",a,b,D[0].M.r[a].rc[b].i,D[0].M.r[a].rc[b].j,D[0].M.r[a].rc[b].k); 
}
}
}
free(cid);
return D;
}
