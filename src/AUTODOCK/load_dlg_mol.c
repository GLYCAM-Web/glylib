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
#include "../inc/load_pdb.h"

//dockinfo *load_dlg_mol(fileset F,atype *AT);

dockinfo *load_dlg_mol(fileset F,types *T){

int a=0,b=0,ndock=0,ai=0,di=0,ri=0,ti=0,distVer=0,startRanking=0;
int runNum[1], dumi[1];
int *AI,nat=0,localdebug=1;
char line[501],dum1[50]; // maybe no lines longer than 500...
char whinetext[501]; char* ptr;
double xl=0,xh=0,yl=0,yh=0,zl=0,zh=0,fullVer=0;
double dumd[1];
atype* AT = T[0].a;
dockinfo *D;
fileslurp sl;
molecule alt;
//dumd=(double*)calloc(1,sizeof(double));
//dumi=(int*)calloc(1,sizeof(int));
//runNum=(int*)calloc(1,sizeof(int));
//First, the function determines which version of Autodock generated the file
while(fgets(line,500,F.F)!=NULL){
	//Find the line with the Autodock version on is
	if(strstr(line,"AutoDock")!=NULL){
		ptr = strstr(line,"AutoDock");
                ptr = ptr+9;
		distVer = sscanf(ptr,"%lf",&fullVer);
		if(distVer!=1){read_eek("location of Autodock version",F.N);}
		distVer = floor(fullVer);
		break;
	}
}
//If this version has not been taken into account, then throw an error
if(distVer != 3 && distVer != 4){
	mywhine("Docking log file does not match known versions.  Must be either version 3.X or 4.X");}

if(localdebug>0){printf("load_dlg_mol: At top.\n");}
//Initialize the dockinfo, and find the input pdb in the file
D=(dockinfo*)calloc(1,sizeof(dockinfo));
D[0].DOCK_PROGRAM = strdup("Autodock");	//Set autodock as the name of the docking program
D[0].VERSION = (char*)calloc(5,sizeof(char));
D[0].numClusters=0;
sprintf(D[0].VERSION,"%.2lf",fullVer);	//Set the version number
rewind(F.F);
sl = slurp_file(F);			//Called from fileslurp.c
sl = isolateInputPDB(sl);		//Called from load_pdb.c
D[0].M = load_pdb_from_slurp(sl);	//Called from load_pdb.c
//Set the boundaries for the box to the values of the first atom
//So comparisons can be made later
xl=xh=D[0].M.r[0].a[0].x.i;
yl=yh=D[0].M.r[0].a[0].x.j;
zl=zh=D[0].M.r[0].a[0].x.k;
for(ri = 0; ri < D[0].M.nr; ri++){
	for(ai = 0; ai < D[0].M.r[ri].na; ai++){
		for(ti = 0; ti < T[0].na; ti++)//For finding the atype that matches
			if(D[0].M.r[ri].a[ai].N[0]==AT[ti].NT[0])//the atom
				break;
		//For finding the boundaries of the box
		if(xl>D[0].M.r[ri].a[ai].x.i){xl=D[0].M.r[ri].a[ai].x.i;}
		if(xh<D[0].M.r[ri].a[ai].x.i){xh=D[0].M.r[ri].a[ai].x.i;}
		if(yl>D[0].M.r[ri].a[ai].x.j){yl=D[0].M.r[ri].a[ai].x.j;}
		if(yh<D[0].M.r[ri].a[ai].x.j){yh=D[0].M.r[ri].a[ai].x.j;}
		if(zl>D[0].M.r[ri].a[ai].x.k){zl=D[0].M.r[ri].a[ai].x.k;}
		if(zh<D[0].M.r[ri].a[ai].x.k){zh=D[0].M.r[ri].a[ai].x.k;}
		D[0].M.r[ri].a[ai].t = ti;
		D[0].M.r[ri].m+=AT[ti].m; // add to residue MW
		D[0].M.r[ri].COM.i+=D[0].M.r[ri].a[ai].x.i*AT[ti].m; // residue COM
		D[0].M.r[ri].COM.j+=D[0].M.r[ri].a[ai].x.j*AT[ti].m;
		D[0].M.r[ri].COM.k+=D[0].M.r[ri].a[ai].x.k*AT[ti].m;
		D[0].M.m+=AT[ti].m; // add to molecule MW
		D[0].M.COM.i+=D[0].M.r[ri].a[ai].x.i*AT[ti].m; // molecule COM
		D[0].M.COM.j+=D[0].M.r[ri].a[ai].x.j*AT[ti].m;
		D[0].M.COM.k+=D[0].M.r[ri].a[ai].x.k*AT[ti].m;
	}//End loop through all atoms in this residue
	nat += D[0].M.r[ri].na;//Generates overall atom total
}//End loop through all residues in this molecule



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
D[0].USE=(double*)calloc(ndock,sizeof(double));// (4) Unbound System's Energy
D[0].M.nrc=ndock;
D[0].M.rc=(coord_3D*)calloc(ndock,sizeof(coord_3D)); // for alt struct mol COM's
D[0].clusterRank=(int*)calloc(ndock,sizeof(int));// Cluster Rank for the corresponding run
for(a=0;a<D[0].M.nr;a++){
	D[0].M.r[a].nrc=ndock;
	D[0].M.r[a].rc=(coord_3D*)calloc(ndock,sizeof(coord_3D)); // for alt struct res COM's
	for(b=0;b<D[0].M.r[a].na;b++){
		D[0].M.r[a].a[b].nalt=ndock;
		D[0].M.r[a].a[b].xa=(coord_3D*)calloc(ndock,sizeof(coord_3D));
		}
	}

//  This part of the function finds the energies for each run, as well as the
//  coordinates for each run, and store it in the dockinfo molecule
di = 0;
while(fgets(line,500,F.F)!=NULL){
	if(strstr(line,"FINAL LAMARCKIAN GENETIC ALGORITHM DOCKED STATE")!=NULL){
		//Call functions, based on which version of Autodock, to get energies 
		if(distVer == 3){findAD3Energies(F,D,di);}
		else if(distVer == 4){findAD4Energies(F,D,di);}
		for(a=0; a<sl.n; a++)	//Free up the memory taken up by the fileslurp
			free(sl.L[a]);
		sl.n = 0;
		//Read in lines from the file, extract any info that is here, then
		//dump all of the lines into the fileslurp
		while(strstr(line,"____________________________________________________________________")==NULL){
			fgets(line,500,F.F);
			//First check for possible data
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
			//Then add the line to the fileslurp
			sl.n++;
			sl.L=(char**)realloc(sl.L,sl.n*sizeof(char*));
			sl.L[sl.n-1]=strdup(line);
		}
		//Now Process the fileslurp lines to remove excess information...
		sl = isolateDockedPDB(sl);
		//...and extract the pdb data from it
		alt = load_pdb_from_slurp(sl);
		for(ri = 0; ri < D[0].M.nr; ri++){
			for(ai = 0; ai < D[0].M.r[ri].na; ai++){
				//For finding the boundaries of the box
				if(xl>alt.r[ri].a[ai].x.i){xl=alt.r[ri].a[ai].x.i;}
				if(xh<alt.r[ri].a[ai].x.i){xh=alt.r[ri].a[ai].x.i;}
				if(yl>alt.r[ri].a[ai].x.j){yl=alt.r[ri].a[ai].x.j;}
				if(yh<alt.r[ri].a[ai].x.j){yh=alt.r[ri].a[ai].x.j;}
				if(zl>alt.r[ri].a[ai].x.k){zl=alt.r[ri].a[ai].x.k;}
				if(zh<alt.r[ri].a[ai].x.k){zh=alt.r[ri].a[ai].x.k;}
				//Assign the alternate data
				D[0].M.r[ri].a[ai].xa[di].i=alt.r[ri].a[ai].x.i;
				D[0].M.r[ri].a[ai].xa[di].j=alt.r[ri].a[ai].x.j;
				D[0].M.r[ri].a[ai].xa[di].k=alt.r[ri].a[ai].x.k;
				ti=D[0].M.r[ri].a[ai].t; // for convenience
				// add to the COM determinations
				D[0].M.rc[di].i+=D[0].M.r[ri].a[ai].xa[di].i*AT[ti].m; // molecule COM
				D[0].M.rc[di].j+=D[0].M.r[ri].a[ai].xa[di].j*AT[ti].m;
				D[0].M.rc[di].k+=D[0].M.r[ri].a[ai].xa[di].k*AT[ti].m;
				D[0].M.r[ri].rc[di].i+=D[0].M.r[ri].a[ai].xa[di].i*AT[ti].m; // residue COM
				D[0].M.r[ri].rc[di].j+=D[0].M.r[ri].a[ai].xa[di].j*AT[ti].m;
				D[0].M.r[ri].rc[di].k+=D[0].M.r[ri].a[ai].xa[di].k*AT[ti].m;

			}//End loop through atoms
		// also for the residues -- plus record atom number correlation
			D[0].M.r[ri].rc[di].i/=D[0].M.r[ri].m;
			D[0].M.r[ri].rc[di].j/=D[0].M.r[ri].m;
			D[0].M.r[ri].rc[di].k/=D[0].M.r[ri].m;
		}//End loop through residues
		// finish COM calculation for molecule
		D[0].M.rc[di].i/=D[0].M.m;
		D[0].M.rc[di].j/=D[0].M.m;
		D[0].M.rc[di].k/=D[0].M.m;
		deallocateMolecule(&alt);//Free the memory allocated in alt
		di++; 			//and incriment to the next alternate set
	}
	if(strstr(line,"CLUSTER ANALYSIS OF CONFORMATIONS") != NULL){
		for(a=0; a<sl.n; a++)	//Free all of the data in the fileslurp
			free(sl.L[a]);	//since it isn't used again
		free(sl.L); sl.n = 0;		
		break;
	}
}

//Now deal with all of the state variables
di = 0;//Reset the alternate place holder
while(fgets(line,500,F.F)!=NULL){
	if(strstr(line,"STATE VARIABLES:") != NULL){
		while(strstr(line,"Torsions") == NULL){
			fgets(line,500,F.F);
			if(strstr(line,"Quaternion nx,ny,nz,angle")!=NULL){
				sscanf(line,"%s %s %s %lf %lf %lf %lf",dum1,dum1,dum1,\
					&D[0].QN[di].i,&D[0].QN[di].j,&D[0].QN[di].k,&D[0].QN[di].d); 
				}
			if(strstr(line,"Quaternion x,y,z,w ")!=NULL){
				sscanf(line,"%s %s %s %lf %lf %lf %lf",dum1,dum1,dum1,\
					&D[0].Q[di].i,&D[0].Q[di].j,&D[0].Q[di].k,&D[0].Q[di].d); 
				}
		}
		di++;
	}
	if(strstr(line,"RMSD TABLE") != NULL){
		while(strstr(line,"_____|_") == NULL){
		fgets(line,500,F.F);
		}
		startRanking=1;
	}
	if(startRanking==1){
		startRanking++;
		fgets(line,500,F.F);
		sscanf(line,"%d %lf %d %lf %lf %lf %s",dumi,dumd,runNum,dumd,dumd,dumd,dum1);
		D[0].clusterRank[runNum[0]-1]=dumi[0];
//printf("load_dlg_mol: 1.) Run %d: Cluster Rank: %d\n",runNum[0],D[0].clusterRank[runNum[0]-1]);
		while(strstr(line,"_______") == NULL){
			fgets(line,500,F.F);
			sscanf(line,"%d %lf %d %lf %lf %lf %s",dumi,dumd,runNum,dumd,dumd,dumd,dum1);
	                D[0].clusterRank[runNum[0]-1]=dumi[0];
//printf("load_dlg_mol: 2.) Run %d: Cluster Rank: %d\n",runNum[0],D[0].clusterRank[runNum[0]-1]);
		}
	}
	}

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

if(localdebug==2){dXprint_molecule(&D[0].M,5);}
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
return D;
}

void findAD3Energies(fileset F,dockinfo* D,int di){
char line[501]; char dum1[50];
fgets(line,500,F.F);
while(strstr(line,"NEWDPF")==NULL){
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
	fgets(line,500,F.F);
}
return ;
}

void findAD4Energies(fileset F,dockinfo* D,int di){
char line[501]; char dum1[50];
fgets(line,500,F.F);
while(strstr(line,"NEWDPF")==NULL){
	if(strstr(line,"Estimated Free Energy of Binding ")!=NULL){
		sscanf(line,"%s %s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,\
			dum1,dum1,dum1,dum1,&D[0].eFEB[di]);
	}
	if(strstr(line,"Estimated Inhibition Constant, Ki")!=NULL){
		sscanf(line,"%s %s %s %s %s %s %s %lf %s %s %lf",dum1,dum1,dum1,dum1,\
			dum1,dum1,dum1,&D[0].eKi[di],dum1,dum1,&D[0].Tmp[di]);
	}
	if(strstr(line,"Final Intermolecular Energy")!=NULL){
		sscanf(line,"%s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,dum1,\
			dum1,dum1,&D[0].fIE[di]);
	}
	if(strstr(line,"Final Total Internal Energy")!=NULL){
		sscanf(line,"%s %s %s %s %s %s %s %s %lf",dum1,dum1,dum1,\
			dum1,dum1,dum1,dum1,dum1,&D[0].fIEL[di]);
	}
	if(strstr(line,"Torsional Free Energy")!=NULL){
		sscanf(line,"%s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,\
			dum1,dum1,dum1,&D[0].TFE[di]);
	}
	if(strstr(line,"Unbound System's Energy")!=NULL){
		sscanf(line,"%s %s %s %s %s %s %s %lf",dum1,dum1,dum1,dum1,dum1,dum1,\
			dum1,&D[0].USE[di]);
	}
	fgets(line,500,F.F);
}
return ;
}

