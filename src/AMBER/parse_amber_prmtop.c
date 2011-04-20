/* File parse_amber_prmtop.c begun on 20080610 by BLFoley
 * Purpose: parse entries from prmtop file (already read in) into
 * 	appropriate data locations within an assembly
 * 	Also adds the prmtop P to the assembly's void pointer
 */
#include "AMBER/amber_prmtop.h"
#include "gly_codeutils.h"

assembly parse_amber_prmtop(amber_prmtop *P){
assembly A;
int *ICO,nICO;
int pa=0,pb=0,pc=0,pd=0,pA1=0,pA2=0,pA3,pA4,pI1=0;
int NextRes=0,nummol=0,*resi; // more utility integers
int sa=0,sr=0,sm=0,ta=0,tr=0,tm=0,ntrue=0,nwhitemid=0;
int IPTRES=0,NSPM=0,NSPSOL=0,*NSP; // solvent/solute info & #atoms per molecule
molbond *MB,*MBTMP;
angle_index *MANG;
torsion_index *MTOR;
molindex *MOLI,*MOLBNDI; // MOLI for here, MOLBNDI for assigning molecules
char **ATNAME,**TREECLASS,*RADTYPE,tmp[80],*PRUNED;
double *R,*SC,*MASS; // radii and screening constants for IS, atom masses

int compare_array_ws(int ai, char **a, char *flag, int bi){ // Function which checks for whitespace in an array and removes it to compare the length to its expected(true) length 
/*	printf("ai is %d\n",ai);
	printf("a[0] is %s\n",a[0]);
	printf("flag is %s\n",flag);
	printf("bi is %d\n",bi);
*/
	int ntrue=0,nwhitemid=0,i=0;// ntrue is the number of non-whitespace fields, nwhitemid is the number of non-terminal whitespace residues (>0 will cause error)
	for(i=0;i<ai;i++){
		PRUNED=strdup(prune_string_whitespace(a[i])); // Remove whitespace from array elements
		if(PRUNED[0]!='\0'){ // Check for only a string terminator
			if(ntrue==bi){ // Check that ntrue never exceeds the expected number of elements (otherwise exit with an error)
	                        printf("ERROR: The declared number of components, %d, in FLAG %s does not equal the number found, %d!\n",(ntrue+1),flag,bi);
				mywhine("\tExiting for an unequal number of elements between expected and found.\n");
				}
			if(i>(ntrue+nwhitemid)){ nwhitemid++; } // Increment the non-terminal whitespace index (note that consecutive non-terminal whitespaces are counted as 1)
			ntrue++;
			}
		free(PRUNED);
		}
	if(ntrue!=ai){ // Test if the expected (ai) value does not equal the whitespace-removed value
		printf("WARNING: Found whitespace at flag %s\n",flag);
		printf("\tFound %d total",(ai-ntrue));
			if(nwhitemid>0){
				printf(" and at least %d non-terminal whitespace(s).\nFATAL WARNING: Non-terminal whitespace(s) can allow for mis-assignment of topology components!\n\tCheck that your topology file is correctly built.\n",nwhitemid);
				mywhine("\tExiting because non-terminal whitespace found.\n");
				}
			else{ printf(" terminal whitespace(s).\n\tIgnoring whitespace(s).\n");}
		}
	return(ntrue);
	}

int count_array_nws(int ai, char **a, char *flag){ // Function to return the number of non-whitespace elements in an array - no comparison to other numbers
        int ntrue=0,nwhitemid=0,i=0;// ntrue is the number of non-whitespace fields, nwhitemid is the number of non-terminal whitespace residues (>0 will cause error)
        for(i=0;i<ai;i++){
                PRUNED=strdup(prune_string_whitespace(a[i])); // Removing whitespace from the read-in data at position pb 
                if(PRUNED[0]!='\0'){ // Not whitespace
                        if(i>(ntrue+nwhitemid)){ nwhitemid++; } // Increment the number of non-terminal whitespaces
                        ntrue++;
                        }
                free(PRUNED);
                }
        if(ntrue!=ai){
                printf("WARNING: Found whitespace at flag %s\n",flag);
                printf("\tFound %d total",(ai-ntrue));
                        if(nwhitemid>0){printf(" and at least %d non-terminal whitespace(s).\nFATAL WARNING: Non-terminal whitespace(s) can allow for mis-assignment of topology components!\n\tCheck that your topology file is correctly built.\n",nwhitemid);}
                        else{printf(" terminal whitespace(s).\n");}
                printf("\tIgnoring whitespace(s).\n");
                }
	return(ntrue);
	}

//fileset F;
//amber_prmtop *aprm;

MOLI=NULL;
ICO=NULL;
MB=NULL;
MANG=NULL;
MTOR=NULL;
MASS=NULL;
NSP=NULL;
R=NULL;
SC=NULL;
MOLBNDI=NULL;
MBTMP=NULL;
resi=NULL;
ATNAME=NULL;
TREECLASS=NULL;

/// Drop the contents of the original file in the void pointer space
A.nVP=1; 
A.VP=P;

// do some initializations
A.nPRM=1;
A.PRM=(parameter_set**)calloc(1,sizeof(parameter_set*));
A.PRM[0]=(parameter_set*)calloc(1,sizeof(parameter_set));

// First, find the pointers and read them into the top structure
// Also add them to the assembly, as needed
// Any that are not added now can be added later -- chances are that 
// 	no one had a use for them previously.  In any case, they will
// 	be present as-is in the void pointer
for(pa=0;pa<P[0].nS;pa++){ // get pointers first...
	if(strcmp(P[0].SN[pa],"POINTERS")==0){ 
		// total number of atoms
		sscanf(P[0].S[pa].D[0],"%d",&P[0].NATOM); 
		A.na = P[0].NATOM; 
		if(A.na>0){
			A.a=(atom**)calloc(A.na,sizeof(atom*));
			MOLI=(molindex*)calloc(A.na,sizeof(molindex)); // to keep up with locations per atom
			} 
		for(pb=0;pb<A.na;pb++){A.a[pb]=(atom*)calloc(1,sizeof(atom));}
		// number of atom types
  		sscanf(P[0].S[pa].D[1],"%d",&P[0].NTYPES); 
  		A.PRM[0][0].nAT=P[0].NTYPES; 
		if(A.PRM[0][0].nAT>0){
			A.PRM[0][0].AT=(atype*)calloc(A.PRM[0][0].nAT,sizeof(atype));
			A.PRM[0][0].nNBT=A.PRM[0][0].nAT*(A.PRM[0][0].nAT+1)/2;
			A.PRM[0][0].NBT=(bond_type*)calloc(A.PRM[0][0].nNBT,sizeof(bond_type));
			nICO=A.PRM[0][0].nAT*A.PRM[0][0].nAT;
			ICO=(int*)calloc(nICO,sizeof(int));
			} 
		// Numbers of actual bonds, angles, dihedrals
		// NOT adding these number-of distinctions to the assembly
		// Instead, adding them all as plain bonds, etc. (see "D" description for distinction)
		// See also NBONA, etc., below (constraints)...
  		sscanf(P[0].S[pa].D[2],"%d",&P[0].NBONH); // number of bonds containing hydrogen
  		sscanf(P[0].S[pa].D[3],"%d",&P[0].MBONA); // number of bonds not containing hydrogen
  		sscanf(P[0].S[pa].D[4],"%d",&P[0].NTHETH); // number of angles containing hydrogen
  		sscanf(P[0].S[pa].D[5],"%d",&P[0].MTHETA); // number of angles not containing hydrogen
  		sscanf(P[0].S[pa].D[6],"%d",&P[0].NPHIH); // number of dihedrals containing hydrogen
  		sscanf(P[0].S[pa].D[7],"%d",&P[0].MPHIA); // number of dihedrals not containing hydrogen
		// Unused parameters and the excluded atoms list (probably not used in glycam)
		// NOT adding these at all
  		sscanf(P[0].S[pa].D[8],"%d",&P[0].NHPARM); // currently not used
  		sscanf(P[0].S[pa].D[9],"%d",&P[0].NPARM); // currently not used
  		sscanf(P[0].S[pa].D[10],"%d",&P[0].NEXT); // number of excluded atoms
		// Number of residues
  		sscanf(P[0].S[pa].D[11],"%d",&P[0].NRES); 
  		A.nr = P[0].NRES; 
		if(A.nr>0){A.r=(residue**)calloc(A.nr,sizeof(residue*));} // will move these later 
		for(pb=0;pb<A.nr;pb++){A.r[pb]=(residue*)calloc(1,sizeof(residue));}
		// The rest of the bonds, angles and torsions (constraints)
  		sscanf(P[0].S[pa].D[12],"%d",&P[0].NBONA); // MBONA + number of constraint bonds
		A.nb = P[0].NBONH + P[0].NBONA; // all the bonds
		if(A.nb>0){A.b=(molbond*)calloc(A.nb,sizeof(molbond));}
		if(A.nb>0){MB=(molbond*)calloc(A.nb,sizeof(molbond));}
  		sscanf(P[0].S[pa].D[13],"%d",&P[0].NTHETA); // MTHETA + number of constraint angles
		A.nANG=P[0].NTHETH+P[0].NTHETA; // all the angles
		if(A.nANG>0){A.ANG=(angle_index*)calloc(A.nANG,sizeof(angle_index));}
		if(A.nANG>0){MANG=(angle_index*)calloc(A.nANG,sizeof(angle_index));}
  		sscanf(P[0].S[pa].D[14],"%d",&P[0].NPHIA); // MPHIA + number of constraint dihedrals
		A.nTOR=P[0].NPHIH+P[0].NPHIA; // all the torsions
		if(A.nTOR>0){A.TOR=(torsion_index*)calloc(A.nTOR,sizeof(torsion_index));}
		if(A.nTOR>0){MTOR=(torsion_index*)calloc(A.nTOR,sizeof(torsion_index));}
		// Numbers of bond, angle and torsion types
  		sscanf(P[0].S[pa].D[15],"%d",&P[0].NUMBND); // number of unique bond types
  		A.PRM[0][0].nBT = P[0].NUMBND; // number of unique bond types
		if(A.PRM[0][0].nBT>0){A.PRM[0][0].BT=(bond_type*)calloc(A.PRM[0][0].nBT,sizeof(bond_type));}
  		sscanf(P[0].S[pa].D[16],"%d",&P[0].NUMANG); // number of unique angle types
  		A.PRM[0][0].nANT = P[0].NUMANG; // number of unique angle types
		if(A.PRM[0][0].nANT>0){A.PRM[0][0].ANT=(angle_type*)calloc(A.PRM[0][0].nANT,sizeof(angle_type));}
  		sscanf(P[0].S[pa].D[17],"%d",&P[0].NPTRA); // number of unique dihedral types
  		A.PRM[0][0].nTRT = P[0].NPTRA; // number of unique dihedral types
		if(A.PRM[0][0].nTRT>0){A.PRM[0][0].TRT=(torsion_type*)calloc(A.PRM[0][0].nTRT,sizeof(torsion_type));}
		// Number of atom types in parameter file
		// NOT adding this at all at present
  		sscanf(P[0].S[pa].D[18],"%d",&P[0].NATYP); // number of atom types in parameter file, see SOLTY below
		// Number of Lennard-Jones 10-12 H-Bond pair types
  		sscanf(P[0].S[pa].D[19],"%d",&P[0].NPHB); // number of distinct 10-12 hydrogen bond pair types
		A.PRM[0][0].nHBT = P[0].NPHB;
		if(A.PRM[0][0].nHBT>0){A.PRM[0][0].HBT=(bond_type*)calloc(A.PRM[0][0].nHBT,sizeof(bond_type));}
		// Perturbation information
		// NOT adding these at all (no longer used as of AMBER 10)
  		sscanf(P[0].S[pa].D[20],"%d",&P[0].IFPERT); // set to 1 if perturbation info is to be read in
  		sscanf(P[0].S[pa].D[21],"%d",&P[0].NBPER); // number of bonds to be perturbed
  		sscanf(P[0].S[pa].D[22],"%d",&P[0].NGPER); // number of angles to be perturbed
  		sscanf(P[0].S[pa].D[23],"%d",&P[0].NDPER); // number of dihedrals to be perturbed
  		sscanf(P[0].S[pa].D[24],"%d",&P[0].MBPER); // number of bonds with atoms completely in perturbed group
  		sscanf(P[0].S[pa].D[25],"%d",&P[0].MGPER); // number of angles with atoms completely in perturbed group
  		sscanf(P[0].S[pa].D[26],"%d",&P[0].MDPER); // number of dihedrals with atoms completely in perturbed groups
		// Box information
  		sscanf(P[0].S[pa].D[27],"%d",&P[0].IFBOX); // set to 1 if standard periodic box, 2 when truncated octahedral
		if(P[0].IFBOX==0){A.nBOX=0;} ///< There are no boxes defined in this prmtop
		else{
			A.nBOX=1;
			A.BOX=(boxinfo*)calloc(A.nBOX,sizeof(boxinfo));
			A.BOX[0].nC=A.BOX[0].nCD=2; ///< One for the "boxang" and the other for the x,y,z values
			A.BOX[0].C=(coord_nD*)calloc(A.BOX[0].nC,sizeof(coord_nD));
			A.BOX[0].CD=(char**)calloc(A.BOX[0].nCD,sizeof(char*));
			A.BOX[0].CD[0]=strdup("Periodic box, angle between the XY and YZ planes in degrees.");
			A.BOX[0].C[0].nD=1; ///< Just an angle ("BETA" in the amber file specs)
			A.BOX[0].C[0].D=(double*)calloc(A.BOX[0].C[0].nD,sizeof(double));
			A.BOX[0].CD[1]=strdup("The periodic box lengths in the X, Y, and Z directions");
			A.BOX[0].C[1].nD=3; ///< 3-D coordinate
			A.BOX[0].C[1].D=(double*)calloc(A.BOX[0].C[1].nD,sizeof(double));
			A.BOX[0].STYPE=strdup("standard periodic");
			if(P[0].IFBOX==1){A.BOX[0].GTYPE=strdup("rectangular");}
			if(P[0].IFBOX==2){A.BOX[0].GTYPE=strdup("truncated octahedral");}
			if(P[0].IFBOX>2){A.BOX[0].GTYPE=strdup("unknown box type");} 
			}

		// NOT adding these at all 
  		sscanf(P[0].S[pa].D[28],"%d",&P[0].NMXRS); // number of atoms in the largest residue
  		sscanf(P[0].S[pa].D[29],"%d",&P[0].IFCAP); // set to 1 if the CAP option from edit was specified
		//
		break; // no need to keep scanning...
		}
	}


for(pa=0;pa<P[0].nS;pa++){// Loop through each of the sections
	P[0].S[pa].is_standard = 1; // Set is_standard to one as default
	// 	check against each standard section
	// 	if there is a match:
	// 		set is_standard to zero
	// 		record to assembly as needed
	//
	//		note that here 1=FALSE and 0=TRUE
	//

//FORMAT(20a4)  (ITITL(i), i=1,20)
if(strcmp(P[0].S[pa].N,"TITLE")==0){ // place in the Assembly description, NOTE: can be whitespace
        P[0].ITITL = pa; // the title section 
	P[0].S[pa].is_standard=0; // set as a standard section
	// copy the title into the assembly's description
	A.D=(char*)calloc(P[0].S[pa].nt*P[0].S[pa].nc+1,sizeof(char));
	for(pb=0;pb<P[0].S[pa].nt;pb++){strcat(A.D,P[0].S[pa].D[pb]);}
}
//FORMAT(12i6) 
if(strcmp(P[0].S[pa].N,"POINTERS")==0){ // these are already added, no whitespace check 
	P[0].POINTERS=pa; // pointer to the section containing the original char-strings
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(20a4)  (IGRAPH(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"ATOM_NAME")==0){ // these eventually go into the molecules/residues/etc
  	P[0].IGRAPH=pa; // IGRAPH : the user atoms names 
	P[0].S[pa].is_standard=0; // set as a standard section
	// for now, add these to the straight list of Assembly atoms
	compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	for(pb=0;pb<A.na;pb++){
		sscanf(P[0].S[pa].D[pb],"%s",tmp);
		A.a[pb][0].N=(char*)calloc((strlen(tmp)+1),sizeof(char));
		strcpy(A.a[pb][0].N,tmp); 
//	        printf("\tP[0].S[pa].D[pb] = %s\n",P[0].S[pa].D[pb]);
		}
}
// FORMAT(5E16.8)  (CHRG(i), i=1,NATOM)
// (Divide by 18.2223 to convert to charge in units of the electron charge) 
if(strcmp(P[0].S[pa].N,"CHARGE")==0){ // these go with each atom
  	P[0].CHRG   =pa; // CHRG   : the atom charges.  
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	for(pb=0;pb<A.na;pb++){
		A.a[pb][0].nch=1;
		A.a[pb][0].ch=(double*)calloc(1,sizeof(double));
		sscanf(P[0].S[pa].D[pb],"%lf",&A.a[pb][0].ch[0]);
		A.a[pb][0].ch[0]/=18.2223;
		}
}
// FORMAT(5E16.8)  (AMASS(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"MASS")==0){ // with each atom
  	P[0].AMASS  =pa; // AMASS  : the atom masses 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	MASS=(double*)calloc(A.na,sizeof(double));
	for(pb=0;pb<A.na;pb++){ sscanf(P[0].S[pa].D[pb],"%lf",&A.a[pb][0].m); }
	// for(pb=0;pb<A.na;pb++){ sscanf(P[0].S[pa].D[pb],"%lf",&MASS[pb]); } // before m in atom struct
}
// FORMAT(12I6)  (IAC(i), i=1,NATOM)
// START HERE -- is this really just an atom type index?  One hopes...
if(strcmp(P[0].S[pa].N,"ATOM_TYPE_INDEX")==0){ // into the atype structure
  	P[0].IAC    =pa; // IAC    : index for the atom types involved in Lennard Jones (6-12) 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	for(pb=0;pb<A.na;pb++){
		sscanf(P[0].S[pa].D[pb],"%d",&A.a[pb][0].t);
		A.a[pb][0].t--;
		if(A.a[pb][0].t<0){mywhine("A.a[pb][0].t<0 in parse_amber_prmtop");}
		if(A.a[pa][0].t>=P[0].NTYPES){mywhine("A.a[pa][0].t>=P[0].NTYPES in parse_amber_prmtop");}
		}
}
           	// interactions.  See ICO below.  
// FORMAT(12I6)  (NUMEX(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"NUMBER_EXCLUDED_ATOMS")==0){ // unused for now (20080612), no whitespace check since unused
  	P[0].NUMEX  =pa; // NUMEX  : total number of excluded atoms for atom "i".  See
	P[0].S[pa].is_standard=0; // set as a standard section
}
           	// NATEX below.  
// FORMAT(12I6)  (ICO(i), i=1,NTYPES*NTYPES)
// START HERE -- figure out how to put this is (efficiently...)
if(strcmp(P[0].S[pa].N,"NONBONDED_PARM_INDEX")==0){ // 
  	P[0].ICO    =pa; // ICO    : provides the index to the nonbon parameter
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,nICO);
	for(pb=0;pb<nICO;pb++){sscanf(P[0].S[pa].D[pb],"%d",&ICO[pb]);}
           	// arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
           	// or 10-12 atoms type interactions are represented.
           	// NOTE: A particular atom type can have either a 10-12
           	// or a 6-12 interaction, but not both.  The index is
           	// calculated as follows:
             	// index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
           	// If index is positive, this is an index into the
           	// 6-12 parameter arrays (CN1 and CN2) otherwise it
           	// is an index into the 10-12 parameter arrays (ASOL and BSOL).  
}
// FORMAT(20A4)  (LABRES(i), i=1,NRES)
if(strcmp(P[0].S[pa].N,"RESIDUE_LABEL")==0){ // names of residues
  	P[0].LABRES =pa; // LABRES : the residue labels 
	P[0].S[pa].is_standard=0; // set as a standard section
	nwhitemid=0; // Total non-terminal whitespace residues
	ntrue=0; // Number of actual (non-whitespace) residues found
	for(pb=0;pb<P[0].S[pa].nt;pb++){
                PRUNED=strdup(prune_string_whitespace(P[0].S[pa].D[pb])); // Removing whitespace from the read-in data at position pb 
//printf("pb eq %d string eq >>%s<< pruned eq >>%s<<\n",pb,P[0].S[pa].D[pb],PRUNED);
                if(PRUNED[0]!='\0'){ // Not whitespace
			if(ntrue==A.nr){ // If the true number of residues equals the total residues in the assembly
//printf("ntrue eq %d and A.nr eq %d\n",ntrue,A.nr);
				mywhine("The declared number of residues does not equal the number found! :-(");
				}
			if(pb>(ntrue+nwhitemid)){ nwhitemid++; } // Increment the number of non-terminal whitespaces, shouldn't this be a "while" loop?
			A.r[ntrue][0].N=(char*)calloc((P[0].S[pa].nc+1),sizeof(char)); // Adjusted for data size, need type-specific functionality
			sscanf(P[0].S[pa].D[pb],"%s",A.r[ntrue][0].N);
			A.r[ntrue][0].N[P[0].S[pa].nc]='\0'; // Adding the 'end of line' character to the end of the string
			ntrue++;
			}
		free(PRUNED);
		}
	if(ntrue!=P[0].S[pa].nt){
		printf("WARNING: Found whitespace at flag %s\n",P[0].S[pa].N);
		printf("\tFound %d total",(P[0].S[pa].nt-ntrue));
			if(nwhitemid>0){printf(" and at least %d non-terminal whitespace(s).\nFATAL WARNING: Non-terminal whitespace(s) can allow for mis-assignment of topology components!\n\tCheck that your topology file is correctly built.\n",nwhitemid);}
			else{printf(" terminal whitespace(s).\n");}
		printf("\tIgnoring whitespace(s).\n");
		}
}
// FORMAT(10a8)
if(strcmp(P[0].S[pa].N,"RESIDUE_ID")==0){// ?? Presumably residue numbers, will store there
	P[0].IRES =pa; // START HERE -- this is probably the wrong name!!!!!
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.nr);
	for(pb=0;pb<A.nr;pb++){sscanf(P[0].S[pa].D[pb],"%d",&A.r[pb][0].n);}
}
// FORMAT(12I6)  (IPRES(i), i=1,NRES)
if(strcmp(P[0].S[pa].N,"RESIDUE_POINTER")==0){ // for adding atoms to residues
  	P[0].IPRES  =pa; // IPRES  : atoms in each residue are listed for atom "i" in
	P[0].S[pa].is_standard=0; // set as a standard section
	pc=0;
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.nr);
	for(pb=0;pb<(A.nr-1);pb++){
		sscanf(P[0].S[pa].D[pb+1],"%d",&NextRes);
		NextRes-=1;
		A.r[pb][0].na=NextRes-pc;
		A.r[pb][0].a=(atom*)calloc(A.r[pb][0].na,sizeof(atom));
		pd=0;
		for(pc=pc;pc<NextRes;pc++){
			A.r[pb][0].a[pd]=A.a[pc][0]; // copy atom into residue
			A.a[pc]=&A.r[pb][0].a[pd]; // set atom pointer to new location
			A.a[pc][0].n=pc+1; // set the amber original number
			MOLI[pc].i=pc; // set absolute atom number in array terms
			MOLI[pc].m=-1; // we don't know this yet
			MOLI[pc].r=pb; // check this later
			MOLI[pc].a=pd; // check this later
			pd++;
			}
		}
		A.r[pb][0].na=P[0].NATOM-pc;
		A.r[pb][0].a=(atom*)calloc(A.r[pb][0].na,sizeof(atom));
		pd=0;
		for(pc=pc;pc<P[0].NATOM;pc++){
			A.r[pb][0].a[pd]=A.a[pc][0]; // copy atom into residue
			A.a[pc]=&A.r[pb][0].a[pd]; // set atom pointer to new location
			A.a[pc][0].n=pc+1; // set the amber original number
			MOLI[pc].i=pc; // set absolute atom number in array terms
			MOLI[pc].m=-1; // we don't know this yet
			MOLI[pc].r=pb; // check this later
			MOLI[pc].a=pd; // check this later
			pd++;
			}

	
}
// FORMAT(5E16.8)  (RK(i), i=1,NUMBND)
if(strcmp(P[0].S[pa].N,"BOND_FORCE_CONSTANT")==0){ // 
  	P[0].RK     =pa; // RK     : force constant for the bonds of each type, kcal/mol 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nBT);
	for(pb=0;pb<A.PRM[0][0].nBT;pb++){
		sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].BT[pb].k);
		}
}
// FORMAT(5E16.8)  (REQ(i), i=1,NUMBND)
if(strcmp(P[0].S[pa].N,"BOND_EQUIL_VALUE")==0){
  	P[0].REQ    =pa; // REQ    : the equilibrium bond length for the bonds of each type, angstroms 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nBT);
	for(pb=0;pb<A.PRM[0][0].nBT;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].BT[pb].l);}
}
// FORMAT(5E16.8)  (TK(i), i=1,NUMANG)
if(strcmp(P[0].S[pa].N,"ANGLE_FORCE_CONSTANT")==0){
  	P[0].TK     =pa; // TK     : force constant for the angles of each type, kcal/mol A**2 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nANT);
	for(pb=0;pb<A.PRM[0][0].nANT;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].ANT[pb].k);}
}
// FORMAT(5E16.8)  (TEQ(i), i=1,NUMANG)
if(strcmp(P[0].S[pa].N,"ANGLE_EQUIL_VALUE")==0){
  	P[0].TEQ    =pa; // TEQ    : the equilibrium angle for the angles of each type, radians 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nANT);
	for(pb=0;pb<A.PRM[0][0].nANT;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].ANT[pb].l);}
}
// FORMAT(5E16.8)  (PK(i), i=1,NPTRA)
if(strcmp(P[0].S[pa].N,"DIHEDRAL_FORCE_CONSTANT")==0){
  	P[0].PK     =pa; // PK     : force constant for the dihedrals of each type, kcal/mol 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nTRT);
	for(pb=0;pb<A.PRM[0][0].nTRT;pb++){
		if(A.PRM[0][0].TRT[pb].n>1){mywhine("A.PRM[0][0].TRT[pb].n>1 in parse_amber_prmtop!");}
		A.PRM[0][0].TRT[pb].n=1;
		A.PRM[0][0].TRT[pb].k=(double*)calloc(1,sizeof(double));
		sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].TRT[pb].k[0]);
		}
}
// FORMAT(5E16.8)  (PN(i), i=1,NPTRA)
if(strcmp(P[0].S[pa].N,"DIHEDRAL_PERIODICITY")==0){
  	P[0].PN     =pa; // PN     : periodicity of the dihedral of a given type 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nTRT);
	for(pb=0;pb<A.PRM[0][0].nTRT;pb++){
		if(A.PRM[0][0].TRT[pb].n>1){mywhine("A.PRM[0][0].TRT[pb].n>1 in parse_amber_prmtop!");}
		A.PRM[0][0].TRT[pb].n=1;
		A.PRM[0][0].TRT[pb].N=(double*)calloc(1,sizeof(double));
		sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].TRT[pb].N[0]);
		}
}
// FORMAT(5E16.8)  (PHASE(i), i=1,NPTRA)
if(strcmp(P[0].S[pa].N,"DIHEDRAL_PHASE")==0){
  	P[0].PHASE  =pa; // PHASE  : phase of the dihedral of a given type, radians 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nTRT);
	for(pb=0;pb<A.PRM[0][0].nTRT;pb++){
		if(A.PRM[0][0].TRT[pb].n>1){mywhine("A.PRM[0][0].TRT[pb].n>1 in parse_amber_prmtop!");}
		A.PRM[0][0].TRT[pb].n=1;
		A.PRM[0][0].TRT[pb].P=(double*)calloc(1,sizeof(double));
		sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].TRT[pb].P[0]);
		}
}
// FORMAT(5E16.8)  (SOLTY(i), i=1,NATYP)
if(strcmp(P[0].S[pa].N,"SOLTY")==0){ // not much to do here at the moment, no whitespace check
  	P[0].SOLTY  =pa; // SOLTY  : currently unused (reserved for future use) 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(5E16.8)  (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
if(strcmp(P[0].S[pa].N,"LENNARD_JONES_ACOEF")==0){
 	P[0].CN1    =pa; // CN1    : Lennard Jones r**12 terms for all possible atom type
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nNBT);
	// record the LJ parameters to the bond type info
	for(pb=0;pb<A.PRM[0][0].nNBT;pb++){ sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].NBT[pb].LJ12_612); }
           	// interactions, indexed by ICO and IAC; for atom i and j
           	// where i < j, the index into this array is as follows
           	// (assuming the value of ICO(index) is positive):
           	// CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).  
}
// FORMAT(5E16.8)  (CN2(i), i=1,NTYPES*(NTYPES+1)/2)
if(strcmp(P[0].S[pa].N,"LENNARD_JONES_BCOEF")==0){
  	P[0].CN2    =pa; // CN2    : Lennard Jones r**6 terms for all possible atom type
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.PRM[0][0].nNBT);
	for(pb=0;pb<A.PRM[0][0].nNBT;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].NBT[pb].LJ6_612);}
           	// interactions.  Indexed like CN1 above.  
}
/* NOTE: the atom numbers in the following arrays that describe bonds, 
	angles, and dihedrals are coordinate array indexes for runtime speed. 
	The true atom number equals the absolute value of the number divided by 
	three, plus one. In the case of the dihedrals, if the fourth atom is negative, 
	this implies that the dihedral is an improper. If the third atom is negative, 
	this implies that the end group interations are to be ignored. End group 
	interactions are ignored, for example, in dihedrals of various ring systems 
	(to prevent double counting of 1-4 interactions) and in multiterm dihedrals.  */
// FORMAT(12I6)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
if(strcmp(P[0].S[pa].N,"BONDS_INC_HYDROGEN")==0){
  	P[0].IBH    =pa; // IBH    : atom involved in bond "i", bond contains hydrogen
  	P[0].JBH    =pa; // JBH    : atom involved in bond "i", bond contains hydrogen
  	P[0].ICBH   =pa; // ICBH   : index into parameter arrays RK and REQ 
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,(3*P[0].NBONH));
	for(pb=0;pb<P[0].NBONH;pb++){
		// read in the values from the section structure
		sscanf(P[0].S[pa].D[3*pb],"%d",&pA1);
		sscanf(P[0].S[pa].D[3*pb+1],"%d",&pA2);
		sscanf(P[0].S[pa].D[3*pb+2],"%d",&pI1); 
		// find the actual atom numbers
		pA1/=3;  // CAREFUL !! these can be negative
		pA2/=3; 
		// set the bond info for these two atoms
		// -- at this point, everything is in the same molecule, so the molecule index
		// 	will be set to -1 (to prevent confusion later).
		// 	Also, since we don't know the ordering in the file, most info will not be set
		MB[pb].i=pI1; // holder for the index into the RK & REQ arrays
		MB[pb].D=strdup("BONDS_INC_HYDROGEN");
		MB[pb].s.m=-1; // we don't know the molecule number yet
		MB[pb].s.r=-1; // we can't assume we know the residue number yet
		MB[pb].s.a=pA1; // 
		MB[pb].s.i=pA1; // 
		MB[pb].t.m=-1; 
		MB[pb].t.r=-1; 
		MB[pb].t.i=pA2; 
		MB[pb].t.a=pA2; 
		}
}
// FORMAT(12I6)  (IB(i),JB(i),ICB(i), i=1,NBONA)
if(strcmp(P[0].S[pa].N,"BONDS_WITHOUT_HYDROGEN")==0){
  	P[0].IB     =pa; // IB     : atom involved in bond "i", bond does not contain hydrogen
  	P[0].JB     =pa; // JB     : atom involved in bond "i", bond does not contain hydrogen
  	P[0].ICB    =pa; // ICB    : index into parameter arrays RK and REQ 
	P[0].S[pa].is_standard=0; // set as a standard section
	//for(pb=P[0].NBONH;pb<(P[0].NBONH+P[0].MBONA);pb++) // save this line for use later...
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,(3*P[0].MBONA));
	for(pb=0;pb<P[0].MBONA;pb++){ // save this line for use later...
		// read in the values from the section structure
		sscanf(P[0].S[pa].D[3*pb],"%d",&pA1);
		sscanf(P[0].S[pa].D[3*pb+1],"%d",&pA2);
		sscanf(P[0].S[pa].D[3*pb+2],"%d",&pI1); 
		// find the actual atom numbers
		if(pA1<0) pA1*=-1;
		if(pA2<0) pA2*=-1;
		pA1/=3; // CAREFUL !! these can be negative
		pA2/=3; 
		// set the bond info for these two atoms
		// -- at this point, everything is in the same molecule, so the molecule index
		// 	will be set to -1 (to prevent confusion later).
		// 	Also, since we don't know the ordering in the file, most info will not be set
		MB[pb+P[0].NBONH].i=pI1; // holder for the index into the RK & REQ arrays
		MB[pb+P[0].NBONH].D=strdup("BONDS_WITHOUT_HYDROGEN");
		MB[pb+P[0].NBONH].s.m=-1; // we don't know the molecule number yet
		MB[pb+P[0].NBONH].s.r=-1; // we can't assume we know the residue number yet
		MB[pb+P[0].NBONH].s.a=pA1; // the best atom number we have at the moment
		MB[pb+P[0].NBONH].s.i=pA1; // for assigning molecules based on bonds, later
		MB[pb+P[0].NBONH].t.m=-1; 
		MB[pb+P[0].NBONH].t.r=-1; 
		MB[pb+P[0].NBONH].t.a=pA2; 
		MB[pb+P[0].NBONH].t.i=pA2; 
		}
}
// FORMAT(12I6)  (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
if(strcmp(P[0].S[pa].N,"ANGLES_INC_HYDROGEN")==0){
  	P[0].ITH    =pa; // ITH    : atom involved in angle "i", angle contains hydrogen
  	P[0].JTH    =pa; // JTH    : atom involved in angle "i", angle contains hydrogen
  	P[0].KTH    =pa; // KTH    : atom involved in angle "i", angle contains hydrogen
  	P[0].ICTH   =pa; // ICTH   : index into parameter arrays TK and TEQ for angle
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,(4*P[0].NTHETH));
	for(pb=0;pb<P[0].NTHETH;pb++){ // angles with hydrogen
		// read in the values from the section structure
		sscanf(P[0].S[pa].D[4*pb],"%d",&pA1);
		sscanf(P[0].S[pa].D[4*pb+1],"%d",&pA2);
		sscanf(P[0].S[pa].D[4*pb+2],"%d",&pA3);
		sscanf(P[0].S[pa].D[4*pb+3],"%d",&pI1); 
		// find the actual atom numbers
		pA1/=3; // CAREFUL !! these can be negative
		pA2/=3; 
		pA3/=3; 
		// set the bond info for these two atoms
		// -- at this point, everything is in the same molecule, so the molecule index
		// 	will be set to -1 (to prevent confusion later).
		// 	Also, since we don't know the ordering in the file, most info will not be set
		MANG[pb].i=pI1; // holder for the index into the RK & REQ arrays
		MANG[pb].D=strdup("ANGLES_INC_HYDROGEN");
		MANG[pb].a.m=-1; // we don't know the molecule number yet
		MANG[pb].a.r=-1; // we can't assume we know the residue number yet
		MANG[pb].a.a=pA1; // we can't assume the atoms are already read in
		MANG[pb].b.m=-1; 
		MANG[pb].b.r=-1;
		MANG[pb].b.a=pA2; 
		MANG[pb].c.m=-1; 
		MANG[pb].c.r=-1; 
		MANG[pb].c.a=pA3; 
		}
}
           	// ITH(i)-JTH(i)-KTH(i) 
// FORMAT(12I6)  (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
if(strcmp(P[0].S[pa].N,"ANGLES_WITHOUT_HYDROGEN")==0){
  	P[0].IT     =pa; // IT     : atom involved in angle "i", angle does not contain hydrogen
  	P[0].JT     =pa; // JT     : atom involved in angle "i", angle does not contain hydrogen
  	P[0].KT     =pa; // KT     : atom involved in angle "i", angle does not contain hydrogen
  	P[0].ICT    =pa; // ICT    : index into parameter arrays TK and TEQ for angle
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,(4*P[0].MTHETA));
	for(pb=0;pb<P[0].MTHETA;pb++){ // angles without hydrogen
		// read in the values from the section structure
		sscanf(P[0].S[pa].D[4*pb],"%d",&pA1);
		sscanf(P[0].S[pa].D[4*pb+1],"%d",&pA2);
		sscanf(P[0].S[pa].D[4*pb+2],"%d",&pA3);
		sscanf(P[0].S[pa].D[4*pb+3],"%d",&pI1); 
		// find the actual atom numbers
		pA1/=3; // CAREFUL !! these can be negative
		pA2/=3; 
		pA3/=3; 
		// set the bond info for these two atoms
		// -- at this point, everything is in the same molecule, so the molecule index
		// 	will be set to -1 (to prevent confusion later).
		// 	Also, since we don't know the ordering in the file, most info will not be set
		MANG[pb+P[0].NTHETH].i=pI1; // holder for the index into the RK & REQ arrays
		MANG[pb+P[0].NTHETH].D=strdup("ANGLES_WITHOUT_HYDROGEN");
		MANG[pb+P[0].NTHETH].a.m=-1; // we don't know the molecule number yet
		MANG[pb+P[0].NTHETH].a.r=-1; // we can't assume we know the residue number yet
		MANG[pb+P[0].NTHETH].a.a=pA1; // we can't assume the atoms are read in already
		MANG[pb+P[0].NTHETH].b.m=-1; 
		MANG[pb+P[0].NTHETH].b.r=-1; 
		MANG[pb+P[0].NTHETH].b.a=pA2; 
		MANG[pb+P[0].NTHETH].c.m=-1; 
		MANG[pb+P[0].NTHETH].c.r=-1; 
		MANG[pb+P[0].NTHETH].c.a=pA3; 
		}
}
           	// IT(i)-JT(i)-KT(i) 
// FORMAT(12I6)  (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
if(strcmp(P[0].S[pa].N,"DIHEDRALS_INC_HYDROGEN")==0){
  	P[0].IPH    =pa; // IPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	P[0].JPH    =pa; // JPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	P[0].KPH    =pa; // KPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	P[0].LPH    =pa; // LPH    : atom involved in dihedral "i", dihedral contains hydrogen
  	P[0].ICPH   =pa; // ICPH   : index into parameter arrays PK, PN, and PHASE for
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,(5*P[0].NPHIH));
	for(pb=0;pb<P[0].NPHIH;pb++){ // dihedrals with hydrogen
		// read in the values from the section structure
		sscanf(P[0].S[pa].D[5*pb],"%d",&pA1);
		sscanf(P[0].S[pa].D[5*pb+1],"%d",&pA2);
		sscanf(P[0].S[pa].D[5*pb+2],"%d",&pA3);
		sscanf(P[0].S[pa].D[5*pb+3],"%d",&pA3);
		sscanf(P[0].S[pa].D[5*pb+4],"%d",&pI1); 
		// find the actual atom numbers
		pA1/=3; // CAREFUL !! these can be negative
		pA2/=3; 
		pA3/=3; 
		pA4/=3; 
		// set the bond info for these two atoms
		// -- at this point, everything is in the same molecule, so the molecule index
		// 	will be set to -1 (to prevent confusion later).
		// 	Also, since we don't know the ordering in the file, most info will not be set
		MTOR[pb].i=pI1; // holder for the index into the RK & REQ arrays
		MTOR[pb].D=strdup("DIHEDRALS_INC_HYDROGEN");
		MTOR[pb].a.m=-1; // we don't know the molecule number yet
		MTOR[pb].a.r=-1; // we can't assume we know the residue number yet
		MTOR[pb].a.a=pA1; // we can't assume the atoms are read in already
		MTOR[pb].b.m=-1; 
		MTOR[pb].b.r=-1; 
		MTOR[pb].b.a=pA2; 
		MTOR[pb].c.m=-1; 
		MTOR[pb].c.r=-1; 
		MTOR[pb].c.a=pA3; 
		MTOR[pb].d.m=-1; 
		MTOR[pb].d.r=-1; 
		MTOR[pb].d.a=pA4; 
		}
}
           	// dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i) 
// FORMAT(12I6)  (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)

if(strcmp(P[0].S[pa].N,"DIHEDRALS_WITHOUT_HYDROGEN")==0){
  	P[0].IP     =pa; // IP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	P[0].JP     =pa; // JP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	P[0].KP     =pa; // KP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	P[0].LP     =pa; // LP     : atom involved in dihedral "i", dihedral does not contain hydrogen
  	P[0].ICP    =pa; // ICP    : index into parameter arrays PK, PN, and PHASE for
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,(5*P[0].MPHIA));
	for(pb=0;pb<P[0].MPHIA;pb++){ // dihedrals with hydrogen
		// read in the values from the section structure
		sscanf(P[0].S[pa].D[5*pb],"%d",&pA1);
		sscanf(P[0].S[pa].D[5*pb+1],"%d",&pA2);
		sscanf(P[0].S[pa].D[5*pb+2],"%d",&pA3);
		sscanf(P[0].S[pa].D[5*pb+3],"%d",&pA3);
		sscanf(P[0].S[pa].D[5*pb+4],"%d",&pI1); 
		// find the actual atom numbers
		pA1/=3; // CAREFUL !! these can be negative
		pA2/=3; 
		pA3/=3; 
		pA4/=3; 
		// set the bond info for these two atoms
		// -- at this point, everything is in the same molecule, so the molecule index
		// 	will be set to -1 (to prevent confusion later).
		// 	Also, since we don't know the ordering in the file, most info will not be set
		MTOR[pb+P[0].NPHIH].i=pI1; // holder for the index into the RK & REQ arrays
		MTOR[pb+P[0].NPHIH].D=strdup("DIHEDRALS_INC_HYDROGEN");
		MTOR[pb+P[0].NPHIH].a.m=-1; // we don't know the molecule number yet
		MTOR[pb+P[0].NPHIH].a.r=-1; // we can't assume we know the residue number yet
		MTOR[pb+P[0].NPHIH].a.a=pA1; // we can't assume the atoms are read in already
		MTOR[pb+P[0].NPHIH].b.m=-1; 
		MTOR[pb+P[0].NPHIH].b.r=-1; 
		MTOR[pb+P[0].NPHIH].b.a=pA2; 
		MTOR[pb+P[0].NPHIH].c.m=-1; 
		MTOR[pb+P[0].NPHIH].c.r=-1; 
		MTOR[pb+P[0].NPHIH].c.a=pA3; 
		MTOR[pb+P[0].NPHIH].d.m=-1; 
		MTOR[pb+P[0].NPHIH].d.r=-1; 
		MTOR[pb+P[0].NPHIH].d.a=pA4; 
		}
}

           	// dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i).  Note, if the
           	// periodicity is negative, this implies the following entry
           	// in the PK, PN, and PHASE arrays is another term in a
           	// multitermed dihedral.  
// FORMAT(12I6)  (NATEX(i), i=1,NEXT)
if(strcmp(P[0].S[pa].N,"EXCLUDED_ATOMS_LIST")==0){  // Not currently stored so no whitespace check
  	P[0].NATEX  =pa; // NATEX  : the excluded atom list.  To get the excluded list for atom
	P[0].S[pa].is_standard=0; // set as a standard section
}
           	// "i" you need to traverse the NUMEX list, adding up all
           	// the previous NUMEX values, since NUMEX(i) holds the number
           	// of excluded atoms for atom "i", not the index into the 
           	// NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
           	// excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).  
// FORMAT(5E16.8)  (ASOL(i), i=1,NPHB)
if(strcmp(P[0].S[pa].N,"HBOND_ACOEF")==0){
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,P[0].NPHB);// NOTE: Not sure if these values match - no files to compare to, MBT 2010-10-26
  	P[0].ASOL   =pa; // ASOL   : the value for the r**12 term for hydrogen bonds of all
	P[0].S[pa].is_standard=0; // set as a standard section
	if(P[0].NPHB>0){
	if(P[0].S[pa].nt>0){
	for(pb=0;pb<P[0].S[pa].nt;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].NBT[pb].LJ12_1012);}
		}
		}
}
           	// possible types.  Index into these arrays is equivalent
           	// to the CN1 and CN2 arrays, however the index is negative.
           	// For example, for atoms i and j, with i < j, the index is
           	// ICO(NTYPES*(IAC(i)-1+IAC(j)).  
// FORMAT(5E16.8)  (BSOL(i), i=1,NPHB)
if(strcmp(P[0].S[pa].N,"HBOND_BCOEF")==0){
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,P[0].NPHB);// NOTE: Not sure if these values match - no files to compare to, MBT 2010-10-26
  	P[0].BSOL   =pa; // BSOL   : the value for the r**10 term for hydrogen bonds of all
	P[0].S[pa].is_standard=0; // set as a standard section
	if(P[0].NPHB>0){
	if(P[0].S[pa].nt>0){
		for(pb=0;pb<P[0].S[pa].nt;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&A.PRM[0][0].NBT[pb].LJ10_1012);}
		}
		}
}
           	// possible types.  Indexed like ASOL.  
// FORMAT(5E16.8)  (HBCUT(i), i=1,NPHB)
if(strcmp(P[0].S[pa].N,"HBCUT")==0){ // no whitespace check used
  	P[0].HBCUT  =pa; // HBCUT  : no longer in use 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(20A4)  (ISYMBL(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"AMBER_ATOM_TYPE")==0){
  	P[0].ISYMBL =pa; // ISYMBL : the AMBER atom types for each atom 
	P[0].S[pa].is_standard=0; // set as a standard section
	// Can't assume we've already read in the type numbers...
	// read one per each atom, just like it is in the file
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	ATNAME=(char**)calloc(A.na,sizeof(char*));
	for(pb=0;pb<A.na;pb++){
		ATNAME[pb]=(char*)calloc(P[0].S[pa].nc+1,sizeof(char));
		sscanf(P[0].S[pa].D[pb],"%s",ATNAME[pb]);
		}
}
// FORMAT(20A4)  (ITREE(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"TREE_CHAIN_CLASSIFICATION")==0){
  	P[0].ITREE  =pa; // ITREE  : the list of tree joining information, classified into five
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	TREECLASS=(char**)calloc(A.na,sizeof(char*));
	for(pb=0;pb<A.na;pb++){
		TREECLASS[pb]=(char*)calloc(P[0].S[pa].nc+1,sizeof(char));
		sscanf(P[0].S[pa].D[pb],"%s",TREECLASS[pb]);
		}
}
           	// types.  M -- main chain, S -- side chain, B -- branch point, 
           	// 3 -- branch into three chains, E -- end of the chain 
// FORMAT(12I6)  (JOIN(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"JOIN_ARRAY")==0){ // no whitespace check used
  	P[0].JOIN   =pa; // JOIN   : tree joining information, potentially used in ancient
           	// analysis programs.  Currently unused in sander or gibbs.  
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  (IROTAT(i), i = 1, NATOM)
if(strcmp(P[0].S[pa].N,"IROTAT")==0){ //no whitespace check used
  	P[0].IROTAT =pa; // IROTAT : apparently the last atom that would move if atom i was
           	// rotated, however the meaning has been lost over time.
           	// Currently unused in sander or gibbs.
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  IPTRES, NSPM, NSPSOL
if(strcmp(P[0].S[pa].N,"SOLVENT_POINTERS")==0){
  	P[0].IPTRES =pa; // IPTRES : final residue that is considered part of the solute,
           	// reset in sander and gibbs
  	P[0].NSPM   =pa; // NSPM   : total number of molecules
  	P[0].NSPSOL =pa; // NSPSOL : the first solvent "molecule" 
	P[0].S[pa].is_standard=0; // set as a standard section
        ntrue=count_array_nws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N);
	if(ntrue==3){
		sscanf(P[0].S[pa].D[0],"%d",&IPTRES);
		sscanf(P[0].S[pa].D[1],"%d",&NSPM);
		sscanf(P[0].S[pa].D[2],"%d",&NSPSOL);
		}
	else{printf("WARNING: Non-standard or non-existent solvent pointer information\n\tIgnoring this section...\n");}
}
// FORMAT(12I6)  (NSP(i), i=1,NSPM)
if(strcmp(P[0].S[pa].N,"ATOMS_PER_MOLECULE")==0){ // Need to add a comparison to A.m[i][0].na (if not set, then add all the A.m[i][0].r[j].na's up) and WHITESPACE CHECK
// Note that a check for this is not currently available until AFTER the A.nm value is defined (near the end of this program)
  	P[0].NSP    =pa; // NSP    : the total number of atoms in each molecule,
           	// necessary to correctly perform the pressure scaling.  
	P[0].S[pa].is_standard=0; // set as a standard section
	NSP=(int*)calloc(P[0].S[pa].nt,sizeof(int));
	for(pb=0;pb<P[0].S[pa].nt;pb++){sscanf(P[0].S[pa].D[pb],"%d",&NSP[pb]);}
}
// FORMAT(5E16.8)  BETA, BOX(1), BOX(2), BOX(3)
if(strcmp(P[0].S[pa].N,"BOX_DIMENSIONS")==0){
  	P[0].BETA   =pa; // BETA   : periodic box, angle between the XY and YZ planes in degrees.
  	P[0].BOX    =pa; // BOX    : the periodic box lengths in the X, Y, and Z directions 
	P[0].S[pa].is_standard=0; // set as a standard section
	ntrue=count_array_nws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N); // Non-comparative whitespace check
// Checks for box definitions
	if(P[0].IFBOX==1 && ntrue==4){ // Cubic box definition from the header
                sscanf(P[0].S[pa].D[0],"%lf",&A.BOX[0].C[0].D[0]); ///< Record value of BETA
                sscanf(P[0].S[pa].D[1],"%lf",&A.BOX[0].C[1].D[0]); ///< Record BOX X
                sscanf(P[0].S[pa].D[2],"%lf",&A.BOX[0].C[1].D[1]); ///< Record BOX Y
                sscanf(P[0].S[pa].D[3],"%lf",&A.BOX[0].C[1].D[2]); ///< Record BOX Z
		}
	else if(P[0].IFBOX==2){ //Octahedral box definition from the header ADD an ntrue comparison for octahedral box definitions
		fprintf(stderr,"\nWARNING: Octahedral box definitions are not currently stored in the assembly.\n\tPlease include a parsing function. Ignoring any contents in this section.\nALERT: Any call to A.BOX or A.nBOX variables may cause memory issues without a segmentation fault.\n");
		}
	else{fprintf(stderr,"\nIn parse_amber_prmtop.c : unexpected number of entries in section BOX_DIMENSIONS.  Ignoring contents.\nALERT: Any call to A.BOX or A.nBOX variables may cause memory issues without a segmentation fault.\n\n");}
}
// The following are only present if IFCAP .gt. 0 
// FORMAT(12I6)  NATCAP
if(strcmp(P[0].S[pa].N,"CAP_INFO")==0){ // no whitespace check
  	P[0].NATCAP =pa; // NATCAP : last atom before the start of the cap of waters placed by edit 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(5E16.8)  CUTCAP, XCAP, YCAP, ZCAP
if(strcmp(P[0].S[pa].N,"CAP_INFO2")==0){ // not putting this in quite yet, no whitespace check
  	P[0].CUTCAP =pa; // CUTCAP : the distance from the center of the cap to the outside
  	P[0].XCAP   =pa; // XCAP   : X coordinate for the center of the cap
  	P[0].YCAP   =pa; // YCAP   : Y coordinate for the center of the cap
  	P[0].ZCAP   =pa; // ZCAP   : Z coordinate for the center of the cap 
	P[0].S[pa].is_standard=0; // set as a standard section
}
if(strcmp(P[0].S[pa].N,"RADIUS_SET")==0){
	P[0].S[pa].is_standard=0; // set as a standard section
	ntrue=count_array_nws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N); // Non-comparative whitespace check
	if(ntrue>0){RADTYPE=strdup(P[0].S[pa].D[0]);}
}
if(strcmp(P[0].S[pa].N,"RADII")==0){
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	R=(double*)calloc(A.na,sizeof(double));
	for(pb=0;pb<A.na;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&R[pb]);}
}
if(strcmp(P[0].S[pa].N,"SCREEN")==0){
	P[0].S[pa].is_standard=0; // set as a standard section
        compare_array_ws(P[0].S[pa].nt,P[0].S[pa].D,P[0].S[pa].N,A.na);
	SC=(double*)calloc(A.na,sizeof(double));
	for(pb=0;pb<A.na;pb++){sscanf(P[0].S[pa].D[pb],"%lf",&SC[pb]);}
}
if(strcmp(P[0].S[pa].N,"LES_NTYP")==0){ // leave these out for now, no whitespace check
	P[0].S[pa].is_standard=0; // set as a standard section
}
if(strcmp(P[0].S[pa].N,"LES_TYPE")==0){ // no whitespace check
	P[0].S[pa].is_standard=0; // set as a standard section
}
if(strcmp(P[0].S[pa].N,"LES_FAC")==0){ // no whitespace check
	P[0].S[pa].is_standard=0; // set as a standard section
}
if(strcmp(P[0].S[pa].N,"LES_CNUM")==0){ // no whitespace check
	P[0].S[pa].is_standard=0; // set as a standard section
}
if(strcmp(P[0].S[pa].N,"LES_ID")==0){ // no whitespace check
	P[0].S[pa].is_standard=0; // set as a standard section
}
// The following is only present if IFPERT .gt. 0
/* Note that the initial state, or equivalently the prep/link/edit state, 
	is represented by lambda=1 and the perturbed state, or final state 
	specified in parm, is the lambda=0 state. */ 
// FORMAT(12I6)  (IBPER(i), JBPER(i), i=1,NBPER)
if(strcmp(P[0].S[pa].N,"PERT_BOND_ATOMS")==0){ // no whitespace check
// The following are only present if IFBOX .gt. 0 
  	P[0].IBPER  =pa; // IBPER  : atoms involved in perturbed bonds
  	P[0].JBPER  =pa; // JBPER  : atoms involved in perturbed bonds 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  (ICBPER(i), i=1,2*NBPER)
if(strcmp(P[0].S[pa].N,"PERT_BOND_PARAMS")==0){ // no whitespace check
  	P[0].ICBPER =pa; // ICBPER : pointer into the bond parameter arrays RK and REQ for the
	P[0].S[pa].is_standard=0; // set as a standard section
}
           	// perturbed bonds.  ICBPER(i) represents lambda=1 and 
           	// ICBPER(i+NBPER) represents lambda=0.  
// FORMAT(12I6)  (ITPER(i), JTPER(i), KTPER(i), i=1,NGPER)
if(strcmp(P[0].S[pa].N,"PERT_ANGLE_ATOMS")==0){ // no whitespace check
  	P[0].IPTER  =pa; // IPTER  : atoms involved in perturbed angles
  	P[0].JTPER  =pa; // JTPER  : atoms involved in perturbed angles
  	P[0].KTPER  =pa; // KTPER  : atoms involved in perturbed angles 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  (ICTPER(i), i=1,2*NGPER)
if(strcmp(P[0].S[pa].N,"PERT_ANGLE_PARAMS")==0){ // no whitespace check
  	P[0].ICTPER =pa; // ICTPER : pointer into the angle parameter arrays TK and TEQ for 
	P[0].S[pa].is_standard=0; // set as a standard section
}
           	// the perturbed angles.  ICTPER(i) represents lambda=0 and 
           	// ICTPER(i+NGPER) represents lambda=1.  
// FORMAT(12I6)  (IPPER(i), JPPER(i), KPPER(i), LPPER(i), i=1,NDPER)
if(strcmp(P[0].S[pa].N,"PERT_DIHEDRAL_ATOMS")==0){ // no whitespace checkf
  	P[0].IPPER  =pa; // IPPER  : atoms involved in perturbed dihedrals
  	P[0].JPPER  =pa; // JPPER  : atoms involved in perturbed dihedrals
  	P[0].KPPER  =pa; // KPPER  : atoms involved in perturbed dihedrals
  	P[0].LPPER  =pa; // LPPER  : atoms involved in pertrubed dihedrals 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  (ICPPER(i), i=1,2*NDPER)
if(strcmp(P[0].S[pa].N,"PERT_DIHEDRAL_PARAMS")==0){ // no whitespace check
  	P[0].ICPPER =pa; // ICPPER : pointer into the dihedral parameter arrays PK, PN and
	P[0].S[pa].is_standard=0; // set as a standard section
}
           	// PHASE for the perturbed dihedrals.  ICPPER(i) represents 
           	// lambda=1 and ICPPER(i+NGPER) represents lambda=0.  
// FORMAT(20A4)  (LABRES(i), i=1,NRES)
if(strcmp(P[0].S[pa].N,"PERT_RESIDUE_NAME")==0){ // no whitespace check
  	P[0].LABRES =pa; // LABRES : residue names at lambda=0 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(20A4)  (IGRPER(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"PERT_ATOM_NAME")==0){ // no whitespace check
  	P[0].IGRPER =pa; // IGRPER : atomic names at lambda=0 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(20A4)  (ISMPER(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"PERT_ATOM_SYMBOL")==0){ // no whitespace check
  	P[0].ISMPER =pa; // ISMPER : atomic symbols at lambda=0 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(5E16.8)  (ALMPER(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"ALMPER")==0){ // no whitespace check
  	P[0].ALMPER =pa; // ALMPER : unused currently in gibbs 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  (IAPER(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"IAPER")==0){ // no whitespace check
  	P[0].IAPER  =pa; // IAPER  : IAPER(i) = 1 if the atom is being perturbed 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// FORMAT(12I6)  (IACPER(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"PERT_ATOM_TYPE_INDEX")==0){ // no whitespace check
  	P[0].IACPER =pa; // IACPER : index for the atom types involved in Lennard Jones
	P[0].S[pa].is_standard=0; // set as a standard section
}
           	// interactions at lambda=0.  Similar to IAC above.  See ICO above.  
// FORMAT(5E16.8)  (CGPER(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"PERT_CHARGE")==0){ // no whitespace check
  	P[0].CGPER  =pa; // CGPER  : atomic charges at lambda=0 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// The following is only present if IPOL .eq. 1 
// FORMAT(5E18.8) (ATPOL(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"POLARIZABILITY")==0){ // leave out for now, // no whitespace check
  	P[0].ATPOL  =pa; // ATPOL  : atomic polarizabilities 
	P[0].S[pa].is_standard=0; // set as a standard section
}
// The following is only present if IPOL .eq. 1 .and. IFPERT .eq. 1 
// FORMAT(5E18.8) (ATPOL1(i), i=1,NATOM)
if(strcmp(P[0].S[pa].N,"PERT_POLARIZABILITY")==0){ // no whitespace check
  	P[0].ATPOL1 =pa; // ATPOL1 : atomic polarizabilities at lambda = 1 (above is at lambda = 0) 
	P[0].S[pa].is_standard=0; // set as a standard section
}
	} // close loop through each section in the prmtop structure

/*  The following code prints out the prmtop file exactly as read in from the unparsed read */
/*
fileset F;
F.N=strdup("test_rewrite_of_prmtop");
F.F=myfopen(F.N,"w");
fprintf(F.F,"%s",P[0].VERSION);
for(pa=0;pa<P[0].nS;pa++){ // for each section found
	fprintf(F.F,"%%FLAG %s\n",P[0].S[pa].N);
	fprintf(F.F,"%%FORMAT(%s)\n",P[0].S[pa].FORMAT);
	for(pb=0;pb<P[0].S[pa].nt;pb++){
		fprintf(F.F,"%s",P[0].S[pa].D[pb]);
		//fprintf(F.F,"D[%d] is >>>%s<<<\n",pb,P[0].S[pa].D[pb]);
		fflush(F.F);
//printf("\npb is %d ; P[0].S[%d].npl is %d ; (pb+1)%%npl is %d\n",pb,pa,P[0].S[pa].npl,);
		if((((pb+1)%P[0].S[pa].npl))==0){fprintf(F.F,"\n");} 
		//if((((pb+1)*P[0].S[pa].nc))%80==0){printf("\n");} 
		}
	if(P[0].S[pa].nt==0){fprintf(F.F,"\n");}
	if(((pb)*P[0].S[pa].nc)%80!=0){fprintf(F.F,"\n");} 
	}
fclose(F.F);
*/

// Set the names of the atom types in the type array and in the atom struct
for(pa=0;pa<A.na;pa++){ 
	A.a[pa][0].T=strdup(ATNAME[pa]); 
	A.PRM[0][0].AT[A.a[pa][0].t].N=strdup(ATNAME[pa]); 
	}
// Now set the names in the non-bonded bond type array
// -- set atom information in the LJ type arrays
//int *ICO,nICO;
if(nICO!=(A.PRM[0][0].nAT*A.PRM[0][0].nAT)){mywhine("nICO!=(A.PRM[0][0].nAT*A.PRM[0][0].nAT) in parse_amber_prmtop");}
// For example, for atoms i and j, with i < j, the index is
// ICO(NTYPES*(IAC(i)-1)+IAC(j)).  
pI1=P[0].NTYPES*(P[0].NTYPES+1)/2;
for(pa=0;pa<A.PRM[0][0].nAT;pa++){ // this is index 'j' in the prmtop file
for(pb=0;pb<(pa+1);pb++){ // this is index 'i' in the prmtop file
	// allocate space for the names of the two atom types
	pc=P[0].NTYPES*pb+pa;
	if(pc>nICO){mywhine("pc>nICO in parse_amber_prmtop");}
	if(ICO[pc]>pI1){mywhine("ICO[pc]>(P[0].NTYPES*(P[0].NTYPES+1)/2) in parse_amber_prmtop");}
	A.PRM[0][0].NBT[ICO[pc]].NT=(char**)calloc(2,sizeof(char*));
	A.PRM[0][0].NBT[ICO[pc]].NT[0]=strdup(A.PRM[0][0].AT[pb].N);
	A.PRM[0][0].NBT[ICO[pc]].NT[1]=strdup(A.PRM[0][0].AT[pa].N);
	}
}

// set molecule information
// 	First, copy MOLI & MB arrays for safety & record
MOLBNDI=(molindex*)calloc(P[0].NATOM,sizeof(molindex));
MBTMP=(molbond*)calloc(A.nb,sizeof(molbond));
for(pa=0;pa<P[0].NATOM;pa++){ MOLBNDI[pa]=MOLI[pa]; }
for(pa=0;pa<A.nb;pa++){ 
	MBTMP[pa]=MB[pa]; 
	}
// call the "find molecules" function
find_molecules_molbond_array(A.nb, MBTMP, P[0].NATOM, MOLBNDI);

// Rearrange the structures to reflect the new molecule information
// 1. Find numbers of molecules and make space
//	The molecules, coming from find_molecules_molbond_array, should count 
//	sequentially starting with zero.
// 2. While we're at it, associate each residue with a molecule in a convenient array (resi)
nummol=0;
resi=(int*)calloc(P[0].NRES,sizeof(int));
for(pa=0;pa<P[0].NRES;pa++){ resi[pa]=-1; }
for(pa=0;pa<P[0].NATOM;pa++){ 
	if(MOLBNDI[pa].m>nummol){ nummol = MOLBNDI[pa].m; }  // find numbers of molecules
	if(resi[MOLBNDI[pa].r]==-1) resi[MOLBNDI[pa].r]=MOLBNDI[pa].m;
	else if(resi[MOLBNDI[pa].r]!=MOLBNDI[pa].m){
		if(resi[MOLBNDI[pa].r]!=-1) mywhine("resi[MOLBNDI[pa].r]!=MOLBNDI[pa].m in parse_amber_prmtop");} 
	}
for(pa=0;pa<P[0].NATOM;pa++){ 
	if(resi[MOLBNDI[pa].r]==-1) { // this one isn't assigned a molecule yet
		nummol++;
		resi[MOLBNDI[pa].r]=MOLBNDI[pa].m=nummol; // Assign a new molecule number
		MOLBNDI[pa].r=0; // since there are no bonds, this is a single atom, so only one residue
		}
	}
nummol++; // because the count starts with zero, but the total number doesn't
A.nm=nummol; // make space for the molecules
A.m=(molecule**)calloc(A.nm,sizeof(molecule*));
for(pa=0;pa<A.nm;pa++){ 
	A.m[pa]=(molecule*)calloc(1,sizeof(molecule)); 
	A.m[pa][0].na=0;
	A.m[pa][0].nr=0;
	A.m[pa][0].a=(atom**)calloc(1,sizeof(atom*));
	A.m[pa][0].r=(residue*)calloc(1,sizeof(residue));
	}
// place residues into molecule structures -- and do a little housekeeping
for(pa=0;pa<P[0].NRES;pa++){
	A.r[pa][0].n=pa+1; // the file-print residue "number"
	A.m[resi[pa]][0].nr++;
	A.m[resi[pa]][0].r=(residue*)realloc(A.m[resi[pa]][0].r,A.m[resi[pa]][0].nr*sizeof(residue));
	A.m[resi[pa]][0].r[A.m[resi[pa]][0].nr-1]=A.r[pa][0]; // copy over residue 
	A.r[pa]=&A.m[resi[pa]][0].r[A.m[resi[pa]][0].nr-1];// reassign pointer
	for(pb=0;pb<A.r[pa][0].na;pb++){
		MOLBNDI[A.r[pa][0].a[pb].n-1].r=A.m[resi[pa]][0].nr-1;
		}
	}


// Use the following to print info for all molecules to screen (for debugging)
//for(pa=0;pa<A.nm;pa++){
//printf("============== there are %d molecules =================\n",A.nm);
//dprint_molecule(A.m[pa],1000);}

// if the NSPM & NSP information is present, check results for consistency
// NSPM = number of molecules
// *NSP = number of atoms in each molecule
// IPTRES = last residue that is solute
// NSPSOL = first solvent residue
//


// Use new molecule information to reset residue numbers and molecule affiliations
//	set the MOLI indices, too...
for(pa=0;pa<A.nm;pa++){
	// set ensemble indices
	A.m[pa][0].nensi=1;
	A.m[pa][0].ensi=(ensindex*)calloc(1,sizeof(ensindex));
	A.m[pa][0].ensi[0].E=A.m[pa][0].ensi[0].r=A.m[pa][0].ensi[0].a=-1;
	A.m[pa][0].ensi[0].A=0;
	A.m[pa][0].ensi[0].m=pa;
	for(pb=0;pb<A.m[pa][0].nr;pb++){
		// set ensemble indices
		A.m[pa][0].r[pb].nensi=1;
		A.m[pa][0].r[pb].ensi=(ensindex*)calloc(1,sizeof(ensindex));
		A.m[pa][0].r[pb].ensi[0].E=A.m[pa][0].r[pb].ensi[0].a=-1;
		A.m[pa][0].r[pb].ensi[0].A=0;
		A.m[pa][0].r[pb].ensi[0].m=pa;
		A.m[pa][0].r[pb].ensi[0].r=pb;
		// set molecule indices
		A.m[pa][0].r[pb].moli.m=pa;
		A.m[pa][0].r[pb].moli.r=pb;
		A.m[pa][0].r[pb].moli.a=-1;
		for(pc=0;pc<A.m[pa][0].r[pb].na;pc++){ 
			// set ensemble indices
			A.m[pa][0].r[pb].a[pc].nensi=1;
			A.m[pa][0].r[pb].a[pc].ensi=(ensindex*)calloc(1,sizeof(ensindex));
			A.m[pa][0].r[pb].a[pc].ensi[0].E=-1;
			A.m[pa][0].r[pb].a[pc].ensi[0].A=0;
			A.m[pa][0].r[pb].a[pc].ensi[0].m=pa;
			A.m[pa][0].r[pb].a[pc].ensi[0].r=pb;
			A.m[pa][0].r[pb].a[pc].ensi[0].a=pc;
			// set molecule indices
			A.m[pa][0].r[pb].a[pc].moli.m=pa;
			A.m[pa][0].r[pb].a[pc].moli.r=pb;
			A.m[pa][0].r[pb].a[pc].moli.a=pc;
			}
		} 
	}

// printf("============== there are %d molecules =================\n",A.nm);

// Set the bond information at the atom level
for(pa=0;pa<A.nm;pa++){
A.m[pa][0].N=strdup("unknown");
A.m[pa][0].D=strdup("Molecule determined by read of Amber prmtop file");
//printf("\tMolecule %d contains %d residues\n",pa,A.m[pa][0].nr);
for(pb=0;pb<A.m[pa][0].nr;pb++){
//printf("\t\tResidue %d contains %d atoms\n",pb,A.m[pa][0].r[pb].na);
for(pc=0;pc<A.m[pa][0].r[pb].na;pc++){
        A.m[pa][0].r[pb].a[pc].mb=(molbond*)calloc(1,sizeof(molbond)); 
//printf("Atom %s (# %d) belongs to residue # %d and molecule # %d\n",A.m[pa][0].r[pb].a[pc].N,pc,pb,pa);
	}
	}
	}

for(pa=0;pa<A.nb;pa++){
	// set global bonds
	A.b[pa].i=MBTMP[pa].i;
	A.b[pa].D=strdup(MBTMP[pa].D);
	A.b[pa].s.m=sm=MOLBNDI[MBTMP[pa].s.i].m;
	A.b[pa].s.r=sr=MOLBNDI[MBTMP[pa].s.i].r;
	A.b[pa].s.a=sa=MOLBNDI[MBTMP[pa].s.i].a;
	A.b[pa].t.m=tm=MOLBNDI[MBTMP[pa].t.i].m;
	A.b[pa].t.r=tr=MOLBNDI[MBTMP[pa].t.i].r;
	A.b[pa].t.a=ta=MOLBNDI[MBTMP[pa].t.i].a;

/* If you're trying to fix the code, the following might be instructive. */
/*
printf("bond number, pa=%d \n",pa);
printf("  MBTMP[pa].s.i is %d \n",MBTMP[pa].s.i);
printf("    MBTMP[pa].t.i is %d\n",MBTMP[pa].t.i);
printf("\tMOLBNDI[(s)].m =%d \n",MOLBNDI[MBTMP[pa].s.i].m);
printf("\t\t .r=%d \n",MOLBNDI[MBTMP[pa].s.i].r);
printf("\t\t    .a=%d\n",MOLBNDI[MBTMP[pa].s.i].a);
printf("\t\tMOLBNDI[(t)].m =%d \n",MOLBNDI[MBTMP[pa].t.i].m);
printf("\t\t\t  .r=%d \n",MOLBNDI[MBTMP[pa].t.i].r);
printf("\t\t\t    .a=%d\n",MOLBNDI[MBTMP[pa].t.i].a);
*/
	// -- set local bonds to molbond structure (later to local bond structure)
	A.m[sm][0].r[sr].a[sa].nmb++;
	A.m[tm][0].r[tr].a[ta].nmb++;
	A.m[sm][0].r[sr].a[sa].mb=(molbond*)realloc(A.m[sm][0].r[sr].a[sa].mb,A.m[sm][0].r[sr].a[sa].nmb*sizeof(molbond));
	A.m[tm][0].r[tr].a[ta].mb=(molbond*)realloc(A.m[tm][0].r[tr].a[ta].mb,A.m[tm][0].r[tr].a[ta].nmb*sizeof(molbond));
	// set both atoms as being bonded to the other
	A.m[sm][0].r[sr].a[sa].mb[A.m[sm][0].r[sr].a[sa].nmb-1]=A.b[pa]; // if already "source", all is ok
	// if the atom is a "target", we have to turn things around a bit
	A.m[tm][0].r[tr].a[ta].mb[A.m[tm][0].r[tr].a[ta].nmb-1].s.a = A.b[pa].t.a;
	A.m[tm][0].r[tr].a[ta].mb[A.m[tm][0].r[tr].a[ta].nmb-1].s.r = A.b[pa].t.r;
	A.m[tm][0].r[tr].a[ta].mb[A.m[tm][0].r[tr].a[ta].nmb-1].s.m = A.b[pa].t.m;
	A.m[tm][0].r[tr].a[ta].mb[A.m[tm][0].r[tr].a[ta].nmb-1].t.a = A.b[pa].s.a ;
	A.m[tm][0].r[tr].a[ta].mb[A.m[tm][0].r[tr].a[ta].nmb-1].t.r = A.b[pa].s.r ;
	A.m[tm][0].r[tr].a[ta].mb[A.m[tm][0].r[tr].a[ta].nmb-1].t.m = A.b[pa].s.m ;
	}

/* uncomment this to check the function's work
for(pa=0;pa<A.nm;pa++){
printf("MOLECULE #%d\n",pa);
for(pb=0;pb<A.m[pa][0].nr;pb++){
printf("\tRESIDUE #%d (name %s)\n",pb,A.m[pa][0].r[pb].N); 
for(pc=0;pc<A.m[pa][0].r[pb].na;pc++){
printf("\t\tATOM #%d (name %s) has %d bonds\n",pc,A.m[pa][0].r[pb].a[pc].N,A.m[pa][0].r[pb].a[pc].nmb); 
for(pd=0;pd<A.m[pa][0].r[pb].a[pc].nmb;pd++){
tm=A.m[pa][0].r[pb].a[pc].mb[pd].t.m;
tr=A.m[pa][0].r[pb].a[pc].mb[pd].t.r;
ta=A.m[pa][0].r[pb].a[pc].mb[pd].t.a;
printf("\t\t\tTo %s (atom number %d)\n",A.m[tm][0].r[tr].a[ta].N,A.m[tm][0].r[tr].a[ta].n);
}
}
}
}
*/


// One day when we know how the connection tree will be structured:
// -- set local torsions and angles
// -- set connection tree after all that...
// 	-- check this tree against the amber tree info for sanity

// free temporary spaces
if(MOLI!=NULL){free(MOLI);}
if(ICO!=NULL){free(ICO);}
if(MB!=NULL){free(MB);}
if(MANG!=NULL){free(MANG);}
if(MTOR!=NULL){free(MTOR);}
if(MASS!=NULL){free(MASS);}
if(NSP!=NULL){free(NSP);}
if(R!=NULL){free(R);}
if(SC!=NULL){free(SC);}
if(MOLBNDI!=NULL){free(MOLBNDI);}
if(MBTMP!=NULL){free(MBTMP);}
if(resi!=NULL){free(resi);}
for(pa=0;pa<A.na;pa++){
	if(ATNAME[pa]!=NULL){free(ATNAME[pa]);}
	if(TREECLASS[pa]!=NULL){free(TREECLASS[pa]);}
}
if(ATNAME!=NULL){free(ATNAME);}
if(TREECLASS!=NULL){free(TREECLASS);}

return A;
}


