/*
 * Oliver Grant 14Mar11
 * Finds clashes between atoms in assemblies in a pairwize fashion. 
 * NOTE: xsA and xsB refer to the co-ord sets to use. This functionality is not implemented yet. so they are dead variables. It could be preferable to create a driver script to call this one for every frame.
*/
#include <glylib.h> // this is clearly stupid and is giving me warnings.
// Bondi definitions of VDW in Angstrom
//#define H 	1.09 // Rowland and Taylor modification
//#define C	1.70
//#define N 	1.55
//#define O	1.52
//#define F 	1.47
//#define P	1.80
//#define S	1.80
//#define Cl	1.75

void find_vdw_clashes_pairwize_between_Assemblies(assembly *A, int xsA, assembly *B, int xsB){ // xs is co-ord set to use. -1 for main co-ord set
int mi=0,mii=0,ri=0,rii=0,ai=0,aii=0;//loop counters for A(i) and B(ii) assemblies
double soa, total_soa=0.0,resid_soa,prev_soa; // Sphere Overlap Area
double rA,rB; // radii
double x,y,z,d; // xyz are holders to make code more readable. d is distance between two atoms
FILE * pFile; // Clash results out file
char filename[25]; // Name of out file
//compare dist between each atom in A assembly and all atoms in B assembly

printf("Entered the FUNction\n");
for(mi=0;mi<(*A).nm;mi++){
        for(ri=0;ri<(*A).m[mi][0].nr;ri++){
                for(ai=0;ai<(*A).m[mi][0].r[ri].na;ai++){// ai is now current atom to check against every atom in assembly B
			for(mii=0;mii<(*B).nm;mii++){
				for(rii=0;rii<(*B).m[mii][0].nr;rii++){
					for(aii=0;aii<(*B).m[mii][0].r[rii].na;aii++){//for current atom check dist to every other atom within the assembly B
						//printf("Comparing Values\n rA=%f\n",rA);
                                		x=((*A).m[mi][0].r[ri].a[ai].x.i-(*B).m[mii][0].r[rii].a[aii].x.i);
                                		y=((*A).m[mi][0].r[ri].a[ai].x.j-(*B).m[mii][0].r[rii].a[aii].x.j);
                                		z=((*A).m[mi][0].r[ri].a[ai].x.k-(*B).m[mii][0].r[rii].a[aii].x.k);
                                		d=sqrt((x*x)+(y*y)+(z*z)); //pythagorus square on hyp in 3D
                                		//printf("\ndist is %f\n",d);
						// set radii of the atoms
						// radii are tuned to get 0 clash for natural substrate
						if ((*A).m[mi][0].r[ri].a[ai].N[0]=='C'){rA=1.70;} 
						if ((*A).m[mi][0].r[ri].a[ai].N[0]=='O'){rA=1.52;} 
						if ((*A).m[mi][0].r[ri].a[ai].N[0]=='N'){rA=1.55;} 
						if ((*A).m[mi][0].r[ri].a[ai].N[0]=='S'){rA=1.80;} 

						if ((*B).m[mii][0].r[rii].a[aii].N[0]=='C'){rB=1.70;}
                                                if ((*B).m[mii][0].r[rii].a[aii].N[0]=='O'){rB=1.52;}
						if ((*B).m[mii][0].r[rii].a[aii].N[0]=='N'){rB=1.55;}
						if ((*B).m[mii][0].r[rii].a[aii].N[0]=='S'){rB=1.80;}
						
						//printf("Dist is %f, rA is %f, rB is %f \n\n",d,rA,rB);
						// if the sum of the radii is greater than the distance between them.
						if (rA + rB > d + 0.4){ // 0.6 overlap is deemed acceptable. (Copying chimera:)
  							soa=2*PI*rA*(rA-d/2-(((rA*rA)-(rB*rB))/(2*d))); 
								// Eqn 1, Rychkov and Petukhov, J. Comput. Chem., 2006, Joint Neighbours...
							total_soa=total_soa + soa;
							strcpy(filename, "Clash_results_atoms.txt");
							pFile = fopen (filename,"a+"); // open file and append
							//if (pFile=NULL){ // this is always true...
							//	printf("File open FAIL!\n");
							//	break ;
							fprintf(pFile,"Clash Area For Resid No.%d Atom_Name:%c%c is %5.1f \n",(*A).m[mi][0].r[ri].n,(*A).m[mi][0].r[ri].a[ai].N[0],(*A).m[mi][0].r[ri].a[ai].N[1],soa);
							fprintf(pFile," ...Clashing with Resid No.%d Atom_Name:%c%c in the receptor.\n",(*B).m[mii][0].r[rii].n,(*B).m[mii][0].r[rii].a[aii].N[0],(*B).m[mii][0].r[rii].a[aii].N[1]);
							fclose (pFile);
						}
					}
				}
                        }
		}
		resid_soa=(total_soa - prev_soa);
		strcpy(filename, "Clash_results_resid.txt");
		pFile = fopen (filename,"a+");
		fprintf(pFile,"Clash Area For Resid No.%d is %5.1f\n",(*A).m[mi][0].r[ri].n,resid_soa);
		fclose (pFile);
		prev_soa=total_soa;
        }
}
total_soa=(total_soa+0.5); // so can round off by truncation
int total_soa_int;
total_soa_int=total_soa;

strcpy(filename, "Clash_results_total.txt");
pFile = fopen (filename,"a+");
fprintf(pFile,"\nTotal_Overlap=%d",total_soa_int);
fclose (pFile);
return;
}
