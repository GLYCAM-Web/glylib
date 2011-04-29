/* The following is a set of functions intended primarily to facilitate 
debugging.  But, they could easily be used for verbose data output within
the program.  They all follow the same basic format:

print_X(X x, int level)

where 	X = one of the structures defined above (e.g., molecule, residue, 
		bondset, atom, bond, atype, plane, vectormag_3D, coord_3D).
	*x = pointer to your structure
	level = the depth to print.  This is only available on structures 
		likely to contain arrays of other structures.  level=0 only
		prints the top-level information.  level=1 prints one level 
		down of other structures, lists, etc. level=2 prints two 
		levels down, etc.   The value of level can be larger than 
		the depth available -- the function will stop when it runs
		out of structures to print.  
  Author: Lachele Foley
*/
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
void dXprint_coord_3D(coord_3D *c){
	printf("coord i j k \t%20.12f \t%20.12f \t%20.12f\n",c[0].i,c[0].j,c[0].k);
return;
}

void dXprint_vectormag_3D(vectormag_3D *v){
	printf("vectormag i j k d \t%20.12f \t%20.12f \t%20.12f \t%20.12f\n",v[0].i,v[0].j,v[0].k,v[0].d); 
return;
}

void dXprint_plane(plane *p){
	printf("plane A B C D \t%20.12f \t%20.12f \t%20.12f \t%20.12f\n",p[0].A,p[0].B,p[0].C,p[0].D); 
return;
}

void dXprint_atype(atype *t,int i){
	int pa=0;
	printf("atype \tn=%d ; N=>>%s<< ; NT=>>%s<< ; m=%20.12f ; nb=%d\n",t[0].n,t[0].N,t[0].NT,t[0].m,t[0].nb);
	printf("\tdesc=>>%s<<\n",t[0].desc);
	if(i>0){ 
		for(pa=0;((pa<i)&&(pa<t[0].nb));pa++){
			printf("\t%20.12f",t[0].bo[pa]);
			}
		printf("\n");
		}
return;
}

void dXprint_molbond(molbond *mb){
	printf("molbond s (m,r,a=%d,%d,%d) t (m,r,a=%d,%d,%d) o %20.12f\n",mb[0].s.m,mb[0].s.r,mb[0].s.a,mb[0].t.m,mb[0].t.r,mb[0].t.a,mb[0].o);
return;
}
void dXprint_bond(bond *b){
	printf("bond s (m,r,a=%d,%d,%d) t (m,r,a=%d,%d,%d) o %20.12f\n",b[0].s.m,b[0].s.r,b[0].s.a,b[0].t.m,b[0].t.r,b[0].t.a,b[0].o);
return;
}



void dXprint_atom(atom *a,int i){
	int pa=0;
	printf("atom \tn=%d ; N=>>%s<< ; t=%d ; nb=%d ; nalt=%d ; nvec=%d\n",a[0].n,a[0].N,\
		a[0].t,a[0].nb,a[0].nalt,a[0].nvec);
	printf("\tcoords x y z \t%20.12f \t%20.12f \t%20.12f\n",a[0].x.i,a[0].x.j,a[0].x.k); 
	if(i>0){
		for(pa=0;((pa<i)&&(pa<a[0].nmb));pa++){
			printf("\tmolbond set %d:  ",pa);
			dXprint_molbond(&(a[0].mb[pa])); 
			}
		for(pa=0;((pa<i)&&(pa<a[0].nalt));pa++){
			printf("\talt coord set %d:  ",pa);
			dXprint_coord_3D(&(a[0].xa[pa])); 
			}
		for(pa=0;((pa<i)&&(pa<a[0].nvec));pa++){
			printf("\tset %d:  ",pa);
			dXprint_vectormag_3D(&(a[0].v[pa])); 
			}
		} 
return;
}

void dXprint_molbondset(molbondset *bs,int i){
	printf("NEED TO WRITE dXprint_molbondset\n");
return;
}
void dXprint_bondset(bondset *bs,int i){
	int pa=0;
	printf("bondset contains n=%d bonds\n",bs[0].n);
	if(i>0){
		for(pa=0;((pa<i)&&(pa<bs[0].n));pa++){
			printf("\tconsec bonds set %d : \t",pa);
			dXprint_bond(&(bs[0].b[pa]));
			}
		}
return;
}

void dXprint_residue(residue *r,int i){
	int pa=0;
	if(r[0].IC==NULL){printf("residue \tn=%d (no IC);",r[0].n);}
	else{printf("residue \tn=%d IC=%s;",r[0].n,r[0].IC);}
	printf(" N=>>%s<< ; na=%d ; nbs=%d ; nring=%d \n",r[0].N,\
		r[0].na,r[0].nbs,r[0].nring); 
	if(i>0){
		for(pa=0;((pa<i)&&(pa<r[0].na));pa++){
			printf("\tatom set %d\n",pa);
			dXprint_atom(&(r[0].a[pa]),i); 
			}
		for(pa=0;((pa<i)&&(pa<r[0].nbs));pa++){
			printf("\tconsec. lin, bondset set %d\n",pa);
			dXprint_molbondset(&(r[0].bs[pa]),i); 
			}
		for(pa=0;((pa<i)&&(pa<r[0].nrc));pa++){
			printf("\tconsec ring coords set %d:  ",pa);
			dXprint_coord_3D(&(r[0].rc[pa])); 
			}
		for(pa=0;((pa<i)&&(pa<r[0].nrp));pa++){
			printf("\tring plane set %d:  ",pa);
			dXprint_plane(&(r[0].rp[pa])); 
			}
		} 
return;
}

void dXprint_molecule(molecule *m,int i){
	int pa=0;
	printf("molecule \ti=%d ; N=>>%s<< ; nr=%d \n",m[0].i,m[0].N,m[0].nr); 
	if(i>0){
		for(pa=0;((pa<i)&&(pa<m[0].nr));pa++){
			printf("\t\tset %d\n",pa);
			dXprint_residue(&(m[0].r[pa]),i); 
			}
		}
return;
}

