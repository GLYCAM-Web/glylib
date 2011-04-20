/* \file set_COM.c 
 * \addtogroup GEOMETRY
 * \brief   Sets COM's and m's (masses) in the standard structures.
 *
 * Replaces assign_COM functions.
 * If ATYPE==NULL or ==0x0, masses will be taken from entries
 * in the atom strucutures.
 *
 * These functions behave differently from the get_COM functions in
 * that they rely, recursively, on nested structures.  In other words,
 * set_molecule_COM, will set each residue's COM, but then will calculate
 * its own COM based on the residues' COM's rather than directly on
 * information in the atoms.  The two approaches, in most circumstances,
 * will give the same result.
 *
 * Begun on 20100310 by BLFoley.
 */

#include <mylib.h>
#include <molecules.h>
/****************** set_residue_COM() *********************/
/* Finds the center of mass of the residue and places it in the
 * COM location.  COM calcualtion based on coords in xs.
 * This function is trivial, but the others aren't so much.
 */ 
void set_residue_COM(residue *r,atype *ATYPE, int xs){
r[0].COM=get_residue_COM(r,ATYPE,xs);
return;
}

/****************** set_molecule_COM() *********************/
/* Finds the center of mass of the molecule and places it in the
 * appropriate location.  Determines it recurively via the residues
 */
void set_molecule_COM(molecule *m,atype *ATYPE, int xs){
int tr=0;
double tmm=0;
m[0].COM.i=0;
m[0].COM.j=0;
m[0].COM.k=0; 
// calculate the center of mass
for(tr=0;tr<m[0].nr;tr++){ // COM = sum(m_i * r_i) / sum(m_i)
	set_residue_COM(&m[0].r[tr],ATYPE,xs);
        m[0].COM.i+=m[0].r[tr].COM.i*m[0].r[tr].m;
        m[0].COM.j+=m[0].r[tr].COM.j*m[0].r[tr].m;
        m[0].COM.k+=m[0].r[tr].COM.k*m[0].r[tr].m;
	tmm+=m[0].r[tr].m;
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
m[0].COM.i/=tmm;
m[0].COM.j/=tmm;
m[0].COM.k/=tmm; 
m[0].m=tmm;
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
}

/****************** set_assembly_molecule_COM() *********************/
/* Finds the center of mass of the assembly and places it in the
 * appropriate location.  This function assumes that the assembly has
 * a traditional structure, with molecules, residues and atoms all in
 * nested structures.  If only the atom or residue pointers contain
 * data, this function won't work.
 */
void set_assembly_molecule_COM(assembly *a,atype *ATYPE,int xs){
int tm=0;
double tmm=0; 

a[0].COM.i=0;
a[0].COM.j=0;
a[0].COM.k=0; 
// calculate the center of mass
for(tm=0;tm<a[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
	set_molecule_COM(&a[0].m[tm][0],ATYPE,xs);
        a[0].COM.i+=a[0].m[tm][0].COM.i*a[0].m[tm][0].m;
        a[0].COM.j+=a[0].m[tm][0].COM.j*a[0].m[tm][0].m;
        a[0].COM.k+=a[0].m[tm][0].COM.k*a[0].m[tm][0].m;
	tmm+=a[0].m[tm][0].m;
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
a[0].COM.i/=tmm;
a[0].COM.j/=tmm;
a[0].COM.k/=tmm; 
a[0].mass=tmm;
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
}

/****************** set_ensemble_COM() *********************/
/* Finds the center of mass of the ensemble and places it in the
 * appropriate location.  
 * NOTE!!  Assemblies within ensembles are assumed to be redundant.
 * 	That is, the molecules in the ensemble are also represented
 * 	within the assemblies (copying the molecule pointers keeps 
 * 	the data from taking extra space).  So, when calculating the
 * 	center of mass for the ensemble, this function ONLY considers
 * 	molecules and not assemblies.
 */
void set_ensemble_COM(ensemble *e,atype *ATYPE, int xs){
int tm=0;
double tmm=0; 
e[0].COM.i=0;
e[0].COM.j=0;
e[0].COM.k=0; 
// calculate the center of mass
for(tm=0;tm<e[0].nm;tm++){ // COM = sum(m_i * r_i) / sum(m_i)
	set_molecule_COM(&e[0].m[tm],ATYPE,xs);
        e[0].COM.i+=e[0].m[tm].COM.i*e[0].m[tm].m;
        e[0].COM.j+=e[0].m[tm].COM.j*e[0].m[tm].m;
        e[0].COM.k+=e[0].m[tm].COM.k*e[0].m[tm].m;
	tmm+=e[0].m[tm].m;
	}
//printf("center of mass is at: %f %f %f \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
e[0].COM.i/=tmm;
e[0].COM.j/=tmm;
e[0].COM.k/=tmm; 
e[0].mass=tmm;
//printf("center of mass is at: %20.12e %20.12e %20.12e \n",m[0].COM.i,m[0].COM.j,m[0].COM.k);
return;
} 

