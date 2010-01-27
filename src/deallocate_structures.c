#include "../inc/mylib.h"
#include "../inc/molecules.h"

void deallocateBondType(bond_type btp){
 //Free up the description
 if(btp.desc != NULL || btp.desc != 0x0){free(btp.desc);}
 //Free up the first and second atom types in bond
 if(btp.NT[0] != NULL || btp.NT[0] != 0x0){free(btp.NT[0]);}
 if(btp.NT[1] != NULL || btp.NT[1] != 0x0){free(btp.NT[1]);}

 return ;
}

void deallocateMolbond(molbond mlb){
 //Free up the Free Form Description
 if(mlb.D != NULL || mlb.D != 0x0){free(mlb.D);}
 //Free up the Type of Bond
 if(mlb.typ != NULL || mlb.typ != 0x0){deallocateBondType(mlb.typ[0]);}

 return ;
}

void deallocateConnectionTree(connection_tree con){
 //Free up the Incoming Bonds
// if((con.i != NULL || con.i != 0x0) && con.ni > 0){free(con.i);}
 //Free up the Outgoing Bonds
// if((con.o != NULL || con.o != 0x0) && con.no > 0){free(con.o);}
 //Free up the Description of Isomerism
// if(con.iso != NULL || con.iso != 0x0){free(con.iso);}
 //Free up the Angles
// if(con.angle != NULL || con.angle != 0x0){free(con.angle);}
 //Free up the Chirality
// if(con.chirot != NULL || con.chirot != 0x0){free(con.chirot);}
 //Free up the Distance
// if(con.distance != NULL || con.distance != 0x0){free(con.distance);}
 
 return ;
}

void deallocateBond(bond bnd){
 //Free up the Free Form Description
 if(bnd.D != NULL || bnd.D != 0x0){free(bnd.D);}
 //Free up the Type of Bond
 if(bnd.typ != NULL || bnd.typ != 0x0){deallocateBondType(bnd.typ[0]);}

 return ;
}

void deallocateBondset(bondset bst){
 int i; 
 //Free up the Bonds
 if((bst.b != NULL || bst.b != 0x0) && bst.n > 0){
  for( i = 0; i < bst.n; i++)
   deallocateBond(bst.b[i]);
  free(bst.b);
 }

 return ; 
}

void deallocateAtom(atom atm){
 int i;
 
 //Free up the Name
 if(atm.N != NULL || atm.N != 0X0){free(atm.N);}
 //Free up the Atom Type Name
 if(atm.T != NULL || atm.T != 0X0){free(atm.T);}
 //Free up the Free-Form Description
 if(atm.D != NULL || atm.D != 0X0){free(atm.D);}
 //Free up the Bond Structure
/* if((atm.b != NULL || atm.b != 0X0) && atm.nb > 0)
  for(i = 0; i < atm.nb; i++)
   deallocateMolbond(atm.b[i]);*/
 //Free up the Alternate Coordinates
 if((atm.xa != NULL || atm.xa != 0X0) && atm.nalt > 0){free(atm.xa);}
 //Free up the Vector Sets
 if((atm.v != NULL || atm.v != 0X0) && atm.nvec > 0){free(atm.v);}
 //Free up the Changes
 if((atm.ch != NULL || atm.ch != 0X0) && atm.nch > 0){free(atm.ch);}
 //Free up the Other Indicies
 if((atm.i != NULL || atm.i != 0X0) && atm.ni > 0){free(atm.i);}
 //Free up the Other Parameters
 if((atm.d != NULL || atm.d != 0X0) && atm.nd > 0){free(atm.d);}
 //Free up the list of Ensemble Indicies
 if((atm.ensi != NULL || atm.ensi != 0X0) && atm.nensi > 0){free(atm.ensi);}
 //Free up the Alternate Residue
 if(atm.sres != NULL || atm.sres != 0X0){free(atm.sres);}
 //Free up the Other Descriptors
 if((atm.OD != NULL || atm.OD != 0X0) && atm.nOD > 0){
  for(i = 0; i < atm.nOD; i++)
   free(atm.OD[i]);
  free(atm.OD);
 }
 return;
}

void deallocateResidue(residue res){
 int i;
 
 //Free up the Name
 if(res.N != NULL || res.N != 0X0){free(res.N);}
 //Free up the Free Form Description
 if(res.D != NULL || res.D != 0X0){free(res.D);}
 //Free up the Atoms
 if((res.a != NULL || res.a != 0X0) && res.na > 0){
  for(i = 0; i < res.na; i++)
   deallocateAtom(res.a[i]);
  free(res.a);
 }
 //Free up the Connection Tree
/* if((res.T != NULL || res.T != 0X0) && res.na > 0){
  for(i = 0; i < res.na; i++)
   deallocateConnectionTree(res.T[i]);
  free(res.T);
 }*/
 //Free up the Bondset
/* if((res.bs != NULL || res.bs != 0X0) && res.nbs > 0){
  for(i = 0; i < res.nbs; i++)
   deallocateBondset(res.bs[i]);
  free(res.bs);
 }*/
 //Free up the Bondset for Rings
/* if((res.rbs != NULL || res.rbs != 0X0) && res.nrbs > 0){
  for(i = 0; i < res.nrbs; i++)
   deallocateBondset(res.rbs[i]);
  free(res.rbs);
 }*/
 //Free up the Rings/Reference Centers 
 if((res.a != NULL || res.a != 0X0) && res.na > 0){free(res.rc);}
 //Free up the Ring Planes  
 if((res.rp != NULL || res.rp != 0X0) && res.nrp > 0){free(res.rp);}
 //Free up the Other Indices 
 if((res.i != NULL || res.i != 0X0) && res.ni > 0){free(res.i);}
 //Free up the Other Parameters 
 if((res.d != NULL || res.d != 0X0) && res.nd > 0){free(res.d);}
 //Free up the list of Ensemble Indices
 if((res.ensi != NULL || res.ensi != 0X0) && res.nensi > 0){free(res.ensi);}
 //Free up the Other Descriptors
 if((res.OD != NULL || res.OD != 0X0) && res.nOD > 0){
  for(i = 0; i < res.nOD; i++)
   free(res.OD[i]);
  free(res.OD);
 }
 return ;
}

void deallocateMolecule(molecule mol){
 int i;

 //Free up the Name
 if(mol.N != NULL || mol.N != 0X0){free(mol.N);}
 //Free up the Free Form Descriptor
 if(mol.D != NULL || mol.D != 0X0){free(mol.D);}
 //The atoms should be removed via their respective residues
 //Free up the Residues
 if((mol.r != NULL || mol.r != 0X0) && mol.nr > 0){
  for(i = 0; i < mol.nr; i++)
   deallocateResidue(mol.r[i]);
  free(mol.r);
 }
 //Free up the Connection Tree
/* if((mol.rT != NULL || mol.rT != 0X0) && res.nr > 0){
  for(i = 0; i < mol.nr; i++)
   deallocateConnectionTree(mol.rT[i]);
  free(mol.rT);
 }*/
 //Free up the Molbond
/* if((mol.rb != NULL || mol.rb != 0X0) && mol.nrb > 0){
  for(i = 0; i < mol.nrb; i++)
   deallocateMolbond(mol.rb[i]);
  free(mol.rb);
 }*/ 
 //Free up the Molbondset
/* if((mol.rbs != NULL || mol.rbs != 0X0) && mol.nrbs > 0){
  for(i = 0; i < mol.nrbs; i++)
   deallocateMolbondset(mol.rbs[i]);
  free(mol.rbs);
 }*/ 
 //Free up the Additional References
 if((mol.rc != NULL || mol.rc != 0X0) && mol.nrc > 0){free(mol.rc);}
 //Free up the Box Info
/* if((mol.BOX != NULL || mol.BOX != 0X0) && mol.nBOX > 0){
  for(i = 0; i < mol.nBOX; i++)
   deallocateBoxinfo(mol.Box[i]);
  free(mol.BOX);
 }*/
 //Free up the Other Indices
 if((mol.oi != NULL || mol.oi != 0X0) && mol.noi > 0){free(mol.oi);}
 //Free up the Other Parameters
 if((mol.d != NULL || mol.d != 0X0) && mol.nd > 0){free(mol.d);}
 //Free up the Ensemble Indices
 if((mol.ensi != NULL || mol.ensi != 0X0) && mol.nensi > 0){free(mol.ensi);}
 //Free up the Other Descriptors
 if((mol.OD != NULL || mol.OD != 0X0) && mol.nOD > 0){
  for(i = 0; i < mol.nOD; i++)
   free(mol.OD[i]);
  free(mol.OD);
 }
 return ;
}
