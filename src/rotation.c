/* Author: Michael Tessier*/
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
//#include <math.h>
void rollMolecule(molecule* mol, double rad)
{
  //Tells the user what's going on
  printf("Rolling molecule approx. %.0lf degrees...\n",(rad*57.2957795));
  int i,j;
  residue* curRes; atom* curAtm;
  coord_3D temp; temp.i = 0; temp.j = 0; temp.k = 0;
  for(i = 0; i < (*mol).nr; i++)
  {
   curRes = ((*mol).r+i);
 
   for(j = 0; j < (*curRes).na; j++)
   {
    curAtm = ((*curRes).a+j);
    //Rotating about the X-axis
    (*curAtm).x.i = (*curAtm).x.i;
    temp.j=((*curAtm).x.j*cos(rad))+((*curAtm).x.k*sin(rad));
    temp.k=((*curAtm).x.j*-sin(rad))+((*curAtm).x.k*cos(rad));
    (*curAtm).x.j = temp.j; (*curAtm).x.k = temp.k;
   }
  }
}

void pitchMolecule(molecule* mol, double rad)
{
  //Tells the user what's going on
  printf("Pitching molecule approx. %.0lf degrees...\n",(rad*57.2957795));
  int i,j;
  residue* curRes; atom* curAtm;
  coord_3D temp; temp.i = 0; temp.j = 0; temp.k = 0;
  for(i = 0; i < (*mol).nr; i++)
  {
   curRes = ((*mol).r+i);
 
   for(j = 0; j < (*curRes).na; j++)
   {
    curAtm = ((*curRes).a+j);
    //Rotating about the Y-axis
    temp.i=((*curAtm).x.i*cos(rad))+((*curAtm).x.k*-sin(rad));
    (*curAtm).x.j = (*curAtm).x.j;
    temp.k=((*curAtm).x.i*sin(rad))+((*curAtm).x.k*cos(rad));
    (*curAtm).x.i = temp.i; (*curAtm).x.k = temp.k;
   }
  }
}

void yawMolecule(molecule* mol, double rad)
{
  //Tells the user what's going on
  printf("Yawing molecule approx. %.0lf degrees...\n",(rad*57.2957795));
  int i,j;
  residue* curRes; atom* curAtm;
  coord_3D temp; temp.i = 0; temp.j = 0; temp.k = 0;
  for(i = 0; i < (*mol).nr; i++)
  {
   curRes = ((*mol).r+i);
 
   for(j = 0; j < (*curRes).na; j++)
   {
    curAtm = ((*curRes).a+j);
    //Rotating about the Z-axis
    temp.i=((*curAtm).x.i*cos(rad))+((*curAtm).x.j*sin(rad));
    temp.j=((*curAtm).x.i*-sin(rad))+((*curAtm).x.j*cos(rad));
    (*curAtm).x.i = temp.i; (*curAtm).x.j = temp.j;
    (*curAtm).x.k = (*curAtm).x.k;

   }
  }
}
void rollAssembly(assembly* asmbl, double rad)
{
  //Tells the user what's going on
  printf("Rolling assembly approx. %.0lf degrees...\n",(rad*57.2957795));
  int i,j,k;
  molecule* mol; residue* curRes; atom* curAtm;
  coord_3D temp; temp.i = 0; temp.j = 0; temp.k = 0;
  for(k = 0; k < (*asmbl).nm; k++)
  {
   mol = ((*asmbl).m+k)[0]; // changed by BLF on 20080622 -- might need revisiting
   for(i = 0; i < (*mol).nr; i++)
   {
    curRes = ((*mol).r+i);
  
    for(j = 0; j < (*curRes).na; j++)
    {
     curAtm = ((*curRes).a+j);
     //Rotating about the X-axis
     (*curAtm).x.i = (*curAtm).x.i;
     temp.j=((*curAtm).x.j*cos(rad))+((*curAtm).x.k*sin(rad));
     temp.k=((*curAtm).x.j*-sin(rad))+((*curAtm).x.k*cos(rad));
     (*curAtm).x.j = temp.j; (*curAtm).x.k = temp.k;
    }
   }
  }
}

void pitchAssembly(assembly* asmbl, double rad)
{
  //Tells the user what's going on
  printf("Pitching assembly approx. %.0lf degrees...\n",(rad*57.2957795));
  int i,j,k;
  molecule* mol; residue* curRes; atom* curAtm;
  coord_3D temp; temp.i = 0; temp.j = 0; temp.k = 0;
  for(k = 0; k < (*asmbl).nm; k++)
  {
   mol = ((*asmbl).m+k)[0]; // changed by BLF on 20080622 -- might need revisiting
   for(i = 0; i < (*mol).nr; i++)
   {
    curRes = ((*mol).r+i);
  
    for(j = 0; j < (*curRes).na; j++)
    {
     curAtm = ((*curRes).a+j);
     //Rotating about the Y-axis
     temp.i=((*curAtm).x.i*cos(rad))+((*curAtm).x.k*-sin(rad));
     (*curAtm).x.j = (*curAtm).x.j;
     temp.k=((*curAtm).x.i*sin(rad))+((*curAtm).x.k*cos(rad));
     (*curAtm).x.i = temp.i; (*curAtm).x.k = temp.k;
    }
   }
  }
}

void yawAssembly(assembly* asmbl, double rad)
{
  //Tells the user what's going on
  printf("Yawing assembly approx. %.0lf degrees...\n",(rad*57.2957795));
  int i,j,k;
  molecule* mol; residue* curRes; atom* curAtm;
  coord_3D temp; temp.i = 0; temp.j = 0; temp.k = 0;
  for(k = 0; k < (*asmbl).nm; k++)
  {
   mol = ((*asmbl).m+k)[0]; // changed by BLF on 20080622 -- might need revisiting
   for(i = 0; i < (*mol).nr; i++)
   {
    curRes = ((*mol).r+i);
  
    for(j = 0; j < (*curRes).na; j++)
    {
     curAtm = ((*curRes).a+j);
     //Rotating about the Z-axis
     temp.i=((*curAtm).x.i*cos(rad))+((*curAtm).x.j*sin(rad));
     temp.j=((*curAtm).x.i*-sin(rad))+((*curAtm).x.j*cos(rad));
     (*curAtm).x.i = temp.i; (*curAtm).x.j = temp.j;
     (*curAtm).x.k = (*curAtm).x.k;
 
    }
   }
  }
}
