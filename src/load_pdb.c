/*Author: Michael Tessier*/
#include <load_pdb.h>
//#include "../inc/load_pdb.h"
assembly* load_pdb(char* file_name)
{
  char* temp;
  WARN=1;SILENT=1;INWC=0;DATE=0;UNCUT=1;UNx=1;UNy=1;UNz=1;
  temp = suf; strcpy(temp,".pdb");
  temp = sufc; strcpy(temp,"_change.txt"); ACT='0';
  DEBUG=-1;LASTRES=-1;LASTOKX=0;LASTOKY=0;LASTOKZ=0;
  UNCTOL=0;CRYX=0;CRYY=0;CRYZ=0;LASTX=0;LASTY=0;LASTZ=0;
  //int resNum;
  int ma;
  //molecule* mol;
  assembly* asmbl;
  IN = myfopen(file_name, "r");
 /*scans pdb file for number of lines
  allocates memory for line structures
  initializes structures that contain line formats*/
  init_struct();
  rewind(IN);
  //Read in lines to ln structure
  printf("Reading in %s...\n",file_name);
  for(ma=0;ma<INWC;ma++){
  if(DEBUG>=1){printf("made it to here main line loop %d ...\n", ma);} 
	rwm_line(ma+1);
	}
  //Determine the number of molecules in the pdb
  printf("There are %d molecule(s)\n",howManyMolecules());
  asmbl = getAssembly();
  printf("PDB Information Successfully Read.\n");
  free(ln);
  return asmbl;
}

int howManyMolecules()
{
 int i; int counter = 0;

 for(i = 0; i < INWC; i++)
  if(endOfMol((ln+i)) == 1)//If it is a TER, LINK, or CONECT card
   counter++;
 if(counter == 0 || isAtom(ln+(INWC-1)) == 1)
  counter ++;
 return counter;
}

assembly* getAssembly()
{
  int molNum = howManyMolecules();
  int i = 0;int j = 0;int k = 0;int l = 0;int atom_num;double x,y,z;
  char* temp;char name[6] = ""; char* atmName = name;
  assembly* asmbl = (assembly*) malloc (1 * sizeof(assembly));
  (*asmbl).nm = molNum;
  molecule* mol = (molecule*) malloc (molNum * sizeof(molecule));
  (*asmbl).m[0] = mol; // changed by BLF on 20080622 -- might need revisiting
  molecule* curMol; residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  curMol = (mol+j);
  (*curMol).nr = findTotalResidue(i);
  (*curMol).r = (residue*) malloc ((*curMol).nr*sizeof(residue));
  getResInfo((*curMol).r,i);
  curRes = ((*curMol).r+k);
  (*curRes).a = (atom*) malloc ((*curRes).na*sizeof(atom));
  curAtm = ((*curRes).a+l);
  for(i = 0; i < INWC; i++)
  {
    //If it is a TER, LINK, or CONECT card
    if(endOfMol((ln+i)) == 1)
    {
      j++;
      curMol = (mol+j);
      (*curMol).nr = findTotalResidue(i+1);
      (*curMol).r = (residue*) malloc ((*curMol).nr*sizeof(residue));
      getResInfo((*curMol).r,i+1);
      curRes = ((*curMol).r+k);
      k = 0;l = 0;
      (*curRes).a = (atom*) malloc ((*curRes).na*sizeof(atom));
      curAtm = ((*curRes).a+l);
    }
    else if(isAtom((ln+i)) == 1)
    {
      if(l == (*curRes).na)
      {
       k++;l = 0;
       curRes = ((*curMol).r+k);
       (*curRes).a = (atom*) malloc ((*curRes).na*sizeof(atom)); 
       //curAtm = ((*curRes).a+l);
      }
      curAtm = ((*curRes).a+l);
      //Getting the atom #
      temp  = (*(ln+i)).f[1].c;
      sscanf(temp,"%d",&atom_num);
      //Getting the atom name
      temp  = (*(ln+i)).f[3].c;
      sscanf(temp,"%s",atmName);
      //Getting the X coordinate
      temp  = (*(ln+i)).f[11].c;
      sscanf(temp,"%lf",&x);
      //Getting the Y coordinate
      temp  = (*(ln+i)).f[12].c;
      sscanf(temp,"%lf",&y);
      //Getting the Z coordinate
      temp  = (*(ln+i)).f[13].c;
      sscanf(temp,"%lf",&z);   

      (*curAtm).n = atom_num;//set the atom #
      strcpy((*curAtm).N,atmName);//set the atom name
      //set the atom's coordinates
      (*curAtm).x.i = x; (*curAtm).x.j = y; (*curAtm).x.k = z;
      l++;
    }
  }
 return asmbl;
}

int findTotalResidue(int start)
{
  int ret = 0;int i = start;int curRes = 0;int thisRes;
  char* temp;
  while(endOfMol((ln+i)) == 0 && i != INWC)
  {
    if(isAtom((ln+i)) == 1)
    {
      temp  = (*(ln+i)).f[8].c;
      sscanf(temp,"%d",&thisRes);
      if(thisRes != curRes)
      { 
        ret ++;
        curRes = thisRes;
      }
    }
   i++;
  }
  return ret;
}

int endOfMol(linedef* line)
{
  if(((*line).a == 4 && (*line).b == 2)||
     ((*line).a == 2 && ((*line).b == 9||(*line).b == 3)))
    return 1;
  else
    return 0;
}

int isAtom(linedef* line)
{
  if(((*line).a == 2 || (*line).a == 3) && (*line).b == 1)
    return 1;
  else
    return 0;
}
void getResInfo(residue* res, int start)
{
  int atmTot = 0;int resNum,thisRes;int i = start;int j = 0;
  char name [10];char* temp;char* resName = name;
  while(isAtom(ln+i) == 0 && i < INWC)
   i++;
  //Getting the residue name
  temp = (*(ln+i)).f[5].c;
  sscanf(temp,"%s",resName);
  strcpy((*(res+j)).N, resName);
  //Getting the residue #
  temp  = (*(ln+i)).f[8].c;
  sscanf(temp,"%d",&thisRes);
  resNum = thisRes;
  (*(res+j)).n = resNum;
  i = start;
  while(endOfMol((ln+i)) == 0 && i != INWC)
  {
    if(isAtom((ln+i)) == 1)
    {
      temp  = (*(ln+i)).f[8].c;
      sscanf(temp,"%d",&thisRes);
      if(thisRes != resNum)
      { 
        //Setting the total # of atoms on the old res
        resNum = thisRes;
        (*(res+j)).na = atmTot;
        atmTot = 0; j++;
        //Getting the residue name for the new res
        temp = (*(ln+i)).f[5].c;
        sscanf(temp,"%s",resName);
        strcpy((*(res+j)).N, resName);
        //Setting the residue # on the new res
        (*(res+j)).n = resNum;
      }
      else
        atmTot++;
    }
   i++;
  }
 (*(res+j)).na = atmTot;
}
