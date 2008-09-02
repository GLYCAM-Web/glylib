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

molecule load_pdb_from_slurp(fileslurp in){
 int i,hldr;
 int minRes = -1;
 char* temp;
 //This should be done outside of this function
 fileslurp pdb = isolateInputPDB(in);
 init_struct_slurp(in);
 for(i = 0; i < pdb.n; i++){
  rwm_line_char(pdb.L[i],i+1);
  if(isAtom((ln+i))){
   temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&hldr);
   if(minRes > hldr || minRes < 0){minRes = hldr;}
  }
 }
 INWC = pdb.n;
 return (*getMolecule(minRes));
}

fileslurp isolateInputPDB(fileslurp S){
 int i,itr;
 char* plc;
 fileslurp O; O.n = 0; itr = 0;
 for(i = 0; i < S.n; i++){
  if(strstr(S.L[i],"INPUT-PDBQ: ")!=NULL ||
     strstr(S.L[i],"INPUT-LIGAND-PDBQT: ")!=NULL)
   O.n++;
 }
 O.L = (char**)calloc(O.n,sizeof(char*));
 for(i = 0; i < S.n; i++){
  if(strstr(S.L[i],"INPUT-PDBQ: ")!=NULL ||
     strstr(S.L[i],"INPUT-LIGAND-PDBQT: ")!=NULL){
   if(strstr(S.L[i],"INPUT-PDBQ: ")!=NULL)
    plc = S.L[i]+12;
   else
    plc = S.L[i]+20;
   O.L[itr] = strdup(plc); itr++;
  }
 }
 return O;
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
  (*asmbl).m = (molecule**)calloc(1,sizeof(molecule*));//added by MNT on 20080806
  (*asmbl).m[0] = mol; // changed by BLF on 20080622 -- might need revisiting
  molecule* curMol; residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  curMol = (mol+j);
  (*curMol).nr = findTotalResidue(i,0);
  (*curMol).r = (residue*) malloc ((*curMol).nr*sizeof(residue));
  getResInfo((*curMol).r,i,1);
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
      (*curMol).nr = findTotalResidue(i+1,0);
      (*curMol).r = (residue*) malloc ((*curMol).nr*sizeof(residue));
      getResInfo((*curMol).r,i+1,1);
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
      (*curAtm).N = strdup(atmName);//set the atom name
      //set the atom's coordinates
      (*curAtm).x.i = x; (*curAtm).x.j = y; (*curAtm).x.k = z;
      l++;
    }
  }
 return asmbl;
}

int findTotalResidue(int start,int minRes)
{
  int ret=0,hRes=0;int i = start;int thisRes;//int curRes = 0;
  char* temp;
  while(endOfMol((ln+i)) == 0 && i != INWC)
  {
    if(isAtom((ln+i)) == 1)
    {
      temp  = (*(ln+i)).f[8].c;
      sscanf(temp,"%d",&thisRes);
      //The New Way
      if(hRes < thisRes)
       hRes = thisRes;       
      //THE OLD WAY
/*      if(thisRes != curRes)
      { 
        ret ++;
        curRes = thisRes;
      }*/
    }
   i++;
  }
  //More of the New Way
  if(minRes < 1){ret = hRes;}
  else{ret = hRes - (minRes-1);}
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
void getResInfo(residue* res, int start, int minRes)
{
  int atmTot = 0;int thisRes,hldr;int i = start;int j = 0;
  int resNum=-1;
  residue* curRes;
  char name [10];char* temp;char* resName = name;
  while(isAtom(ln+i) == 0 && i < INWC)
   i++;
  while(endOfMol((ln+i)) == 0 && i != INWC){
   temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&hldr);
   if(resNum != hldr){	//If this is a different residue than the one in the previous line
    resNum = hldr;	//Set resNum as the actual residue number
    thisRes=hldr-minRes;//Set thisRes as the index of the residue
    curRes = (res+thisRes);
    temp = (*(ln+i)).f[5].c; sscanf(temp,"%s",resName);
   }
   //If this residue has not been found yet
   if((*curRes).N == 0X0 || strcmp(resName,(*curRes).N) != 0){
    (*curRes).N = strdup(resName);	//Set the residue name
    (*curRes).n = resNum;		//Set the actual residue number
    (*curRes).na = 0;			//...and make sure the total # of atoms is 0
   }//Otherwise we can assume all of these have already been set
   (*curRes).na++;
   //..and then incriment to the next line
   i++;
  }

/*  //Getting the residue name
  temp = (*(ln+i)).f[5].c;
  sscanf(temp,"%s",resName);
  (*(res+j)).N = strdup(resName);
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
      
      //THE OLD WAY
      if(thisRes != resNum)
      { 
        //Setting the total # of atoms on the old res
        resNum = thisRes;
        (*(res+j)).na = atmTot;
        atmTot = 0; j++;
        //Getting the residue name for the new res
        temp = (*(ln+i)).f[5].c;
        sscanf(temp,"%s",resName);
        (*(res+j)).N = strdup(resName);
        //Setting the residue # on the new res
        (*(res+j)).n = resNum;
      }
      else
        atmTot++;
    }
   i++;
  }
 (*(res+j)).na = atmTot;*/
}

molecule* getMolecule(int minRes){
  int molNum = 1;
  int i = 0,rI = -1;//int j = 0;int k = 0;
  int ntX;
  //double x,y,z;
  double dblY;
  char* temp;char name[6] = ""; //char* atmName = name;
  molecule* mol = (molecule*)calloc(molNum,sizeof(molecule));
  residue* curRes;// residue* res;
  atom* curAtm;// atom* atm;
  //Get all the information for the first molecule
  (*mol).nr = findTotalResidue(i,minRes);
  (*mol).r = (residue*)calloc((*mol).nr,sizeof(residue));
  getResInfo((*mol).r,i,minRes);
  int l[(*mol).nr];
  for(i = 0; i < (*mol).nr; i++){
   (*mol).r[i].a = (atom*)calloc((*mol).r[i].na,sizeof(atom));
   l[i] = 0;
  }
  //resIdx = -1;
  for(i = 0; i < INWC; i++)
  {
    if(isAtom((ln+i)) == 1)
    {
      temp  = (*(ln+i)).f[8].c; sscanf(temp,"%d",&ntX);
      if(rI != ntX - minRes)
      {
       rI = ntX -minRes;
       curRes = ((*mol).r+rI);
      }
      curAtm = ((*curRes).a+l[rI]);
      //Getting the record name
      temp = (*(ln+i)).f[0].c; sscanf(temp,"%s",name);
      (*curAtm).D = strdup(name);
      //Getting the atom #
      temp  = (*(ln+i)).f[1].c; sscanf(temp,"%d",&ntX);
      (*curAtm).n = ntX;
      //Getting the atom name
      temp  = (*(ln+i)).f[3].c; sscanf(temp,"%s",name);
      (*curAtm).N = strdup(name);
      //Getting the chain identifier
      (*curAtm).cID = (*(ln+i)).f[7].c[0];

      //Getting the X coordinate
      temp  = (*(ln+i)).f[11].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.i = dblY;
      //Getting the Y coordinate
      temp  = (*(ln+i)).f[12].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.j = dblY;
      //Getting the Z coordinate
      temp  = (*(ln+i)).f[13].c; sscanf(temp,"%lf",&dblY); (*curAtm).x.k = dblY;   

      l[rI]++;
    }
  }
 return mol;
}

