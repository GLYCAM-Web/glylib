/* File read_amber_prmtop.c begun on 20080425 by BLFoley
	Purpose: read amber prmtop file info into an assembly 
	See also: amber_prmtop_structs.h */
//#include "molecules.h"
//#include "general.h"
#include "AMBER/amber_prmtop.h"
//#include "amber_prmtop_structs.h"
//#include "declarations.h"
//void read_amber_prmtop_asis(fileset F,amber_prmtop *P);
// the amber_prmtop *P should be -already- allocated and initialized


/********** assembly *read_amber_prmtop(fileset F) **********/
/* Main function for reading the prmtop file 
 * the amber_prmtop *P should be -already- allocated and initialized 
 * BUT, it should be otherwise empty  */
void read_amber_prmtop_asis(fileset F,amber_prmtop *P){
int aa=0,anl=0,sec=0,ent=0,fld=0; // generic, line, section, entry & field counters
int newflag=0,al=0,newline=0,iseof=1; // utility flags
char line[501],over[501],tmp[81],*s1,*s2; 

// Open files and initialize some variables
//F.F=myfreopen(F.N,"r",F.F); // don't depend on calling function to have opened
F.F=myfopen(F.N,"r"); // 
P[0].FN=strdup(F.N); // the name of the file from which these data were taken
P[0].nS=0; // the prmtop structure should be empty
P[0].S=(amber_prmtop_section*)calloc(1,sizeof(amber_prmtop_section));
P[0].SN=(char**)calloc(1,sizeof(char*));

//
// Get the first line, which should be the "VERSION"
//
if(fgets(line,500,F.F)==NULL){read_fneek("Problem in read_amber_prmtop",1,0,F.N);}
// See if the first 500 chars includes a newline...
if(strchr(line,'\n')==NULL){ // if we didn't get the whole line
	fld=500;
	fprintf(stderr,"Unusually long VERSION line.\n");
       	fprintf(stderr,"-- Ignoring characters beyond 500.\n");
	if(fgets(over,500,F.F)==NULL){// try to get rest of line
		read_fneek("Problem in read_amber_prmtop.",1,fld,F.N);}
	while(strchr(over,'\n')==NULL){// try to get rest of line
		fld+=500;
		if(fgets(over,500,F.F)==NULL){
			read_fneek("Problem in read_amber_prmtop.",1,fld,F.N);}
		}
	} // close lengthy line check...
if(strncmp("\%VERSION",line,8)!=0){
	read_eek("First line not \%VERSION.",F.N);}
P[0].VERSION=strdup(line); 
anl=2;
fld=0;
//
// grab next line, which should be a flag line
//
if(fgets(line,500,F.F)==NULL){ 
	read_fneek("Problem in read_amber_prmtop",anl,fld,F.N);} 
// Check the line -- it should be a short-ish FLAG line
if(strchr(line,'\n')==NULL){// die if line too long
read_eek("FLAG line over 500 chars... \n",F.N);}
//
// Read the FLAG entry
//
if(strncmp("\%FLAG",line,5)!=0){ // must start with %FLAG
read_fneek("Entry not FLAG in read_amber_prmtop",anl,fld,F.N);}
// Parse the FLAG line
if(sscanf(line,"%s %s",over,tmp)!=2){mywhine("Problem in read_amber_prmtop with 1st FLAG read.");}
newflag=0;
//
// Start main read of topology file
// 
while(newflag==0){// while we have found another flag entry
	// We have a FLAG entry, so increment the number of sections
//printf("line (FLAG) is >>%s<<\n",line);
	sec=P[0].nS;
	P[0].nS++; 
	// Allocate memory 
	P[0].S=(amber_prmtop_section*)realloc(P[0].S,P[0].nS*sizeof(amber_prmtop_section));
	P[0].SN=(char**)realloc(P[0].SN,P[0].nS*sizeof(char*));
	P[0].SN[sec]=strdup(tmp); // record FLAG entry in
	P[0].S[sec].N=strdup(tmp); // both locations 
	//
	// Read the FORMAT entry
	// 
	anl++;
	if(fgets(line,500,F.F)==NULL){ // if we can't read the line, die
		read_fneek("Problem in read_amber_prmtop.",anl,fld,F.N);}
	if(strchr(line,'\n')==NULL){ // FORMAT lines should be short-ish
		read_fneek("FORMAT line over 500 chars... \n",anl,fld,F.N);}
	if(strncmp("\%FORMAT",line,7)!=0){// must start with %FORMAT
		read_fneek("Entry not FORMAT in read_amber_prmtop",anl,fld,F.N);} 
//printf("line (FORMAT) is >>%s<<\n",line); 
	// split line by parentheses as tokens
	s1=strdup(strtok(line,"()"));// toss out the "%FORMAT" part
	P[0].S[sec].FORMAT=strdup(strtok(NULL,"()")); // save next
	s1=strdup(P[0].S[sec].FORMAT); // save a char string to tear apart
	// parse the format string
	s2=strdup(strtok(s1,"aAiIfFeE")); // up to the format identifier
	sscanf(s2,"%d",&P[0].S[sec].npl); // max entries expected per line
	s2=strdup(strtok(NULL,".")); // up to any periods
	if((strchr(s2,'S')!=NULL)||(strchr(s2,'N')!=NULL)){over[0]='0';} // EN or ES
	sscanf(s2,"%d",&P[0].S[sec].nc); // number of chars in format 
//printf("parsed format:  npl is %d ; nc is %d \n",P[0].S[sec].npl,P[0].S[sec].nc);
	//
	// Read the section
	// 
	P[0].S[sec].D=(char**)calloc(1,sizeof(char*)); // initially allocate **D
	P[0].S[sec].nt=0; // total pieces of data found 
	ent=0;
	while(fgets(line,6,F.F)!=NULL){
//printf("the fgets's line is >>%s<<\n",line);
//printf("in ascii decimal: %d %d %d %d %d %d\n",line[0],line[1],line[2],line[3],line[4],line[5]);
		if(strcspn(line,"\0")==0){
//printf("found empty line\n");
			continue;
			}
		//ae=fgetc(F.F);
//printf("ae is >>%c<<\n",ae);
		//if(ae==EOF){
//printf("Found end of file, stopping.\n");
			//break;}
		//else{ungetc(ae,F.F);}
		anl++; // we have a new line
		fld=0; // start field count over
		// check to see if first 5 chars match "%FLAG"
		if(strncmp(line,"\%FLAG",5)==0){ 
//printf("found a FLAG line\n");
			break;} // if so, stop this
		if(strcspn(line,"\0")==0){
			continue;
			}
		if(strcspn(line,"\n")==0){
			fprintf(stderr,"Unexpected blank line in read_amber_prmtop\n");
			fprintf(stderr,"   at line number %d in file %s.  Ignoring.\n",anl,F.N);
			fprintf(stderr,"The section (%d) FLAG is:  %s\n",sec,P[0].SN[sec]);
			continue;
			}
// 
// Since the "flag check" field length might not be five...
//
if(P[0].S[sec].nc<5){ // the flag read contains more than one entry
	// find out how much memory we need
//printf("P[0].S[sec].nc<5\n");
	// copy extra bits to the new strings
	// being long here because it's easier to write...
	ent=P[0].S[sec].nt;
	switch(P[0].S[sec].nc){
		case 0:
			mywhine("format character length reports as zero in read_amber_prmtop\n");
			break;
		case 1:
			P[0].S[sec].nt+=5;
			P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
			al=0;
			for(aa=0;aa<5;aa++){
				P[0].S[sec].D[ent]=(char*)calloc(2,sizeof(char));
				P[0].S[sec].D[ent][0]=line[al];
				P[0].S[sec].D[ent][1]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
				al++;
				ent++;
				}
			break;
		case 2:
			P[0].S[sec].nt+=3;
			P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
			P[0].S[sec].D[ent]=(char*)calloc(3,sizeof(char));
			P[0].S[sec].D[ent][0]=line[0];
			P[0].S[sec].D[ent][1]=line[1];
			P[0].S[sec].D[ent][2]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
			P[0].S[sec].D[ent+1]=(char*)calloc(3,sizeof(char));
			P[0].S[sec].D[ent+1][0]=line[2];
			P[0].S[sec].D[ent+1][1]=line[3];
			P[0].S[sec].D[ent+1][2]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent+1,P[0].S[sec].D[ent+1]);
			P[0].S[sec].D[ent+2]=(char*)calloc(3,sizeof(char));
			P[0].S[sec].D[ent+2][0]=line[4];
			P[0].S[sec].D[ent+2][1]=fgetc(F.F);
			P[0].S[sec].D[ent+2][2]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent+2,P[0].S[sec].D[ent+1]);
			ent+=3;
			break;
		case 3:
			P[0].S[sec].nt+=2;
			P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
			P[0].S[sec].D[ent]=(char*)calloc(4,sizeof(char));
			P[0].S[sec].D[ent][0]=line[0];
			P[0].S[sec].D[ent][1]=line[1];
			P[0].S[sec].D[ent][2]=line[2];
			P[0].S[sec].D[ent][3]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
			P[0].S[sec].D[ent+1]=(char*)calloc(4,sizeof(char));
			P[0].S[sec].D[ent+1][0]=line[3];
			P[0].S[sec].D[ent+1][1]=line[4];
			P[0].S[sec].D[ent+1][2]=fgetc(F.F);
			P[0].S[sec].D[ent+1][3]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent+1,P[0].S[sec].D[ent+1]);
			ent+=2;
			break;
		case 4:
			P[0].S[sec].nt+=2;
			P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
			P[0].S[sec].D[ent]=(char*)calloc(5,sizeof(char));
			P[0].S[sec].D[ent][0]=line[0];
			P[0].S[sec].D[ent][1]=line[1];
			P[0].S[sec].D[ent][2]=line[2];
			P[0].S[sec].D[ent][3]=line[3];
			P[0].S[sec].D[ent][4]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
			P[0].S[sec].D[ent+1]=(char*)calloc(5,sizeof(char));
			P[0].S[sec].D[ent+1][0]=line[4];
			P[0].S[sec].D[ent+1][1]=fgetc(F.F);
			P[0].S[sec].D[ent+1][2]=fgetc(F.F);
			P[0].S[sec].D[ent+1][3]=fgetc(F.F);
			P[0].S[sec].D[ent][4]='\0';
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent+1,P[0].S[sec].D[ent+1]);
			ent+=2;
			break;
		default:
			mywhine("nc should be less than five in read_amber_prmtop, but seems like more than five...\n");
		}
	}
if(P[0].S[sec].nc==5){ // the read contains exactly one entry
//printf("P[0].S[sec].nc==5\n");
	ent=P[0].S[sec].nt; // for convenience
	// allocate space for the next entry
	P[0].S[sec].nt++;
	P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
	P[0].S[sec].D[ent]=(char*)calloc(6,sizeof(char));
	strncpy(P[0].S[sec].D[ent],line,6);
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
	}
if(P[0].S[sec].nc>5){ // the flag read contains less than one entry
//printf("P[0].S[sec].nc>5\n"); 
	ent=P[0].S[sec].nt; // for convenience
	P[0].S[sec].nt++; // we will only need one more character space
	P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
//printf("line is >>>%s<<<\n",line);
	// fill out the "line" string to the proper length
	for(aa=5;aa<(P[0].S[sec].nc);aa++){ 
		line[aa]=fgetc(F.F); 
//printf("the new character is >>>%c<<<\n",line[aa]);
		}
	line[aa]='\0'; // terminate the string 
//printf("line is -now- >>>%s<<<\n",line);
	P[0].S[sec].D[ent]=(char*)calloc(P[0].S[sec].nc+1,sizeof(char));
	strncpy(P[0].S[sec].D[ent],line,P[0].S[sec].nc+1);
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
	}
// 
// Now read the rest of the line...
// 
//printf("reading the rest of the line.\n");
newline=1;
while(newline!=0){
	if(fgets(line,(P[0].S[sec].nc+1),F.F)==NULL){
		read_fneek("problem in read_amber_prmtop",anl,ent,F.N);}
//printf("the line is >>>%s<<<\n",line);
//printf("in ascii decimal: ");
//for(ae=0;ae<P[0].S[sec].nc+1;ae++){printf(" %d ",line[ae]);}
//printf("\n");
		if(strcspn(line,"\0")==0){
//printf("found empty string\n");
			newline=0;
			break;
			}
	if(strcspn(line,"\n")==0){ // at end of line, no new data
		newline=0;
		break;}
	if(strchr(line,'\n')!=NULL){ // this is the end of a line, but with data 
		for(aa=0;aa<strlen(line);aa++){ //check for all whitespace
			if(isspace(line[aa])==0) newline=-1; // something isn't white
			}
		if(newline==-1){
			line[strcspn(line,"\n")]='\0'; 
			ent=P[0].S[sec].nt; // for convenience
			// allocate space for the next entry
			P[0].S[sec].nt++;
			P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
			P[0].S[sec].D[ent]=(char*)calloc(P[0].S[sec].nc+1,sizeof(char));
			strncpy(P[0].S[sec].D[ent],line,P[0].S[sec].nc+1);
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
			}
		newline=0;
		break;
		}
	ent=P[0].S[sec].nt; // for convenience
	// allocate space for the next entry
	P[0].S[sec].nt++;
	P[0].S[sec].D=(char**)realloc(P[0].S[sec].D,P[0].S[sec].nt*sizeof(char*));
	P[0].S[sec].D[ent]=(char*)calloc(P[0].S[sec].nc+1,sizeof(char));
	strncpy(P[0].S[sec].D[ent],line,P[0].S[sec].nc+1); 
//printf("the entry for P[0].S[%d].D[%d] is -now- >>>%s<<<\n",sec,ent,P[0].S[sec].D[ent]);
	}// close while not a new line
	}// close fgets within a section 

	if(newflag==-1){break;}
	if(fgets(line,500,F.F)==NULL){ 
		if(fgetc(F.F)==EOF){
			iseof=1;
			newflag=-1;
			break;}
		else{read_fneek("Problem in FLAG line read read_amber_prmtop",anl,0,F.N);}
		} 
	if(newflag==-1){break;}
	// If here, then we found a FLAG, read it in, exit loop and start new
	// Check the line -- it should be a short-ish FLAG line
	if(strchr(line,'\n')==NULL){// die if line too long
	read_eek("FLAG line over 500 chars... \n",F.N);}
	//
	// Read the FLAG entry
	//
	// Parse the FLAG line
	if(sscanf(line,"%s",tmp)!=1){
		read_fneek("Problem FLAG entry scan, read_amber_prmtop.",anl,1,F.N);}
	} // close fgets for reading file

fclose(F.F);

/*F.F=myfopen("test_rewrite_of_prmtop","w");
fprintf(F.F,"%s",P[0].VERSION);
for(sec=0;sec<P[0].nS;sec++){ // for each section found
	fprintf(F.F,"%%FLAG %s\n",P[0].S[sec].N);
	fprintf(F.F,"%%FORMAT(%s)\n",P[0].S[sec].FORMAT);
	for(aa=0;aa<P[0].S[sec].nt;aa++){
		fprintf(F.F,"%s",P[0].S[sec].D[aa]);
		//fprintf(F.F,"D[%d] is >>>%s<<<\n",aa,P[0].S[sec].D[aa]);
		fflush(F.F);
//printf("\naa is %d ; P[0].S[%d].npl is %d ; (aa+1)%%npl is %d\n",aa,sec,P[0].S[sec].npl,);
		if((((aa+1)%P[0].S[sec].npl))==0){fprintf(F.F,"\n");} 
		//if((((aa+1)*P[0].S[sec].nc))%80==0){printf("\n");} 
		}
	if(((aa)*P[0].S[sec].nc)%80!=0){fprintf(F.F,"\n");} 
	}
fclose(F.F);
*/

return;
} 
