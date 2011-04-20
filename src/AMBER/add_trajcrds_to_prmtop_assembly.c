/** \file add_trajcrds_to_prmtop_assembly.c 
Function for adding the contents of an AMBER coordinate trajectory to an 
existing assembly.

Please read the following carefully if this is your first use.

Usage: 
void add_trajcrds_to_prmtop_assembly(
	fileset F, ///< The file containing the data *Must* be open
	assembly *A, ///< Assembly to which the crds should be added
	char ftype, ///< 'c'=coordinate; 'v'=velocity; 'r'=restart
	int offset ///< read in the offset-th trajectory -- zero=current
	)


*** Storage of Coordinates ***

This function only adds *one* set of coordinates at a time.  If the file
pointer is not altered between calls, each time this function is called it
will read in the next set of coordinates.

The coordinates are stored via the assembly's atom double-pointers.  The 
**atom pointers must be set correctly because the the coordinates will be 
added in the order dictated by the **atom locations.  

The coordinates are added into an 'xa' location in the incoming assembly. 
Regarding the 'xa' locations:  The program will check through the "nalt" for 
each atom.  It will reallocate them all to the same size if they differ.  
The new coordinates will be added starting in the lowest xa array location 
that is empty for all atoms.

You can figure out where they were added because the function will only
add one set of coordinates, and the value stored in "nalt" tells you where.
Alternately, you can set nalt=0 and free(xa) between each call.

If the file is a velocity trajectory, substitute velocity locations for
alternate coordinate locations in the last two paragraphs.

If the file is specified as a restart file, and if velocities are found, 
velocities will be added into array locations to match the xa locations.  
Any unused velocity locations will be cleared to equal zero.


*** The input file pointer status ***

The calling function *must* open the file.  In most cases, the file pointer 
should not be altered between calls to this function.  This is because the
function will begin reading coordinates from the current position pointed to
by the file pointer.  At the end of each read, it will be pointing to the
next set of coordinates.  Unless the calling program wants to change this
behavior, it should not alter the file pointer.


*** Regarding BOX information ***

This assembly *MUST* have nBOX set correctly.  If nBOX is less than or equal
to zero, this function will assume there is no box information contained in 
the coordinate file.  If nBOX is greater than or equal to one, it will assume
there is box information.  This distinction can have a significant impact on
the sanity of structures defined by the coordinate sets.  

LOCATION OF BOX INFO:  It is assumed most folks will ignore the box info, so
it will be put in a slightly annoying location.  In the assembly's nBOX-long
BOX array, the box info will be offset in the BOX array by an index of +1 from 
the locations of the coordinate and velocity information.  The zeroth and first
locations will be the box info originally contained in the TOP file.

To reiterate:  BOX INFO IS CONTAINED STARTING AT INDEX +1

...but, the three coordinates are in the first coordinate slot, unlike the 
prmtop info, which also contains an angle... confused yet?


*** An example for coordinate and box information location ***

For example, if there are 5 sets of alternative coordinates saved in the 
incoming assembly, and this function is reading in a coordinate trajectory
containing box info, the function will save new coordinates into the sixth
coordinate location (.xa[5]), and box info into the seventh (.BOX[6]).

Begun on 20080811 by BLFoley.  Significantly altered starting 20090803 BLF.
*/
#include <mylib.h>
#include <molecules.h>
#include <AMBER/amber.h>

void add_trajcrds_to_prmtop_assembly(
	fileset F, ///< The trajectory file containing the data *Must* be open
	assembly *A, ///< Assembly ("incoming") to which the crds should be added
	char ftype, ///< 'c' if coordinate ; 'v' if velocity ; 'r' if restart (might have both)
	int offset ///< read in the offset-th trajectory -- starts with zero=current
	){
char tmpc,readnum[13]; ///< temporary char holders
int allocated_xa=0, ///< alternate coords or vectors allocated outgoing
    nalt_here=0, ///< alternate coords or vectors allocated incoming
    scan_tst=0, ///< for testing [f,s]scanf results
    ai=0, ///< atom index
    xi=0, ///< dummy coordinate index
    newline=0; ///< counter to test if there should be a newline
long foffset;
fpos_t here,here2;

for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';

/// DO NOT RE/OPEN the file 
// But do check to see that it is already open
if(F.F==NULL){mywhine("add_trajcrds_to_prmtop_assembly: called file is not already open.");}
// check other sanity
if(A[0].na<=0){mywhine("add_trajcrds_to_prmtop_assembly: Incoming assembly contains no atom pointers.");}
if((ftype!='c')&&(ftype!='v')&&(ftype!='r')){mywhine("add_trajcrds_to_prmtop_assembly: unknown input file type.");}
if(offset<0){mywhine("add_trajcrds_to_prmtop_assembly: Negative coordinate offset/skip values are not supported.");}
if((ftype=='r')&&(offset>0)){
	printf("add_trajcrds_to_prmtop_assembly: WARNING! asked for non-zero offset in a restart file.\n");
	printf("\t\tA restart file should contain precisely one data set.  Ignoring the offset request.\n");}
if((A[0].nBOX>0)&&(A[0].BOX[0].nC!=2)&&(A[0].BOX[0].C[1].nD!=3)){
	mywhine("add_trajcrds_to_prmtop_assembly: non-prmtop BOXINFO not supported.");}

/// See if we are at the start of the file or not
foffset=ftell(F.F);
if(foffset==0){ // we are at the start of the file
// Ignore any titles -- let the calling program read them if it wants to
	tmpc=fgetc(F.F); // read and ignore the first line
	while((tmpc!='\n')&&(tmpc!='\r')&&(tmpc!=EOF)){ 
		tmpc=fgetc(F.F); 
		//printf("%c",tmpc);
		}
	//printf("\n");
	if (ftype=='r'){ // if restart, also read and ignore the second line
		tmpc=fgetc(F.F);
		while((tmpc!='\n')&&(tmpc!='\r')&&(tmpc!=EOF)){ tmpc=fgetc(F.F);
		//printf("%c",tmpc);
	       	}
	//printf("\n");
		}
	if(tmpc==EOF){printf("add_trajcrds_to_prmtop_assembly: WARNING: EOF found at top of file.\n");}
	}

/// Don't make lots of silly checks -- just add coordinates to atom structures...

/// First, allocate space as needed, and clear any unused middle-space
if(ftype=='v'){
	for(ai=0;ai<A[0].na;ai++){ if(A[0].a[ai][0].nvec) allocated_xa=A[0].a[ai][0].nvec; }
	if(allocated_xa==0){
		for(ai=0;ai<A[0].na;ai++){
			A[0].a[ai][0].nvec=allocated_xa=1;
			A[0].a[ai][0].v=(vectormag_3D*)calloc(1,sizeof(vectormag_3D));
			}
		}
	else{	
		allocated_xa++;
		for(ai=0;ai<A[0].na;ai++){
			if(A[0].a[ai][0].nvec==0){
				A[0].a[ai][0].nvec=allocated_xa;
				A[0].a[ai][0].v=(vectormag_3D*)calloc(allocated_xa,sizeof(vectormag_3D));
				}
			else{
				nalt_here=A[0].a[ai][0].nvec;
				A[0].a[ai][0].nvec=allocated_xa;
				A[0].a[ai][0].v=(vectormag_3D*)realloc(A[0].a[ai][0].v,allocated_xa*sizeof(vectormag_3D));
				for(xi=nalt_here;xi<allocated_xa;xi++){
					A[0].a[ai][0].v[xi].i=A[0].a[ai][0].v[xi].j=0;
					A[0].a[ai][0].v[xi].k=A[0].a[ai][0].v[xi].d=0;
					}
				}
			}
		}
	}
else{
	for(ai=0;ai<A[0].na;ai++){ if(A[0].a[ai][0].nalt>allocated_xa) allocated_xa=A[0].a[ai][0].nalt; }
	if(allocated_xa==0){
		for(ai=0;ai<A[0].na;ai++){
			A[0].a[ai][0].nalt=allocated_xa=1;
			A[0].a[ai][0].xa=(coord_3D*)calloc(1,sizeof(coord_3D));
			}
		}
	else{	
		allocated_xa++;
		for(ai=0;ai<A[0].na;ai++){
			if(A[0].a[ai][0].nalt==0){
				A[0].a[ai][0].nalt=allocated_xa;
				A[0].a[ai][0].xa=(coord_3D*)calloc(allocated_xa,sizeof(coord_3D));
				}
			else{
				nalt_here=A[0].a[ai][0].nalt;
				A[0].a[ai][0].nalt=allocated_xa;
				A[0].a[ai][0].xa=(coord_3D*)realloc(A[0].a[ai][0].xa,allocated_xa*sizeof(coord_3D));
				for(xi=nalt_here;xi<allocated_xa;xi++){
					A[0].a[ai][0].xa[xi].i=A[0].a[ai][0].xa[xi].j=A[0].a[ai][0].xa[xi].k=0;
					}
				}
			}
		}
	}

/// Now scan through and add info to the structure
xi=allocated_xa-1;
if(ftype=='r') newline=1;
else newline=0;
for(ai=0;ai<A[0].na;ai++){
	//if(ftype=='v') scan_tst=fscanf(F.F,"%lf",&A[0].a[ai][0].v[xi].i);
	//else scan_tst=fscanf(F.F,"%lf",&A[0].a[ai][0].xa[xi].i);
	//if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass 1-i.");}
	//if(ftype=='v') scan_tst=fscanf(F.F,"%lf",&A[0].a[ai][0].v[xi].j);
	//else scan_tst=fscanf(F.F,"%lf",&A[0].a[ai][0].xa[xi].j);
	//if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass 1-j.");}
	//if(ftype=='v') scan_tst=fscanf(F.F,"%lf",&A[0].a[ai][0].v[xi].k);
	//else scan_tst=fscanf(F.F,"%lf",&A[0].a[ai][0].xa[xi].k);
	
	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	if(ftype=='r') scan_tst=fscanf(F.F,"%12c",readnum);
	else {
		scan_tst=fscanf(F.F,"%8c",readnum);
		newline++;
		}
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass 1-i.");}

//printf("readnum (i) is >>>%s<<< newline is %d\n",readnum,newline);
	if( ((ftype!='r')&&(newline==10)) ) 
		{
		scan_tst=fgetc(F.F);
//printf("  newline is %d (10, we hope), scan_tst is >>>%c<<<\n",newline, scan_tst);
		if(scan_tst!='\n') {mywhine("i: newline expected but not found during first read-set");}
		if(ftype!='r') newline=0;	
		} 

	if(ftype=='v') scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].v[xi].i);
	else scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].xa[xi].i);

//printf("A[0].a[ai][0].xa[xi].i is %f\n",A[0].a[ai][0].xa[xi].i);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	if(ftype=='r') scan_tst=fscanf(F.F,"%12c",readnum);
	else {
		scan_tst=fscanf(F.F,"%8c",readnum);
		newline++;
		}
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass 1-j.");}

//printf("readnum (j) is >>>%s<<< newline is %d\n",readnum,newline);
	if( ((ftype!='r')&&(newline==10)) ) 
		{
		scan_tst=fgetc(F.F);
//printf("  newline is %d (10, we hope), scan_tst is >>>%c<<<\n",newline, scan_tst);
		if(scan_tst!='\n') {mywhine("j: newline expected but not found during first read-set");}
		if(ftype!='r') newline=0;	
		} 

	if(ftype=='v') scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].v[xi].j);
	else scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].xa[xi].j);

//printf("A[0].a[ai][0].xa[xi].j is %f\n",A[0].a[ai][0].xa[xi].j);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	if(ftype=='r') {
		scan_tst=fscanf(F.F,"%12c",readnum);
		newline*=-1;
		}
	else {
		scan_tst=fscanf(F.F,"%8c",readnum);
		newline++;
		}
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass 1-k.");}

//printf("readnum (k) is >>>%s<<< newline is %d\n",readnum,newline);

	if(ftype=='v') scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].v[xi].k);
	else scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].xa[xi].k);

//printf("A[0].a[ai][0].xa[xi].k is %f\n",A[0].a[ai][0].xa[xi].k);

	if( ((ftype=='r')&&(newline==1)) || ((ftype!='r')&&(newline==10)) ) 
		{
		scan_tst=fgetc(F.F);
//printf("  newline is %d (10, we hope), scan_tst is >>>%c<<<\n",newline, scan_tst);
		if(scan_tst!='\n') {mywhine("1. newline not found during first read-set");}
		if(ftype!='r') newline=0;	
		} 

	if(ftype=='v'){
		A[0].a[ai][0].v[xi].d=sqrt(A[0].a[ai][0].v[xi].i*A[0].a[ai][0].v[xi].i+\
			A[0].a[ai][0].v[xi].j*A[0].a[ai][0].v[xi].j+\
			A[0].a[ai][0].v[xi].k*A[0].a[ai][0].v[xi].k);}
	}


if( ((ftype=='r')&&(newline!=1)) || ((ftype!='r')&&(newline!=0)) ) { 
	scan_tst=fgetc(F.F);
//printf("newline is %d, scan_tst is >>>%c<<<\n",newline, scan_tst);
	if(scan_tst!='\n') {mywhine("2. newline not found during first read-set");}
	}

if(fgetpos(F.F,&here)!=0){mywhine("problem getting file position -- before first BOX read");};
//printf("before first BOX read and xallocted_xa=%d\n",allocated_xa);
// Do BOX info, too, if appropriate:
if(A[0].nBOX>0){
	nalt_here=A[0].nBOX;
	A[0].nBOX=allocated_xa+1;
	A[0].BOX=(boxinfo*)realloc(A[0].BOX,(allocated_xa+1)*sizeof(boxinfo));
	for(xi=nalt_here;xi<(allocated_xa+1);xi++){
		A[0].BOX[xi].nC=1;
		A[0].BOX[xi].C=(coord_nD*)calloc(1,sizeof(coord_nD));
		A[0].BOX[xi].C[0].nD=3;
		A[0].BOX[xi].C[0].D=(double*)calloc(3,sizeof(double));
		}
	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	if(ftype=='r') scan_tst=fscanf(F.F,"%12c",readnum);
	else scan_tst=fscanf(F.F,"%8c",readnum); 
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass BOX 1-0.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].BOX[allocated_xa].C[0].D[0]);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	if(ftype=='r') scan_tst=fscanf(F.F,"%12c",readnum);
	else scan_tst=fscanf(F.F,"%8c",readnum); 
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass BOX 1-1.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].BOX[allocated_xa].C[0].D[1]);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	if(ftype=='r') scan_tst=fscanf(F.F,"%12c",readnum);
	else scan_tst=fscanf(F.F,"%8c",readnum); 
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass BOX 1-2.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].BOX[allocated_xa].C[0].D[2]);
	if(ftype=='r'){
		scan_tst=fscanf(F.F,"%12c",readnum);
		//printf("readnum is >>>%s<<<\n",readnum);
		scan_tst=fscanf(F.F,"%12c",readnum);
		//printf("readnum is >>>%s<<<\n",readnum);
		scan_tst=fscanf(F.F,"%12c",readnum);
		//printf("readnum is >>>%s<<<\n",readnum);
		}
	scan_tst=fgetc(F.F);
	if(scan_tst!='\n') {ungetc(scan_tst,F.F);}
	}
// START HERE -- one day there might be non-square box info in these files...
// as of now, though, there isn't...
if(ftype!='r'){return;}

/// If this is type 'r', see if there is any more info in the file
if(fgetpos(F.F,&here2)!=0){mywhine("problem getting file position -- before check for extra-r read");};
scan_tst=fscanf(F.F,"%12s",readnum);
//printf("readnum is >>>%s<<< scan_tst is %d\n",readnum,scan_tst);
if(scan_tst!=1){return;} // no more info in file

// Still here? Assume there are velocities in the file, but don't do any more fancy checks --
//	just assume the calling program has a clue...  :-)
if(fsetpos(F.F,&here)!=0){mywhine("problem setting file position -- before restart velocities read");};

// allocate velocity space 
if(allocated_xa==1){
	for(ai=0;ai<A[0].na;ai++){
		A[0].a[ai][0].nvec=allocated_xa;
		A[0].a[ai][0].v=(vectormag_3D*)calloc(1,sizeof(vectormag_3D));
		}
	}
else{	
	for(ai=0;ai<A[0].na;ai++){
		if(A[0].a[ai][0].nvec==0){
			A[0].a[ai][0].nvec=allocated_xa;
			A[0].a[ai][0].v=(vectormag_3D*)calloc(allocated_xa,sizeof(vectormag_3D));
			}
		else{
			nalt_here=A[0].a[ai][0].nvec;
			A[0].a[ai][0].nvec=allocated_xa;
			A[0].a[ai][0].v=(vectormag_3D*)realloc(A[0].a[ai][0].v,allocated_xa*sizeof(vectormag_3D));
			for(xi=nalt_here;xi<allocated_xa;xi++){
				A[0].a[ai][0].v[xi].i=A[0].a[ai][0].v[xi].j=0;
				A[0].a[ai][0].v[xi].k=A[0].a[ai][0].v[xi].d=0;
				}
			}
		}
	}

// read in velocities 
xi=allocated_xa-1;
newline=1;
for(ai=0;ai<A[0].na;ai++){
	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	scan_tst=fscanf(F.F,"%12c",readnum);
//printf("readnum, v-i is %12s\n",readnum);
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass r-v-i.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].v[xi].i);
//printf("the value saved is  %12.7f\n",A[0].a[ai][0].v[xi].i);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	scan_tst=fscanf(F.F,"%12c",readnum);
//printf("readnum, v-i is %12s\n",readnum);
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass r-v-i.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].v[xi].j);
//printf("the value saved is  %12.7f\n",A[0].a[ai][0].v[xi].j);
	
	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	scan_tst=fscanf(F.F,"%12c",readnum);
//printf("readnum, v-i is %12s\n",readnum);
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass r-v-i.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].a[ai][0].v[xi].k);
//printf("the value saved is  %12.7f\n",A[0].a[ai][0].v[xi].k);

	newline*=-1;
	if(newline==1)
		{
		scan_tst=fgetc(F.F);
		if(scan_tst!='\n')
			{
			printf("ai is %d and newline is %d\n",ai,newline);
			printf("the velocities are:\n");
			printf("%12.7f%12.7f%12.7f\n",A[0].a[ai][0].v[xi].i,A[0].a[ai][0].v[xi].j,A[0].a[ai][0].v[xi].k);
			mywhine("expected newline not found during velocity read in restart file");
			}
		}

	A[0].a[ai][0].v[xi].d=sqrt(A[0].a[ai][0].v[xi].i*A[0].a[ai][0].v[xi].i+\
		A[0].a[ai][0].v[xi].j*A[0].a[ai][0].v[xi].j+\
		A[0].a[ai][0].v[xi].k*A[0].a[ai][0].v[xi].k);
	}

if(newline==-1)
	{
	scan_tst=fgetc(F.F);
	if(scan_tst!='\n')
		{
		printf("ai is %d and newline is %d\n",ai,newline);
		printf("the velocities are:\n");
		printf("%12.7f%12.7f%12.7f\n",A[0].a[ai][0].v[xi].i,A[0].a[ai][0].v[xi].j,A[0].a[ai][0].v[xi].k);
		mywhine("expected newline not found during velocity read in restart file");
		}
	}
//for(ai=0;ai<A[0].na;ai++){
//printf("The coords are: %12.7f%12.7f%12.7f\n",A[0].a[ai][0].v[xi].i,A[0].a[ai][0].v[xi].j,A[0].a[ai][0].v[xi].k); }

// Do BOX info, too, if appropriate:
// Was already allocated -- just need to read it in
//printf("before rst box read and allocated_xa is %d\n",allocated_xa);
if(A[0].nBOX>0){
	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	scan_tst=fscanf(F.F,"%12c",readnum);
//printf("readnum, v-i is %12s\n",readnum);
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass BOX r-v-1.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].BOX[allocated_xa].C[0].D[0]);
//printf("the value saved is  %12.7f\n",A[0].BOX[allocated_xa].C[0].D[0]);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	scan_tst=fscanf(F.F,"%12c",readnum);
//printf("readnum, v-i is %12s\n",readnum);
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass BOX r-v-1.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].BOX[allocated_xa].C[0].D[1]);
//printf("the value saved is  %12.7f\n",A[0].BOX[allocated_xa].C[0].D[1]);

	for(scan_tst=0;scan_tst<13;scan_tst++) readnum[scan_tst]='\0';
	scan_tst=fscanf(F.F,"%12c",readnum);
//printf("readnum, v-i is %12s\n",readnum);
	if(scan_tst!=1){mywhine("add_trajcrds_to_prmtop_assembly: File read error, pass BOX r-v-1.");}
	scan_tst=sscanf(readnum,"%lf",&A[0].BOX[allocated_xa].C[0].D[2]);
//printf("the value saved is  %12.7f\n",A[0].BOX[allocated_xa].C[0].D[2]);
	}
//printf("After rst box read and BOX is %12.7f%12.7f%12.7f\n",A[0].BOX[allocated_xa].C[0].D[0],A[0].BOX[allocated_xa].C[0].D[1],A[0].BOX[allocated_xa].C[0].D[2]);
//for(ai=0;ai<A[0].na;ai++){
//printf("The coords are now: %12.7f%12.7f%12.7f\n",A[0].a[ai][0].v[xi].i,A[0].a[ai][0].v[xi].j,A[0].a[ai][0].v[xi].k); }
return;
} 
