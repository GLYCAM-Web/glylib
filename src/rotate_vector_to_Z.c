// Function written by B. Lachele Foley, 2007
#include <mylib.h>
#include <molecules.h>
//#include "../inc/mylib.h"
//#include "../inc/molecules.h"
/* This function rotates the coordinate system so that the
z-axis is aligned with the specified vector.  The x and y axes
are in the proper plane, but their exact position within that
plane is determined somewhat arbitrarily.   

There are two versions of the function.  Note that the meaning 
of the input integer is different in the two functions.  In each
function, the program will check for memory available and exit 
with complaint if it doesn't find enough.

The first works on a list of coord_3D coordinate sets. 

	void rotate_vector_to_Z_list(coord_3D* XS,int n,vectormag_3D v)

	*XS = pointer to list of coordinates to rotate
	  n = number of coordinates in list XS to rotate
	  v = the vector into which the z axis should be rotated.

Note that it will overwrite the coordinates you send to it.  So, if
you want to keep a copy of the originals, do that first.  The
int, here, is the number of coordinates pointed to by the 
coord_3D pointer.  

The second works on all the coordinate sets for a given molecule. 

	void rotate_vector_to_Z_M(molecule *m,int xs, int xl,int vs,int vl,vectormag_3D v)

	  m = a molecule structure (see molecules.h)
	 xs = the position of the initial coordinates in the molecule structure 
		if xs==-1 then the coords are in m.x.  If xs>=0, then
		the coords are in m.xa[xs].
	 xl = the position in the structure where you want the new
		coords to go; same rules as for xs
	 vs = the initial set of vectors to rotate.  If you don't need
		vectors rotated, set this to -1.
	 vl = the position in the array (m.v[vl]) where the rotated
		vectors should go.  Set to -1 for no vector rotation.
	  v = the vector into which the z axis should be rotated.

*/
/******************* rotate_vector_to_Z_list() *****************/
void rotate_vector_to_Z_list(coord_3D* XS,int n,vectormag_3D v){
int ra=0;
vectormag_3D cX,cY,cZ; // direction cosines for X Y and Z
coord_3D RC; // temporary coords

// If the requested direction already points along z
if((v.i==0)&&(v.j==0)){ return; }

// get direction cosines for unit vector along new Z axis (the rotation axis)
cZ=normalize_vec(v);
if(cZ.d==0){mywhine("zero-length vector passed to rotate_vector_to_Z (can't rotate to that)");}

if(cZ.j!=0){ // if the j component of the target direction isn't zero
	// set abitrary point in new XY plane for positive Y position
	cY.i=cZ.i; // set xY=iZ
	cY.k=cZ.k; // set zY=kZ
	cY.j=-(cZ.k*cZ.k+cZ.i*cZ.i)/cZ.j;  ; // solve for yY in X-Y plane
	// turn that into a unit vector (make into direction cosines)
	cY=normalize_vec(cY); 
	// get direction cosines for unit vector along new X axis
	// (this is the X-prod of the new Y and Z axes)
	cX=get_crossprod(cY,cZ); // note ordering looks backwards...
	}
else{ // the i component isn't zero (or we would have exited)
	// set abitrary point in new XY plane for positive X position
	cX.j=cZ.j; // set yX=jZ, which is zero or we wouldn't be here...
	cX.k=cZ.k; // set zX=kZ
	cX.i=-(cZ.k*cZ.k)/cZ.i;  ; // solve for yY in X-Y plane
	// turn that into a unit vector (make into direction cosines)
	cX=normalize_vec(cX); 
	// get direction cosines for unit vector along new X axis
	// (this is the X-prod of the new Y and Z axes)
	cY=get_crossprod(cZ,cX);
	}

// rotate molecule 
for(ra=0;ra<n;ra++){ // for each residue
	RC.i=cX.i*XS[ra].i + cX.j*XS[ra].j + cX.k*XS[ra].k;
	RC.j=cY.i*XS[ra].i + cY.j*XS[ra].j + cY.k*XS[ra].k;
	RC.k=cZ.i*XS[ra].i + cZ.j*XS[ra].j + cZ.k*XS[ra].k;
	XS[ra]=RC;
	}

return;
}

/******************* rotate_vector_to_Z_M() *****************/
void rotate_vector_to_Z_M(molecule *m,int xs, int xl,int vs,int vl,vectormag_3D v){
int ra=0,rb=0,localdebug=0;;
vectormag_3D cX,cY,cZ,RV; // direction cosines for X Y and Z
residue *R; // residue pointer for convenience
coord_3D RC;

// get direction cosines for unit vector along new Z axis (the rotation axis)
if(localdebug>=2){
printf("in rotate_vector_to_Z_M xs is %d ; xl is %d ; vs is %d ; vl is %d \n",xs,xl,vs,vl);
printf("normalizing rotation vector.  Before (v) \n");
dprint_vectormag_3D(&v);
}


// If the requested direction already points along z
if((v.i==0)&&(v.j==0)){ return; }

// get direction cosines for unit vector along new Z axis (the rotation axis)
cZ=normalize_vec(v);
if(localdebug>=2){
printf("normalizing rotation vector.  After (cZ)\n");
dprint_vectormag_3D(&cZ);
}
if(cZ.d==0){mywhine("zero-length vector passed to rotate_vector_to_Z (can't rotate to that)");}

if(cZ.j!=0){ // if the j component of the target direction isn't zero
	// set abitrary point in new XY plane for positive Y position
	cY.i=cZ.i; // set xY=iZ
	cY.k=cZ.k; // set zY=kZ
	cY.j=-(cZ.k*cZ.k+cZ.i*cZ.i)/cZ.j;  ; // solve for yY in X-Y plane
	// turn that into a unit vector (make into direction cosines)
if(localdebug>=2){
printf("normalizing Y axis vector. Before:\n");
dprint_vectormag_3D(&cY);
}
	cY=normalize_vec(cY); 
if(localdebug>=2){
printf("normalized Y axis vector:\n");
dprint_vectormag_3D(&cY);
}
	// get direction cosines for unit vector along new X axis
	// (this is the X-prod of the new Y and Z axes)
	cX=get_crossprod(cY,cZ); // note ordering looks backwards...
if(localdebug>=2){
printf("Just got X axis vector:\n");
dprint_vectormag_3D(&cX);
}
	}
else{ // the i component isn't zero (or we would have exited)
	// set abitrary point in new XY plane for positive X position
	cX.j=cZ.j; // set yX=jZ, which is zero or we wouldn't be here...
	cX.k=cZ.k; // set zX=kZ
	cX.i=-(cZ.k*cZ.k)/cZ.i;  ; // solve for yY in X-Y plane
	// turn that into a unit vector (make into direction cosines)
if(localdebug>=2){
printf("normalizing X axis vector (cZ.j=0).  Before:\n");
dprint_vectormag_3D(&cX);
}
	cX=normalize_vec(cX); 
if(localdebug>=2){
printf("normalized X axis vector:\n");
dprint_vectormag_3D(&cX);
}
	// get direction cosines for unit vector along new X axis
	// (this is the X-prod of the new Y and Z axes)
	cY=get_crossprod(cZ,cX);
if(localdebug>=2){
printf("Just got Y axis vector:\n");
dprint_vectormag_3D(&cY);
}
	}


//dprint_molecule(&m[0],24);
//dprint_residue(&m[0].r[0],24);
//printf("m[0].nr is %d\n",m[0].nr);

// rotate molecule 
if((xs==-1)&&(xl==-1)){ //this makes code harder to read, but execution faster
if(localdebug>=2){
printf("(xs==-1)&&(xl==-1)\n");
printf("m[0].nr is %d\n",m[0].nr);
}
for(ra=0;ra<m[0].nr;ra++){ // for each residue
	R=&m[0].r[ra]; // this isn't necessary, and certainly slows the execution 
		//by a few nanoseconds, bit it makes reading the code a bit easier
	for(rb=0;rb<R[0].na;rb++){ // for each atom
		RC.i=cX.i*R[0].a[rb].x.i + cX.j*R[0].a[rb].x.j + cX.k*R[0].a[rb].x.k;
		RC.j=cY.i*R[0].a[rb].x.i + cY.j*R[0].a[rb].x.j + cY.k*R[0].a[rb].x.k;
		RC.k=cZ.i*R[0].a[rb].x.i + cZ.j*R[0].a[rb].x.j + cZ.k*R[0].a[rb].x.k;
		R[0].a[rb].x=RC;
if(localdebug>=2){dprint_coord_3D(&R[0].a[rb].x); }
		} // for each atom
	} // for each residue
	}
if((xs!=-1)&&(xl==-1)){ //this makes code harder to read, but execution faster
for(ra=0;ra<m[0].nr;ra++){ // for each residue
	R=&m[0].r[ra]; // this isn't necessary, and certainly slows the execution 
		//by a few nanoseconds, bit it makes reading the code a bit easier
	for(rb=0;rb<R[0].na;rb++){ // for each atom
		RC.i=cX.i*R[0].a[rb].xa[xs].i + cX.j*R[0].a[rb].xa[xs].j + cX.k*R[0].a[rb].xa[xs].k;
		RC.j=cY.i*R[0].a[rb].xa[xs].i + cY.j*R[0].a[rb].xa[xs].j + cY.k*R[0].a[rb].xa[xs].k;
		RC.k=cZ.i*R[0].a[rb].xa[xs].i + cZ.j*R[0].a[rb].xa[xs].j + cZ.k*R[0].a[rb].xa[xs].k; 
		R[0].a[rb].x=RC;
if(localdebug>=2){dprint_coord_3D(&R[0].a[rb].x); }
		} // for each atom
	} // for each residue
	}
if((xs==-1)&&(xl!=-1)){ //this makes code harder to read, but execution faster
for(ra=0;ra<m[0].nr;ra++){ // for each residue
	R=&m[0].r[ra]; // this isn't necessary, and certainly slows the execution 
		//by a few nanoseconds, bit it makes reading the code a bit easier
	for(rb=0;rb<R[0].na;rb++){ // for each atom
		RC.i=cX.i*R[0].a[rb].x.i + cX.j*R[0].a[rb].x.j + cX.k*R[0].a[rb].x.k;
		RC.j=cY.i*R[0].a[rb].x.i + cY.j*R[0].a[rb].x.j + cY.k*R[0].a[rb].x.k;
		RC.k=cZ.i*R[0].a[rb].x.i + cZ.j*R[0].a[rb].x.j + cZ.k*R[0].a[rb].x.k; 
		R[0].a[rb].xa[xl]=RC;
if(localdebug>=2){dprint_coord_3D(&R[0].a[rb].x); }
		} // for each atom
	} // for each residue
	}
if((xs!=-1)&&(xl!=-1)){ //this makes code harder to read, but execution faster
for(ra=0;ra<m[0].nr;ra++){ // for each residue
	R=&m[0].r[ra]; // this isn't necessary, and certainly slows the execution 
		//by a few nanoseconds, bit it makes reading the code a bit easier
	for(rb=0;rb<R[0].na;rb++){ // for each atom

		RC.i=cX.i*R[0].a[rb].xa[xs].i + cX.j*R[0].a[rb].xa[xs].j + cX.k*R[0].a[rb].xa[xs].k;
		RC.j=cY.i*R[0].a[rb].xa[xs].i + cY.j*R[0].a[rb].xa[xs].j + cY.k*R[0].a[rb].xa[xs].k;
		RC.k=cZ.i*R[0].a[rb].xa[xs].i + cZ.j*R[0].a[rb].xa[xs].j + cZ.k*R[0].a[rb].xa[xs].k; 
		R[0].a[rb].xa[xl]=RC; 
if(localdebug>=2){dprint_coord_3D(&R[0].a[rb].x); }
		} // for each atom
	} // for each residue
	}


// if vector rotation is specified
if((vs>=0)&&(vl>=0)){
for(ra=0;ra<m[0].nr;ra++){ // for each residue
	R=&m[0].r[ra]; // this isn't necessary, and certainly slows the execution 
		//by a few nanoseconds, bit it makes reading the code a bit easier
	for(rb=0;rb<R[0].na;rb++){ // for each atom
		RV.i=cX.i*R[0].a[rb].v[vs].i + cX.j*R[0].a[rb].v[vs].j + cX.k*R[0].a[rb].v[vs].k;
		RV.j=cY.i*R[0].a[rb].v[vs].i + cY.j*R[0].a[rb].v[vs].j + cY.k*R[0].a[rb].v[vs].k;
		RV.k=cZ.i*R[0].a[rb].v[vs].i + cZ.j*R[0].a[rb].v[vs].j + cZ.k*R[0].a[rb].v[vs].k; 
		RV.d=get_magnitude(RV); // not that it should change...
		R[0].a[rb].v[vl]=RV; 
if(localdebug>=2){
	printf("Source (%d), intermediate and rotated (%d) vectors for set %d:\n",vs,vl,rb);
	dprint_vectormag_3D(&R[0].a[rb].v[vs]); 
	dprint_vectormag_3D(&RV); 
	dprint_vectormag_3D(&R[0].a[rb].v[vl]); 
	}
		} // for each atom
	} // for each residue
	}

return;
} 
