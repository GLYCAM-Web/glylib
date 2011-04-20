/** \file 0u.check_fd_status.c  Utilities for checking the existence, 
	size and/or other attributes of a file or directory. 
*/

#include <mylib.h>
#include <gly_fileutils.h>


/****************** check_file_existence() ******************/
/** Figures out whether a file exists, but that's all. 

Name is the name, including path if needed, of the file
to check. 
*/

///  START HERE
// uhm.... this won't work...

int check_file_existence(const char *Name){
//int status=0;
///return status;
fprintf(stderr,"ERROR!  The file check_file_existence doesn't work.  Rewrite it or don't use it.\n");
return -1;
}

/****************** check_directory_existence() ******************/
/** Figures out whether a directory exists, but that's all. 

Name is the name, including path if needed, of the directory 
to check. 
*/

///  START HERE
// uhm.... this won't work...

int check_directory_existence(const char *Name){
//int status=0;
//return status;
fprintf(stderr,"ERROR!  The file check_directory_existence doesn't work.  Rewrite it or don't use it.\n");
return -1;
}

/*************** gly_get_current_working_directory() **************/
/** Returns a pointer to a string corresponding to the current working 
directory.  

Typical usage:

char *myCWD;

myCWD=strdup(gly_get_current_working_directory());
*/
const char *gly_get_current_working_directory(void){
char *dirstr;
int sz=250;
dirstr=(char*)calloc(sz,sizeof(char));
while(getcwd(dirstr, 250)==NULL){
	if(sz>5000) {mywhine("Unbelievably long (> 5000 chars) current working directory string\n\tfound in gly_get_current_working_directory.");}
	sz+=50;
	dirstr=realloc(dirstr,sz*(sizeof(char))); } 
return dirstr;
}
