/* 
\file slurp_directory.c 
\addtogroup FILE_UTILS
\brief  Read an entire directory into a structure to be processed
        later within the calling function.

        Argument: A directory path.

        Returns: Fileslurp structure containing each entry in the
		directory except for "." and "..".

	Date started: 24th February 2010.

\file glyopendir.c 
\addtogroup FILE_UTILS
\brief  Verifies if the path to a directory structure is valid,  and if the directory itself is a valid directory.

        Argument: A directory path.

        Returns: A pointer to the directory structure of interest. 

	Date started: 7th March 2010.


Author: Anita K. Nivedha
	minor modifications by BLFoley on 20100225

*/

#include<mylib.h>
//char glyopendir(char *pathToDir);
//fileslurp slurp_directory(char *pathToDir);

char glyopendir(char *pathToDir)
{
DIR *__restrict__ dir;
int errsv=0;
dir=opendir(pathToDir);
        if(dir==NULL){
                errsv = errno;
                switch(errsv) {
                        case EACCES:
                                printf("Search or Read Permission denied.\n");
                                exit(1);
                                break;
                        case ELOOP:
                                printf("A loop exists in symbolic links encountered during resolution of the dirname argument.\n");
                                exit(1);
                                break;
                        case ENAMETOOLONG:
                                printf("Length of directory name too long.\n");
                                exit(1);
                                break;
                        case ENOENT:
                                printf("Invalid directory name in input.\n");
                                exit(1);
                                break;
                        case ENOTDIR:
                                printf("A component of the input directory name or path does not exist.\n");
                                exit(1);
                                break;
                        case EMFILE:
                                printf("File descriptors are currently open in the calling process.");
                                exit(1);
                                break;
                        case ENFILE:
                                printf("Too many files are currently open in the system.\n")
                                exit(1);
                                break;
               	             }
			}
	else
		{
		return dir;
		}	
}


fileslurp slurp_directory(char *pathToDir)
{
DIR *__restrict__ dir;
struct dirent *__restrict__ entry;
fileslurp S;
int errsv=0;
dir=glyopendir(pathToDir);  
entry=readdir(dir);     //reading contents of the directory            
S.n=0;
S.L=(char**)calloc(1,sizeof(char*));
while(entry!=NULL)
	{
  	if((strcmp(entry[0].d_name,".")!=0)&&(strcmp(entry[0].d_name,"..")!=0))
		{
		S.n++;
		S.L=(char**)realloc(S.L,S.n*sizeof(char*));
		S.L[S.n-1]=strdup(entry[0].d_name);			
		}
	entry=readdir(dir);
	}
return S;
}

