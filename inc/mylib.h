/* mylib.h
	-- generic header file that loads utilities 
		often needed by programmers and that
		are needed by the library */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#define _GNU_SOURCE
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
//#include <general.h>

FILE *myfopen(const char*,const char*);
FILE *myfreopen(const char*,const char*,FILE*);
void  myusage(const char *);
void  mywhine(const char *);
void  mywhineusage(const char *,const char *);
void  read_eek(const char *, const char *);
void read_neek(const char*,int,int); 
void read_fneek(const char *NEEK, int x, int y, const char *FILENAME);
int   crdsnaps(const char *, const char *, int);
const char * mytime(const char*);
int sscan_file(FILE*,const char*);
int cscan_file(FILE*,char);

