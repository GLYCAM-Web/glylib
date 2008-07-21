/* general.h -- 
	convenience structures used by various programs 
*/

// Header file written by Lachele Foley, summer 2007

// for convenient arrays of strings ("Text")
// get rid of these ASAP...  Can they go now?
typedef struct {
	char *T;
} llcharset;
typedef struct {
	char T[501];
} lcharset;
typedef struct {
	char T[251];
} mcharset;
typedef struct {
	char T[51];
} scharset;
typedef struct {
	char T[11];
} sscharset;

// for ease in passing sets of file info between functions
typedef struct {
	char *N; // the name of the file -- use strdup to allocate/copy
	FILE *F; // the file pointer
} fileset;
// similar to the perl notion of "slurping" a file
typedef struct {
	int n; // number of lines
	char **L; // each line
} fileslurp; 
