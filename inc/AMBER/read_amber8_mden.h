/* header file read_amber8_mden.h for reading in mden files
output by amber programs starting at around version 8.
Begun 20070910 by BLFoley as a different file
Changed to current form starting on 20080528
*/

#include <mylib.h>
#include <general.h>
#include <molecules.h>

double **read_amber8_mden(int NENT,char **ENT,int NDATA,fileset F);

