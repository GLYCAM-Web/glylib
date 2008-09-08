/* declarations.h begun 20071016 by BLFoley
 * for functions requiring multiple header files...
 */
#if !defined(GLYLIB_DECLARES)
#define GLYLIB_DECLARES
dockinfo *load_dlg_mol(fileset F,types *T);
void load_atypes(fileset FT, types *T);
void findAD3Energies(fileset F,dockinfo* D,int di);
void findAD4Energies(fileset F,dockinfo* D,int di);
#endif
