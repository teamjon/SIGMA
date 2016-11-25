#ifndef _ASSYM_CRYSTAL_H
#define _ASSYM_CRYSTAL_H
/*Karin. Stuff specific to the assymetric Gretina crystals
 (included when SYMMETRIC_CRYSTAL is not defined)
*/

#define N_CRYSTAL_TYPES 2
enum CTYPE{CRYSTAL_A=0, CRYSTAL_B=1};//numbers have to be consecutive, start at 0
#define NCORNERS 12 /*because I assume hexagonal cylindrical symmetry*/
#define DEFAULT_CRYSTAL_TYPE CRYSTAL_A

/* set crystal type. Can be either CRYSTAL_A or CRYSTAL_B
   invalid type => nothing happens. Returns new value
*/
int set_crystal_geometry(int type);

int get_crystal_geometry(void);


#endif /*#ifndef _ASSYM_CRYSTAL_H*/

