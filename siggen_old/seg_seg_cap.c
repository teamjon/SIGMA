#include <stdio.h>
#include <stdlib.h>

#include "signal_calc_util.h"
#include "detector_geometry.h"
#include "fields.h"

#define GEOMETRY_FILE "geometry_setup_ppc.dat"
#define FIELDS_DIR_FILE "field_setup_ppc.dat"

int seg_caps(void);


int main(void) {
  int nsegments;

  printf("Reading detector geometry from file: %s\n", GEOMETRY_FILE);
  if ((nsegments = geometry_init(GEOMETRY_FILE)) <= 0) {
    error("setup of detector geometry failed\n");
    return -1;
  }
  printf("Reading field configuration data from file: %s\n", FIELDS_DIR_FILE);
  if (field_setup(FIELDS_DIR_FILE, nsegments) != 0) {
    error("Field setup failed\n");
    return -1;
  }
  printf("Setup done\n");

  return seg_caps();
}
