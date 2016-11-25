/* program to histogram chi-squared values from
   GRETINA / GRETA decomposition routine.

   Author:  D.C. Radford    Sept 2004
*/

#include <stdio.h>
#include <string.h>

int main(int argc, char **argv)
{
  float  energy, c2a, c2;
  int    i, j, k, n = 0, nhits, mat[4096];
  char   in_fn[256], line[256];
  FILE   *in_file = 0, *spn_file = 0;

  /* process command line options:
     fn = file name of input file
  */
  if (argc < 2) {
    printf("Usage: his_c2 input_file_name\n");
    return 1;
  }
  strncpy(in_fn, argv[1], 256);
  if (!(in_file=fopen(in_fn, "r"))) {
    printf("ERROR -- Cannot open file %s for reading.\n", in_fn);
    return 1;
  }
  if (strlen(in_fn) > 250) in_fn[250] = 0;
  strcat(in_fn, ".spn");
  if (!(spn_file=fopen(in_fn, "w+"))) {
    printf("ERROR -- Cannot open file %s for writing.\n", in_fn);
    return 1;
  }

  for (i=0; i<4096; i++) mat[i] = 0;

  while (fgets(line, 256, in_file) &&
	 6 == sscanf(line, "%d %d %f %d %f %f",
		     &nhits, &j, &energy, &k, &c2a, &c2)) {
    j = (100.0 * c2) + 0.5;
    if (j >= 0 && j < 4096) mat[j]++;
    for (i=0; i<nhits; i++) {
      fgets(line, 256, in_file);
    }
    n++;
  }
  printf("%d events processed.\n", n);

  fwrite(mat, 4096*4, 1, spn_file);
  return 0;
}
