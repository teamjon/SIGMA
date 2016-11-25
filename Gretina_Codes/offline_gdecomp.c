/*** offl_decomp.c: process events for decomposition offline using decompLib */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "gdecomp.h"
#include "wrap.h"

int quiet;               /* set to 1 to prevent output of diagnostic info */

int main(int argc, char *argv[]) {

  struct decomp_state *state;
  Event_Signal *e;
  struct crys_intpts *x;
  FILE *fin, *fou;
  int num;
  char *s;

  e = (Event_Signal*) Calloc(1, sizeof(Event_Signal));

  fin = fopen(argv[1], "r");
  assert(fin != 0);
  fou = fopen(argv[2], "w");
  assert(fou != 0);

  state = dl_decomp_init("/Users/Mario/d0128/xtalk_basis.dat", 1);
  assert(state != 0);

  while ( (num = fread(e, sizeof(Event_Signal), 1, fin)) == 1) {
    x = dl_decomp(state, e);
    if (x != 0) {
      s = dl_crys_intpts_2s(x);
      num = fwrite(s, strlen(s), 1, fou);
      assert(num == 1);
      free(x);
      free(s);
    }
  }
  printf("nonet = %d\n", state->err->nonet);
  printf("bad chisq = %d\n", state->err->badchisq);

  printf("# decomp = %d\n", state->cnt);
  return 0;
}
