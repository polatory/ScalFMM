
#include "../src/Container/instrset.h"


int main() {
  int iset = instrset_detect();
  bool hasAVX = iset >= 7;
  printf("Has AVX? %s (iset = %d)\n", hasAVX ? "Yes!" : "No.", iset);
  return iset;
}
