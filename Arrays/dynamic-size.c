#include <stdio.h>
#include <stdlib.h>

#define NUM_BITS 12
#define SIZE (1<<NUM_BITS)
//int main() {
//  double **array;
//  int ii, jj;
//  printf("Size = %d\n", SIZE);
//  printf("Array memory used: %ld\n", (long)SIZE*(long)SIZE*(long)8);
//  array = (double**)calloc((size_t)SIZE, sizeof(double*));
//  for (ii = 0; ii < SIZE; ii++) {
//    array[ii] = (double*)calloc(SIZE, sizeof(double));
//    for (jj = 0; jj < SIZE; jj++) {
//      array[ii][jj] = (ii*SIZE + jj) * 0.1;
//      /* array2[ii][jj] = (ii*SIZE + jj) * 0.1; */
//      /* array3[ii][jj] = (ii*SIZE + jj) * 0.1; */
//      /* array4[ii][jj] = (ii*SIZE + jj) * 0.1; */
//    }
//  }
//  /* Do some other stuff with the data. */
//  /* for (ii = 0; ii < SIZE; ii++) { */
//  /*   free(array[ii]); */
//  /* } */
//  /* free(array); */
//  return 0;
//}