#include <stdio.h>
#include <stdlib.h>

#define NUM_BITS 9
#define SIZE (1<<NUM_BITS)

static void staticOverrun()
{
  double array[SIZE], array2[SIZE], array3[SIZE];
  int ii, jj;
  printf("Size = %d\n", SIZE);
  printf("Array memory used: %d\n", SIZE*8);
  for (ii = 0; ii < SIZE+10; ii++) {
    array[ii] = array2[ii] = array3[ii] = 0;
  }
  for (ii = 0; ii < SIZE+10; ii++) {
    array3[ii] = ii;
    array2[ii] = (ii) * 0.1;
  }

  for (ii = 0; ii < SIZE; ii++) {
    if (array[ii] != 0 || array3[ii] != ii) {
      printf("%d: %12f %12f %12f\n", ii, array[ii], array2[ii], array3[ii]);
    }
  }
}

static void dynamicOverrun()
{
  double *array, *array2, *array3;
  int ii, jj;

  array = (double*)calloc(SIZE, sizeof(double));
  array2 = (double*)calloc(SIZE, sizeof(double));
  array3 = (double*)calloc(SIZE, sizeof(double));

  printf("Size = %d\n", SIZE);
  printf("Array memory used: %d\n", SIZE*8);
  for (ii = 0; ii < SIZE; ii++) {
    array[ii] = array2[ii] = array3[ii] = 0;
  }
  for (ii = 0; ii < SIZE+10; ii++) {
    array3[ii] = ii;
    array2[ii] = (ii) * 0.1;
  }

  for (ii = 0; ii < SIZE; ii++) {
    if (array[ii] != 0 || array3[ii] != ii) {
      printf("%d: %12f %12f %12f\n", ii, array[ii], array2[ii], array3[ii]);
    }
  }
  free(array);
  free(array2);
  free(array3);
}

int main() {
  staticOverrun();
  dynamicOverrun();
  exit(0);
}