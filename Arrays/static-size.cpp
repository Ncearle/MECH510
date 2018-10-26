#include <stdio.h>

constexpr int NUM_BITS = 6;
constexpr int SIZE = (1 << NUM_BITS);
//int main() {
//  double array[SIZE][SIZE], array2[SIZE][SIZE], array3[SIZE][SIZE], array4[SIZE][SIZE];
//  int ii, jj;
//  printf("Size = %d\n", SIZE);
//  printf("Array memory used: %d\n", SIZE*SIZE*8);
//  for (ii = 0; ii < SIZE; ii++) {
//    for (jj = 0; jj < SIZE; jj++) {
//      array[ii][jj] = (ii*SIZE + jj) * 0.1;
//      /* array2[ii][jj] = (ii*SIZE + jj) * 0.1; */
//      /* array3[ii][jj] = (ii*SIZE + jj) * 0.1; */
//      /* array4[ii][jj] = (ii*SIZE + jj) * 0.1; */
//    }
//  }
//  return 0;
//}