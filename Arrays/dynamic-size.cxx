#include <iostream>

#define NUM_BITS 12
#define SIZE (1<<NUM_BITS)

//int main() {
//  double **array;
//  int ii, jj;
//  std::cout << "Size = " << SIZE << std::endl;
//  std::cout << "Array memory used: " << (long)SIZE*SIZE*8 << std::endl;
//  array = new double*[SIZE];
//  for (ii = 0; ii < SIZE; ii++) {
//    array[ii] = new double[SIZE];
//    for (jj = 0; jj < SIZE; jj++) {
//      array[ii][jj] = (ii*SIZE) * 0.1;
//      /* array2[ii][jj] = (ii*SIZE + jj) * 0.1; */
//      /* array3[ii][jj] = (ii*SIZE + jj) * 0.1; */
//      /* array4[ii][jj] = (ii*SIZE + jj) * 0.1; */
//    }
//  }
//  /* Do some other stuff with the data. */
//  /* for (ii = 0; ii < SIZE; ii++) { */
//  /*   delete [] array[ii]; */
//  /* } */
//  delete [] array;
//  return 0;
//}