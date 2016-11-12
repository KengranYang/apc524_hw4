// hw4

#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char const *argv[]) {

  // initialize the grid
  const int grid_size = atof(argv[1]);

  printf("grid_size: %d\n", grid_size);
  typedef vector<double> Row;
  typedef vector<Row> Matrix;

  Matrix my_matrix(grid_size,Row(1));

  printf("number of column: %lu\n", my_matrix.size());
  printf("number of row: %lu\n", my_matrix[0].size());
  my_matrix[1][3] = 1;
  printf("%f\n", my_matrix[1][3]);

  
  return 0;
}
