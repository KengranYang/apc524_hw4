// hw4

#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
#define _USE_MATH_DEFINES // get value of pi
#include <cmath>

// double step(Matrix mat){
//
// }

int main(int argc, char const *argv[]) {

  // initialize the grid
  const int grid_size = atof(argv[1]);
  printf("grid_size: %d\n", grid_size);

  typedef vector<double> Row;
  typedef vector<Row> Matrix;

  Matrix my_matrix(grid_size,Row(grid_size));
  printf("number of column: %lu\n", my_matrix.size());
  printf("number of row: %lu\n", my_matrix[0].size());

  // step size
  const double delta_x = M_PI/(grid_size-1);
  // const double delta_y = M_PI/grid_size;

  const double kappa = 1;
  const double total_time  = 0.5*M_PI*M_PI/kappa;
  const double time_step = total_time/1;//delta_x*delta_x/8/kappa; // make sure it is small enough

  // add boundary conditions here
  for (int i = 0; i < grid_size; ++i) {
    my_matrix[0][i] = cos(i*delta_x)*cos(i*delta_x);
    my_matrix[grid_size-1][i] = sin(i*delta_x)*sin(i*delta_x);
    // printf("i = %d\n", i);
    // printf("%f\n", my_matrix[grid_size-1][i]);
  }
  printf("0,0; pi,0: %f %f\n", my_matrix[0][0], my_matrix[0][grid_size-1]);
  printf("0,pi; pi,pi: %f %f\n", my_matrix[grid_size-1][0], my_matrix[grid_size-1][grid_size-1]);

  for (int i = 0; i < grid_size; ++i) {
    my_matrix[i][grid_size-1] = my_matrix[i][0];
    // printf("i = %d\n", i);
    // printf("%f\n", my_matrix[i][0]);
  }

  // add time loop
  int t=0;
  while (t*time_step < total_time) {
    printf("t = %d\n", t);
    // make it loop over the grid
    for (int i = 1; i < grid_size-1; ++i) {
      for (int j = 1; j < grid_size-1; ++j) {
        my_matrix[i][j] += time_step*kappa/delta_x/delta_x*(my_matrix[i-1][j]+my_matrix[i+1][j]+my_matrix[i][j-1]+my_matrix[i][j+1]-4*my_matrix[i][j]);
        // printf("i=%f, j=%f\n",i,j );
        // printf("%f\n", my_matrix[i][j]);
      }
    }
    ++t;
  }

  // output data
  
  // for (int i = 0; i < grid_size; ++i) {
  //   for (int j = 0; j < grid_size; ++j) {
  //     printf("i=%d, j=%d\n",i,j );
  //     printf("%f\n", my_matrix[i][j]);
  //   }
  // }

  return 0;
}
