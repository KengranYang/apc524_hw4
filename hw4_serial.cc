// hw4

#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
#define _USE_MATH_DEFINES // get value of pi
#include <cmath>
#include <fstream>

int main(int argc, char const *argv[]) {
  clock_t cpu_t;
  cpu_t = clock();

  // initialize the grid
  const int grid_size = atof(argv[1]) +2;
  // printf("grid_size: %d\n", grid_size);

  typedef vector<double> Row;
  typedef vector<Row> Matrix;

  Matrix T(grid_size,Row(grid_size));
  // printf("number of rows: %lu\n", T.size());
  // printf("number of columns: %lu\n", T[0].size());

  // step size
  const double delta_x = M_PI/(grid_size-1);

  const double kappa = 1;
  const double total_time  = 0.5*M_PI*M_PI/kappa;
  const double time_step = delta_x*delta_x/kappa/5.0; // make sure it is smaller than delta_x*delta_x/kappa/4


  // add boundary conditions here
  for (int i = 0; i < grid_size; ++i) {
    T[0][i] = cos(i*delta_x)*cos(i*delta_x);
    T[grid_size-1][i] = sin(i*delta_x)*sin(i*delta_x);
    // printf("i = %d\n", i);
    // printf("%f\n", T[grid_size-1][i]);
  }
  // printf("0,0; pi,0: %f %f\n", T[0][0], T[0][grid_size-1]);
  // printf("0,pi; pi,pi: %f %f\n", T[grid_size-1][0], T[grid_size-1][grid_size-1]);


  // add time loop
  int t=0;
  while (t*time_step < total_time) {
    // printf("t = %d\n", t);

    // make it loop over the grid
    for (int i = 1; i < grid_size-1; ++i) {

      for (int j = 1; j < grid_size-1; ++j) {
        if (j==1) {
          T[i][0] = T[i][grid_size-2];
        }else if (j == grid_size-2) {
          T[i][grid_size-1] = T[i][1];
        }
        T[i][j] += time_step*kappa/delta_x/delta_x*(T[i-1][j]+T[i+1][j]+T[i][j-1]+T[i][j+1]-4*T[i][j]);
        // printf("i=%f, j=%f\n",i,j );
        // printf("%f\n", T[i][j]);
      }
    }
    ++t;
  }//end while loop

  // output data
  ofstream myfile;
  myfile.open ("example.txt");
  if (myfile.is_open()){
    myfile << "i j T\n";
    for (int i = 1; i < grid_size-1; ++i) {
      for (int j = 1; j < grid_size-1; ++j) {
        // myfile <<  i << " "<<j<<" "<<" "<<T[i][j]<<"\n";
        myfile <<T[i][j]<<" ";
        // printf("%f\n", T[i][j]);
      }
      myfile <<"\n";
    }
    myfile.close();
  }else{
    printf("unable to open file.\n");
  }

  // record the time used
  t = clock() - t;
  printf("CPU time: %d (%f seconds).\n",t , ((float)t)/CLOCKS_PER_SEC);


  return 0;
}
