// hw4

#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
#define _USE_MATH_DEFINES // get value of pi
#include <cmath>
#include <fstream>
#include <assert.h>
#include <stdlib.h> /*atof*/

int main(int argc, char const *argv[]) {

  // clock_t cpu_t = clock();
  time_t start = time(0);

  //obtain user inputs
  if (argc!=2) {
    printf("Incorrect input! Please follow the example: ./heat_omp 128 \n" );
    assert(argc==2);
  }
  const int grid_size = atof(argv[1]) +2;

  // initialize the grid
  typedef vector<double> Row;
  typedef vector<Row> Matrix;

  Matrix T(grid_size,Row(grid_size));

  // step size
  const double delta_x = M_PI/(grid_size-1);

  const double kappa = 1;
  const double total_time  = 0.5*M_PI*M_PI/kappa;
  const double time_step =  delta_x*delta_x/kappa/4.0*0.99; // make sure it is smaller than delta_x*delta_x/kappa/4


  // add boundary conditions here
  for (int i = 0; i < grid_size; ++i) {
    T[0][i] = cos(i*delta_x)*cos(i*delta_x);
    T[grid_size-1][i] = sin(i*delta_x)*sin(i*delta_x);

  }

  // add time loop
  int t=0;
  while (t*time_step < total_time) {

    // make it loop over the grid
    for (int i = 1; i < grid_size-1; ++i) {

      for (int j = 1; j < grid_size-1; ++j) {
        if (j==1) {
          T[i][0] = T[i][grid_size-2];
        }else if (j == grid_size-2) {
          T[i][grid_size-1] = T[i][1];
        }
        T[i][j] += time_step*kappa/delta_x/delta_x*(T[i-1][j]+T[i+1][j]+T[i][j-1]+T[i][j+1]-4*T[i][j]);

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

        myfile <<T[i][j]<<" ";

      }
      myfile <<"\n";
    }
    myfile.close();
  }else{
    printf("unable to open file.\n");
  }

  // record the time used
  // cpu_t = clock() - cpu_t;
  time_t end = time(0);
  double time = difftime(end, start) * 1000.0;
  // printf("CPU time: %lu (%f seconds).\n",cpu_t , ((float)cpu_t)/CLOCKS_PER_SEC);
  printf("time: %f\n", time);

  return 0;
}
