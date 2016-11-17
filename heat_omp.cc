// hw4

#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <vector>
using namespace std;
#define _USE_MATH_DEFINES // get value of pi
#include <cmath>
#include <fstream>
#include <assert.h>
#include <stdlib.h> /*atof*/
#include <chrono>

int main(int argc, char const *argv[]) {
  //record time elapsed
  chrono::time_point<chrono::system_clock> start, end;
  start = chrono::system_clock::now();

  //obtain user inputs
  if (argc!=3) {
    printf("Incorrect input! Please follow the example: ./heat_omp 128 8 \n" );
    assert(argc==3);
  }
  const int grid_size = atof(argv[1]) +2;
  const int nthr = atof(argv[2]) ;

  // initialize the grid
  typedef vector<double> Row;
  typedef vector<Row> Matrix;
  Matrix T(grid_size,Row(grid_size));

  // step size
  const double delta_x = M_PI/(grid_size-1);

  const double kappa = 1;
  const double total_time  = 0.5*M_PI*M_PI/kappa;
  const double time_step = delta_x*delta_x/kappa/4.0*0.99; // make sure it is smaller than delta_x*delta_x/kappa/4


  // add boundary conditions here
  for (int i = 0; i < grid_size; ++i) {
    T[0][i] = cos(i*delta_x)*cos(i*delta_x);
    T[grid_size-1][i] = sin(i*delta_x)*sin(i*delta_x);

  }

  // add time loop
  int t=0;
  while (t*time_step < total_time) {

    // add openMP directive here
    int part_rows= (grid_size-2)/nthr;
    int th_id;  //th_id holds the thread number for each thread

    omp_set_num_threads(nthr); //set the number of threads
    #pragma omp parallel shared(T, part_rows) private(th_id)
    {
      th_id = omp_get_thread_num();
      //Split the first for loop among the threads
      #pragma omp for schedule(guided,part_rows)
      // make it loop over the grid
      for (int i = 1; i < grid_size-1; ++i) {
        // printf("Thread #%d is doing row %d.\n",th_id,i);
        for (int j = 1; j < grid_size-1; ++j) {
          if (j==1) {
            T[i][0] = T[i][grid_size-2];
          }else if (j == grid_size-2) {
            T[i][grid_size-1] = T[i][1];
          }
          T[i][j] += time_step*kappa/delta_x/delta_x*(T[i-1][j]+T[i+1][j]+T[i][j-1]+T[i][j+1]-4*T[i][j]);

        }
      }
    }/* end of parallel section */

    ++t;
  }//end while loop

  // output data
  ofstream myfile;
  char file_name [100];
  sprintf (file_name, "heat_omp_%d_%d.out", grid_size-2, nthr);
  myfile.open (file_name);
  if (myfile.is_open()){
    myfile << "\n";
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
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("Time used: %fs.\n", elapsed_seconds.count());

  return 0;
}
