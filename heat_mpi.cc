#include "mpi.h"
#include <cstdio> //printf
#include <stdlib.h> /*atof*/
// #include <iostream>
#include <vector>
using namespace std;
#define _USE_MATH_DEFINES // get value of pi
#include <cmath>
// #include <fstream>
#include <assert.h>
// # include <cstdlib>
// # include <ctime>
// # include <iomanip>
// # include <iostream>

int main(int argc, char *argv[]) {

  // clock_t cpu_t = clock();
  int id;
  int ierr;
  int p;
  // double wtime;



//
//  Initialize MPI.
//
  ierr = MPI_Init(&argc,&argv);
  if (ierr != MPI_SUCCESS) {
    printf ("Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, ierr);
  }
//
//  Get the number of processes.
//
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &p );
//
//  Get the individual process ID.
//
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &id );
  // printf ("MPI init success.\n");
  cout << "  Process " << id << " says:\n";

  //obtain user inputs
  if (argc!=2) {

    printf("Incorrect input! Please follow the example: srun ./heat_mpi 128 \n" );
    assert(argc==2);
  }
  const int grid_size = atoi(argv[1]) +2;
  // printf("grid_size: %d\n", grid_size);


  // initialize the grid
  typedef vector<double> Row;
  typedef vector<Row> Matrix;

  const int num_rows = (grid_size-2)/p+2;
  Matrix T((grid_size),Row(num_rows));
  // printf("number of rows: %lu\n", T.size());
  // printf("number of columns: %lu\n", T[0].size());

  // step size
  const double delta_x = M_PI/(grid_size-1);
  printf("delta_x: %f\n", delta_x);

  const double kappa = 1;
  const double total_time  = 0.5*M_PI*M_PI/kappa;
  const double time_step = total_time; //delta_x*delta_x/kappa/4.0*0.99; // make sure it is smaller than delta_x*delta_x/kappa/4
  // printf("total_time: %f\n", total_time);

  // add boundary conditions here
  for (int j = 0; j < p; j++) { //loop over the number of processors
    if (id == j) {

        for (int i = 0; i < num_rows; ++i) {
          const int k = i + j*(num_rows-2);
          T[0][i] = pow(cos(k*delta_x),2);
          T[num_rows-1][i] = pow(sin(k*delta_x),2);
      }
    }
  }//end of for loop
  printf("0,0; pi,0: %f %f\n", T[0][0], T[0][num_rows-1]);
  printf("0,pi; pi,pi: %f %f\n", T[num_rows-1][0], T[num_rows-1][num_rows-1]);



  //
  //  Terminate MPI.
  //
  MPI_Finalize ( );
  // //
  // //  Terminate.
  // //
  // if ( id == 0 )
  // {
  //   cout << "\n";
  //   cout << "HELLO_MPI:\n";
  //   cout << "  Normal end of execution.\n";
  //   cout << "\n";
  // }
  return 0;
}//end of main
