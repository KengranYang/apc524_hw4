#include "mpi.h"
#include <cstdio> //printf
#include <stdlib.h> /*atof*/
#include <vector>
using namespace std;
#define _USE_MATH_DEFINES // get value of pi
#include <cmath>
#include <fstream>
#include <assert.h>
#include <string.h>
// # include <cstdlib>
// # include <iomanip>
// # include <iostream>

int main(int argc, char *argv[]) {

  double start_time = MPI_Wtime();//we're timing this run

  int id;
  int ierr;
  int p;
  int tag;

  MPI_File file;
  MPI_Status status;


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
  const int num_rows = (grid_size-2)/p+2;


  // initialize the grid
  typedef vector<double> Row;
  typedef vector<Row> Matrix;

  Matrix T((grid_size),Row(num_rows));

  // step size
  const double delta_x = M_PI/(grid_size-1);
  // printf("delta_x: %f\n", delta_x);

  const double kappa = 1;
  const double total_time  = 0.5*M_PI*M_PI/kappa;
  const double time_step =  delta_x*delta_x/kappa/4.0*0.99; // make sure it is smaller than delta_x*delta_x/kappa/4

  // add boundary conditions here
  for (int j = 0; j < p; j++) { //loop over the number of processors
    if (id == j) {

        for (int i = 0; i < num_rows; ++i) {
          const int k = i + j*(num_rows-2);
          T[0][i] = pow(cos(k*delta_x),2);
          T[grid_size-1][i] = pow(sin(k*delta_x),2);
          // printf("T[][i]: %f\n", T[grid_size-1][i] );
      }
    }
  }//end of for loop

  // add time loop
  int t=0;
  while (t*time_step < total_time) {

  // ------- send and receive boundaries -------------
    //send column(1) to column(num_rows-1) of processor id-1
  if (id > 0) {
    Row temp_send_1(grid_size);  //temp vector to store the boundary column for send and receive
    for (int i = 0; i < grid_size; i++) {

      temp_send_1[i] = T[i][1];

    }
    tag = 1;
    MPI_Send(&temp_send_1[0], grid_size, MPI_DOUBLE, id-1, tag, MPI_COMM_WORLD);

  }
  if (id < p-1) {
    Row col_right(grid_size);
    tag = 1;
    MPI_Recv ( &col_right[0], grid_size,  MPI_DOUBLE, id+1, tag, MPI_COMM_WORLD, &status );

    for (int i = 0; i < grid_size; i++) {
      T[i][num_rows-1] = col_right[i] ;

    }
  }


    //send column(num_rows-2) to column(0) of processor id+1
  if (id<p-1) {
    Row temp_send_num_rows_2(grid_size);
    for (int i = 0; i < grid_size; i++) {
      temp_send_num_rows_2[i] = T[i][num_rows-2];

    }
    tag = 2;
    MPI_Send(&temp_send_num_rows_2[0], grid_size, MPI_DOUBLE, id+1, tag, MPI_COMM_WORLD);
  }
  if (id > 0 ) {
    Row col_left(grid_size);
    tag = 2;
    MPI_Recv ( &col_left[0], grid_size,  MPI_DOUBLE, id-1, tag, MPI_COMM_WORLD, &status );

    for (int i = 0; i < grid_size; i++) {
      T[i][0] = col_left[i] ;

    }
  }



    //send colum(num_rows-2) of proc.(p-1) to colum(0) of proc. 0
    //send column(1) of proc. 0 to column(num_rows-1) of proc. p-1
    if (id==p-1) {
      Row temp_send(grid_size);
      Row col_right(grid_size);
      for (int i = 0; i < grid_size; i++) {
        temp_send[i] = T[i][num_rows-2];

      }
      tag = 3;
      MPI_Send(&temp_send[0], grid_size, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
      MPI_Recv ( &col_right[0], grid_size,  MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status );


      for (int i = 0; i < grid_size; i++) {
        T[i][num_rows-1] = col_right[i] ;

      }
    }
    if (id==0) {
      Row col_left(grid_size);
      Row temp_send(grid_size);
      for (int i = 0; i < grid_size; i++) {
        temp_send[i] = T[i][1];

      }
      tag = 3;
      MPI_Recv ( &col_left[0], grid_size,  MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD, &status );
      MPI_Send(&temp_send[0], grid_size, MPI_DOUBLE, p-1, tag, MPI_COMM_WORLD);

      for (int i = 0; i < grid_size; i++) {
        T[i][0] = col_left[i] ;

      }
    }

    // make it loop over the grid
    for (int i = 1; i < grid_size-1; ++i) {

      for (int j = 1; j < num_rows-1; ++j) {

        T[i][j] += time_step*kappa/delta_x/delta_x*(T[i-1][j]+T[i+1][j]+T[i][j-1]+T[i][j+1]-4*T[i][j]);

        // printf("T[%d][%d]: %f  ",i,j, T[i][j]);
      }
      // printf("\n");
    }

    ++t;
  }//end of while loop

  //calculate average temperature
  double T_sum_processor=0;
  double T_sum_all = 0;
  for (int i = 1; i < grid_size-1; ++i) {

    for (int j = 1; j < num_rows-1; ++j) {
      T_sum_processor += T[i][j];

    }

  }
  MPI_Allreduce (&T_sum_processor, &T_sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (id==0) {
    printf("Average temperature: %f\n", T_sum_all/pow((grid_size-2),2));
  }


//
//   ierr = MPI_File_open(MPI_COMM_WORLD, "final", MPI_MODE_CREATE|MPI_MODE_WRONLY,
//                           MPI_INFO_NULL, &file);
//   if(ierr != MPI_SUCCESS){
//       printf ("Error starting MPI program. Terminating.\n");
//       MPI_Abort(MPI_COMM_WORLD, ierr);
//     }
//
// // Writing the result of each local array to the shared final file:
//   // specifying Offset for writing to file
//   MPI_Offset offset = sizeof(double)*id*num_rows;
//
//   ierr = MPI_File_set_view (file, offset, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);
//   if(ierr != MPI_SUCCESS){
//         cerr << "Problem setting process view" << endl;
//         return EXIT_FAILURE;
//     }
//
//   MPI_File_write_all(file, &T[0], grid_size, MPI_DOUBLE, &status);
//   MPI_File_close(&file);


  // output data
  ofstream myfile;
  char file_name [100];
  sprintf (file_name, "heat_mpi_%d_%d_%d.dat", grid_size-2, p, id);

  myfile.open (file_name);
  if (myfile.is_open()){
    myfile << "id: " <<id << "\n" ;
    for (int i = 1; i < grid_size-1; ++i) {
      for (int j = 1; j < num_rows-1; ++j) {
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

  //
  //  Terminate MPI.
  //
  MPI_Finalize ( );

  // Report time consumed by each processor
  printf("Time used: %fs.\n", MPI_Wtime()-start_time);

  return 0;
}//end of main
