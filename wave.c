/* -*- c-basic-offset:2; tab-width:2; indent-tabs-mode:nil -*-
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>

//#define WRITE_TO_FILE
//#define VERIFY

double timer();
double initialize(double x, double y, double t);
void save_solution(double *u, int Ny, int Nx, int n);

int main(int argc, char *argv[]) {
  
  // Start up MPI
  MPI_Init(&argc, &argv);
  
  // Find out number of processes
  int N_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &N_procs);
  
  // Define variables
  int Nx, Ny, Nt;
  double dt, dx, lambda_sq;
  double *u;
  double *u_old;
  double *u_new;
  double begin, end;
  int p1, p2;
  
  // Number of processes
  p1 = 1;
  if (argc > 1)
    p1 = atoi(argv[1]);
  p2 = 1;
  if (argc > 2)
    p2 = atoi(argv[2]);
  
  // Create the cartesian topology
  int ndims = 2;
  int dims[2] = {p1, p2};
  int periods[2] = {0, 0};
  int reorder = 1;
  int rank;
  int coords[2];
  MPI_Comm proc_grid, proc_x, proc_y;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &proc_grid);
  
  // Divide into rows and get ranks
  MPI_Comm_rank(proc_grid, &rank);
  MPI_Cart_coords(proc_grid, rank, ndims, coords);
  MPI_Comm_split(proc_grid, coords[0], coords[1], &proc_x);
  MPI_Comm_split(proc_grid, coords[1], coords[0], &proc_y);
  int px = coords[0];
  int py = coords[1];
  //MPI_Comm_rank(proc_x, &rank_x);
  //MPI_Comm_rank(proc_y, &rank_y);
  
  // Partition the finite differences grid
  int N_points = 128;
  if (argc > 3)
    N_points = atoi(argv[3]);
  int N_tasks1 = (double)N_points/p1;
  int N_tasks2 = (double)N_points/p2;
  int mod1 = N_points%p1;
  int mod2 = N_points%p2;
  Nx = N_tasks1 + 1;
  Ny = N_tasks2 + 1;
  if (px < mod1)
    Nx++;
  if (py < mod2)
    Ny++;
  if (px != 0 && px != (p1-1))
    Nx++;
  if (py != 0 && py != (p2-1))
    Ny++;
  
  // Create a block type
  MPI_Datatype block_type;
  MPI_Type_vector(Ny, 1, Ny, MPI_DOUBLE, &block_type);
  MPI_Type_commit(&block_type);
  
  // Time step and step sizes
  Nt = N_points;
  if (argc > 4)
    Nt = atoi(argv[4]);
  dx = 1.0/(Nx - 1);
  dt = 0.50*dx;
  lambda_sq = (dt/dx)*(dt/dx);
  
  // Allocate solution matrices
  u = malloc(Nx*Ny*sizeof(double));
  u_old = malloc(Nx*Ny*sizeof(double));
  u_new = malloc(Nx*Ny*sizeof(double));

  /* Setup IC */
  memset(u, 0, Nx*Ny*sizeof(double));
  memset(u_old, 0, Nx*Ny*sizeof(double));
  memset(u_new, 0, Nx*Ny*sizeof(double));
  for (int i = 0; i < Ny; ++i) {
    for (int j = 0; j < Nx; ++j) {
      double x, y;
      
      // x-value
      if (px == 0)
        x = j*dx;
      else if (px == (p1-1))
        x = j*dx + 1 - (Nx-1)*dx;
      else if (px < mod1)
        x = j*dx + px*(Nx-2)*dx - dx;
      else
        x = j*dx + mod1*(Nx-1)*dx + (px-mod1)*(Nx-2)*dx - dx;
      
      // y-value
      if (py == 0)
        y = i*dx;
      else if (py == (p2-1))
        y = i*dx + 1 - (Ny-1)*dx;
      else if (py < mod2)
        y = i*dx + py*(Ny-2)*dx - dx;
      else
        y = i*dx + mod2*(Ny-1)*dx + (py-mod2)*(Ny-2)*dx - dx;
      
      // Compute initial value
      u[i*Nx + j] = initialize(x, y, 0);
      u_new[i*Nx + j] = initialize(x, y, dt);
    }
  }
  
  /* Setup BC */
  if (px == 0) {
    for (int i = 0; i < Ny; i++) {
      u[i*Nx] = 0;
      u_new[i*Nx] = 0;
    }
  }
  if (py == 0) {
    for (int i = 0; i < Nx; i++) {
      u[i] = 0;
      u_new[i] = 0;
    }
  }
  if (px == (p1-1)) {
    for (int i = 0; i < Ny; i++) {
      u[i*Nx + (Nx-1)] = 0;
      u_new[i*Nx + (Nx-1)] = 0;
    }
  }
  if (py == (p2-1)) {
    for (int i = 0; i < Nx; i++) {
      u[(Ny-1)*Nx + i] = 0;
      u_new[(Ny-1)*Nx + i] = 0;
    }
  }
  
  #ifdef WRITE_TO_FILE
  save_solution(u_new, Ny, Nx, 1);
  #endif
  
  #ifdef VERIFY
  double max_error = 0.0;
  #endif

  /* Integrate */
  begin = timer();
  for (int n = 2; n < Nt; ++n) {
    
    /* Swap ptrs */
    double *tmp = u_old;
    u_old = u;
    u = u_new;
    u_new = tmp;
    
    // Send data in the x-direction
    if (p1 > 1) {
      if (px == 0) {
        MPI_Isend(&(u[Nx-1]), 1, block_type, 1, 0, proc_y, NULL);    // send to the right
      } else if (px == (p1-1)) {
        MPI_Isend(&(u[0]), 1, block_type, p1-2, 0, proc_y, NULL);    // send to the left
      } else {
        MPI_Isend(&(u[Nx-1]), 1, block_type, px+1, 0, proc_y, NULL); // send to the right
        MPI_Isend(&(u[0]), 1, block_type, px-1, 0, proc_y, NULL);    // send to the left
      }
    }
    
    // Send data in the y-direction
    if (p2 > 1) {
      if (py == 0) {
        MPI_Isend(&(u[(Ny-1)*Nx]), Nx, MPI_DOUBLE, 1, 0, proc_x, NULL);    // send upwards
      } else if (py == (p2-1)) {
        MPI_Isend(&(u[0]), Nx, MPI_DOUBLE, p2-2, 0, proc_x, NULL);         // send downwards
      } else {
        MPI_Isend(&(u[(Ny-1)*Nx]), Nx, MPI_DOUBLE, py+1, 0, proc_x, NULL); // send upwards
        MPI_Isend(&(u[0]), Nx, MPI_DOUBLE, py-1, 0, proc_x, NULL);         // send downwards
      }
    }
    
    // Receive data in the x-direction
    if (p1 > 1) {
      if (px == 0) {
        MPI_Recv(&(u[Nx-1]), 1, block_type, 1, 0, proc_y, NULL);    // receive from the right
      } else if (px == (p1-1)) {
        MPI_Recv(&(u[0]), 1, block_type, p1-2, 0, proc_y, NULL);    // receive from the left
      } else {
        MPI_Recv(&(u[Nx-1]), 1, block_type, px+1, 0, proc_y, NULL); // receive the right
        MPI_Recv(&(u[0]), 1, block_type, px-1, 0, proc_y, NULL);    // receive the left
      }
    }
    
    // Receive data in the y-direction
    if (p2 > 1) {
      if (py == 0) {
        MPI_Recv(&(u[(Ny-1)*Nx]), Nx, MPI_DOUBLE, 1, 0, proc_x, NULL);    // receive from above
      } else if (py == (p2-1)) {
        MPI_Recv(&(u[0]), Nx, MPI_DOUBLE, p2-2, 0, proc_x, NULL);         // receive from below
      } else {
        MPI_Recv(&(u[(Ny-1)*Nx]), Nx, MPI_DOUBLE, py+1, 0, proc_x, NULL); // receive from above
        MPI_Recv(&(u[0]), Nx, MPI_DOUBLE, py-1, 0, proc_x, NULL);         // receive from below
      }
    }
    
    /* Apply stencil */
    for (int i = 1; i < (Ny-1); ++i) {
      for (int j = 1; j < (Nx-1); ++j) {
        u_new[i*Nx + j] = 2*u[i*Nx + j] - u_old[i*Nx + j] + lambda_sq*
          (u[(i+1)*Nx + j] + u[(i-1)*Nx + j] + u[i*Nx + j+1] + u[i*Nx + j-1] - 4*u[i*Nx + j]);
      }
    }

    #ifdef VERIFY
    double error = 0.0;
    for (int i = 0; i < Ny; ++i) {
      for (int j = 0; j < Nx; ++j) {
        double e = fabs(u_new[i*Nx + j] - initialize(j*dx, i*dx, n*dt));
        if (e > error)
          error = e;
      }
    }
    if (error > max_error)
      max_error = error;
    #endif

    #ifdef WRITE_TO_FILE
    save_solution(u_new, Ny, Nx, n);
    #endif
  }
  end = timer();

  printf("Time elapsed: %g s\n", (end - begin));

  #ifdef VERIFY
  printf("Maximum error: %g\n", max_error);
  #endif

  free(u);
  free(u_old);
  free(u_new);
  
  // Shut down MPI
  MPI_Type_free(&block_type);
  MPI_Barrier(proc_grid);
  MPI_Finalize();
  
  return 0;
}

double timer() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

double initialize(double x, double y, double t) {
  double value = 0;
  
  /* standing wave */
  #ifdef VERIFY
  value = sin(3*M_PI*x)*sin(4*M_PI*y)*cos(5*M_PI*t);
  
  /* squared-cosine hump */
  #else
  const double width = 0.1;
  double centerx = 0.25;
  double centery = 0.5;
  double dist = sqrt((x-centerx)*(x-centerx) + (y-centery)*(y-centery));
  if (dist < width) {
    double cs = cos(M_PI_2*dist/width);
    value = cs*cs;
  }
  #endif
  
  return value;
}

void save_solution(double *u, int Ny, int Nx, int n) {
  char fname[50];
  sprintf(fname,"solution-%d.dat", n);
  FILE *fp = fopen(fname,"w");

  fprintf(fp, "%d %d\n", Nx, Ny);

  for (int j = 0; j < Ny; ++j) {
    for (int k = 0; k < Nx; ++k) {
      fprintf(fp, "%e\n", u[j*Nx + k]);
    }
  }

  fclose(fp);
}

