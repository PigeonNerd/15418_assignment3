/* Copyright 2013 Harry Q Bovik (hbovik) */

#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdlib.h>

#include "./wsp.h"

/* The number of processors in use. */
extern size_t procs;

/* Our process id. */
extern int procId;

/* The number of cities in the graph. */
extern size_t ncities;

/* The adjacency matrix of the graph. adj[i][j] is the distance from i to j. */
extern int adj[MAX_N][MAX_N];

/*
 * Solves WSP from the current city with the current distance traveled,
 * and the provided set of unvisited nodes.
 *
 * Will use best_solution to prune unnecessary work and will modify
 * best_solution so that it contains the best solution possible with the given
 * starting parameters.
 *
 * Note: this should be initially called with unvisited == current_route + 1.
 * There is a probably-too-clever trick here where it modifies unvisited so
 * that it leaves behind an array that is sorted in the order of traversal.
 *
 * The easiest way to use this is:
 *    unsigned char* current_route = new unsigned char[ncities];
 *    for (size_t i = 0; i < ncities; i++) { current_route[i] = i; }
 *    current_route[0] = start_city;
 *    current_route[start_city] = 0;
 *    solve_wsp_serial(start_city, 0, current_route,
 *                     &current_route[1], ncities-1,
 *                     best_solution);
 */
void solve_wsp_serial(size_t current_city, int current_distance,
                      unsigned char *current_route,
                      unsigned char *unvisited, size_t num_unvisited,
                      solution_t *best_solution, MPI_Datatype mpi_solution_type);

/* Approximates a solution to the wsp and loads it into solution. */
void approx_wsp_greedy(solution_t *solution);

MPI_Datatype mpi_solution_type = MPI_DATATYPE_NULL;

static void init_mpi_solution_type() {
  MPI_Datatype field_types[] = { MPI_UNSIGNED_CHAR, MPI_INT };
  int field_lengths[] = { MAX_N, 1 };
  MPI_Aint field_offsets[] = { 0, MAX_N };
  size_t num_fields = 2;

  int err = MPI_Type_create_struct(num_fields,
                                   field_lengths,
                                   field_offsets,
                                   field_types,
                                   &mpi_solution_type);

  assert(err == MPI_SUCCESS);
  err = MPI_Type_commit(&mpi_solution_type);
  assert(err == MPI_SUCCESS);

  int mpi_struct_size;

  assert(MPI_SUCCESS == MPI_Type_size(mpi_solution_type, &mpi_struct_size));
  assert(mpi_struct_size == sizeof(solution_t));
}

void run_master(){


}

void run_worker(){


}

void solve_wsp(solution_t *solution) {
  init_mpi_solution_type();
  /* Make sure all cores initialize the type before proceeding. */
  int err = MPI_Barrier(MPI_COMM_WORLD);

  assert(err == MPI_SUCCESS);

  //if (procId == 0) {
    /*
     * Approximate with a greedy solution first so we start with a reasonably
     * tight bound.
     */
    approx_wsp_greedy(solution);
    int u;
    for(u = 0; u < 12; u ++){
        printf("Thread %d got iteration %d\n", procId, u);
    }
    
    
    unsigned char unvisited[MAX_N];
    size_t start_city;
    /*
     * We never start at city 0 - any path starting with 0 will be equivalent to
     * another path passing through 0.
     */
    for (start_city = 1; start_city + procId < ncities; start_city += procs) {
      printf("Thread %d got iteration %d\n", procId,(int)start_city + procId);
        for (size_t i = 0; i < ncities; i++) {
        unvisited[i] = i;
      }
      unvisited[0] = start_city + procId;
      unvisited[start_city] = 0;
      solve_wsp_serial(start_city, 0, unvisited,
                       &unvisited[1], ncities - 1, solution, mpi_solution_type);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /* Tell all the other processors the answer. */
   // size_t i;
   // for (i = 1; i < procs; i++) {
   //   err = MPI_Send(solution, 1, mpi_solution_type, i, 0, MPI_COMM_WORLD);
   //   assert(err == MPI_SUCCESS);
   // }
  //} else {
  //  MPI_Status status;
  //  err =
  //    MPI_Recv(solution, 1, mpi_solution_type, 0, 0, MPI_COMM_WORLD, &status);
  //  assert(err == MPI_SUCCESS);
    /* 
     * The status object contains information about the request, including the
     * sender. In this case, it really should be the master since nobody else
     * is sending anything.
     */
    //assert(status.MPI_SOURCE == 0);
  //}
}
