/* Copyright 2013 Harry Q Bovik (hbovik) */

#include <assert.h>
#include <limits.h>
#include <mpi.h>
#include <stdlib.h>

#include "./wsp.h"

#define GET_TREE_TAG 1
#define PUT_BEST_SOLUTION_TAG 2
#define REPLY_TREE_TAG 3
#define DIE_TAG 4

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

void run_master(solution_t* best_solution){
  unsigned char unvisited[MAX_N];
  solution_t solution;
  MPI_Status status;
  int numWorkers = procs - 1;
  int start_city = 1;
  int secondLevel = 0;
  //printf("master is running\n");
  while(1){
    MPI_Recv (&solution, 1, mpi_solution_type, MPI_ANY_SOURCE, MPI_ANY_TAG,
	      MPI_COMM_WORLD, &status);
    switch (status.MPI_TAG) {
    case GET_TREE_TAG:
    //printf("master recieves request from %d\n", status.MPI_SOURCE);
     if(start_city < (int)ncities){
        for (size_t i = 0; i < ncities; i++) {
             unvisited[i] = i;
        }
        unvisited[0] = start_city;
        unvisited[start_city] = 0;
        if(secondLevel != start_city){
           int tmp = unvisited[1];
           unvisited[1] = unvisited[secondLevel];
           unvisited[secondLevel] = tmp;
        }
        MPI_Send (unvisited, MAX_N, MPI_CHARACTER, status.MPI_SOURCE,
		        REPLY_TREE_TAG, MPI_COMM_WORLD);
        secondLevel ++;
        if (secondLevel == (int)ncities){
             secondLevel = 0;
             start_city++;
           }
        }
     else{
        // here we run out of job
        // so we let this worker die
        //printf("worker %d done its job\n", status.MPI_SOURCE);
        numWorkers --;
        MPI_Send (unvisited, 1, MPI_INT, status.MPI_SOURCE,
		                DIE_TAG, MPI_COMM_WORLD);
        if( numWorkers == 0){
            size_t i;
            for(i = 1; i < procs ; i++){
                MPI_Send (best_solution, 1, mpi_solution_type, i,
		                DIE_TAG, MPI_COMM_WORLD);
            }
            return;
        }
     }
      break;
    case PUT_BEST_SOLUTION_TAG:
      //printf("master recieves update best solution\n");
      if (solution.distance < best_solution->distance) {
	        //update best solution, put some lock
            best_solution->distance = solution.distance;
            memcpy(best_solution->path, solution.path, MAX_N);
      }
      break;
    }
  }
}

void run_worker(solution_t* best_solution){
  unsigned char unvisited[MAX_N];
  MPI_Status status;
  int start_city;
  //printf("worker %d starts to run\n", procId);
  while (1) {
    //printf("worker %d sends GET_TREE_MESSAGE\n", procId);
    MPI_Send (best_solution, 1, mpi_solution_type, 0,
	      GET_TREE_TAG, MPI_COMM_WORLD);
    MPI_Recv(unvisited, MAX_N, MPI_CHARACTER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);   
    switch(status.MPI_TAG){
        case REPLY_TREE_TAG:
            //printf("worker %d receives work [%d, %d, %d, %d]\n", procId, unvisited[0],unvisited[1],unvisited[2],unvisited[4]);
            //unvisited[0] = start_city;
            //unvisited[start_city] = 0;
            solve_wsp_serial(unvisited[1], adj[unvisited[0]][unvisited[1]], unvisited,
                       &unvisited[2], ncities - 2, best_solution, mpi_solution_type);
            break;
        case DIE_TAG:
            //printf("worker %d recieves die and wait\n", procId);
            MPI_Recv(best_solution, 1, mpi_solution_type, 0, DIE_TAG, MPI_COMM_WORLD, &status);   
            //printf("Thread %d returns\n", procId);
            return;
    }
  }
}

void solve_wsp(solution_t *solution) {
  approx_wsp_greedy(solution);
  init_mpi_solution_type();
  /* Make sure all cores initialize the type before proceeding. */
  int err = MPI_Barrier(MPI_COMM_WORLD);
  assert(err == MPI_SUCCESS);

  if (procId == 0) {
    run_master(solution);
  }else{
    run_worker(solution);
  }
  //printf("THERAD %d finishes\n", procId);
}
