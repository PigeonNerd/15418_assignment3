/* Copyright 2013 Harry Q Bovik (hbovik) */

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>
#include "./wsp.h"

/* The number of processors in use. */
extern size_t procs;

/* Our process id. */
extern int procId;

/* The number of cities in the graph. */
extern size_t ncities;

/* The adjacency matrix of the graph. adj[i][j] is the distance from i to j. */
extern int adj[MAX_N][MAX_N];
extern int shortestEdge;
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
                      solution_t *general_solution);

/* Approximates a solution to the wsp and loads it into solution. */
void approx_wsp_greedy(solution_t *solution);

void solve_wsp(solution_t *solution) {
  /*
   * Approximate with a greedy solution first so we start with a reasonably
   * tight bound.
OA   */
    approx_wsp_greedy(solution);
    /* The iterations of the for loop will be split up accross all threads. */
    omp_set_nested(1);
#pragma omp parallel for default(shared) schedule(dynamic, 1)
    for (size_t start_city = 1; start_city < ncities; start_city++) {

      unsigned char unvisited[MAX_N];
      size_t j;
      for (j = 0; j < ncities; j ++) {
        unvisited[j] = j;
      }
      unvisited[0] = start_city;
      unvisited[start_city] = 0;
      //printf("Thread %d got iteration %d ------1 start as [%d, %d, %d, %d]\n",omp_get_thread_num(), start_city, unvisited[0], unvisited[1], unvisited[2], unvisited[3]);
#pragma omp parallel for default(shared) firstprivate(unvisited) schedule(dynamic, 1)
      for(size_t i = 1; i < ncities; i ++){
	int tmp = unvisited[1];
	unvisited[1] = unvisited[i];
	unvisited[i] = tmp;
	//printf("               Thread %d got iterarion %d ------2 start as [%d, %d, %d, %d]\n", omp_get_thread_num(), i,unvisited[0], unvisited[1], unvisited[2], unvisited[3]);

	//#pragma omp parallel for firstprivate(unvisited) default(shared) schedule(dynamic, 1)

	for(size_t k = 2; k < ncities; k++) {
	  int tmp2 = unvisited[2];
	  unvisited[2] = unvisited[k];
	  unvisited[k] = tmp2;

	  //printf("                                   Thread %d got iterarion %d ------2 start as [%d, %d, %d, %d]\n", omp_get_thread_num(), k, unvisited[0], unvisited[1], unvisited[2], unvisited[3]);
	  //#pragma omp parallel for firstprivate(unvisited) default(shared) schedule(dynamic, 1)
	  /*
	    for(size_t p = 3; p < ncities; p++){
	    int tmp = unvisited[3];
	    unvisited[3] = unvisited[p];
	    unvisited[p] = tmp;
	  */
	  //	printf("Thread %d got iterarion %d ------ [%d, %d, %d, %d]\n", omp_get_thread_num(), p,unvisited[0], unvisited[1], unvisited[2], unvisited[3]);
	int current_dist = adj[start_city][unvisited[1]] + adj[unvisited[1]][unvisited[2]];
	if(current_dist + shortestEdge * (ncities - 3) < solution->distance){
	  solve_wsp_serial(unvisited[2], current_dist, unvisited,
			   &unvisited[3], ncities - 3, solution);
	}
	}
      }
    }
    // printf("Thread %d got iteration %lu with distance %d\n", omp_get_thread_num(), start_city, localSolution.distance);

    /*
      #pragma omp critical
      {
      // printf("Thread %d got iteration %lu with distance %d\n", omp_get_thread_num(), start_city, localSolution.distance);
      if ( localSolution.distance < solution->distance ) {
      solution->distance = localSolution.distance;
      memcpy(solution->path, localSolution.path, ncities);
      }
      }
    */
}
