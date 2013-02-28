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
                      solution_t *general_solution, solution_t* local_solution);

/* Approximates a solution to the wsp and loads it into solution. */
void approx_wsp_greedy(solution_t *solution);

void solve_wsp(solution_t *solution) {
  /*
   * Approximate with a greedy solution first so we start with a reasonably
   * tight bound.
OA */
    approx_wsp_greedy(solution);
    /* The iterations of the for loop will be split up accross all threads. */
    int level = 3;
    size_t num_cities = ncities - 1;
    size_t totalTasks = (ncities - 1) * (ncities - 1) * (ncities - 2) * (ncities - 3);
    unsigned char unvisited[MAX_N];
    unsigned char init[MAX_N];
#pragma omp parallel for  schedule(dynamic)
    for (size_t i = 0; i < ncities; i++) {
      init[i] = i;
    }

#pragma omp parallel for private(unvisited) schedule(dynamic, 1)
    for (size_t taskId = 0; taskId < totalTasks; taskId++) {

      //printf("Task ID = %zd ........\n", taskId);
      size_t parentIndex = (taskId / ((num_cities)*(num_cities -1) * (num_cities -2))) + 1;
      size_t childIndex = (taskId  % num_cities) + 1;
      size_t thirdLevel = (taskId  % (num_cities - 1)) + 2;
      size_t fourthLevel = (taskId % (num_cities - 2)) + 3;
      //printf("Parent Index is %zd, child index is %zd, third index is %zd \n", parentIndex, childIndex, thirdLevel);

      memcpy(unvisited, init, ncities);

      //size_t tmpP = unvisited[parentIndex];
      //unvisited[parentIndex] = unvisited[0];
      //unvisited[0] = tmpP;
      unvisited[0] = parentIndex;
      unvisited[parentIndex] = 0;

      size_t tmpC = unvisited[childIndex];
      unvisited[childIndex] = unvisited[1];
      unvisited[1] = tmpC;

      size_t tmpT = unvisited[thirdLevel];
      unvisited[thirdLevel] = unvisited[2];
      unvisited[2] = tmpT;

      size_t tmpF = unvisited[fourthLevel];
      unvisited[fourthLevel] = unvisited[3];
      unvisited[3] = tmpF;


      int curr_dist = adj[unvisited[0]][unvisited[1]] + adj[unvisited[1]][unvisited[2]] + adj[unvisited[2]][unvisited[3]];
      //printf("Task %d: [%d, %d, %d]\n", taskId, unvisited[0], unvisited[1], unvisited[2]);
      solution_t local_solution;
      local_solution.distance = solution->distance;
      if (( curr_dist + shortestEdge * (ncities-2)) < solution->distance) {
	solve_wsp_serial(unvisited[3], curr_dist, unvisited,
			 &unvisited[4], ncities - 4, solution, &local_solution);
/*#pragma omp critical
	    {
        if( local_solution.distance < solution->distance){
	      // printf("Thread %d got iteration %lu with distance %d\n", omp_get_thread_num(), start_city, localSolution.distance);
	      solution->distance = local_solution.distance;
	      memcpy(solution->path, local_solution.path, ncities);
        }
        }*/
      }
    }
}





    /*
#pragma omp parallel for default(shared) schedule(dynamic,1)
    for (size_t start_city = 1; start_city < ncities; start_city++) {

      unsigned char unvisited[MAX_N];
      size_t j;

#pragma omp parallel for private(j) schedule(dynamic)
      for (j = 0; j < ncities; j ++) {
	unvisited[j] = j;
      }

      unvisited[0] = start_city;
      unvisited[start_city] = 0;

      for(size_t i = 1; i < ncities; i ++){
	int tmp = unvisited[1];
	unvisited[1] = unvisited[i];
	unvisited[i] = tmp;

#pragma omp parallel for firstprivate(unvisited) default(shared) schedule(dynamic, 1)

	for(size_t k = 2; k < ncities; k++) {
	  int tmp2 = unvisited[2];
	  unvisited[2] = unvisited[k];
	  unvisited[k] = tmp2;

	  int current_dist = adj[start_city][unvisited[1]] + adj[unvisited[1]][unvisited[2]];
	  if(current_dist + shortestEdge * (ncities - 3) < (int)solution->distance){
	    solve_wsp_serial(unvisited[2], current_dist, unvisited,
			     &unvisited[3], ncities - 3, solution);
	  }
    */
	  /*
	  if (localSolution.distance < solution->distance) {
#pragma omp critical
	    {
	      // printf("Thread %d got iteration %lu with distance %d\n", omp_get_thread_num(), start_city, localSolution.distance);
	      solution->distance = localSolution.distance;
	      memcpy(solution->path, localSolution.path, ncities);
	    }
	  }
	  */
