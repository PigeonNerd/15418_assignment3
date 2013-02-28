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

    /* Traversing three levels of the tree before splitting into subtress. */
    int level = 3;
    size_t num_cities = ncities - 1;

    /* Total number of subtrees we are creating */
    size_t totalTasks = (ncities - 1) * (ncities - 1) * (ncities - 2) * (ncities - 3);
    unsigned char unvisited[MAX_N];
    unsigned char init[MAX_N];

    /* init is the intialized array that we will set
       unvisited to at the start of every new subtree. */
#pragma omp parallel for  schedule(dynamic)
    for (size_t i = 0; i < ncities; i++) {
      init[i] = i;
    }

    /* Single for loop that assigns threads to all subtrees. */
#pragma omp parallel for private(unvisited) schedule(dynamic, 1)
    for (size_t taskId = 0; taskId < totalTasks; taskId++) {

      /* Based on the current task ID, grab the indices of the parent
	 city, first child city, second child city, and third child
	 city. These indices will be the index of each city in the
	 initialized unvisited array. */
      size_t parentIndex = (taskId / ((num_cities)*(num_cities -1) * (num_cities -2))) + 1;
      size_t childIndex = (taskId  % num_cities) + 1;
      size_t thirdLevel = (taskId  % (num_cities - 1)) + 2;
      size_t fourthLevel = (taskId % (num_cities - 2)) + 3;

      /* Reset unvisited to the initialized values. */
      memcpy(unvisited, init, ncities);

      /* Swap the parent index with the starting city. */
      unvisited[0] = parentIndex;
      unvisited[parentIndex] = 0;

      /* Swap the first child index with unvisited[1],
	 or the first child of the parent in the path. */
      size_t tmpC = unvisited[childIndex];
      unvisited[childIndex] = unvisited[1];
      unvisited[1] = tmpC;

      /* Swap the second child index with unvisited[2],
	 or the second child of the parent in the path. */
      size_t tmpT = unvisited[thirdLevel];
      unvisited[thirdLevel] = unvisited[2];
      unvisited[2] = tmpT;

      /* Swap the third child index with unvisited[3],
	 or the third child fo the parent in the path. */
      size_t tmpF = unvisited[fourthLevel];
      unvisited[fourthLevel] = unvisited[3];
      unvisited[3] = tmpF;

      /* Calculate the current distance of the path from the
	 parent to the third child. */
      int curr_dist = adj[unvisited[0]][unvisited[1]] + adj[unvisited[1]][unvisited[2]] + adj[unvisited[2]][unvisited[3]];

      solution_t local_solution;
      local_solution.distance = solution->distance;

      /* If the distance of the path so far is less than the best general solution,
	 run the serial code on the remaining subtree. */
      if (( curr_dist + shortestEdge * (ncities-2)) < solution->distance) {
	solve_wsp_serial(unvisited[3], curr_dist, unvisited,
			 &unvisited[4], ncities - 4, solution, &local_solution);
      }
    }
}
