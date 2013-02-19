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
                      solution_t *best_solution);

/* Approximates a solution to the wsp and loads it into solution. */
void approx_wsp_greedy(solution_t *solution);

void solve_wsp(solution_t *solution) {
  /*
   * Approximate with a greedy solution first so we start with a reasonably
   * tight bound.
   */
  approx_wsp_greedy(solution);
  solution_t localSolution;
  printf("The initial distance is: %d\n", localSolution.distance);
  /* The iterations of the for loop will be split up accross all threads. */
  #pragma omp parallel for private(localSolution)
  for (size_t start_city = 0; start_city < ncities; start_city++) {
    /* This block will be executed by at most one thread at a time. */
    unsigned char unvisited[MAX_N];
    localSolution.distance = solution->distance;
    memcpy(localSolution.path, solution->path, MAX_N);
    size_t j;
    for (j = 0; j < ncities; j ++) {
        unvisited[j] = j;
    }
    unvisited[0] = start_city;
    unvisited[start_city] = 0;
    printf("before call, best solution is: %d\n", localSolution.distance);
    solve_wsp_serial(start_city, 0, unvisited,
                     &unvisited[1], ncities - 1, &localSolution);
    #pragma omp critical
    {
      printf("Thread %d got iteration %lu with distance %d\n", omp_get_thread_num(), start_city, localSolution.distance);
        if ( localSolution.distance < solution->distance ) {
            solution->distance = localSolution.distance;
            memcpy(solution->path, localSolution.path, ncities);
        }
    }
  }
  // size_t start_city;
  /*
   * We never start at city 0 - any path starting with 0 will be equivalent to
   * another path passing through 0.
   */
  // for (start_city = 1; start_city < ncities; start_city++) {
  //  unsigned char unvisited[MAX_N];
  //  for (size_t i = 0; i < ncities; i++) {
  //    unvisited[i] = i;
  //  }
  //  unvisited[0] = start_city;
  //  unvisited[start_city] = 0;
  //  solve_wsp_serial(start_city, 0, unvisited,
  //                   &unvisited[1], ncities - 1, solution);
  // }
}
