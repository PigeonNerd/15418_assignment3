/* Copyright 2013 15418 Course Staff */
//This is the serial code
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "./wsp.h"
#include <omp.h>
/* The number of cities in the graph. */
extern size_t ncities;

/* The adjacency matrix of the graph. adj[i][j] is the distance from i to j. */
extern int adj[MAX_N][MAX_N];

/* The global shortest edge in the graph. Used for a quick pruning heuristic. */
extern int shortestEdge;

/*
 * Starting from the specified start city and unvisited nodes, computes the
 * remaining distance and modifies unvisited so that it is sorted in the order
 * we traverse them for the greedy solution.
 *
 * This doesn't need to be implemented recursively, but is more understandable
 * because it is.
 */
static int approx_wsp_greedy_unvisited(size_t start,
                                       unsigned char *unvisited,
                                       size_t num_unvisited) {
  if (num_unvisited == 0) {
    return 0;
  }

  /* Search for the closest city to the start city. */
  int local_min = INT_MAX;
  size_t local_min_i = 0;
  size_t i;

  for (i = 0; i < num_unvisited; i++) {
    size_t next_city = unvisited[i];

    if (adj[start][next_city] < local_min) {
      local_min_i = i;
      local_min = adj[start][next_city];
    }
  }

  /*
   * Move the closest city to the front of the list by swapping it with
   * whatever is already there.
   */
  unsigned char temp = unvisited[0];
  unvisited[0] = unvisited[local_min_i];
  unvisited[local_min_i] = temp;

  /* Continue on our journey with the remaining cities. */
  return local_min + approx_wsp_greedy_unvisited(unvisited[0],
                                                 &unvisited[1],
                                                 num_unvisited - 1);
}

void approx_wsp_greedy(solution_t *solution) {
  size_t i;

  for (i = 0; i < ncities; i++) {
    solution->path[i] = i;
  }

  solution->distance =
    approx_wsp_greedy_unvisited(0, &solution->path[1], ncities - 1);
}

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
                      solution_t *general_solution, solution_t* local_solution) {

  if (num_unvisited == 0) {

    if (current_distance < general_solution->distance) {
      #pragma omp critical
      {
	    general_solution->distance = current_distance;
	    memcpy(local_solution->path, current_route, ncities);
      }
    }
    return;
  }

  size_t i;
  int num_unvisited_i = static_cast < int >(num_unvisited);
  for (i = 0; i < num_unvisited; i++) {
    size_t next_city = unvisited[i];
    int next_dist = current_distance + adj[current_city][next_city];
    /*
     * Break early if there isn't a way for us to beat the best solution.
     *
     * In this case, we speculate that our total distance will be at least
     * the distance so far plus a hop of size at least the shortest distance
     * betweeen any pair of cities for each unvisited city.
     */
    if (next_dist + shortestEdge * (num_unvisited_i - 1) >=
        general_solution->distance) {
      continue;
    }
    /*
     * Mark next_city as unvisited by moving it to the front of the list and
     * advancing unvisited. This also sorts unvisited in the order we have
     * traversed the cities.
     */
    unvisited[i] = unvisited[0];
    unvisited[0] = next_city;
    solve_wsp_serial(next_city, next_dist, current_route,
                     &unvisited[1], num_unvisited - 1, general_solution, local_solution);
    /* Restore unvisited to its former state. */
    unvisited[0] = unvisited[i];
    unvisited[i] = next_city;
    }
}
