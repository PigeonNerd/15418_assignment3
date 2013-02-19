/* Copyright 2013 15418 Staff */

#include <fcntl.h>
#include <getopt.h>
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // For STDIN_FILENO.

#include "./wsp.h"

/* The number of processors in use. */
int procs = 1;

/* Our process id. */
int procId = 1;

/* The number of cities in the graph. */
size_t ncities = 1;

/* The adjacency matrix of the graph. adj[i][j] is the distance from i to j. */
int adj[MAX_N][MAX_N];

/* The global shortest edge in the graph. Used for a quick pruning heuristic. */
int shortestEdge = INT_MAX;

// ============================================================================
// Timer: returns time in seconds
// ============================================================================

void usage(char *program) {
  printf("Usage: %s [options]\n"
         "  Solves the Wandering Salesman Problem for the distances provided \n"
         "  on stdin.\n"
         "\n"
         "Program Options:\n"
         "  -b  --bench        Run benchmarking\n"
         "  -i <file>          Use <file> instead of standard in\n"
         "  -?  --help         This message\n", program);
}

void parse_args(int argc, char **argv) {
  int opt;

  static struct option long_opts[] = {
    {"benchmark", 0, 0, 'b'},
    {"input", 1, 0, 'i'},
    {"help", 0, 0, '?'},
    {0, 0, 0, 0}
  };

  while ((opt = getopt_long(argc, argv, "i:?h", long_opts, NULL)) != EOF) {
    switch (opt) {
    case 'i':
      if (dup2(open(optarg, O_RDONLY), STDIN_FILENO) < 0) {
        perror(optarg);
        exit(EXIT_FAILURE);
      }
      break;
    case 'h':                  /* Explicit fall through */
    case '?':
      usage(argv[0]);
      exit(EXIT_SUCCESS);
    default:
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
}

/* Parses standard in for an adjaceny matrix, loading into the global array. */
void parse_adj_matrix() {
  scanf("%lu", &ncities);
  if (ncities > MAX_N) {
    printf("Invalid number of cities: saw %lu but max is %d\n", ncities, MAX_N);
    exit(EXIT_FAILURE);
  }

  size_t i, j;

  for (i = 1; i < ncities; i++) {
    for (j = 0; j < i; j++) {
      scanf("%d", &adj[i][j]);
      adj[j][i] = adj[i][j];
      /* Update shortest edge. */
      if (adj[i][j] < shortestEdge) {
        shortestEdge = adj[i][j];
      }
    }
  }
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  parse_args(argc, argv);
  parse_adj_matrix();

  /* Find out how many processes we are using and which one we are. */
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);

  solution_t solution;

  double startTime = MPI_Wtime();

  solve_wsp(&solution);
  double endTime = MPI_Wtime();

  /*
   * Normally, we would do this for procId == 0, but we want to test that our
   * solution gets propogated to all processes.
   */
  if (procId == procs - 1) {
    size_t i;

    printf("Solution found: (%d", solution.path[0]);
    for (i = 1; i < ncities; i++) {
      printf(" -> %d", solution.path[i]);
    }
    printf(") of distance %d\n", solution.distance);
    printf("Solution took %.4fs on %d processors\n",
           endTime - startTime, procs);
  }

  MPI_Finalize();
  return 0;
}
