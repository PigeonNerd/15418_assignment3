#ifndef _WSP_H_
#define _WSP_H_

/* The maximum number of possible cities. */
#define MAX_N 20

typedef struct solution_s {
  unsigned char path[MAX_N];
  int distance;
} solution_t;

typedef struct message_s {
   int taskId;
   int distance;
} message_t;

void solve_wsp(solution_t* solution);

#endif
