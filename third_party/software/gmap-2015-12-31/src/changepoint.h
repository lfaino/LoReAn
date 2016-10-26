#ifndef CHANGEPOINT_INCLUDED
#define CHANGEPOINT_INCLUDED

extern int
Changepoint_left (int *nmatches_left, int *ntotal_left, int *matchscores, int length);
extern int
Changepoint_right (int *nmatches_right, int *ntotal_right, int *matchscores, int length);

#endif


