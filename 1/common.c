#include <math.h>
#include <stdlib.h>

// Calculates the infinity norm-induced distance of n-sized vectors x and y
double distance(double *x, double *y, size_t n) {

  double diff = 0;
  double entry_diff = 0;

  for (int j = 0; j < n; j++) {
    entry_diff = fabs(y[j]-x[j]);
    diff = diff>entry_diff ? diff : entry_diff;
  }

  return diff;

}
