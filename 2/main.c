#include <math.h>
#include <stdio.h>

typedef struct v2 {
    double x;
    double y;
} v2;

double f(double x, double y) {
  return +1*x*x*x*x +
    +2*x*x*x*y +
    +3*x*x*y*y +
    +2*x*y*y*y +
    +2*y*y*y*y +
    -2.6*x*x*x +
    -3.2*x*x*y +
    -2.6*x*y*y +
    -3.2*y*y*y +
    -7.22*x*x +
    -16*x*y +
    -15.22*y*y +
    +20.8*x +
    +25.6*y +
    -5.94;
}

v2 grad_f(double x, double y) {
  return (v2) {

    +4*x*x*x +
    +6*x*x*y +
    +6*x*y*y +
    +2*y*y*y +
    -7.8*x*x +
    -6.4*x*y +
    -2.6*y*y +
    -14.44*x +
    -16*y +
    +20.8,

    +2*x*x*x +
    +6*x*x*y +
    +6*x*y*y +
    +8*y*y*y +
    -3.2*x*x +
    -5.2*x*y +
    -9.6*y*y +
    -16*x +
    -30.44*y +
    +25.6

  };
}
