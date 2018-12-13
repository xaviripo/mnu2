#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct v2 {
    double x;
    double y;
} v2;

double norm(v2 a) {
    return pow(pow(a.x, 2.) + pow(a.y, 2.), .5);
}

double distance(v2 a, v2 b) {
    return pow(pow(a.x - b.x, 2.) + pow(a.y - b.y, 2.), .5);
}

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

v2 get_initial_point(double tol) {

    // We use Newton's method to find an initial solution on y=0
    double x = 0.; // x_0
    double x_next = 0.;

    do {
        x = x_next;
        x_next = x - f(x, 0.)/grad_f(x, 0.).x;
    } while (fabs(f(x_next, 0.)) > tol);

    return (v2) {x_next, 0.};

}

v2 get_next_point(v2 point, v2 vect, double dis, double tol) {

    v2 point_next;
    point_next.x = point.x + vect.x;
    point_next.y = point.y + vect.y;

    v2 point_init;
    point_init.x = point.x;
    point_init.y = point.y;

    double dissqr = dis*dis;
    double
        f_xy,   // f(x,y)
        g_xy,   // g(x,y) = (x-x0)^2 + (y-y0)^2 - dis^2
        df_x,   // (df/dx)(x,y)
        df_y,   // (df/dy)(x,y)
        dg_x,   // (dg/dx)(x,y)
        dg_y,   // (dg/dy)(x,y)
        scalar; // -1/det(DF(x,y))
    v2 df_xy;  // (grad f)(x,y)

    do {

        printf("\tNext point: (%f, %f)\n", point_next.x, point_next.y);

        point = point_next;

        // F(x,y)
        f_xy = f(point.x, point.y);
        g_xy = pow(point.x - point_init.x, 2.) + pow(point.y - point_init.y, 2.) - dissqr;

        // (grad F)(x,y)
        df_xy = grad_f(point.x, point.y);
        df_x = df_xy.x;
        df_y = df_xy.y;
        dg_x = 2*(point.x - point_init.x);
        dg_y = 2*(point.y - point_init.y);
        scalar = -1./(df_x*dg_y - df_y*dg_x);
        printf("\tdf_x: %f\n", df_x);
        printf("\tdg_y: %f\n", dg_y);
        printf("\tdf_y: %f\n", df_y);
        printf("\tdg_x: %f\n", dg_x);
        printf("\tScalar: %f\n", scalar);

        point_next.x = point.x - scalar*(dg_y*f_xy - df_y*g_xy);
        point_next.y = point.y - scalar*(-dg_x*f_xy + df_x*g_xy);

    } while (f(point_next.x, point_next.y) > tol);

    return point_next;

}

void get_list_points() {

    FILE *f = fopen("output", "w");
    if (f == NULL) {
        printf("Error opening file\n");
        exit(1);
    }

    const double tol = 1e-10;
    const double dis = .01;

    const v2 point_init = get_initial_point(tol);

    v2 point = point_init;
    v2 point_next = point;

    v2 grad, vect, vect_next;

    fprintf(f, "%f %f\n", point_next.x, point_next.y);

    double vect_norm;

    int i = 0;
    do {
        printf("-------- ITERATION %d --------\n", ++i);

        // Update point and vector from previous iteration
        printf("Update point and vector from previous iteration\n");
        point = point_next;
        vect = vect_next;
        printf("\tPoint: (%f, %f)\n", point.x, point.y);
        printf("\tVector: (%f, %f)\n", vect.x, vect.y);

        // Calculate next initial approximation
        printf("Calculate next initial approximation\n");
        grad = grad_f(point.x, point.y);
        printf("\tGradient: (%f, %f)\n", grad.x, grad.y);
        vect_next = (v2) {-grad.y, grad.x};

        // Change it to length dis
        printf("Change it to length dis\n");
        printf("\tPre-normalized vector: (%f, %f)\n", vect_next.x, vect_next.y);
        vect_norm = norm(vect_next);
        vect_next.x = dis*vect_next.x/vect_norm;
        vect_next.y = dis*vect_next.y/vect_norm;

        // Make sure it's got the correct direction
        printf("Make sure it's got the correct direction\n");
        printf("\tOld direction: (%f, %f)\n", vect_next.x, vect_next.y);
        if (vect.x*vect_next.x + vect.y*vect_next.y < 0) {
            vect_next.x = -vect_next.x;
            vect_next.y = -vect_next.y;
        }
        printf("\tNew direction: (%f, %f)\n", vect_next.x, vect_next.y);

        // Get the next point
        printf("Get the next point\n");
        point_next = get_next_point(point, vect_next, dis, tol);
        printf("\tPoint next: (%f, %f)\n", point_next.x, point_next.y);

        fprintf(f, "%f %f\n", point_next.x, point_next.y);

    } while (distance(point_init, point_next) > dis);

    fclose(f);

}

int main() {
  printf("\x1B[36m");
  get_list_points();
  printf("\x1B[0m");
}
