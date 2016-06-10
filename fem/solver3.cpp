#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "consts.h"
#include "timer.h"
#include "utils.h"
#include "common.h"

int TYPE_EXACT = 5;

inline void print_3D(const char *filename, int ox_len, int oy_len,
                          double hx, double hy, int iter, double a, double c, double *data) {
    char name[150];
    sprintf(name, "%s_nx=%d_h=%5.3f_iter=%d_a=%g_c=%g.dat",
            filename, ox_len + 1, hx, iter, a, c);
    FILE *file = fopen(name, "w");
    fprintf(file, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'\nVARIABLES = 'x' 'y' %s\nZONE T='SubZone'",filename);
    fprintf(file, "\nI=%d J=%d K=%d ZONETYPE=Ordered", oy_len + 1, ox_len + 1, 1);
    fprintf(file, "\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE)");
    for (int i = 0; i < ox_len + 1; i++)
        for (int j = 0; j < oy_len + 1; j++)
            fprintf(file, "\n%-30.20g  %-30.20g %-30.20g", i * hx, j * hy,
                    data[(oy_len + 1) * i + j]);

    fclose(file);
}

inline void print_data_to_files(const char *filename,
                                double *u, double *u_xx, double *u_yy, int iter) {
    char name[50];
    print_3D(filename, OX_LEN, OY_LEN, HX, HY, iter, A, C, u);
    sprintf(name, "%s_xx", filename);
    print_3D(name, OX_LEN, OY_LEN, HX, HY, iter, A, C, u_xx);
    sprintf(name, "%s_yy", filename);
    print_3D(name, OX_LEN, OY_LEN, HX, HY, iter, A, C, u_yy);
}


inline static double analytical_slv(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, xy; A=0, B=1, C=0, D=1!
            return x * y;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return x * (1. - x) * y;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return x * (1. - y) * y;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return x * (1. - x);
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return y * (1. - y);
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return x * (1. - x) * y * (1. - y);
    }

}

inline static double analytical_slv_x(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, xy; A=0, B=1, C=0, D=1!
            return y;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return (1. - 2. * x) * y;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return (1. - y) * y;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 1. - 2. * x;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 0.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return (1. - 2. * x) * y * (1. - y);
    }

}

inline static double analytical_slv_y(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, xy; A=0, B=1, C=0, D=1!
            return x;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return x * (1. - x);
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return x * (1. - 2. * y);
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 0.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 1. - 2. * y;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return x * (1. - x) *  (1. - 2. * y);
    }

}

inline static double analytical_slv_xx(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return -2. * y;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return 0.;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return -2.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 0.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return -2. * y * (1. - y);
    }

}

inline static double analytical_slv_yy(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return 0.;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return -2. * x;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 0.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return -2.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return -2. * x * (1. - x);
    }

}

inline static double f_rp(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return 2. * y;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return 2. * x;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 2.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 2.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return 2. * (x * (1. - x) + y * (1. - y));
    }

}

inline static double f_rp_x(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return 0.;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return 2.;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 0.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 0.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return 2. * (1. - 2. * x);
    }

}

inline static double f_rp_y(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return 2.;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return 0.;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 0.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 0.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return 2. * (1. - 2. * y);
    }

}

inline static double f_rp_xx(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return 0.;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return 0.;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 0.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 0.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return -4.;
    }

}

inline static double f_rp_yy(double x, double y) {

    switch (TYPE_EXACT) {
        case 0:
            // BiLinear, A=0, B=1, C=0, D=1!
            return 0.;
        case 1:
            // from P3, xxy; A=0, B=1, C=0, D=1!
            return 0.;
        case 2:
            // from P3, xyy; A=0, B=1, C=0, D=1!
            return 0.;
        case 3:
            // from P2, xx; A=0, B=1, C=0, D=1!
            return 0.;
        case 4:
            // from P2, yy; A=0, B=1, C=0, D=1!
            return 0.;
        case 5:
            // BiQuadrat, A=0, B=1, C=0, D=1!
            return -4.;
    }
}

double *calc_error_u(double hx, double hy, double *solution) {
    double *err = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            err[i * OY_LEN_1 + j] = solution[i * OY_LEN_1 + j]
                                    - analytical_slv(A + hx * i, C + hy * j);
    return err;
}

double *calc_error_u_xx(double hx, double hy, double *solution) {
    double *err = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            err[i * OY_LEN_1 + j] = solution[i * OY_LEN_1 + j]
                                    - analytical_slv_xx(A + hx * i, C + hy * j);
    return err;
}

double *calc_error_u_yy(double hx, double hy, double *solution) {
    double *err = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            err[i * OY_LEN_1 + j] = solution[i * OY_LEN_1 + j]
                                    - analytical_slv_yy(A + hx * i, C + hy * j);
    return err;
}

inline static double b_phi0(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    return (1. - xx / HX) * (1. - yy / HY);
}

inline static double b_phi0_x(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    if ((x-xi)<0)
        return (1. - yy / HY) / HX;
    else return -(1. - yy / HY) / HX;
}

inline static double b_phi0_y(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    if ((y-yj)<0)
        return (1. - xx / HX) / HY;
    else return -(1. - xx / HX) / HY;
}

inline static double b_phi1(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    return HX * xx * (1. - xx / HX) * (xx / HX - 2.) * (1. - yy / HY) / 6.;
}

inline static double b_phi1_x(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    if ((x-xi)>0)
        return HX * (1. - yy / HY) * (6. * xx / HX - 3. * (xx / HX) * (xx / HX) - 2.) / 6.;
    else return - HX * (1. - yy / HY) * (6. * xx / HX - 3. * (xx / HX) * (xx / HX) - 2.) / 6.;
}

inline static double b_phi1_y(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);
    double H = HX / HY / 6.;

    if ((y-yj)>0)
        return - H * xx * (1. - xx / HX) * (xx / HX - 2.);
    else return H * xx * (1. - xx / HX) * (xx / HX - 2.);
}

inline static double b_phi2(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    return HY * yy * (1. - yy / HY) * (yy / HY - 2.) * (1. - xx / HX) / 6.;
}

inline static double b_phi2_x(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);
    double H = HY / HX / 6.;

    if ((x-xi)>0)
        return - H * yy * (1. - yy / HY) * (yy / HY - 2.);
    else return H * yy * (1. - yy / HY) * (yy / HY - 2.);
}

inline static double b_phi2_y(double x, double y, int ii, int jj) {

    double xi = A + ii * HX;
    double yj = C + jj * HY;
    double xx = fabs(x - xi);
    double yy = fabs(y - yj);

    if ((y-yj)>0)
        return HY * (1. - xx / HX) * (6. * yy / HY - 3. * (yy / HY) * (yy / HY) - 2.) / 6.;
    else return - HY * (1. - xx / HX) * (6. * yy / HY - 3. * (yy / HY) * (yy / HY) - 2.) / 6.;
}

inline static double u_h_x(double *u_h, double *u_h_xx, double *u_h_yy,
                           double x, double y, int ii, int jj) {

    double val = u_h[OY_LEN_1 * ii + jj] * b_phi0_x(x, y, ii, jj)
                 + u_h[OY_LEN_1 * (ii + 1) + jj] * b_phi0_x(x, y, ii + 1, jj)
                 + u_h[OY_LEN_1 * (ii + 1) + jj + 1] * b_phi0_x(x, y, ii + 1, jj + 1)
                 + u_h[OY_LEN_1 * ii + jj + 1] * b_phi0_x(x, y, ii, jj + 1)
                 + u_h_xx[OY_LEN_1 * ii + jj] * b_phi1_x(x, y, ii, jj)
                 + u_h_xx[OY_LEN_1 * (ii + 1) + jj] * b_phi1_x(x, y, ii + 1, jj)
                 + u_h_xx[OY_LEN_1 * (ii + 1) + jj + 1] * b_phi1_x(x, y, ii + 1, jj + 1)
                 + u_h_xx[OY_LEN_1 * ii + jj + 1] * b_phi1_x(x, y, ii, jj + 1)
                 + u_h_yy[OY_LEN_1 * ii + jj] * b_phi2_x(x, y, ii, jj)
                 + u_h_yy[OY_LEN_1 * (ii + 1) + jj] * b_phi2_x(x, y, ii + 1, jj)
                 + u_h_yy[OY_LEN_1 * (ii + 1) + jj + 1] * b_phi2_x(x, y, ii + 1, jj + 1)
                 + u_h_yy[OY_LEN_1 * ii + jj + 1] * b_phi2_x(x, y, ii, jj + 1);
    return val;
}

inline static double u_h_y(double *u_h, double *u_h_xx, double *u_h_yy,
                           double x, double y, int ii, int jj) {

    double val = u_h[OY_LEN_1 * ii + jj] * b_phi0_y(x, y, ii, jj)
                 + u_h[OY_LEN_1 * (ii + 1) + jj] * b_phi0_y(x, y, ii + 1, jj)
                 + u_h[OY_LEN_1 * (ii + 1) + jj + 1] * b_phi0_y(x, y, ii + 1, jj + 1)
                 + u_h[OY_LEN_1 * ii + jj + 1] * b_phi0_y(x, y, ii, jj + 1)
                 + u_h_xx[OY_LEN_1 * ii + jj] * b_phi1_y(x, y, ii, jj)
                 + u_h_xx[OY_LEN_1 * (ii + 1) + jj] * b_phi1_y(x, y, ii + 1, jj)
                 + u_h_xx[OY_LEN_1 * (ii + 1) + jj + 1] * b_phi1_y(x, y, ii + 1, jj + 1)
                 + u_h_xx[OY_LEN_1 * ii + jj + 1] * b_phi1_y(x, y, ii, jj + 1)
                 + u_h_yy[OY_LEN_1 * ii + jj] * b_phi2_y(x, y, ii, jj)
                 + u_h_yy[OY_LEN_1 * (ii + 1) + jj] * b_phi2_y(x, y, ii + 1, jj)
                 + u_h_yy[OY_LEN_1 * (ii + 1) + jj + 1] * b_phi2_y(x, y, ii + 1, jj + 1)
                 + u_h_yy[OY_LEN_1 * ii + jj + 1] * b_phi2_y(x, y, ii, jj + 1);
    return val;
}

inline static double Simpson_1D_x(double *u_h, double *u_h_xx, double *u_h_yy,
                           int i, int j, int ii, int jj) {

    double x1 = A + i * HX;
    double y1 = C + j * HY;
    double y2 = C + j * HY + HY / 2.;
    double y3 = C + (j + 1) * HY;

    double val = HX * (u_h_x(u_h, u_h_xx, u_h_yy, x1, y1, ii, jj) * b_phi1(x1, y1, i, j)
                 + 4. * u_h_x(u_h, u_h_xx, u_h_yy, x1, y2, ii, jj) * b_phi1(x1, y2, i, j)
                 + u_h_x(u_h, u_h_xx, u_h_yy, x1, y3, ii, jj) * b_phi1(x1, y3, i, j)) / 6;
    return val;
}

inline static double Simpson_1D_y(double *u_h, double *u_h_xx, double *u_h_yy,
                                  int i, int j, int ii, int jj) {

    double x1 = A + i * HX;
    double y1 = C + j * HY;
    double x2 = A + i * HX + HX / 2.;
    double x3 = A + (i + 1) * HX;

    double val = HY * (u_h_y(u_h, u_h_xx, u_h_yy, x1, y1, ii, jj) * b_phi2(x1, y1, i, j)
                 + 4. * u_h_y(u_h, u_h_xx, u_h_yy, x2, y1, ii, jj) * b_phi2(x2, y1, i, j)
                 + u_h_y(u_h, u_h_xx, u_h_yy, x3, y1, ii, jj) * b_phi2(x3, y1, i, j))/ 6.;
    return val;
}



// jj+1 -------------
//      |     |     |
//      |  2  |  1  |
//      |     |     |
// jj   ------X------
//      |     |     |
//      |  3  |  4  |
//      |     |     |
// jj-1 |------------               X: (i,j)
//     ii-1   ii    ii+1


inline static double ff_qrule_Simpson(int ii, int jj, int i, int j) {

    double x1 = A + ii * HX;
    double x2 = A + ii * HX + HX/2.;
    double x3 = A + (ii+1) * HX;
    double y1 = C + jj * HY;
    double y2 = C + jj * HY + HY/2.;
    double y3 = C + (jj + 1) * HY;
    double x = A + i * HX;
    double y = C + j * HY;


    double simpson = HX * HY * (  f_rp(x1, y1)*b_phi0(x1, y1, i, j)
                                + f_rp(x3, y1)*b_phi0(x3, y1, i, j)
                                + f_rp(x3, y3)*b_phi0(x3, y3, i, j)
                                + f_rp(x1, y3)*b_phi0(x1, y3, i, j)
                                + 4. * ( f_rp(x2, y1)*b_phi0(x2, y1, i, j)
                                       + f_rp(x3, y2)*b_phi0(x3, y2, i, j)
                                       + f_rp(x2, y3)*b_phi0(x2, y3, i, j)
                                       + f_rp(x1, y2)*b_phi0(x1, y2, i, j))
                                + 16.  * f_rp(x2, y2)*b_phi0(x2, y2, i, j)) / 36.;

    return simpson;
}

inline static double ff_xx_qrule_Simpson(int ii, int jj, int i, int j) {

    double x1 = A + ii * HX;
    double x2 = A + ii * HX + HX/2.;
    double x3 = A + (ii+1) * HX;
    double y1 = C + jj * HY;
    double y2 = C + jj * HY + HY/2.;
    double y3 = C + (jj + 1) * HY;
    double x = A + i * HX;
    double y = C + j * HY;


    double simpson = HX * HY * (  f_rp(x1, y1)*b_phi1(x1, y1, i, j)
                                  + f_rp(x3, y1)*b_phi1(x3, y1, i, j)
                                  + f_rp(x3, y3)*b_phi1(x3, y3, i, j)
                                  + f_rp(x1, y3)*b_phi1(x1, y3, i, j)
                                  + 4. * ( f_rp(x2, y1)*b_phi1(x2, y1, i, j)
                                           + f_rp(x3, y2)*b_phi1(x3, y2, i, j)
                                           + f_rp(x2, y3)*b_phi1(x2, y3, i, j)
                                           + f_rp(x1, y2)*b_phi1(x1, y2, i, j))
                                  + 16.  * f_rp(x2, y2)*b_phi1(x2, y2, i, j)) / 36.;

    return simpson;
}

inline static double ff_yy_qrule_Simpson(int ii, int jj, int i, int j) {

    double x1 = A + ii * HX;
    double x2 = A + ii * HX + HX/2.;
    double x3 = A + (ii+1) * HX;
    double y1 = C + jj * HY;
    double y2 = C + jj * HY + HY/2.;
    double y3 = C + (jj + 1) * HY;
    double x = A + i * HX;
    double y = C + j * HY;


    double simpson = HX * HY * (  f_rp(x1, y1)*b_phi2(x1, y1, i, j)
                                  + f_rp(x3, y1)*b_phi2(x3, y1, i, j)
                                  + f_rp(x3, y3)*b_phi2(x3, y3, i, j)
                                  + f_rp(x1, y3 )*b_phi2(x1, y3, i, j)
                                  + 4. * ( f_rp(x2, y1)*b_phi2(x2, y1, i, j)
                                           + f_rp(x3, y2)*b_phi2(x3, y2, i, j)
                                           + f_rp(x2, y3)*b_phi2(x2, y3, i, j)
                                           + f_rp(x1, y2)*b_phi2(x1, y2, i, j))
                                  + 16.  * f_rp(x2, y2)*b_phi2(x2, y2, i, j)) / 36.;

    return simpson;
}



void fill_stencils_u_uxx(double H_SQ, double *stnl_1, double *stnl_2, double *stnl_3) {

    // stansil u --> u for u_xx
    stnl_1[0] = 16. / 3.;
    stnl_1[1] = -8. / 3.;
    stnl_1[2] = -2. / 3.;
    stnl_1[3] = 4. / 3.;
    stnl_1[4] = -2. / 3.;
    stnl_1[5] = -8. / 3.;
    stnl_1[6] = -2. / 3.;
    stnl_1[7] = 4. / 3.;
    stnl_1[8] = -2. / 3.;

    // stansil u --> u_xx for u_xx
    stnl_2[0] = -2. * H_SQ / 9.;
    stnl_2[1] = H_SQ / 9.;
    stnl_2[2] = - H_SQ / 18.;
    stnl_2[3] = H_SQ / 9.;
    stnl_2[4] = - H_SQ / 18.;
    stnl_2[5] = H_SQ / 9.;
    stnl_2[6] = - H_SQ / 18.;
    stnl_2[7] = H_SQ / 9.;
    stnl_2[8] = - H_SQ / 18.;


    // stansil u --> u_yy for u_xx
    stnl_3[0] = -8. * H_SQ / 9.;
    stnl_3[1] = 4. * H_SQ / 9.;
    stnl_3[2] = - H_SQ / 18.;
    stnl_3[3] = H_SQ / 9.;
    stnl_3[4] = - H_SQ / 18.;
    stnl_3[5] = 4. * H_SQ / 9.;
    stnl_3[6] = - H_SQ / 18.;
    stnl_3[7] = H_SQ / 9.;
    stnl_3[8] = - H_SQ / 18.;

}

void fill_stencils_uxx_uxx(double H_SQ, double *stnl_1, double *stnl_2, double *stnl_3) {

    double H_SQ_SQ=H_SQ * H_SQ;


    // stansil u_xx --> u for u_xx
    stnl_1[0] = -2. * H_SQ / 9.;
    stnl_1[1] = H_SQ / 9.;
    stnl_1[2] = - H_SQ / 18.;
    stnl_1[3] = H_SQ / 9.;
    stnl_1[4] = - H_SQ / 18.;
    stnl_1[5] = H_SQ / 9.;
    stnl_1[6] = - H_SQ / 18.;
    stnl_1[7] = H_SQ / 9.;
    stnl_1[8] = - H_SQ / 18.;

    // stansil u_xx --> u_xx for u_xx
    stnl_2[0] = -11. * H_SQ_SQ / 27.;
    stnl_2[1] = 11. * H_SQ_SQ / 54.;
    stnl_2[2] = H_SQ_SQ / 27.;
    stnl_2[3] = 5. * H_SQ_SQ / 54.;
    stnl_2[4] = H_SQ_SQ / 27.;
    stnl_2[5] = 11. * H_SQ_SQ / 54.;
    stnl_2[6] = H_SQ_SQ / 27.;
    stnl_2[7] = 5. * H_SQ_SQ / 54.;
    stnl_2[8] = H_SQ_SQ / 27.;

    // stansil u_xx --> u_yy for u_xx
    stnl_3[0] = -2. * H_SQ_SQ / 27.;
    stnl_3[1] = H_SQ_SQ / 27.;
    stnl_3[2] = - H_SQ_SQ / 216.;
    stnl_3[3] = H_SQ_SQ / 108.;
    stnl_3[4] = - H_SQ_SQ / 216.;
    stnl_3[5] = H_SQ_SQ / 27.;
    stnl_3[6] = - H_SQ_SQ / 216.;
    stnl_3[7] = H_SQ_SQ / 108.;
    stnl_3[8] = - H_SQ_SQ / 216.;

}

void fill_stencils_uyy_uxx(double H_SQ, double *stnl_1, double *stnl_2, double *stnl_3) {

    double H_SQ_SQ=H_SQ * H_SQ;

    // stansil u_yy --> u for u_xx
    stnl_1[0] = -8. * H_SQ / 9.;
    stnl_1[1] = 4. * H_SQ / 9.;
    stnl_1[2] = - H_SQ / 18.;
    stnl_1[3] = H_SQ / 9.;
    stnl_1[4] = - H_SQ / 18.;
    stnl_1[5] = 4. * H_SQ / 9.;
    stnl_1[6] = - H_SQ / 18.;
    stnl_1[7] = H_SQ / 9.;
    stnl_1[8] = - H_SQ / 18.;

    // stansil u_yy --> u_xx for u_xx
    stnl_2[0] = -2. * H_SQ_SQ / 27.;
    stnl_2[1] = H_SQ_SQ / 27.;
    stnl_2[2] = - H_SQ_SQ / 216.;
    stnl_2[3] = H_SQ_SQ / 108.;
    stnl_2[4] = - H_SQ_SQ / 216.;
    stnl_2[5] = H_SQ_SQ / 27.;
    stnl_2[6] = - H_SQ_SQ / 216.;
    stnl_2[7] = H_SQ_SQ / 108.;
    stnl_2[8] = - H_SQ_SQ / 216.;

    // stansil u_yy --> u_yy for u_xx
    stnl_3[0] = -5. * H_SQ_SQ / 27.;
    stnl_3[1] = 5. * H_SQ_SQ / 54.;
    stnl_3[2] = H_SQ_SQ / 27.;
    stnl_3[3] = - 2. * H_SQ_SQ / 27.;
    stnl_3[4] = H_SQ_SQ / 27.;
    stnl_3[5] = 5. * H_SQ_SQ / 54.;
    stnl_3[6] = H_SQ_SQ / 27.;
    stnl_3[7] = - 2. * H_SQ_SQ / 27.;
    stnl_3[8] = H_SQ_SQ / 27.;

}
 // !!!!!!!!!!!!!!!!!

void fill_stencils_u_uyy(double H_SQ, double *stnl_1, double *stnl_2, double *stnl_3) {

    // stansil u --> u for u_yy
    stnl_1[0] = 16. / 3.;
    stnl_1[1] = 4. / 3.;
    stnl_1[2] = -2. / 3.;
    stnl_1[3] = -8. / 3.;
    stnl_1[4] = -2. / 3.;
    stnl_1[5] = 4. / 3.;
    stnl_1[6] = -2. / 3.;
    stnl_1[7] = -8. / 3.;
    stnl_1[8] = -2. / 3.;

    // stansil u --> u_xx for u_yy
    stnl_2[0] = -8. * H_SQ / 9.;
    stnl_2[1] = H_SQ / 9.;
    stnl_2[2] = - H_SQ / 18.;
    stnl_2[3] = 4. * H_SQ / 9.;
    stnl_2[4] = - H_SQ / 18.;
    stnl_2[5] = H_SQ / 9.;
    stnl_2[6] = - H_SQ / 18.;
    stnl_2[7] = 4. * H_SQ / 9.;
    stnl_2[8] = - H_SQ / 18.;

    // stansil u --> u_yy for u_yy
    stnl_3[0] = -2. * H_SQ / 9.;
    stnl_3[1] = H_SQ / 9.;
    stnl_3[2] = - H_SQ / 18.;
    stnl_3[3] = H_SQ / 9.;
    stnl_3[4] = - H_SQ / 18.;
    stnl_3[5] = H_SQ / 9.;
    stnl_3[6] = - H_SQ / 18.;
    stnl_3[7] = H_SQ / 9.;
    stnl_3[8] = - H_SQ / 18.;



}

void fill_stencils_uxx_uyy(double H_SQ, double *stnl_1, double *stnl_2, double *stnl_3) {

    double H_SQ_SQ=H_SQ * H_SQ;

    // stansil u_xx --> u for u_yy
    stnl_1[0] = -8. * H_SQ / 9.;
    stnl_1[1] = H_SQ / 9.;
    stnl_1[2] = - H_SQ / 18.;
    stnl_1[3] = 4. * H_SQ / 9.;
    stnl_1[4] = - H_SQ / 18.;
    stnl_1[5] = H_SQ / 9.;
    stnl_1[6] = - H_SQ / 18.;
    stnl_1[7] = 4. * H_SQ / 9.;
    stnl_1[8] = - H_SQ / 18.;

    // stansil u_xx --> u_xx for u_yy
    stnl_2[0] = -5. * H_SQ_SQ / 27.;
    stnl_2[1] = -2. * H_SQ_SQ / 27.;
    stnl_2[2] = H_SQ_SQ / 27.;
    stnl_2[3] = 5. * H_SQ_SQ / 54.;
    stnl_2[4] = H_SQ_SQ / 27.;
    stnl_2[5] = -2. * H_SQ_SQ / 27.;
    stnl_2[6] = H_SQ_SQ / 27.;
    stnl_2[7] = 5. * H_SQ_SQ / 54.;
    stnl_2[8] = H_SQ_SQ / 27.;

    // stansil u_xx --> u_yy for u_yy
    stnl_3[0] = -2. * H_SQ_SQ / 27.;
    stnl_3[1] = H_SQ_SQ / 108.;
    stnl_3[2] = - H_SQ_SQ / 216.;
    stnl_3[3] = H_SQ_SQ / 27.;
    stnl_3[4] = - H_SQ_SQ / 216.;
    stnl_3[5] = H_SQ_SQ / 108.;
    stnl_3[6] = - H_SQ_SQ / 216.;
    stnl_3[7] = H_SQ_SQ / 27.;
    stnl_3[8] = - H_SQ_SQ / 216.;

}

void fill_stencils_uyy_uyy(double H_SQ, double *stnl_1, double *stnl_2, double *stnl_3) {

    double H_SQ_SQ=H_SQ * H_SQ;

    // stansil u_yy --> u for u_xx
    stnl_1[0] = -2. * H_SQ / 9.;
    stnl_1[1] = H_SQ / 9.;
    stnl_1[2] = - H_SQ / 18.;
    stnl_1[3] = H_SQ / 9.;
    stnl_1[4] = - H_SQ / 18.;
    stnl_1[5] = H_SQ / 9.;
    stnl_1[6] = - H_SQ / 18.;
    stnl_1[7] = H_SQ / 9.;
    stnl_1[8] = - H_SQ / 18.;

    // stansil u_yy --> u_xx for u_xx
    stnl_2[0] = -2. * H_SQ_SQ / 27.;
    stnl_2[1] = H_SQ_SQ / 108.;
    stnl_2[2] = - H_SQ_SQ / 216.;
    stnl_2[3] = H_SQ_SQ / 27.;
    stnl_2[4] = - H_SQ_SQ / 216.;
    stnl_2[5] = H_SQ_SQ / 108.;
    stnl_2[6] = - H_SQ_SQ / 216.;
    stnl_2[7] = H_SQ_SQ / 27.;
    stnl_2[8] = - H_SQ_SQ / 216.;

    // stansil u_yy --> u_yy for u_xx
    stnl_3[0] = -11. * H_SQ_SQ / 27.;
    stnl_3[1] = 5. * H_SQ_SQ / 54.;
    stnl_3[2] = H_SQ_SQ / 27.;
    stnl_3[3] = 11. * H_SQ_SQ / 54.;
    stnl_3[4] = H_SQ_SQ / 27.;
    stnl_3[5] = 5. * H_SQ_SQ / 54.;
    stnl_3[6] = H_SQ_SQ / 27.;
    stnl_3[7] = 11. * H_SQ_SQ / 54.;
    stnl_3[8] = H_SQ_SQ / 27.;

}


double Residual_u(double *res_u, double *u, double *u_xx, double *u_yy, double *rp, int flag)
{
    double *st_11_u_xx = new double[9]; // the 1st equation, u for u_xx
    double *st_12_u_xx = new double[9]; // the 1st equation, u_xx for u_xx
    double *st_13_u_xx = new double[9]; // the 1st equation, u_yy for u_xx

    double *st_11_u_yy = new double[9]; // the 1st equation, u for u_yy
    double *st_12_u_yy = new double[9]; // the 1st equation, u_xx for u_yy
    double *st_13_u_yy = new double[9]; // the 1st equation, u_yy for u_yy

    double H = HX; // ! let HX = HY
    double H_SQ = H * H;
    double H_coef = 1.;
    if (flag) H_coef = 1. / H;

    fill_stencils_u_uxx(H_SQ, st_11_u_xx, st_12_u_xx, st_13_u_xx);
    fill_stencils_u_uyy(H_SQ, st_11_u_yy, st_12_u_yy, st_13_u_yy);

    double Res = FLT_MIN;

    // inner points
    for (int i = 1; i < OX_LEN; ++i)
        for (int j = 1; j < OY_LEN; ++j) {
            double val1 = st_11_u_xx[0] * u[OY_LEN_1 * i + j]
                          + st_11_u_xx[1] * u[OY_LEN_1 * (i + 1) + j]
                          + st_11_u_xx[2] * u[OY_LEN_1 * (i + 1) + j + 1]
                          + st_11_u_xx[3] * u[OY_LEN_1 * i + j + 1]
                          + st_11_u_xx[4] * u[OY_LEN_1 * (i - 1) + j + 1]
                          + st_11_u_xx[5] * u[OY_LEN_1 * (i -1) + j]
                          + st_11_u_xx[6] * u[OY_LEN_1 * (i - 1) + j - 1]
                          + st_11_u_xx[7] * u[OY_LEN_1 * i + j - 1]
                          + st_11_u_xx[8] * u[OY_LEN_1 * (i + 1) + j - 1];
            double val2 = st_12_u_xx[0] * u_xx[OY_LEN_1 * i + j]
                          + st_12_u_xx[1] * u_xx[OY_LEN_1 * (i + 1) + j]
                          + st_12_u_xx[2] * u_xx[OY_LEN_1 * (i + 1) + j + 1]
                          + st_12_u_xx[3] * u_xx[OY_LEN_1 * i + j + 1]
                          + st_12_u_xx[4] * u_xx[OY_LEN_1 * (i - 1) + j + 1]
                          + st_12_u_xx[5] * u_xx[OY_LEN_1 * (i -1) + j]
                          + st_12_u_xx[6] * u_xx[OY_LEN_1 * (i - 1) + j - 1]
                          + st_12_u_xx[7] * u_xx[OY_LEN_1 * i + j - 1]
                          + st_12_u_xx[8] * u_xx[OY_LEN_1 * (i + 1) + j - 1];
            double val3 = st_13_u_xx[0] * u_yy[OY_LEN_1 * i + j]
                          + st_13_u_xx[1] * u_yy[OY_LEN_1 * (i + 1) + j]
                          + st_13_u_xx[2] * u_yy[OY_LEN_1 * (i + 1) + j + 1]
                          + st_13_u_xx[3] * u_yy[OY_LEN_1 * i + j + 1]
                          + st_13_u_xx[4] * u_yy[OY_LEN_1 * (i - 1) + j + 1]
                          + st_13_u_xx[5] * u_yy[OY_LEN_1 * (i -1) + j]
                          + st_13_u_xx[6] * u_yy[OY_LEN_1 * (i - 1) + j - 1]
                          + st_13_u_xx[7] * u_yy[OY_LEN_1 * i + j - 1]
                          + st_13_u_xx[8] * u_yy[OY_LEN_1 * (i + 1) + j - 1];
            double val4 = st_11_u_yy[0] * u[OY_LEN_1 * i + j]
                          + st_11_u_yy[1] * u[OY_LEN_1 * (i + 1) + j]
                          + st_11_u_yy[2] * u[OY_LEN_1 * (i + 1) + j + 1]
                          + st_11_u_yy[3] * u[OY_LEN_1 * i + j + 1]
                          + st_11_u_yy[4] * u[OY_LEN_1 * (i - 1) + j + 1]
                          + st_11_u_yy[5] * u[OY_LEN_1 * (i -1) + j]
                          + st_11_u_yy[6] * u[OY_LEN_1 * (i - 1) + j - 1]
                          + st_11_u_yy[7] * u[OY_LEN_1 * i + j - 1]
                          + st_11_u_yy[8] * u[OY_LEN_1 * (i + 1) + j - 1];
            double val5 = st_12_u_yy[0] * u_xx[OY_LEN_1 * i + j]
                          + st_12_u_yy[1] * u_xx[OY_LEN_1 * (i + 1) + j]
                          + st_12_u_yy[2] * u_xx[OY_LEN_1 * (i + 1) + j + 1]
                          + st_12_u_yy[3] * u_xx[OY_LEN_1 * i + j + 1]
                          + st_12_u_yy[4] * u_xx[OY_LEN_1 * (i - 1) + j + 1]
                          + st_12_u_yy[5] * u_xx[OY_LEN_1 * (i -1) + j]
                          + st_12_u_yy[6] * u_xx[OY_LEN_1 * (i - 1) + j - 1]
                          + st_12_u_yy[7] * u_xx[OY_LEN_1 * i + j - 1]
                          + st_12_u_yy[8] * u_xx[OY_LEN_1 * (i + 1) + j - 1];
            double val6 = st_13_u_yy[0] * u_yy[OY_LEN_1 * i + j]
                          + st_13_u_yy[1] * u_yy[OY_LEN_1 * (i + 1) + j]
                          + st_13_u_yy[2] * u_yy[OY_LEN_1 * (i + 1) + j + 1]
                          + st_13_u_yy[3] * u_yy[OY_LEN_1 * i + j + 1]
                          + st_13_u_yy[4] * u_yy[OY_LEN_1 * (i - 1) + j + 1]
                          + st_13_u_yy[5] * u_yy[OY_LEN_1 * (i -1) + j]
                          + st_13_u_yy[6] * u_yy[OY_LEN_1 * (i - 1) + j - 1]
                          + st_13_u_yy[7] * u_yy[OY_LEN_1 * i + j - 1]
                          + st_13_u_yy[8] * u_yy[OY_LEN_1 * (i + 1) + j - 1];

            res_u[OY_LEN_1 * i + j] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6)
                                      - rp[OY_LEN_1 * i + j];

            double val = fabs(res_u[i * OY_LEN_1 + j]);
            if (val > Res) Res = val;
        }

    return Res;
}

double Residual_u_xx(double *res_u_xx, double *u, double *u_xx, double *u_yy, double *rp, int flag)
{
    double *st_11_u_xx = new double[9]; // the 1st equation, u for u_xx
    double *st_12_u_xx = new double[9]; // the 1st equation, u_xx for u_xx
    double *st_13_u_xx = new double[9]; // the 1st equation, u_yy for u_xx

    double *st_11_u_yy = new double[9]; // the 1st equation, u for u_yy
    double *st_12_u_yy = new double[9]; // the 1st equation, u_xx for u_yy
    double *st_13_u_yy = new double[9]; // the 1st equation, u_yy for u_yy

    double H = HX; // ! let HX = HY
    double H_SQ = H * H;
    double H_coef = 1.;
    if (flag) H_coef = 1. / H;

    fill_stencils_uxx_uxx(H_SQ, st_11_u_xx, st_12_u_xx, st_13_u_xx);
    fill_stencils_uxx_uyy(H_SQ, st_11_u_yy, st_12_u_yy, st_13_u_yy);

    double Res = FLT_MIN;

    // inner points
    for (int i = 1; i < OX_LEN; ++i)
        for (int j = 1; j < OY_LEN; ++j) {
            double val1 = st_11_u_xx[0] * u[OY_LEN_1 * i + j]
                          + st_11_u_xx[1] * u[OY_LEN_1 * (i + 1) + j]
                          + st_11_u_xx[2] * u[OY_LEN_1 * (i + 1) + j + 1]
                          + st_11_u_xx[3] * u[OY_LEN_1 * i + j + 1]
                          + st_11_u_xx[4] * u[OY_LEN_1 * (i - 1) + j + 1]
                          + st_11_u_xx[5] * u[OY_LEN_1 * (i -1) + j]
                          + st_11_u_xx[6] * u[OY_LEN_1 * (i - 1) + j - 1]
                          + st_11_u_xx[7] * u[OY_LEN_1 * i + j - 1]
                          + st_11_u_xx[8] * u[OY_LEN_1 * (i + 1) + j - 1];
            double val2 = st_12_u_xx[0] * u_xx[OY_LEN_1 * i + j]
                          + st_12_u_xx[1] * u_xx[OY_LEN_1 * (i + 1) + j]
                          + st_12_u_xx[2] * u_xx[OY_LEN_1 * (i + 1) + j + 1]
                          + st_12_u_xx[3] * u_xx[OY_LEN_1 * i + j + 1]
                          + st_12_u_xx[4] * u_xx[OY_LEN_1 * (i - 1) + j + 1]
                          + st_12_u_xx[5] * u_xx[OY_LEN_1 * (i -1) + j]
                          + st_12_u_xx[6] * u_xx[OY_LEN_1 * (i - 1) + j - 1]
                          + st_12_u_xx[7] * u_xx[OY_LEN_1 * i + j - 1]
                          + st_12_u_xx[8] * u_xx[OY_LEN_1 * (i + 1) + j - 1];
            double val3 = st_13_u_xx[0] * u_yy[OY_LEN_1 * i + j]
                          + st_13_u_xx[1] * u_yy[OY_LEN_1 * (i + 1) + j]
                          + st_13_u_xx[2] * u_yy[OY_LEN_1 * (i + 1) + j + 1]
                          + st_13_u_xx[3] * u_yy[OY_LEN_1 * i + j + 1]
                          + st_13_u_xx[4] * u_yy[OY_LEN_1 * (i - 1) + j + 1]
                          + st_13_u_xx[5] * u_yy[OY_LEN_1 * (i -1) + j]
                          + st_13_u_xx[6] * u_yy[OY_LEN_1 * (i - 1) + j - 1]
                          + st_13_u_xx[7] * u_yy[OY_LEN_1 * i + j - 1]
                          + st_13_u_xx[8] * u_yy[OY_LEN_1 * (i + 1) + j - 1];
            double val4 = st_11_u_yy[0] * u[OY_LEN_1 * i + j]
                          + st_11_u_yy[1] * u[OY_LEN_1 * (i + 1) + j]
                          + st_11_u_yy[2] * u[OY_LEN_1 * (i + 1) + j + 1]
                          + st_11_u_yy[3] * u[OY_LEN_1 * i + j + 1]
                          + st_11_u_yy[4] * u[OY_LEN_1 * (i - 1) + j + 1]
                          + st_11_u_yy[5] * u[OY_LEN_1 * (i -1) + j]
                          + st_11_u_yy[6] * u[OY_LEN_1 * (i - 1) + j - 1]
                          + st_11_u_yy[7] * u[OY_LEN_1 * i + j - 1]
                          + st_11_u_yy[8] * u[OY_LEN_1 * (i + 1) + j - 1];
            double val5 = st_12_u_yy[0] * u_xx[OY_LEN_1 * i + j]
                          + st_12_u_yy[1] * u_xx[OY_LEN_1 * (i + 1) + j]
                          + st_12_u_yy[2] * u_xx[OY_LEN_1 * (i + 1) + j + 1]
                          + st_12_u_yy[3] * u_xx[OY_LEN_1 * i + j + 1]
                          + st_12_u_yy[4] * u_xx[OY_LEN_1 * (i - 1) + j + 1]
                          + st_12_u_yy[5] * u_xx[OY_LEN_1 * (i -1) + j]
                          + st_12_u_yy[6] * u_xx[OY_LEN_1 * (i - 1) + j - 1]
                          + st_12_u_yy[7] * u_xx[OY_LEN_1 * i + j - 1]
                          + st_12_u_yy[8] * u_xx[OY_LEN_1 * (i + 1) + j - 1];
            double val6 = st_13_u_yy[0] * u_yy[OY_LEN_1 * i + j]
                          + st_13_u_yy[1] * u_yy[OY_LEN_1 * (i + 1) + j]
                          + st_13_u_yy[2] * u_yy[OY_LEN_1 * (i + 1) + j + 1]
                          + st_13_u_yy[3] * u_yy[OY_LEN_1 * i + j + 1]
                          + st_13_u_yy[4] * u_yy[OY_LEN_1 * (i - 1) + j + 1]
                          + st_13_u_yy[5] * u_yy[OY_LEN_1 * (i -1) + j]
                          + st_13_u_yy[6] * u_yy[OY_LEN_1 * (i - 1) + j - 1]
                          + st_13_u_yy[7] * u_yy[OY_LEN_1 * i + j - 1]
                          + st_13_u_yy[8] * u_yy[OY_LEN_1 * (i + 1) + j - 1];

            res_u_xx[OY_LEN_1 * i + j] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6)
                                         - rp[OY_LEN_1 * i + j];

            double val = fabs(res_u_xx[i * OY_LEN_1 + j]);
            if (val > Res) Res = val;
        }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 1; j < OY_LEN; ++j) {

        double val1 = st_11_u_xx[0] * u[OY_LEN_1 * OX_LEN + j] / 2.
                      + st_11_u_xx[3] * u[OY_LEN_1 * OX_LEN + j + 1] / 2.
                      + st_11_u_xx[4] * u[OY_LEN_1 * (OX_LEN - 1) + j + 1]
                      + st_11_u_xx[5] * u[OY_LEN_1 * (OX_LEN -1) + j]
                      + st_11_u_xx[6] * u[OY_LEN_1 * (OX_LEN - 1) + j - 1]
                      + st_11_u_xx[7] * u[OY_LEN_1 * OX_LEN + j - 1] / 2.;
        double val2 = st_12_u_xx[0] * u_xx[OY_LEN_1 * OX_LEN + j] / 2.
                      + st_12_u_xx[3] * u_xx[OY_LEN_1 * OX_LEN + j + 1] / 2.
                      + st_12_u_xx[4] * u_xx[OY_LEN_1 * (OX_LEN - 1) + j + 1]
                      + st_12_u_xx[5] * u_xx[OY_LEN_1 * (OX_LEN -1) + j]
                      + st_12_u_xx[6] * u_xx[OY_LEN_1 * (OX_LEN - 1) + j - 1]
                      + st_12_u_xx[7] * u_xx[OY_LEN_1 * OX_LEN + j - 1] / 2.;
        double val3 = st_13_u_xx[0] * u_yy[OY_LEN_1 * OX_LEN + j] / 2.
                      + st_13_u_xx[3] * u_yy[OY_LEN_1 * OX_LEN + j + 1] / 2.
                      + st_13_u_xx[4] * u_yy[OY_LEN_1 * (OX_LEN - 1) + j + 1]
                      + st_13_u_xx[5] * u_yy[OY_LEN_1 * (OX_LEN -1) + j]
                      + st_13_u_xx[6] * u_yy[OY_LEN_1 * (OX_LEN - 1) + j - 1]
                      + st_13_u_xx[7] * u_yy[OY_LEN_1 * OX_LEN + j - 1] / 2.;
        double val4 = st_11_u_yy[0] * u[OY_LEN_1 * OX_LEN + j] / 2.
                      + st_11_u_yy[3] * u[OY_LEN_1 * OX_LEN + j + 1] / 2.
                      + st_11_u_yy[4] * u[OY_LEN_1 * (OX_LEN - 1) + j + 1]
                      + st_11_u_yy[5] * u[OY_LEN_1 * (OX_LEN -1) + j]
                      + st_11_u_yy[6] * u[OY_LEN_1 * (OX_LEN - 1) + j - 1]
                      + st_11_u_yy[7] * u[OY_LEN_1 * OX_LEN + j - 1] / 2.;
        double val5 = st_12_u_yy[0] * u_xx[OY_LEN_1 * OX_LEN + j] / 2.
                      + st_12_u_yy[3] * u_xx[OY_LEN_1 * OX_LEN + j + 1] / 2.
                      + st_12_u_yy[4] * u_xx[OY_LEN_1 * (OX_LEN - 1) + j + 1]
                      + st_12_u_yy[5] * u_xx[OY_LEN_1 * (OX_LEN -1) + j]
                      + st_12_u_yy[6] * u_xx[OY_LEN_1 * (OX_LEN - 1) + j - 1]
                      + st_12_u_yy[7] * u_xx[OY_LEN_1 * OX_LEN + j - 1] / 2.;
        double val6 = st_13_u_yy[0] * u_yy[OY_LEN_1 * OX_LEN + j] / 2.
                      + st_13_u_yy[3] * u_yy[OY_LEN_1 * OX_LEN + j + 1] / 2.
                      + st_13_u_yy[4] * u_yy[OY_LEN_1 * (OX_LEN - 1) + j + 1]
                      + st_13_u_yy[5] * u_yy[OY_LEN_1 * (OX_LEN -1) + j]
                      + st_13_u_yy[6] * u_yy[OY_LEN_1 * (OX_LEN - 1) + j - 1]
                      + st_13_u_yy[7] * u_yy[OY_LEN_1 * OX_LEN + j - 1] / 2.;

        double val7 = Simpson_1D_x(u, u_xx, u_yy, OX_LEN, j, OX_LEN - 1, j)
                      + Simpson_1D_x(u, u_xx, u_yy, OX_LEN, j, OX_LEN - 1, j - 1);

        res_u_xx[OY_LEN_1 * OX_LEN + j] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6) - val7
                                          - rp[OY_LEN_1 * OX_LEN + j];

        double val = fabs(res_u_xx[OY_LEN_1 * OX_LEN + j]);
        if (val > Res) Res = val;


        val1 = st_11_u_xx[0] * u[j] / 2.
                      + st_11_u_xx[1] * u[OY_LEN_1 + j]
                      + st_11_u_xx[2] * u[OY_LEN_1 + j + 1]
                      + st_11_u_xx[3] * u[j + 1] / 2.
                      + st_11_u_xx[7] * u[j - 1] / 2.
                      + st_11_u_xx[8] * u[OY_LEN_1 + j - 1];
        val2 = st_12_u_xx[0] * u_xx[j] / 2.
                      + st_12_u_xx[1] * u_xx[OY_LEN_1 + j]
                      + st_12_u_xx[2] * u_xx[OY_LEN_1 + j + 1]
                      + st_12_u_xx[3] * u_xx[j + 1] / 2.
                      + st_12_u_xx[7] * u_xx[j - 1] / 2.
                      + st_12_u_xx[8] * u_xx[OY_LEN_1 + j - 1];
        val3 = st_13_u_xx[0] * u_yy[j] / 2.
                      + st_13_u_xx[1] * u_yy[OY_LEN_1 + j]
                      + st_13_u_xx[2] * u_yy[OY_LEN_1 + j + 1]
                      + st_13_u_xx[3] * u_yy[j + 1] / 2.
                      + st_13_u_xx[7] * u_yy[j - 1] / 2.
                      + st_13_u_xx[8] * u_yy[OY_LEN_1 + j - 1];
        val4 = st_11_u_yy[0] * u[j] / 2.
                      + st_11_u_yy[1] * u[OY_LEN_1 + j]
                      + st_11_u_yy[2] * u[OY_LEN_1 + j + 1]
                      + st_11_u_yy[3] * u[j + 1] / 2.
                      + st_11_u_yy[7] * u[j - 1] / 2.
                      + st_11_u_yy[8] * u[OY_LEN_1 + j - 1];
        val5 = st_12_u_yy[0] * u_xx[j] / 2.
                      + st_12_u_yy[1] * u_xx[OY_LEN_1 + j]
                      + st_12_u_yy[2] * u_xx[OY_LEN_1 + j + 1]
                      + st_12_u_yy[3] * u_xx[j + 1] / 2.
                      + st_12_u_yy[7] * u_xx[j - 1] / 2.
                      + st_12_u_yy[8] * u_xx[OY_LEN_1 + j - 1];
        val6 = st_13_u_yy[0] * u_yy[j] / 2.
                      + st_13_u_yy[1] * u_yy[OY_LEN_1 + j]
                      + st_13_u_yy[2] * u_yy[OY_LEN_1 + j + 1]
                      + st_13_u_yy[3] * u_yy[j + 1] / 2.
                      + st_13_u_yy[7] * u_yy[j - 1] / 2.
                      + st_13_u_yy[8] * u_yy[OY_LEN_1 + j - 1];

        val7 = Simpson_1D_x(u, u_xx, u_yy, 0, j, 0, j)
               + Simpson_1D_x(u, u_xx, u_yy, 0, j, 0, j-1);

        res_u_xx[j] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6) + val7 - rp[j];

        val = fabs(res_u_xx[j]);
        if (val > Res) Res = val;
    }

    return Res;
}

double Residual_u_yy(double *res_u_yy, double *u, double *u_xx, double *u_yy, double *rp, int flag)
{
    double *st_11_u_xx = new double[9]; // the 1st equation, u for u_xx
    double *st_12_u_xx = new double[9]; // the 1st equation, u_xx for u_xx
    double *st_13_u_xx = new double[9]; // the 1st equation, u_yy for u_xx

    double *st_11_u_yy = new double[9]; // the 1st equation, u for u_yy
    double *st_12_u_yy = new double[9]; // the 1st equation, u_xx for u_yy
    double *st_13_u_yy = new double[9]; // the 1st equation, u_yy for u_yy

    double H = HX; // ! let HX = HY
    double H_SQ = H * H;
    double H_coef = 1.;
    if (flag) H_coef = 1. / H;

    fill_stencils_uyy_uxx(H_SQ, st_11_u_xx, st_12_u_xx, st_13_u_xx);
    fill_stencils_uyy_uyy(H_SQ, st_11_u_yy, st_12_u_yy, st_13_u_yy);

    double Res = FLT_MIN;

    // inner points
    for (int i = 1; i < OX_LEN; ++i)
        for (int j = 1; j < OY_LEN; ++j) {
            double val1 = st_11_u_xx[0] * u[OY_LEN_1 * i + j]
                          + st_11_u_xx[1] * u[OY_LEN_1 * (i + 1) + j]
                          + st_11_u_xx[2] * u[OY_LEN_1 * (i + 1) + j + 1]
                          + st_11_u_xx[3] * u[OY_LEN_1 * i + j + 1]
                          + st_11_u_xx[4] * u[OY_LEN_1 * (i - 1) + j + 1]
                          + st_11_u_xx[5] * u[OY_LEN_1 * (i -1) + j]
                          + st_11_u_xx[6] * u[OY_LEN_1 * (i - 1) + j - 1]
                          + st_11_u_xx[7] * u[OY_LEN_1 * i + j - 1]
                          + st_11_u_xx[8] * u[OY_LEN_1 * (i + 1) + j - 1];
            double val2 = st_12_u_xx[0] * u_xx[OY_LEN_1 * i + j]
                          + st_12_u_xx[1] * u_xx[OY_LEN_1 * (i + 1) + j]
                          + st_12_u_xx[2] * u_xx[OY_LEN_1 * (i + 1) + j + 1]
                          + st_12_u_xx[3] * u_xx[OY_LEN_1 * i + j + 1]
                          + st_12_u_xx[4] * u_xx[OY_LEN_1 * (i - 1) + j + 1]
                          + st_12_u_xx[5] * u_xx[OY_LEN_1 * (i -1) + j]
                          + st_12_u_xx[6] * u_xx[OY_LEN_1 * (i - 1) + j - 1]
                          + st_12_u_xx[7] * u_xx[OY_LEN_1 * i + j - 1]
                          + st_12_u_xx[8] * u_xx[OY_LEN_1 * (i + 1) + j - 1];
            double val3 = st_13_u_xx[0] * u_yy[OY_LEN_1 * i + j]
                          + st_13_u_xx[1] * u_yy[OY_LEN_1 * (i + 1) + j]
                          + st_13_u_xx[2] * u_yy[OY_LEN_1 * (i + 1) + j + 1]
                          + st_13_u_xx[3] * u_yy[OY_LEN_1 * i + j + 1]
                          + st_13_u_xx[4] * u_yy[OY_LEN_1 * (i - 1) + j + 1]
                          + st_13_u_xx[5] * u_yy[OY_LEN_1 * (i -1) + j]
                          + st_13_u_xx[6] * u_yy[OY_LEN_1 * (i - 1) + j - 1]
                          + st_13_u_xx[7] * u_yy[OY_LEN_1 * i + j - 1]
                          + st_13_u_xx[8] * u_yy[OY_LEN_1 * (i + 1) + j - 1];
            double val4 = st_11_u_yy[0] * u[OY_LEN_1 * i + j]
                          + st_11_u_yy[1] * u[OY_LEN_1 * (i + 1) + j]
                          + st_11_u_yy[2] * u[OY_LEN_1 * (i + 1) + j + 1]
                          + st_11_u_yy[3] * u[OY_LEN_1 * i + j + 1]
                          + st_11_u_yy[4] * u[OY_LEN_1 * (i - 1) + j + 1]
                          + st_11_u_yy[5] * u[OY_LEN_1 * (i -1) + j]
                          + st_11_u_yy[6] * u[OY_LEN_1 * (i - 1) + j - 1]
                          + st_11_u_yy[7] * u[OY_LEN_1 * i + j - 1]
                          + st_11_u_yy[8] * u[OY_LEN_1 * (i + 1) + j - 1];
            double val5 = st_12_u_yy[0] * u_xx[OY_LEN_1 * i + j]
                          + st_12_u_yy[1] * u_xx[OY_LEN_1 * (i + 1) + j]
                          + st_12_u_yy[2] * u_xx[OY_LEN_1 * (i + 1) + j + 1]
                          + st_12_u_yy[3] * u_xx[OY_LEN_1 * i + j + 1]
                          + st_12_u_yy[4] * u_xx[OY_LEN_1 * (i - 1) + j + 1]
                          + st_12_u_yy[5] * u_xx[OY_LEN_1 * (i -1) + j]
                          + st_12_u_yy[6] * u_xx[OY_LEN_1 * (i - 1) + j - 1]
                          + st_12_u_yy[7] * u_xx[OY_LEN_1 * i + j - 1]
                          + st_12_u_yy[8] * u_xx[OY_LEN_1 * (i + 1) + j - 1];
            double val6 = st_13_u_yy[0] * u_yy[OY_LEN_1 * i + j]
                          + st_13_u_yy[1] * u_yy[OY_LEN_1 * (i + 1) + j]
                          + st_13_u_yy[2] * u_yy[OY_LEN_1 * (i + 1) + j + 1]
                          + st_13_u_yy[3] * u_yy[OY_LEN_1 * i + j + 1]
                          + st_13_u_yy[4] * u_yy[OY_LEN_1 * (i - 1) + j + 1]
                          + st_13_u_yy[5] * u_yy[OY_LEN_1 * (i -1) + j]
                          + st_13_u_yy[6] * u_yy[OY_LEN_1 * (i - 1) + j - 1]
                          + st_13_u_yy[7] * u_yy[OY_LEN_1 * i + j - 1]
                          + st_13_u_yy[8] * u_yy[OY_LEN_1 * (i + 1) + j - 1];

            res_u_yy[OY_LEN_1 * i + j] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6)
                                         - rp[OY_LEN_1 * i + j];


            double val = fabs(res_u_yy[i * OY_LEN_1 + j]);
            if (val > Res) Res = val;
        }

    // G1 -- (x_i, 0=C) -- bottom boundary
    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 1; i < OX_LEN; ++i) {
        double val1 = st_11_u_xx[0] * u[OY_LEN_1 * i] / 2.
                      + st_11_u_xx[1] * u[OY_LEN_1 * (i + 1)] / 2.
                      + st_11_u_xx[2] * u[OY_LEN_1 * (i + 1) + 1]
                      + st_11_u_xx[3] * u[OY_LEN_1 * i + 1]
                      + st_11_u_xx[4] * u[OY_LEN_1 * (i - 1) + 1]
                      + st_11_u_xx[5] * u[OY_LEN_1 * (i -1)] / 2.;
        double val2 = st_12_u_xx[0] * u_xx[OY_LEN_1 * i] / 2.
                      + st_12_u_xx[1] * u_xx[OY_LEN_1 * (i + 1)] / 2.
                      + st_12_u_xx[2] * u_xx[OY_LEN_1 * (i + 1) + 1]
                      + st_12_u_xx[3] * u_xx[OY_LEN_1 * i + 1]
                      + st_12_u_xx[4] * u_xx[OY_LEN_1 * (i - 1) + 1]
                      + st_12_u_xx[5] * u_xx[OY_LEN_1 * (i -1)] / 2.;
        double val3 = st_13_u_xx[0] * u_yy[OY_LEN_1 * i] / 2.
                      + st_13_u_xx[1] * u_yy[OY_LEN_1 * (i + 1)] / 2.
                      + st_13_u_xx[2] * u_yy[OY_LEN_1 * (i + 1) + 1]
                      + st_13_u_xx[3] * u_yy[OY_LEN_1 * i + 1]
                      + st_13_u_xx[4] * u_yy[OY_LEN_1 * (i - 1) + 1]
                      + st_13_u_xx[5] * u_yy[OY_LEN_1 * (i -1)] / 2.;
        double val4 = st_11_u_yy[0] * u[OY_LEN_1 * i] / 2.
                      + st_11_u_yy[1] * u[OY_LEN_1 * (i + 1)] / 2.
                      + st_11_u_yy[2] * u[OY_LEN_1 * (i + 1) + 1]
                      + st_11_u_yy[3] * u[OY_LEN_1 * i + 1]
                      + st_11_u_yy[4] * u[OY_LEN_1 * (i - 1) + 1]
                      + st_11_u_yy[5] * u[OY_LEN_1 * (i -1)] / 2.;
        double val5 = st_12_u_yy[0] * u_xx[OY_LEN_1 * i] / 2.
                      + st_12_u_yy[1] * u_xx[OY_LEN_1 * (i + 1)] / 2.
                      + st_12_u_yy[2] * u_xx[OY_LEN_1 * (i + 1) + 1]
                      + st_12_u_yy[3] * u_xx[OY_LEN_1 * i + 1]
                      + st_12_u_yy[4] * u_xx[OY_LEN_1 * (i - 1) + 1]
                      + st_12_u_yy[5] * u_xx[OY_LEN_1 * (i -1)] / 2.;
        double val6 = st_13_u_yy[0] * u_yy[OY_LEN_1 * i] / 2.
                      + st_13_u_yy[1] * u_yy[OY_LEN_1 * (i + 1)] / 2.
                      + st_13_u_yy[2] * u_yy[OY_LEN_1 * (i + 1) + 1]
                      + st_13_u_yy[3] * u_yy[OY_LEN_1 * i + 1]
                      + st_13_u_yy[4] * u_yy[OY_LEN_1 * (i - 1) + 1]
                      + st_13_u_yy[5] * u_yy[OY_LEN_1 * (i - 1)] / 2.;

        double val7 = Simpson_1D_y(u, u_xx, u_yy, i, 0, i, 0)
                      + Simpson_1D_y(u, u_xx, u_yy, i-1, 0, i, 0);

        res_u_yy[OY_LEN_1 * i] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6) + val7
                                 - rp[OY_LEN_1 * i];

        double val = fabs(res_u_yy[i * OY_LEN_1]);
        if (val > Res) Res = val;


        val1 = st_11_u_xx[0] * u[OY_LEN_1 * i + OY_LEN] / 2.
               + st_11_u_xx[1] * u[OY_LEN_1 * (i + 1) + OY_LEN] / 2.
               + st_11_u_xx[5] * u[OY_LEN_1 * (i -1) + OY_LEN] / 2.
               + st_11_u_xx[6] * u[OY_LEN_1 * (i - 1) + OY_LEN - 1]
               + st_11_u_xx[7] * u[OY_LEN_1 * i + OY_LEN - 1]
               + st_11_u_xx[8] * u[OY_LEN_1 * (i + 1) + OY_LEN - 1];
        val2 = st_12_u_xx[0] * u_xx[OY_LEN_1 * i + OY_LEN] / 2.
               + st_12_u_xx[1] * u_xx[OY_LEN_1 * (i + 1) + OY_LEN] / 2.
               + st_12_u_xx[5] * u_xx[OY_LEN_1 * (i -1) + OY_LEN] / 2.
               + st_12_u_xx[6] * u_xx[OY_LEN_1 * (i - 1) + OY_LEN - 1]
               + st_12_u_xx[7] * u_xx[OY_LEN_1 * i + OY_LEN - 1]
               + st_12_u_xx[8] * u_xx[OY_LEN_1 * (i + 1) + OY_LEN - 1];
        val3 = st_13_u_xx[0] * u_yy[OY_LEN_1 * i + OY_LEN] / 2.
               + st_13_u_xx[1] * u_yy[OY_LEN_1 * (i + 1) + OY_LEN] / 2.
               + st_13_u_xx[5] * u_yy[OY_LEN_1 * (i -1) + OY_LEN] / 2.
               + st_13_u_xx[6] * u_yy[OY_LEN_1 * (i - 1) + OY_LEN - 1]
               + st_13_u_xx[7] * u_yy[OY_LEN_1 * i + OY_LEN - 1]
               + st_13_u_xx[8] * u_yy[OY_LEN_1 * (i + 1) + OY_LEN - 1];
        val4 = st_11_u_yy[0] * u[OY_LEN_1 * i + OY_LEN] / 2.
               + st_11_u_yy[1] * u[OY_LEN_1 * (i + 1) + OY_LEN] / 2.
               + st_11_u_yy[5] * u[OY_LEN_1 * (i -1) + OY_LEN] / 2.
               + st_11_u_yy[6] * u[OY_LEN_1 * (i - 1) + OY_LEN - 1]
               + st_11_u_yy[7] * u[OY_LEN_1 * i + OY_LEN - 1]
               + st_11_u_yy[8] * u[OY_LEN_1 * (i + 1) + OY_LEN - 1];
        val5 = st_12_u_yy[0] * u_xx[OY_LEN_1 * i + OY_LEN] / 2.
               + st_12_u_yy[1] * u_xx[OY_LEN_1 * (i + 1) + OY_LEN] / 2.
               + st_12_u_yy[5] * u_xx[OY_LEN_1 * (i -1) + OY_LEN] / 2.
               + st_12_u_yy[6] * u_xx[OY_LEN_1 * (i - 1) + OY_LEN - 1]
               + st_12_u_yy[7] * u_xx[OY_LEN_1 * i + OY_LEN - 1]
               + st_12_u_yy[8] * u_xx[OY_LEN_1 * (i + 1) + OY_LEN - 1];
        val6 = st_13_u_yy[0] * u_yy[OY_LEN_1 * i + OY_LEN] / 2.
               + st_13_u_yy[1] * u_yy[OY_LEN_1 * (i + 1) + OY_LEN] / 2.
               + st_13_u_yy[5] * u_yy[OY_LEN_1 * (i -1) + OY_LEN] / 2.
               + st_13_u_yy[6] * u_yy[OY_LEN_1 * (i - 1) + OY_LEN - 1]
               + st_13_u_yy[7] * u_yy[OY_LEN_1 * i + OY_LEN - 1]
               + st_13_u_yy[8] * u_yy[OY_LEN_1 * (i + 1) + OY_LEN - 1];

        val7 = Simpson_1D_y(u, u_xx, u_yy, i, OY_LEN - 1, i, OY_LEN)
               + Simpson_1D_y(u, u_xx, u_yy, i - 1, OY_LEN - 1, i, OY_LEN);

        res_u_yy[OY_LEN_1 * i + OY_LEN] = 0.25 *(val1 + val2 + val3 + val4 + val5 + val6) - val7
                                          - rp[OY_LEN_1 * i + OY_LEN];

        val = fabs(res_u_yy[i * OY_LEN_1 + OY_LEN]);
        if (val > Res) Res = val;
    }

    return Res;
}

double *solve_3(double &tme) {
    double H_coef = 1. / HX;
    StartTimer();

    fflush(stdout);

    int ic = 0;
    double *rp_u = new double[XY_LEN];
    double *rp_u_xx = new double[XY_LEN];
    double *rp_u_yy = new double[XY_LEN];
    double *prev_u = new double[XY_LEN];
    double *u = new double[XY_LEN];
    double *prev_u_xx = new double[XY_LEN];
    double *u_xx = new double[XY_LEN];
    double *prev_u_yy = new double[XY_LEN];
    double *u_yy = new double[XY_LEN];
    double *residual_u = new double[XY_LEN];
    double *residual_u_xx = new double[XY_LEN];
    double *residual_u_yy = new double[XY_LEN];
    double *dif_u = new double[XY_LEN];
    double *dif_u_xx = new double[XY_LEN];
    double *dif_u_yy = new double[XY_LEN];

    double *u_exact = new double[XY_LEN];
    double *u_exact_xx = new double[XY_LEN];
    double *u_exact_yy = new double[XY_LEN];
    double *f_u_check = new double[XY_LEN];
    double *f_u_xx_check = new double[XY_LEN];
    double *f_u_yy_check = new double[XY_LEN];


    //<editor-fold desc="Fill initial data">

    for (int i = 0; i < OX_LEN_1; ++i) {
        for (int j = 0; j < OY_LEN_1; ++j) {
            u[OY_LEN_1 * i + j] = 0.;
            u_xx[OY_LEN_1 * i + j] = 0.;
            u_yy[OY_LEN_1 * i + j] = 0.;
            prev_u[OY_LEN_1 * i + j] = 0.;
            prev_u_xx[OY_LEN_1 * i + j] = 0.;
            prev_u_yy[OY_LEN_1 * i + j] = 0.;
            residual_u[OY_LEN_1 * i + j] = 0.;
            residual_u_xx[OY_LEN_1 * i + j] = 0.;
            residual_u_yy[OY_LEN_1 * i + j] = 0.;
            rp_u[OY_LEN_1 * i + j] = 0.;
            rp_u_xx[OY_LEN_1 * i + j] = 0.;
            rp_u_yy[OY_LEN_1 * i + j] = 0.;
            dif_u[OY_LEN_1 * i + j] = 0.;
            dif_u_xx[OY_LEN_1 * i + j] = 0.;
            dif_u_yy[OY_LEN_1 * i + j] = 0.;
            u_exact[OY_LEN_1 * i + j] = 0.;
            u_exact_xx[OY_LEN_1 * i + j] = 0.;
            u_exact_yy[OY_LEN_1 * i + j] = 0.;
            f_u_check[OY_LEN_1 * i + j] = 0.;
            f_u_xx_check[OY_LEN_1 * i + j] = 0.;
            f_u_yy_check[OY_LEN_1 * i + j] = 0.;
        }
    }

    // G1 -- (x_i, 0=C) -- bottom boundary
    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 0; i < OX_LEN_1; ++i) {
        u_exact[OY_LEN_1 * i] = analytical_slv(A + HX * i, C);
        u_exact[OY_LEN_1 * i + OY_LEN] = analytical_slv(A + HX * i, C + HY * OY_LEN);
        u_exact_xx[OY_LEN_1 * i] = analytical_slv_xx(A + HX * i, C);
        u_exact_xx[OY_LEN_1 * i + OY_LEN] = analytical_slv_xx(A + HX * i, C + HY * OY_LEN);
        u_exact_yy[OY_LEN_1 * i] = analytical_slv_yy(A + HX * i, C);
        u_exact_yy[OY_LEN_1 * i + OY_LEN] = analytical_slv_yy(A + HX * i, C + HY * OY_LEN);

        prev_u[OY_LEN_1 * i] = analytical_slv(A + HX * i, C);
        prev_u[OY_LEN_1 * i + OY_LEN] = analytical_slv(A + HX * i, C + HY * OY_LEN);
        prev_u_xx[OY_LEN_1 * i] = analytical_slv_xx(A + HX * i, C);
        prev_u_xx[OY_LEN_1 * i + OY_LEN] = analytical_slv_xx(A + HX * i, C + HY * OY_LEN);
    }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 0; j < OY_LEN_1; ++j) {
        u_exact[OY_LEN_1 * OX_LEN + j] = analytical_slv(A + HX * OX_LEN, C + HY * j);
        u_exact[j] = analytical_slv(A, C + HY * j);
        u_exact_xx[OY_LEN_1 * OX_LEN + j] = analytical_slv_xx(A + HX * OX_LEN, C + HY * j);
        u_exact_xx[j] = analytical_slv_xx(A, C + HY * j);
        u_exact_yy[OY_LEN_1 * OX_LEN + j] = analytical_slv_yy(A + HX * OX_LEN, C + HY * j);
        u_exact_yy[j] = analytical_slv_yy(A, C + HY * j);

        prev_u[OY_LEN_1 * OX_LEN + j] = analytical_slv(A + HX * OX_LEN, C + HY * j);
        prev_u[j] = analytical_slv(A, C + HY * j);
        prev_u_yy[OY_LEN_1 * OX_LEN + j] = analytical_slv_yy(A + HX * OX_LEN, C + HY * j);
        prev_u_yy[j] = analytical_slv_yy(A, C + HY * j);
    }


    // to check a residual, inner points are calculated
    for (int i = 1; i < OX_LEN; ++i)
        for (int j = 1; j < OY_LEN; ++j) {
            u_exact[OY_LEN_1 * i + j] = analytical_slv(A + HX * i, C + HY * j);
            u_exact_xx[OY_LEN_1 * i + j] = analytical_slv_xx(A + HX * i, C + HY * j);
            u_exact_yy[OY_LEN_1 * i + j] = analytical_slv_yy(A + HX * i, C + HY * j);
        }

    memcpy(u, prev_u, XY_LEN * sizeof(double));
    memcpy(u_xx, prev_u_xx, XY_LEN * sizeof(double));
    memcpy(u_yy, prev_u_yy, XY_LEN * sizeof(double));



    //</editor-fold>

    //<editor-fold desc="Calculate right-hand side rp">

    // G1 -- (x_i, 0=C) -- bottom boundary
    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 1; i < OX_LEN; ++i) {
        rp_u_yy[OY_LEN_1 * i] = ff_yy_qrule_Simpson(i, 0, i, 0)
                                + ff_yy_qrule_Simpson(i - 1, 0, i, 0);
        rp_u_yy[OY_LEN_1 * i + OY_LEN] = ff_yy_qrule_Simpson(i - 1, OY_LEN - 1, i, OY_LEN)
                                         + ff_yy_qrule_Simpson(i, OY_LEN - 1, i, OY_LEN);

        f_u_check[OY_LEN_1 * i] = ff_qrule_Simpson(i, 0, i, 0)
                                  + ff_qrule_Simpson(i - 1, 0, i, 0);
        f_u_xx_check[OY_LEN_1 * i] = ff_xx_qrule_Simpson(i, 0, i, 0)
                                  + ff_xx_qrule_Simpson(i - 1, 0, i, 0);
        f_u_yy_check[OY_LEN_1 * i] = ff_yy_qrule_Simpson(i, 0, i, 0)
                                  + ff_yy_qrule_Simpson(i - 1, 0, i, 0);

        f_u_check[OY_LEN_1 * i + OY_LEN] = ff_qrule_Simpson(i - 1, OY_LEN - 1, i, OY_LEN)
                                         + ff_qrule_Simpson(i, OY_LEN - 1, i, OY_LEN);
        f_u_xx_check[OY_LEN_1 * i + OY_LEN] = ff_xx_qrule_Simpson(i - 1, OY_LEN - 1, i, OY_LEN)
                                           + ff_xx_qrule_Simpson(i, OY_LEN - 1, i, OY_LEN);
        f_u_yy_check[OY_LEN_1 * i + OY_LEN] = ff_yy_qrule_Simpson(i - 1, OY_LEN - 1, i, OY_LEN)
                                           + ff_yy_qrule_Simpson(i, OY_LEN - 1, i, OY_LEN);
    }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 1; j < OY_LEN; ++j) {
        rp_u_xx[OY_LEN_1 * OX_LEN + j] = ff_xx_qrule_Simpson(OX_LEN - 1, j, OX_LEN, j)
                                         + ff_xx_qrule_Simpson(OX_LEN - 1, j - 1, OX_LEN, j);
        rp_u_xx[j] = ff_xx_qrule_Simpson(0, j, 0, j)
                     + ff_xx_qrule_Simpson(0, j - 1, 0, j);

        f_u_check[OY_LEN_1 * OX_LEN + j] = ff_qrule_Simpson(OX_LEN - 1, j, OX_LEN, j)
                                         + ff_qrule_Simpson(OX_LEN - 1, j - 1, OX_LEN, j);
        f_u_xx_check[OY_LEN_1 * OX_LEN + j] = ff_xx_qrule_Simpson(OX_LEN - 1, j, OX_LEN, j)
                                           + ff_xx_qrule_Simpson(OX_LEN - 1, j - 1, OX_LEN, j);
        f_u_yy_check[OY_LEN_1 * OX_LEN + j] = ff_yy_qrule_Simpson(OX_LEN - 1, j, OX_LEN, j)
                                           + ff_yy_qrule_Simpson(OX_LEN - 1, j - 1, OX_LEN, j);

        f_u_check[j] = ff_qrule_Simpson(0, j, 0, j)
                     + ff_qrule_Simpson(0, j - 1, 0, j);
        f_u_xx_check[j] = ff_xx_qrule_Simpson(0, j, 0, j)
                     + ff_xx_qrule_Simpson(0, j - 1, 0, j);
        f_u_yy_check[j] = ff_yy_qrule_Simpson(0, j, 0, j)
                     + ff_yy_qrule_Simpson(0, j - 1, 0, j);

    }

    // inner points
    for (int i = 1; i < OX_LEN; ++i)
        for (int j = 1; j < OY_LEN; ++j) {
            rp_u[OY_LEN_1 * i + j] = ff_qrule_Simpson(i, j, i, j)
                                     + ff_qrule_Simpson(i - 1, j, i, j)
                                     + ff_qrule_Simpson(i - 1, j - 1, i, j)
                                     + ff_qrule_Simpson(i, j - 1, i, j);
            rp_u_xx[OY_LEN_1 * i + j] = ff_xx_qrule_Simpson(i, j, i, j)
                                        + ff_xx_qrule_Simpson(i - 1, j, i, j)
                                        + ff_xx_qrule_Simpson(i - 1, j - 1, i, j)
                                        + ff_xx_qrule_Simpson(i, j - 1, i, j);
            rp_u_yy[OY_LEN_1 * i + j] = ff_yy_qrule_Simpson(i, j, i, j)
                                        + ff_yy_qrule_Simpson(i - 1, j, i, j)
                                        + ff_yy_qrule_Simpson(i - 1, j - 1, i, j)
                                        + ff_yy_qrule_Simpson(i, j - 1, i, j);

            f_u_check[OY_LEN_1 * i + j] = ff_qrule_Simpson(i, j, i, j)
                                     + ff_qrule_Simpson(i - 1, j, i, j)
                                     + ff_qrule_Simpson(i - 1, j - 1, i, j)
                                     + ff_qrule_Simpson(i, j - 1, i, j);
            f_u_xx_check[OY_LEN_1 * i + j] = ff_xx_qrule_Simpson(i, j, i, j)
                                        + ff_xx_qrule_Simpson(i - 1, j, i, j)
                                        + ff_xx_qrule_Simpson(i - 1, j - 1, i, j)
                                        + ff_xx_qrule_Simpson(i, j - 1, i, j);
            f_u_yy_check[OY_LEN_1 * i + j] = ff_yy_qrule_Simpson(i, j, i, j)
                                        + ff_yy_qrule_Simpson(i - 1, j, i, j)
                                        + ff_yy_qrule_Simpson(i - 1, j - 1, i, j)
                                        + ff_yy_qrule_Simpson(i, j - 1, i, j);
        }

    //(0,0)
    f_u_check[0] = ff_qrule_Simpson(0, 0, 0, 0);
    f_u_xx_check[0] = ff_xx_qrule_Simpson(0, 0, 0, 0);
    f_u_yy_check[0] = ff_yy_qrule_Simpson(0, 0, 0, 0);

    //(1,0) = (OX_LEN,0)
    f_u_check[OY_LEN_1 * OX_LEN] = ff_qrule_Simpson(OX_LEN-1, 0, OX_LEN, 0);
    f_u_xx_check[OY_LEN_1 * OX_LEN] = ff_xx_qrule_Simpson(OX_LEN-1, 0, OX_LEN, 0);
    f_u_yy_check[OY_LEN_1 * OX_LEN] = ff_yy_qrule_Simpson(OX_LEN-1, 0, OX_LEN, 0);

    //(1,1) = (OX_LEN,OY_LEN)
    f_u_check[OY_LEN_1 * OX_LEN + OY_LEN] = ff_qrule_Simpson(OX_LEN-1, OY_LEN-1, OX_LEN, OY_LEN);
    f_u_xx_check[OY_LEN_1 * OX_LEN + OY_LEN] = ff_xx_qrule_Simpson(OX_LEN-1, OY_LEN-1, OX_LEN, OY_LEN);
    f_u_yy_check[OY_LEN_1 * OX_LEN + OY_LEN] = ff_yy_qrule_Simpson(OX_LEN-1, OY_LEN-1, OX_LEN, OY_LEN);

    //(0,1) = (0,OY_LEN)
    f_u_check[OY_LEN] = ff_qrule_Simpson(0, OY_LEN-1, 0, OY_LEN);
    f_u_xx_check[OY_LEN] = ff_xx_qrule_Simpson(0, OY_LEN-1, 0, OY_LEN);
    f_u_yy_check[OY_LEN] = ff_yy_qrule_Simpson(0, OY_LEN-1, 0, OY_LEN);

    //</editor-fold>

    fflush(stdout);

    double *maxRes = new double[4];
    maxRes[0] = maxRes[1] = maxRes[2] = maxRes[3] = FLT_MIN;
    double *maxErr = new double[4];
    maxErr[0] = maxErr[1] = maxErr[2] = maxErr[3] = FLT_MIN;
    double *maxErr_l2 = new double[4];
    maxErr_l2[0] = maxErr_l2[1] = maxErr_l2[2] = maxErr_l2[3] = FLT_MIN;
    double *maxDif = new double[4];
    maxDif[0] = maxDif[1] = maxDif[2] = maxDif[3] = FLT_MIN;

    int iter = 0;

    double gamma = 1.e-1;


    maxRes[1] = Residual_u(residual_u, u_exact, u_exact_xx, u_exact_yy, f_u_check, 1);
    maxRes[2] = Residual_u_xx(residual_u_xx, u_exact, u_exact_xx, u_exact_yy, f_u_xx_check, 1);
    maxRes[3] = Residual_u_yy(residual_u_yy, u_exact, u_exact_xx, u_exact_yy, f_u_yy_check, 1);

    printf("%d :  maxRes u = %le   maxRes u_xx = %le   maxRes u_yy = %le\n",
           iter, maxRes[1], maxRes[2], maxRes[3]);
    fflush(stdout);

    print_data_to_files("exact_u", u_exact, u_exact_xx, u_exact_yy, iter);
    print_data_to_files("f_u", f_u_check, f_u_xx_check, f_u_yy_check, iter);
    print_data_to_files("residual_u", residual_u, residual_u_xx, residual_u_yy, iter);

    for (int i = 0; i < OX_LEN_1; ++i) {
        for (int j = 0; j < OY_LEN_1; ++j) {
            residual_u[OY_LEN_1 * i + j] = 0.;
            residual_u_xx[OY_LEN_1 * i + j] = 0.;
            residual_u_yy[OY_LEN_1 * i + j] = 0.;
        }
    }

    maxRes[1] = Residual_u(residual_u, prev_u, prev_u_xx, prev_u_yy, rp_u, 1);
    maxRes[2] = Residual_u_xx(residual_u_xx, prev_u, prev_u_xx, prev_u_yy, rp_u_xx, 1);
    maxRes[3] = Residual_u_yy(residual_u_yy, prev_u, prev_u_xx, prev_u_yy, rp_u_yy, 1);

    while ((maxErr[1] > EPS || maxRes[1] > RES_EPS) && iter < 10 * OX_LEN_1 * OY_LEN_1) {

        iter++;

        // G1 -- (x_i, 0=C) -- bottom boundary
        // G3 -- (x_i, OY_LEN=D) -- top boundary
        for (int i = 1; i < OX_LEN; ++i) {
            u_yy[OY_LEN_1 * i] = prev_u_yy[OY_LEN_1 * i]
                                 - gamma * residual_u_yy[OY_LEN_1 * i];
            u_yy[OY_LEN_1 * i + OY_LEN] = prev_u_yy[OY_LEN_1 * i + OY_LEN]
                                          - gamma * residual_u_yy[OY_LEN_1 * i + OY_LEN];
        }

        // G2 -- (OX_LEN=B, y_j) -- right boundary
        // G4 -- (0=A, y_j) -- left boundary
        for (int j = 1; j < OY_LEN; ++j) {
            u_xx[OY_LEN_1 * OX_LEN + j] = prev_u_xx[OY_LEN_1 * OX_LEN + j]
                                          - gamma * residual_u_xx[OY_LEN_1 * OX_LEN + j];
            u_xx[j] = prev_u_xx[j]
                      - gamma * residual_u_xx[j];
        }

        // inner points
        for (int i = 1; i < OX_LEN; ++i)
            for (int j = 1; j < OY_LEN; ++j) {
                u[OY_LEN_1 * i + j] = prev_u[OY_LEN_1 * i + j]
                                      - gamma * residual_u[OY_LEN_1 * i + j];
                u_xx[OY_LEN_1 * i + j] = prev_u_xx[OY_LEN_1 * i + j]
                                         - gamma * residual_u_xx[OY_LEN_1 * i + j];
                u_yy[OY_LEN_1 * i + j] = prev_u_yy[OY_LEN_1 * i + j]
                                         - gamma * residual_u_yy[OY_LEN_1 * i + j];
            }

        maxDif[0] = maxDif[1] = maxDif[2] = maxDif[3] = FLT_MIN;
        for (int i = 0; i < OX_LEN_1; ++i)
            for (int j = 0; j < OY_LEN_1; ++j) {

                double val = fabs(u[OY_LEN_1 * i + j] - prev_u[OY_LEN_1 * i + j]);
                dif_u[OY_LEN_1 * i + j] = val;
                if (val > maxDif[1]) maxDif[1] = val;

                val = fabs(u_xx[OY_LEN_1 * i + j] - prev_u_xx[OY_LEN_1 * i + j]);
                dif_u_xx[OY_LEN_1 * i + j] = val;
                if (val > maxDif[2]) maxDif[2] = val;

                val = fabs(u_yy[OY_LEN_1 * i + j] - prev_u_yy[OY_LEN_1 * i + j]);
                dif_u_yy[OY_LEN_1 * i + j] = val;
                if (val > maxDif[3]) maxDif[3] = val;
            }

        memcpy(prev_u, u, XY_LEN * sizeof(double));
        memcpy(prev_u_xx, u_xx, XY_LEN * sizeof(double));
        memcpy(prev_u_yy, u_yy, XY_LEN * sizeof(double));

        maxRes[1] = Residual_u(residual_u, prev_u, prev_u_xx, prev_u_yy, rp_u, 1);
        maxRes[2] = Residual_u_xx(residual_u_xx, prev_u, prev_u_xx, prev_u_yy, rp_u_xx, 1);
        maxRes[3] = Residual_u_yy(residual_u_yy, prev_u, prev_u_xx, prev_u_yy, rp_u_yy, 1);

        if (iter % 10 == 0) {
            printf("%d :  maxDif u = %le   maxDif u_xx = %le   maxDif u_yy = %le   "
                           "maxRes u = %le   maxRes u_xx = %le   maxRes u_yy = %le\n",
                   iter, maxDif[1], maxDif[2], maxDif[3], maxRes[1], maxRes[2], maxRes[3]);
            fflush(stdout);
        }
    }



    double *err_u = calc_error_u(HX, HY, u);
    double *err_u_xx = calc_error_u_xx(HX, HY, u_xx);
    double *err_u_yy = calc_error_u_yy(HX, HY, u_yy);

    maxErr[1] = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err_u);
    maxErr[2] = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err_u_xx);
    maxErr[3] = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err_u_yy);

    maxErr_l2[1] = get_l2_norm_vec(OX_LEN_1, OY_LEN_1, err_u);
    maxErr_l2[2] = get_l2_norm_vec(OX_LEN_1, OY_LEN_1, err_u_xx);
    maxErr_l2[3] = get_l2_norm_vec(OX_LEN_1, OY_LEN_1, err_u_yy);

    printf("\n***************************************\n");
           printf("%d :  maxErr u = %le   maxErr u_xx = %le   maxErr u_yy = %le\n",
           iter, maxErr[1], maxErr[2], maxErr[3]);
    printf("%d :  maxRes u = %le   maxRes u_xx = %le   maxRes u_yy = %le\n",
           iter, maxRes[1], maxRes[2], maxRes[3]);
    printf("%d :  maxErr_L2 u = %le   maxErr_L2 u_xx = %le   maxErr_L2 u_yy = %le\n",
           iter, maxErr_l2[1], maxErr_l2[2], maxErr_l2[3]);

    print_data_to_files("u", u, u_xx, u_yy, iter);
    print_data_to_files("err_u", err_u, err_u_xx, err_u_yy, iter);


/*
    printf("tl = %d IterCount = %d Max(Residual) = %le Sum(Rho) = %le Sum(absRho) = %le\n",
           tl, ic, maxRes, calc_array_sum(density, OX_LEN_1, OY_LEN_1, 0),
           calc_array_sum(density, OX_LEN_1, OY_LEN_1, 1));
    fflush(stdout);

    if (tl % 5 == 0) {
        print_data_to_files(phi, density, residual, tl);
        int fixed_x = (int) (get_center_x() / HX);
        int fixed_y = (int) (get_center_y() / HY);
        print_line_along_x("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U_VELOCITY, V_VELOCITY, density, fixed_y);
        print_line_along_y("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                           U_VELOCITY, V_VELOCITY, density, fixed_x);
    }



    double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
    double l1_err_vec = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
    double l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, OX_LEN, OY_LEN, err); // note! a loop boundary
//    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err_vec, l1_err_tr, maxRes, TIME_STEP_CNT);
    extrems = calc_array_extrems(density, OX_LEN_1, OY_LEN_1);
    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err_vec, l1_err_tr, maxRes, extrems,
                      extrems, TIME_STEP_CNT); // !!!!!!!! tmp stab
*/

    delete[] rp_u;
    delete[] rp_u_xx;
    delete[] rp_u_yy;
    delete[] prev_u;
    //delete[] u;
    delete[] prev_u_xx;
    delete[] u_xx;
    delete[] prev_u_yy;
    delete[] u_yy;
    delete[] residual_u;
    delete[] residual_u_xx;
    delete[] residual_u_yy;
    delete[] err_u;
    delete[] err_u_xx;
    delete[] err_u_yy;

    tme = GetTimer() / 1000;
    return u;
}



double *get_exact_solution_3(double hx, double hy) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(analytical_slv(A + hx * i, C + hy * j));
    return res;
}