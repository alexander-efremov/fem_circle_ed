#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "consts.h"
#include "timer.h"
#include "utils.h"
#include "common.h"

inline void print_data_to_files(double *phi, double *density, double *residual, int tl) {
    print_surface("phi", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                  U_VELOCITY, V_VELOCITY, phi);
    print_surface("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                  U_VELOCITY, V_VELOCITY, density);
    print_surface("res", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                  U_VELOCITY, V_VELOCITY, residual);
    double *err_lock = calc_error_4(HX, HY, tl * TAU, density);
    print_surface("err-l", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(),
                  TAU, U_VELOCITY, V_VELOCITY, err_lock);
    delete[] err_lock;
}

inline static double func_u(double t, double x, double y) { return U_VELOCITY; }

inline static double func_v(double t, double x, double y) { return  V_VELOCITY; }

inline static double analytical_solution_circle(double t, double x, double y) {
    double x0 = get_center_x() + t * func_u(t, x, y);
    double y0 = get_center_y() + t * func_v(t, x, y);
    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    if (value <= R_SQ) return INN_DENSITY;
    return OUT_DENSITY;
}

static double get_phi_integ_midpoint(double x0, double y0, double x1, double y1,
                                     double x2, double y2, double x3, double y3,
                                     double *density, double time_value, int ind_omega, int ii, int jj) {
    double phi=0.;

    double u = func_u(time_value, x0, y0);
    double v = func_v(time_value, x0, y0);
    x0 = x0 - TAU * u;
    y0 = y0 - TAU * v;

    u = func_u(time_value, x1, y1);
    v = func_v(time_value, x1, y1);
    x1 = x1 - TAU * u;
    y1 = y1 - TAU * v;

    u = func_u(time_value, x2, y2);
    v = func_v(time_value, x2, y2);
    x2 = x2 - TAU * u;
    y2 = y2 - TAU * v;

    u = func_u(time_value, x3, y3);
    v = func_v(time_value, x3, y3);
    x3 = x3 - TAU * u;
    y3 = y3 - TAU * v;

    if (x0 < A || x0 > B || x1 < A || x1 > B || x2 < A || x2 > B || x3 < A || x3 > B
            || y0 < A || y0 > B || y1 < C || y1 > D || y2 < C || y2 > D || y3 < C || y3 > D) {
        phi = 0.; // !!! must work up Gamma_in
        printf("PREV TL %.8le POINT ERROR: (x0=%.8le %.8le; y0=%.8le %.8le) ** (x1=%.8le %.8le; y1=%.8le %.8le) ** "
                       "(x2=%.8le %.8le; y2=%.8le %.8le) ** (x3=%.8le %.8le; y3=%.8le %.8le)\n ",
                   time_value, x0 + TAU * u, x0, y0 + TAU * v, y0,
               x1 + TAU * u, x1, y1 + TAU * v, y1,
               x2 + TAU * u, x2, y2 + TAU * v, y2,
               x3 + TAU * u, x3, y3 + TAU * v, y3);
        printf("PREV TL %d, %d POINT ERROR! ind_omega= %d\n", ii, jj, ind_omega);
        return phi;
    }

    int nx = IDEAL_SQ_SIZE_X;
    int ny = IDEAL_SQ_SIZE_Y;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

    double mes = x_step * y_step;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            double ideal_x = i * x_step + x_step / 2.;
            double ideal_y = j * y_step + y_step / 2.;

            double Psi=0;
            switch (ind_omega) {
                case 0:
                    Psi = (1 - ideal_x) * (1 - ideal_y);
                    break;
                case 1:
                    Psi = ideal_x * (1 - ideal_y);
                    break;
                case 2:
                    Psi = ideal_x * ideal_y;
                    break;
                case 3:
                    Psi = (1 - ideal_x) * ideal_y;
                    break;
            }

            double a11 = (x1 - x0) + (x0 + x2 - x1 - x3) * ideal_y;
            double a12 = (x3 - x0) + (x0 + x2 - x1 - x3) * ideal_x;
            double a21 = (y1 - y0) + (y0 + y2 - y1 - y3) * ideal_y;
            double a22 = (y3 - y0) + (y0 + y2 - y1 - y3) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            double real_x = x0 + (x1 - x0) * ideal_x + (x3 - x0) * ideal_y
                                + (x0 + x2 - x1 - x3) * ideal_x * ideal_y;
            double real_y = y0 + (y1 - y0) * ideal_x + (y3 - y0) * ideal_y
                                + (y0 + y2 - y1 - y3) * ideal_x * ideal_y;

            // find out in which square omega_{ij} (i=sq_i, j=sq_j) real point (real_x, real_y) was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            // find density at the point (real_x, real_y)
            double dens = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            phi += mes * dens * jakob * Psi;
        }
    }

    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0.;
    return phi;

}

static double get_entire_phi(int ii, int jj, double *density, double time_value, int ind_QR) {
    double phi=0.;
    double x0, y0, x1, y1, x2, y2, x3, y3;

    // it is not calculated at Gamma_in
    if (ii==0 && jj>0 && jj<OY_LEN && G4[jj]==0)  {
        return phi; }
    if (ii==OX_LEN && jj>0 && jj<OY_LEN && G2[jj]==0)  {
        return phi; }
    if (jj==0 && ii>0 && ii<OX_LEN && G1[ii]==0)  {
        return phi; }
    if (jj==OY_LEN && ii>0 && ii<OX_LEN && G3[ii]==0)  {
        return phi; }
    if (ii==0 && jj==0 && CP00==0)  {
        return phi; }
    if (ii==OX_LEN && jj==0 && CP10==0)  {
        return phi; }
    if (ii==OX_LEN && jj==OY_LEN && CP11==0)  {
        return phi; }
    if (ii==0 && jj==OY_LEN && CP01==0)  {
        return phi; }

    double val0, val1, val2, val3;

    //inner point
    if (ii>0 && ii<OX_LEN && jj>0 && jj<OY_LEN) {
        // omega_ij = [x_i, x_{i+1}] x [y_j, y_{j+1}]; ind_omega=0
        if (
                ( (jj+1==OY_LEN) && (ii+1<OX_LEN) && (G3[ii]==0 || G3[ii+1]==0) ) ||
                ( (ii+1==OX_LEN) && (jj+1<OY_LEN) && (G2[jj]==0 || G2[jj+1]==0) ) ||
                ( (ii+1==OX_LEN) && (jj+1==OY_LEN) && (CP11==0))
           ) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_i, y_j)
            x0 = A + ii * HX;
            y0 = C + jj * HY;
            // (x_{i+1}, y_j)
            x1 = A + (ii + 1) * HX;
            y1 = C + jj * HY;
            // (x_{i+1}, y_{j+1})
            x2 = A + (ii + 1) * HX;
            y2 = C + (jj + 1) * HY;
            // (x_i, y_{j+1})
            x3 = A + ii * HX;
            y3 = C + (jj + 1) * HY;

//          printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val0 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 0, ii, jj);
        }

        // omega_{i-1,j} = [x_{i-1}, x_i] x [y_j, y_{j+1}]; ind_omega=1
        if (
                ( (jj+1==OY_LEN) && (ii-1>0) && (G3[ii]==0 || G3[ii-1]==0) ) ||
                ( (ii-1==0) && (jj+1<OY_LEN) && (G4[jj]==0 || G4[jj+1]==0) ) ||
                ( (ii-1==0) && (jj+1==OY_LEN) && (CP01==0))
                ) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_{i-1}, y_j)
            x0 = A + (ii - 1) * HX;
            y0 = C + jj * HY;
            // (x_i, y_j)
            x1 = A + ii * HX;
            y1 = C + jj * HY;
            // (x_i, y_{j+1})
            x2 = A + ii * HX;
            y2 = C + (jj + 1) * HY;
            // (x_{i-1), y_{j+1})
            x3 = A + (ii - 1) * HX;
            y3 = C + (jj + 1) * HY;

//          printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val1 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 1, ii, jj);
        }

        // omega_{i-1,j-1} = [x_{i-1}, x_i] x [y_{j-1}, y_j]; ind_omega=2
        if (
                ( (jj-1==0) && (ii-1>0) && (G1[ii]==0 || G1[ii-1]==0) ) ||
                ( (ii-1==0) && (jj-1>0) && (G4[jj]==0 || G4[jj-1]==0) ) ||
                ( (ii-1==0) && (jj-1==0) && (CP00==0))
                ) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_{i-1}, y_{j-1})
            x0 = A + (ii - 1) * HX;
            y0 = C + (jj - 1) * HY;
            // (x_i, y_{j-1})
            x1 = A + ii * HX;
            y1 = C + (jj - 1) * HY;
            // (x_i, y_j)
            x2 = A + ii * HX;
            y2 = C + jj * HY;
            // (x_{i-1), y_j)
            x3 = A + (ii - 1) * HX;
            y3 = C + jj * HY;

//          printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val2 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 2, ii, jj);
        }

        // omega_{i,j-1} = [x_i, x_{i+1}] x [y_{j-1}, y_j]; ind_omega=3
        if (
                ( (jj-1==0) && (ii+1<OX_LEN) && (G1[ii]==0 || G1[ii+1]==0) ) ||
                ( (ii+1==OX_LEN) && (jj-1>0) && (G2[jj]==0 || G2[jj-1]==0) ) ||
                ( (ii+1==OX_LEN) && (jj-1==0) && (CP10==0))
                ) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_i, y_{j-1})
            x0 = A + ii * HX;
            y0 = C + (jj - 1) * HY;
            // (x_{i+1}, y_{j-1})
            x1 = A + (ii + 1) * HX;
            y1 = C + (jj - 1) * HY;
            // (x_{i+1), y_j)
            x2 = A + (ii + 1) * HX;
            y2 = C + jj * HY;
            // (x_i, y_j)
            x3 = A + ii * HX;
            y3 = C + jj * HY;


//          printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val3 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 3, ii, jj);

            phi = val0 + val1 + val2 + val3;
            return phi;
        }
    }

    else if (ii==0 && jj>0 && jj<OY_LEN && G4[jj]==1)  { // Gamma_4 \ Gamma_in
        // omega_ij = [x_i, x_{i+1}] x [y_j, y_{j+1}]; ind_omega=0
        if ((jj+1==OY_LEN) && (CP01==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_i, y_j)
            x0 = A + ii * HX;
            y0 = C + jj * HY;
            // (x_{i+1}, y_j)
            x1 = A + (ii + 1) * HX;
            y1 = C + jj * HY;
            // (x_{i+1}, y_{j+1})
            x2 = A + (ii + 1) * HX;
            y2 = C + (jj + 1) * HY;
            // (x_i, y_{j+1})
            x3 = A + ii * HX;
            y3 = C + (jj + 1) * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val0 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 0, ii, jj);
        }

        // omega_{i,j-1} = [x_i, x_{i+1}] x [y_{j-1}, y_j]; ind_omega=3
        if ((jj-1==0) && (CP00==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_i, y_{j-1})
            x0 = A + ii * HX;
            y0 = C + (jj - 1) * HY;
            // (x_{i+1}, y_{j-1})
            x1 = A + (ii + 1) * HX;
            y1 = C + (jj - 1) * HY;
            // (x_{i+1), y_j)
            x2 = A + (ii + 1) * HX;
            y2 = C + jj * HY;
            // (x_i, y_j)
            x3 = A + ii * HX;
            y3 = C + jj * HY;


//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val3 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 3, ii, jj);
        }

        phi = val0 + val3;
        return phi;

    }
    else if (jj==0 && ii>0 && ii<OX_LEN && G1[ii]==1) { // Gamma_1 \ Gamma_in
        // omega_ij = [x_i, x_{i+1}] x [y_j, y_{j+1}]; ind_omega=0
        if ((ii+1==OX_LEN) && (CP10==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_i, y_j)
            x0 = A + ii * HX;
            y0 = C + jj * HY;
            // (x_{i+1}, y_j)
            x1 = A + (ii + 1) * HX;
            y1 = C + jj * HY;
            // (x_{i+1}, y_{j+1})
            x2 = A + (ii + 1) * HX;
            y2 = C + (jj + 1) * HY;
            // (x_i, y_{j+1})
            x3 = A + ii * HX;
            y3 = C + (jj + 1) * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val0 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 0, ii, jj);
        }

        // omega_{i-1,j} = [x_{i-1}, x_i] x [y_j, y_{j+1}]; ind_omega=1
        if ((ii-1==0) && (CP00==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_{i-1}, y_j)
            x0 = A + (ii - 1) * HX;
            y0 = C + jj * HY;
            // (x_i, y_j)
            x1 = A + ii * HX;
            y1 = C + jj * HY;
            // (x_i, y_{j+1})
            x2 = A + ii * HX;
            y2 = C + (jj + 1) * HY;
            // (x_{i-1), y_{j+1})
            x3 = A + (ii - 1) * HX;
            y3 = C + (jj + 1) * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val1 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 1, ii, jj);
        }

        phi = val0 + val1;
        return phi;
    }
    else if (ii==OX_LEN && jj>0 && jj<OY_LEN && G2[jj]==1) { // Gamma_2 \ Gamma_in
        // omega_{i-1,j} = [x_{i-1}, x_i] x [y_j, y_{j+1}]; ind_omega=1
        if ((jj+1==OY_LEN) && (CP11==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_{i-1}, y_j)
            x0 = A + (ii - 1) * HX;
            y0 = C + jj * HY;
            // (x_i, y_j)
            x1 = A + ii * HX;
            y1 = C + jj * HY;
            // (x_i, y_{j+1})
            x2 = A + ii * HX;
            y2 = C + (jj + 1) * HY;
            // (x_{i-1), y_{j+1})
            x3 = A + (ii - 1) * HX;
            y3 = C + (jj + 1) * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val1 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 1, ii, jj);
        }

        // omega_{i-1,j-1} = [x_{i-1}, x_i] x [y_{j-1}, y_j]; ind_omega=2
        if ((jj-1==0) && (CP10==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_{i-1}, y_{j-1})
            x0 = A + (ii - 1) * HX;
            y0 = C + (jj - 1) * HY;
            // (x_i, y_{j-1})
            x1 = A + ii * HX;
            y1 = C + (jj - 1) * HY;
            // (x_i, y_j)
            x2 = A + ii * HX;
            y2 = C + jj * HY;
            // (x_{i-1), y_j)
            x3 = A + (ii - 1) * HX;
            y3 = C + jj * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val2 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 2, ii, jj);
        }

        phi = val1 + val2;
        return phi;
    }
    else if (jj==OY_LEN && ii>0 && ii<OX_LEN && G3[ii]==1) { // Gamma_3 \ Gamma_in
        // omega_{i-1,j-1} = [x_{i-1}, x_i] x [y_{j-1}, y_j]; ind_omega=2
        if ((ii-1==0) && (CP01==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_{i-1}, y_{j-1})
            x0 = A + (ii - 1) * HX;
            y0 = C + (jj - 1) * HY;
            // (x_i, y_{j-1})
            x1 = A + ii * HX;
            y1 = C + (jj - 1) * HY;
            // (x_i, y_j)
            x2 = A + ii * HX;
            y2 = C + jj * HY;
            // (x_{i-1), y_j)
            x3 = A + (ii - 1) * HX;
            y3 = C + jj * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val2 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 2, ii, jj);
        }

        // omega_{i,j-1} = [x_i, x_{i+1}] x [y_{j-1}, y_j]; ind_omega=3
        if ((ii+1==OX_LEN) && (CP11==0)) {
            return phi; // stubs at Gamma_in
        }
        else {
            // (x_i, y_{j-1})
            x0 = A + ii * HX;
            y0 = C + (jj - 1) * HY;
            // (x_{i+1}, y_{j-1})
            x1 = A + (ii + 1) * HX;
            y1 = C + (jj - 1) * HY;
            // (x_{i+1), y_j)
            x2 = A + (ii + 1) * HX;
            y2 = C + jj * HY;
            // (x_i, y_j)
            x3 = A + ii * HX;
            y3 = C + jj * HY;


//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
            val3 = get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 3, ii, jj);
        }

        phi = val2 + val3;
        return phi;
    }
    else if (ii==0 && jj==0 && CP00==1) { // (0,0) \ Gamma_in
        // omega_ij = [x_i, x_{i+1}] x [y_j, y_{j+1}]; ind_omega=0
        // (x_i, y_j)
        x0 = A + ii * HX;
        y0 = C + jj * HY;
        // (x_{i+1}, y_j)
        x1 = A + (ii + 1) * HX;
        y1 = C + jj * HY;
        // (x_{i+1}, y_{j+1})
        x2 = A + (ii + 1) * HX;
        y2 = C + (jj + 1) * HY;
        // (x_i, y_{j+1})
        x3 = A + ii * HX;
        y3 = C + (jj + 1) * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
        val0=get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 0, ii, jj);

        phi = val0;
        return phi;
    }
    else if (ii==OX_LEN && jj==0 && CP10==1) { // (1,0) \ Gamma_in
        // omega_{i-1,j} = [x_{i-1}, x_i] x [y_j, y_{j+1}]; ind_omega=1
        // (x_{i-1}, y_j)
        x0 = A + (ii - 1) * HX;
        y0 = C + jj * HY;
        // (x_i, y_j)
        x1 = A + ii * HX;
        y1 = C + jj * HY;
        // (x_i, y_{j+1})
        x2 = A + ii * HX;
        y2 = C + (jj + 1) * HY;
        // (x_{i-1), y_{j+1})
        x3 = A + (ii - 1) * HX;
        y3 = C + (jj + 1) * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
        val1=get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 1, ii, jj);

        phi = val1;
        return phi;
    }
    else if (ii==OX_LEN && jj==OY_LEN && CP11==1) { // (1,1) \ Gamma_in
        // omega_{i-1,j-1} = [x_{i-1}, x_i] x [y_{j-1}, y_j]; ind_omega=2
        // (x_{i-1}, y_{j-1})
        x0 = A + (ii - 1) * HX;
        y0 = C + (jj - 1) * HY;
        // (x_i, y_{j-1})
        x1 = A + ii * HX;
        y1 = C + (jj - 1) * HY;
        // (x_i, y_j)
        x2 = A + ii * HX;
        y2 = C + jj * HY;
        // (x_{i-1), y_j)
        x3 = A + (ii - 1) * HX;
        y3 = C + jj * HY;

//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
        val2=get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 2, ii, jj);

        phi = val2;
        return phi;
    }
    else if (ii==0 && jj==OY_LEN && CP01==1) { // (0,1) \ Gamma_in
        // omega_{i,j-1} = [x_i, x_{i+1}] x [y_{j-1}, y_j]; ind_omega=3
        // (x_i, y_{j-1})
        x0 = A + ii * HX;
        y0 = C + (jj - 1) * HY;
        // (x_{i+1}, y_{j-1})
        x1 = A + (ii + 1) * HX;
        y1 = C + (jj - 1) * HY;
        // (x_{i+1), y_j)
        x2 = A + (ii + 1) * HX;
        y2 = C + jj * HY;
        // (x_i, y_j)
        x3 = A + ii * HX;
        y3 = C + jj * HY;


//    printf("POINT: (x0=%.8le ; y0=%.8le) * (x1=%.8le ; y1=%.8le) * (x2=%.8le ; y2=%.8le) * "
//                   "(x3=%.8le ; y3=%.8le)\n", x0,y0, x1,y1, x2,y2, x3,y3);
        val3=get_phi_integ_midpoint(x0, y0, x1, y1, x2, y2, x3, y3, density, time_value, 3, ii, jj);

        phi = val3;
        return phi;
    }

    else {
        printf("INDEX ERROR in PHI! i=%d   j=%d\n", ii, jj);
        return phi;
    }
}

double Residual(double *density, double *phi, double *residual) {

    double rpCoef = HX * HY / 36.;

    // G1 -- (x_i, 0=C) -- bottom boundary
    for (int i = 1; i < OX_LEN; ++i) {
        if (G1[i] == 1) {
            residual[OY_LEN_1 * i] = rpCoef * (
                    8. * density[OY_LEN_1 * i] +
                    4. * density[OY_LEN_1 * i + 1] +
                    2. * (
                            density[OY_LEN_1 * (i - 1)] +
                            density[OY_LEN_1 * (i + 1)]
                    ) +
                    density[OY_LEN_1 * (i - 1) + 1] +
                    density[OY_LEN_1 * (i + 1) + 1]
            ) - phi[OY_LEN_1 * i];
        }
    }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    for (int j = 1; j < OY_LEN; ++j) {
        if (G2[j] == 1) {
            residual[OY_LEN_1 * OX_LEN + j] = rpCoef * (
                    8. * density[OY_LEN_1 * OX_LEN + j] +
                    4. * density[OY_LEN_1 * (OX_LEN - 1) + j] +
                    2. * (
                            density[OY_LEN_1 * OX_LEN + j - 1] +
                            density[OY_LEN_1 * OX_LEN + j + 1]
                    ) +
                    density[OY_LEN_1 * (OX_LEN - 1) + j - 1] +
                    density[OY_LEN_1 * (OX_LEN - 1) + j + 1]
            ) - phi[OY_LEN_1 * OX_LEN + j];
        }
    }

    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 1; i < OX_LEN; ++i) {
        if (G3[i] == 1) {
            residual[OY_LEN_1 * i + OY_LEN] = rpCoef * (
                    8. * density[OY_LEN_1 * i + OY_LEN] +
                    4. * density[OY_LEN_1 * i + OY_LEN - 1] +
                    2. * (
                            density[OY_LEN_1 * (i + 1) + OY_LEN] +
                            density[OY_LEN_1 * (i - 1) + OY_LEN]
                    ) +
                    density[OY_LEN_1 * (i + 1) + OY_LEN - 1] +
                    density[OY_LEN_1 * (i - 1) + OY_LEN - 1]
            ) - phi[OY_LEN_1 * i + OY_LEN];
        }
    }

    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 1; j < OY_LEN; ++j) {
        if (G4[j] == 1) {
            residual[j] = rpCoef * (
                    8. * density[j] +
                    4. * density[OY_LEN_1 + j] +
                    2. * (density[j + 1] + density[j - 1]) +
                    density[OY_LEN_1 + j + 1] + density[OY_LEN_1 + j - 1]
            ) - phi[j];
        }
    }

    // point (0,0)
    if (CP00 == 1) {
        residual[0] = rpCoef * (
                4. * density[0] +
                2. * (density[1] + density[OY_LEN_1]) +
                density[OY_LEN_1 + 1]
        ) - phi[0];
    }

    // point (1,0)
    if (CP10 == 1) {
        residual[OY_LEN_1 * OX_LEN] = rpCoef * (
                4. * density[OY_LEN_1 * OX_LEN] +
                2. * (
                        density[OY_LEN_1 * (OX_LEN - 1)] +
                        density[OY_LEN_1 * OX_LEN + 1]
                ) +
                density[OY_LEN_1 * (OX_LEN - 1) + 1]
        ) - phi[OY_LEN_1 * OX_LEN];
    }

    // point (0,1)
    if (CP01 == 1) {
        residual[OY_LEN] = rpCoef * (
                4. * density[OY_LEN] +
                2. * (density[OY_LEN - 1] + density[OY_LEN_1 + OY_LEN]) +
                density[OY_LEN_1 + OY_LEN - 1]
        ) - phi[OY_LEN];
    }

    // point (1,1)
    if (CP11 == 1) {
        residual[OY_LEN_1 * OX_LEN + OY_LEN] = rpCoef * (
                4. * density[OY_LEN_1 * OX_LEN + OY_LEN] +
                2. * (
                        density[OY_LEN_1 * OX_LEN + OY_LEN - 1] +
                        density[OY_LEN_1 * (OX_LEN - 1) + OY_LEN]
                ) +
                density[OY_LEN_1 * (OX_LEN - 1) + OY_LEN - 1]
        ) - phi[OY_LEN_1 * OX_LEN + OY_LEN];
    }

    // inner points
    for (int i = 1; i < OX_LEN; ++i) {
        for (int j = 1; j < OY_LEN; ++j) {
            residual[OY_LEN_1 * i + j] = rpCoef * (
                    16. * density[OY_LEN_1 * i + j] +
                    4. * (
                            density[OY_LEN_1 * i + j - 1] + // left
                            density[OY_LEN_1 * (i - 1) + j] + // upper
                            density[OY_LEN_1 * i + j + 1] + // right
                            density[OY_LEN_1 * (i + 1) + j] // bottom
                    ) +
                    density[OY_LEN_1 * (i + 1) + j + 1] + // bottom right
                    density[OY_LEN_1 * (i + 1) + j - 1] + // bottom left
                    density[OY_LEN_1 * (i - 1) + j + 1] + // upper right
                    density[OY_LEN_1 * (i - 1) + j - 1]  // upper left
            ) - phi[OY_LEN_1 * i + j];
        }
    }

    double maxRes = FLT_MIN;
    for (int i = 0; i < OX_LEN_1; ++i) {
        for (int j = 0; j < OY_LEN_1; ++j) {
            double val = fabs(residual[i * OY_LEN_1 + j]);
            if (val > maxRes) maxRes = val;
        }
    }

    return maxRes;
}

double *solve_4(double &tme) {
    StartTimer();

    fflush(stdout);

    double gamma = 1./(36*HX*HY);

    int ic = 0;
    double *phi = new double[XY_LEN];
    double *prev_density = new double[XY_LEN];
    double *density = new double[XY_LEN];
    double *residual = new double[XY_LEN];

    //<editor-fold desc="Fill initial data">

    for (int i = 0; i < OX_LEN_1; ++i) {
        for (int j = 0; j < OY_LEN_1; ++j) {
            density[OY_LEN_1 * i + j] = 0.;
            prev_density[OY_LEN_1 * i + j] = 0.;
            residual[OY_LEN_1 * i + j] = 0.;
            phi[OY_LEN_1 * i + j] = 0.;
        }
    }

    // G1 -- (x_i, 0=C) -- bottom boundary
    for (int i = 0; i < OX_LEN_1; ++i) {
        prev_density[OY_LEN_1 * i] = analytical_solution_circle(0., A + HX * i, C);
        if (fabs(prev_density[OY_LEN_1 * i]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i] = 0.;

    }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    for (int j = 1; j < OY_LEN; ++j) {
        prev_density[OY_LEN_1 * OX_LEN + j] = analytical_solution_circle(0., A + HX * OX_LEN, C + HY * j);
        if (fabs(prev_density[OY_LEN_1 * OX_LEN + j]) < fabs(DBL_MIN_TRIM))
            prev_density[OY_LEN_1 * OX_LEN + j] = 0.;

    }

    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 0; i < OX_LEN_1; ++i) {
        prev_density[OY_LEN_1 * i + OY_LEN] = analytical_solution_circle(0., A + HX * i, C + HY * OY_LEN);
        if (fabs(prev_density[OY_LEN_1 * i + OY_LEN]) < fabs(DBL_MIN_TRIM))
            prev_density[OY_LEN_1 * i + OY_LEN] = 0.;

    }

    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 1; j < OY_LEN; ++j) {
        prev_density[j] = analytical_solution_circle(0., A, C + HY * j);
        if (fabs(prev_density[j]) < fabs(DBL_MIN_TRIM)) prev_density[j] = 0.;
    }

    memcpy(density, prev_density, XY_LEN * sizeof(double));

    // inner points
    for (int i = 1; i < OX_LEN; ++i) {
        for (int j = 1; j < OY_LEN; ++j) {
            prev_density[OY_LEN_1 * i + j] = analytical_solution_circle(0., A + HX * i, C + HY * j);
            if (fabs(prev_density[OY_LEN_1 * i + j]) < fabs(DBL_MIN_TRIM))
                prev_density[OY_LEN_1 * i + j] = 0.;
        }
    }

    //</editor-fold>

    printf("SUM RHO INIT = %le\n", calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 0));
    printf("SUM ABS RHO INIT= %le\n", calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 1));
    fflush(stdout);

    double *extrems = new double[2];
    double maxRes = FLT_MIN;

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {

        //<editor-fold desc="Calculate phi">

        // with usage of prev_density we calculate phi function values

        for (int i = 0; i < OX_LEN_1; ++i)
            for (int j = 0; j < OY_LEN_1; ++j) {
                double value = 0.;
                if (INTEGR_TYPE == 1) {
                    value = get_entire_phi(i, j, prev_density, TAU * tl, 1);
                }
                else if (INTEGR_TYPE == 2) {
                   value = get_entire_phi(i, j, prev_density, TAU * tl, 2);
                }
                phi[OY_LEN_1 * i + j] = value;
            }


        //</editor-fold>

        ic = 0;
        double maxDiff = FLT_MIN;

        maxRes = Residual(prev_density, phi, residual);

        printf("ITER STEP = %d Max(Residual) = %le Sum(Rho) = %le Sum(absRho) = %le\n",
               ic, maxRes, calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 0),
               calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 1));

        print_data_to_files(phi, prev_density, residual, 0);

        while ((maxDiff > EPS || maxRes > RES_EPS) && ic < 1000) {

            //<editor-fold desc="Calculate density">

            for (int i = 0; i < OX_LEN_1; ++i) {
                for (int j = 0; j < OY_LEN_1; ++j) {
                    density[i * OY_LEN_1 + j]=prev_density[i * OY_LEN_1 + j]
                            - gamma * residual[i * OY_LEN_1 + j];
                    if (fabs(density[i * OY_LEN_1 + j]) < fabs(DBL_MIN_TRIM))
                        density[i * OY_LEN_1 + j] = 0.;
                }
            }

            maxDiff = FLT_MIN;
            for (int i = 0; i < OX_LEN_1; ++i)
                for (int j = 0; j < OY_LEN_1; ++j) {
                    double val = fabs(density[i * OY_LEN_1 + j] - prev_density[i * OY_LEN_1 + j]);
                    if (val > maxDiff) maxDiff = val;
                }

            maxRes = Residual(density, phi, residual);

            ++ic;

            memcpy(prev_density, density, XY_LEN * sizeof(double));

            printf("ITER STEP = %d Max(Diff) = %le Max(Residual) = %le Sum(Rho) = %le Sum(absRho) = %le\n",
                   ic, maxDiff, maxRes, calc_array_sum(density, OX_LEN_1, OY_LEN_1, 0),
                   calc_array_sum(density, OX_LEN_1, OY_LEN_1, 1));
        }

        printf("tl = %d IterCount = %d Max(Residual) = %le Sum(Rho) = %le Sum(absRho) = %le\n",
               tl, ic, maxRes, calc_array_sum(density, OX_LEN_1, OY_LEN_1, 0),
               calc_array_sum(density, OX_LEN_1, OY_LEN_1, 1));
        fflush(stdout);

        if (tl % 1 == 0) {
            print_data_to_files(phi, density, residual, tl);
            /*int fixed_x = (int) (get_center_x() / HX);
            int fixed_y = (int) (get_center_y() / HY);
            print_line_along_x("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                               U_VELOCITY, V_VELOCITY, density, fixed_y);
            print_line_along_y("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                               U_VELOCITY, V_VELOCITY, density, fixed_x);*/
        }
    }

    double *err = calc_error_4(HX, HY, TAU * TIME_STEP_CNT, density);
    double l1_err_vec = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
    double l1_err_tr = get_l1_norm_int_trapezoidal(HX, HY, OX_LEN, OY_LEN, err); // note! a loop boundary
//    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err_vec, l1_err_tr, maxRes, TIME_STEP_CNT);
    extrems = calc_array_extrems(density, OX_LEN_1, OY_LEN_1);
    append_statistics(OX_LEN_1, OY_LEN_1, TAU, ic, l1_err_vec, l1_err_tr, maxRes, extrems,
                      extrems, TIME_STEP_CNT); // !!!!!!!! tmp stab


    delete[] prev_density;
    delete[] phi;
    delete[] err;
    delete[] residual;
    delete[] extrems;
    tme = GetTimer() / 1000;
    return density;
}


double *calc_error_4(double hx, double hy, double tt, double *solution) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(solution[i * OY_LEN_1 + j]
                                         - analytical_solution_circle(tt, A + hx * i, C + hy * j));
    return res;
}

double *get_exact_solution_4(double hx, double hy, double t) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = analytical_solution_circle(t, A + hx * i, C + hy * j);
    return res;
}