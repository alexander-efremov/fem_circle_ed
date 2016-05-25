#include <math.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "consts.h"
#include "timer.h"
#include "utils.h"
#include "common.h"

inline void print_data_to_files(double *phi, double *density, double *residual, int tl) {
//    print_surface("phi", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
//                  U_VELOCITY, V_VELOCITY, phi);
    print_surface("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                  U_VELOCITY, V_VELOCITY, density);
//    print_surface("res", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
//                  U_VELOCITY, V_VELOCITY, residual);
    double *err_lock = calc_error_3(HX, HY, tl * TAU, density);
    print_surface("err-l", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(),
                  TAU, U_VELOCITY, V_VELOCITY, err_lock);
    delete[] err_lock;
}

inline static double func_u(double t, double x, double y) { return (-y+CENTER_OFFSET_Y)*U_VELOCITY; }

inline static double func_v(double t, double x, double y) { return (x-CENTER_OFFSET_X)*V_VELOCITY; }

inline static double analytical_solution_circle(double t, double x, double y) {
    double r = 0.2;
    double da = OMEGA;

    double x0 = get_center_x() - r*sin(t*da);
    double y0 = get_center_y() + r*cos(t*da);
    double value = (x - x0) * (x - x0) + (y - y0) * (y - y0);
    if (value < R_SQ)
        return INN_DENSITY;
    return OUT_DENSITY;
}

static double get_phi_integ_trapezium(int ii, int jj, double *density, double time_value) {
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    double x3 = 0.;
    double y3 = 0.;
    double x4 = 0.;
    double y4 = 0.;

    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);

    double u = func_u(time_value, x1, y1);
    double v = func_v(time_value, x1, y1);
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
    u = func_u(time_value, x4, y4);
    v = func_v(time_value, x4, y4);
    x4 = x4 - TAU * u;
    y4 = y4 - TAU * v;
    if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);

    int nx = IDEAL_SQ_SIZE_X;
    int ny = IDEAL_SQ_SIZE_Y;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

    // get right part for jakoby
    double phi = 0.;
    double mes = x_step * y_step;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            double ideal_x = i * x_step + x_step / 2.;
            double ideal_y = j * y_step + y_step / 2.;

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            // point (x_i,y_j)
            ideal_x = i * x_step;
            ideal_y = j * y_step;

            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;

            if (real_x < A || real_y < C || real_x > B || real_y > D) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : REAL x=%.8le * y=%.8le ** IDEAL x=%.8le * y=%.8le \n",
                       time_value, ii, jj, real_x, real_y, ideal_x, ideal_y);
                printf("1: %.8le * %.8le ** 2: %.8le * %.8le ** 3: %.8le * %.8le ** 4: %.8le * %.8le\n",
                       x1, y1, x2, y2, x3, y3, x4, y4);
            }

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            if (sq_i < 0 || sq_j < 0 || sq_i > OX_LEN - 1 || sq_j > OY_LEN - 1) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : i*=%d j*=%d\n", time_value,
                       ii, jj, sq_i, sq_j);
            }
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            // formula 4
            double dens_1 = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                            + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            // point (x_{i+1},y_j)
            ideal_x = (i + 1) * x_step;
            ideal_y = j * y_step;
            real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                     + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                     + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;
            if (real_x < A || real_y < C || real_x > B || real_y > D) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : REAL x=%.8le * y=%.8le ** IDEAL x=%.8le * y=%.8le \n",
                       time_value, ii, jj, real_x, real_y, ideal_x, ideal_y);
                printf("1: %.8le * %.8le ** 2: %.8le * %.8le ** 3: %.8le * %.8le ** 4: %.8le * %.8le\n",
                       x1, y1, x2, y2, x3, y3, x4, y4);
            }

            // find out in which square real point was placed
            sq_i = (int) ((real_x - A) / HX);
            sq_j = (int) ((real_y - C) / HY);
            if (sq_i < 0 || sq_j < 0 || sq_i > OX_LEN - 1 || sq_j > OY_LEN - 1) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : i*=%d j*=%d\n", time_value,
                       ii, jj, sq_i, sq_j);
            }
            x = A + sq_i * HX;
            y = C + sq_j * HY;

            // formula 4
            double dens_2 = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                            + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);


            // point (x_{i+1},y_{j+1})
            ideal_x = (i + 1) * x_step;
            ideal_y = (j + 1) * y_step;
            real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                     + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                     + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;
            if (real_x < A || real_y < C || real_x > B || real_y > D) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : REAL x=%.8le * y=%.8le ** IDEAL x=%.8le * y=%.8le \n",
                       time_value, ii, jj, real_x, real_y, ideal_x, ideal_y);
                printf("1: %.8le * %.8le ** 2: %.8le * %.8le ** 3: %.8le * %.8le ** 4: %.8le * %.8le\n",
                       x1, y1, x2, y2, x3, y3, x4, y4);
            }

            // find out in which square real point was placed
            sq_i = (int) ((real_x - A) / HX);
            sq_j = (int) ((real_y - C) / HY);
            if (sq_i < 0 || sq_j < 0 || sq_i > OX_LEN - 1 || sq_j > OY_LEN - 1) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : i*=%d j*=%d\n", time_value,
                       ii, jj, sq_i, sq_j);
            }
            x = A + sq_i * HX;
            y = C + sq_j * HY;

            // formula 4
            double dens_3 = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                            + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            // point (x_i,y_{j+1})
            ideal_x = i * x_step;
            ideal_y = (j + 1) * y_step;
            real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                     + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                     + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;
            if (real_x < A || real_y < C || real_x > B || real_y > D) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : REAL x=%.8le * y=%.8le ** IDEAL x=%.8le * y=%.8le \n",
                       time_value, ii, jj, real_x, real_y, ideal_x, ideal_y);
                printf("1: %.8le * %.8le ** 2: %.8le * %.8le ** 3: %.8le * %.8le ** 4: %.8le * %.8le\n",
                       x1, y1, x2, y2, x3, y3, x4, y4);
            }

            // find out in which square real point was placed
            sq_i = (int) ((real_x - A) / HX);
            sq_j = (int) ((real_y - C) / HY);
            if (sq_i < 0 || sq_j < 0 || sq_i > OX_LEN - 1 || sq_j > OY_LEN - 1) {
                printf("Time level %.8le! ERROR INDEX i=%d j=%d : i*=%d j*=%d\n", time_value,
                       ii, jj, sq_i, sq_j);
            }
            x = A + sq_i * HX;
            y = C + sq_j * HY;

            // formula 4
            double dens_4 = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                            + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                            + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            phi += (dens_1 + dens_2 + dens_3 + dens_4) * jakob;
        }
    }

    phi = 0.25 * phi * mes;
    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0;

    return phi;
}

static double get_phi_integ_midpoint(int ii, int jj, double *density, double time_value) {
    double x1 = 0.;
    double y1 = 0.;
    double x2 = 0.;
    double y2 = 0.;
    double x3 = 0.;
    double y3 = 0.;
    double x4 = 0.;
    double y4 = 0.;

    get_coordinates_on_curr(ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);

//    printf("POINT: %d   %d :  x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le **
//                   x4=%.8le * y4=%.8le\n", ii,jj, x1,y1, x2,y2, x3,y3, x4,y4);

    double u = func_u(time_value, x1, y1);
    double v = func_v(time_value, x1, y1);
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
    u = func_u(time_value, x4, y4);
    v = func_v(time_value, x4, y4);
    x4 = x4 - TAU * u;
    y4 = y4 - TAU * v;
    /*
     if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
        || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
        printf("PREV Time level %.8le! ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                       "x4=%.8le * y4=%.8le\n ", time_value, ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
*/
    int nx = IDEAL_SQ_SIZE_X;
    int ny = IDEAL_SQ_SIZE_Y;

    double x_step = 1. / nx;
    double y_step = 1. / ny;

    // get right part for jakoby
    double phi = 0.;
    double mes = x_step * y_step;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {

            double ideal_x = i * x_step + x_step / 2.;
            double ideal_y = j * y_step + y_step / 2.;

            double real_x = x1 + (x2 - x1) * ideal_x + (x4 - x1) * ideal_y
                            + (x1 + x3 - x2 - x4) * ideal_x * ideal_y;
            double real_y = y1 + (y2 - y1) * ideal_x + (y4 - y1) * ideal_y
                            + (y1 + y3 - y2 - y4) * ideal_x * ideal_y;

            // find out in which square real point was placed
            int sq_i = (int) ((real_x - A) / HX);
            int sq_j = (int) ((real_y - C) / HY);
            double x = A + sq_i * HX;
            double y = C + sq_j * HY;

            double a11 = (x2 - x1) + (x1 + x3 - x2 - x4) * ideal_y;
            double a12 = (x4 - x1) + (x1 + x3 - x2 - x4) * ideal_x;
            double a21 = (y2 - y1) + (y1 + y3 - y2 - y4) * ideal_y;
            double a22 = (y4 - y1) + (y1 + y3 - y2 - y4) * ideal_x;
            double jakob = a11 * a22 - a21 * a12;

            // formula 4
            double dens = density[sq_i * OY_LEN_1 + sq_j] * (1 - (real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j] * ((real_x - x) / HX) * (1 - (real_y - y) / HY)
                          + density[(sq_i + 1) * OY_LEN_1 + sq_j + 1] * ((real_x - x) / HX) * ((real_y - y) / HY)
                          + density[sq_i * OY_LEN_1 + sq_j + 1] * (1 - (real_x - x) / HX) * ((real_y - y) / HY);

            phi += mes * dens * jakob;
        }
    }

    if (fabs(phi) < fabs(DBL_MIN_TRIM)) phi = 0;
    return phi;
}

double *solve_3(double &tme) {
    StartTimer();

    fflush(stdout);

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
        if (fabs(prev_density[OY_LEN_1 * i]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i] = 0;

    }

    // G2 -- (OX_LEN=B, y_j) -- right boundary
    for (int j = 1; j < OY_LEN; ++j) {
        prev_density[OY_LEN_1 * OX_LEN + j] = analytical_solution_circle(0., A + HX * OX_LEN, C + HY * j);
        if (fabs(prev_density[OY_LEN_1 * OX_LEN + j]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * OX_LEN + j] = 0;

    }

    // G3 -- (x_i, OY_LEN=D) -- top boundary
    for (int i = 0; i < OX_LEN_1; ++i) {
        prev_density[OY_LEN_1 * i + OY_LEN] = analytical_solution_circle(0., A + HX * i, C + HY * OY_LEN);
        if (fabs(prev_density[OY_LEN_1 * i + OY_LEN]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i + OY_LEN] = 0;

    }

    // G4 -- (0=A, y_j) -- left boundary
    for (int j = 1; j < OY_LEN; ++j) {
        prev_density[j] = analytical_solution_circle(0., A, C + HY * j);
        if (fabs(prev_density[j]) < fabs(DBL_MIN_TRIM)) prev_density[j] = 0;
    }

    memcpy(density, prev_density, XY_LEN * sizeof(double));

    // inner points
    for (int i = 1; i < OX_LEN; ++i) {
        for (int j = 1; j < OY_LEN; ++j) {
            prev_density[OY_LEN_1 * i + j] = analytical_solution_circle(0., A + HX * i, C + HY * j);
            //if (fabs(prev_density[OY_LEN_1 * i + j]) < fabs(DBL_MIN_TRIM)) prev_density[OY_LEN_1 * i + j] = 0;
        }
    }

    //</editor-fold>

    printf("SUM RHO INIT = %le\n", calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 0));
    printf("SUM ABS RHO INIT= %le\n", calc_array_sum(prev_density, OX_LEN_1, OY_LEN_1, 1));
    fflush(stdout);

    double maxRes = FLT_MAX;
    double *extrems = new double[2];

    for (int tl = 1; tl <= TIME_STEP_CNT; tl++) {

        //<editor-fold desc="Calculate phi">

        // with usage of prev_density we calculate phi function values

        // G1 -- (x_i, 0=C) -- bottom boundary
        for (int i = 1; i < OX_LEN; ++i) {
            if (G1[i] == 1) {
                double value = 0.;
                if (INTEGR_TYPE == 1) {
                    value = get_phi_integ_midpoint(i, 0, prev_density, TAU * tl);
                }
                else if (INTEGR_TYPE == 2) {
                    value = get_phi_integ_trapezium(i, 0, prev_density, TAU * tl);
                }
                phi[OY_LEN_1 * i] = value;
            }
        }

        // G2 -- (OX_LEN=B, y_j) -- right boundary
        for (int j = 1; j < OY_LEN; ++j) {
            if (G2[j] == 1) {
                double value = 0.;
                if (INTEGR_TYPE == 1) {
                    value = get_phi_integ_midpoint(OX_LEN, j, prev_density, TAU * tl);
                }
                else if (INTEGR_TYPE == 2) {
                    value = get_phi_integ_trapezium(OX_LEN, j, prev_density, TAU * tl);
                }
                phi[OY_LEN_1 * OX_LEN + j] = value;
            }
        }

        // G3 -- (x_i, OY_LEN=D) -- top boundary
        for (int i = 1; i < OX_LEN; ++i) {
            if (G3[i] == 1) {
                double value = 0.;
                if (INTEGR_TYPE == 1) {
                    value = get_phi_integ_midpoint(i, OY_LEN, prev_density, TAU * tl);
                }
                else if (INTEGR_TYPE == 2) {
                    value = get_phi_integ_trapezium(i, OY_LEN, prev_density, TAU * tl);
                }
                phi[OY_LEN_1 * i + OY_LEN] = value;
            }
        }

        // G4 -- (0=A, y_j) -- left boundary
        for (int j = 1; j < OY_LEN; ++j) {
            if (G4[j] == 1) {
                double value = 0.;
                if (INTEGR_TYPE == 1) {
                    value = get_phi_integ_midpoint(0, j, prev_density, TAU * tl);
                }
                else if (INTEGR_TYPE == 2) {
                    value = get_phi_integ_trapezium(0, j, prev_density, TAU * tl);
                }
                phi[j] = value;
            }
        }

        // point (0.0)
        if (CP00 == 1) {
            double value = 0.;
            if (INTEGR_TYPE == 1) {
                value = get_phi_integ_midpoint(0, 0, prev_density, TAU * tl);
            }
            else if (INTEGR_TYPE == 2) {
                value = get_phi_integ_trapezium(0, 0, prev_density, TAU * tl);
            }
            phi[0] = value;
        }

        // point (1.0)
        if (CP10 == 1) {
            double value = 0.;
            if (INTEGR_TYPE == 1) {
                value = get_phi_integ_midpoint(OX_LEN, 0, prev_density, TAU * tl);
            }
            else if (INTEGR_TYPE == 2) {
                value = get_phi_integ_trapezium(OX_LEN, 0, prev_density, TAU * tl);
            }
            phi[OY_LEN_1 * OX_LEN] = value;
        }

        // point (0.1)
        if (CP01 == 1) {
            double value = 0.;
            if (INTEGR_TYPE == 1) {
                value = get_phi_integ_midpoint(0, OY_LEN, prev_density, TAU * tl);
            }
            else if (INTEGR_TYPE == 2) {
                value = get_phi_integ_trapezium(0, OY_LEN, prev_density, TAU * tl);
            }
            phi[OY_LEN] = value;
        }

        // point (1,1)
        if (CP11 == 1) {
            double value = 0.;
            if (INTEGR_TYPE == 1) {
                value = get_phi_integ_midpoint(OX_LEN, OY_LEN, prev_density, TAU * tl);
            }
            else if (INTEGR_TYPE == 2) {
                value = get_phi_integ_trapezium(OX_LEN, OY_LEN, prev_density, TAU * tl);
            }
            phi[OY_LEN_1 * OX_LEN + OY_LEN] = value;
        }

        // inner points
        for (int i = 1; i < OX_LEN; ++i)
            for (int j = 1; j < OY_LEN; ++j) {
                double value = 0.;
                if (INTEGR_TYPE == 1) {
                    value = get_phi_integ_midpoint(i, j, prev_density, TAU * tl);
                }
                else if (INTEGR_TYPE == 2) {
                    value = get_phi_integ_trapezium(i, j, prev_density, TAU * tl);
                }
                phi[OY_LEN_1 * i + j] = value;
            }

        //</editor-fold>

        ic = 0;
        double maxErr = FLT_MAX;

        while ((maxErr > EPS || maxRes > RES_EPS) && ic < 5 * OX_LEN_1) {

            //<editor-fold desc="Calculate density">

            // G1 -- (x_i, 0=C) -- bottom boundary
            double rpCoef = 32. / (9. * HX * HY);
            for (int i = 1; i < OX_LEN; ++i) {
                if (G1[i] == 1) {
                    density[OY_LEN_1 * i] = -1. / 3. * prev_density[OY_LEN_1 * i + 1]
                                            - 1. / 6. * (prev_density[OY_LEN_1 * (i + 1)] +
                                                         prev_density[OY_LEN_1 * (i - 1)])
                                            - 1. / 18. * (prev_density[OY_LEN_1 * (i + 1) + 1] +
                                                          prev_density[OY_LEN_1 * (i - 1) + 1])
                                            + rpCoef * phi[OY_LEN_1 * i];
                    if (fabs(density[OY_LEN_1 * i]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * i] = 0;
                }
            }

            // G2 -- (OX_LEN=B, y_j) -- right boundary
            for (int j = 1; j < OY_LEN; ++j) {
                if (G2[j] == 1) {
                    density[OY_LEN_1 * OX_LEN + j] = -1. / 3. * prev_density[OY_LEN_1 * (OX_LEN - 1) + j]
                                                     - 1. / 6. * (prev_density[OY_LEN_1 * OX_LEN + j - 1] +
                                                                  prev_density[OY_LEN_1 * OX_LEN + j + 1])
                                                     - 1. / 18. * (prev_density[OY_LEN_1 * (OX_LEN - 1) + j - 1] +
                                                                   prev_density[OY_LEN_1 * (OX_LEN - 1) + j + 1])
                                                     + rpCoef * phi[OY_LEN_1 * OX_LEN + j];
                    if (fabs(density[OY_LEN_1 * OX_LEN + j]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * OX_LEN + j] = 0;
                }
            }

            // G3 -- (x_i, OY_LEN=D) -- top boundary
            for (int i = 1; i < OX_LEN; ++i) {
                if (G3[i] == 1) {
                    density[OY_LEN_1 * i + OY_LEN] = -1. / 3. * prev_density[OY_LEN_1 * i + OY_LEN - 1]
                                                     - 1. / 6. * (prev_density[OY_LEN_1 * (i + 1) + OY_LEN] +
                                                                  prev_density[OY_LEN_1 * (i - 1) + OY_LEN])
                                                     - 1. / 18. * (prev_density[OY_LEN_1 * (i + 1) + OY_LEN - 1] +
                                                                   prev_density[OY_LEN_1 * (i - 1) + OY_LEN - 1])
                                                     + rpCoef * phi[OY_LEN_1 * i + OY_LEN];
                    if (fabs(density[OY_LEN_1 * i + OY_LEN]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * i + OY_LEN] = 0;
                }
            }

            // G4 -- (0=A, y_j) -- left boundary
            for (int j = 1; j < OY_LEN; ++j) {
                if (G4[j] == 1) { // проверить коэф-ты
                    density[j] = -1. / 3. * prev_density[OY_LEN_1 + j]
                                 - 1. / 6. * (prev_density[j + 1] +
                                              prev_density[j - 1])
                                 - 1. / 18. * (prev_density[OY_LEN_1 + j + 1] +
                                               prev_density[OY_LEN_1 + j - 1])
                                 + rpCoef * phi[j];
                    if (fabs(density[j]) < fabs(DBL_MIN_TRIM)) density[j] = 0;
                }
            }

            rpCoef = 64. / (9. * HX * HY);

            // point (0,0)
            if (CP00 == 1) {
                density[0] = -1. / 3. * (prev_density[1]
                                         + prev_density[OY_LEN_1])
                             - 1. / 9. * prev_density[OY_LEN_1 + 1]
                             + rpCoef * phi[0];
                if (fabs(density[0]) < fabs(DBL_MIN_TRIM))
                    density[0] = 0;
            }

            // point (1,0)
            if (CP10 == 1) {
                density[OY_LEN_1 * OX_LEN] = -1. / 3. * (prev_density[OY_LEN_1 * OX_LEN + 1]
                                                         + prev_density[OY_LEN_1 * (OX_LEN - 1)])
                                             - 1. / 9. * prev_density[OY_LEN_1 * (OX_LEN - 1) + 1]
                                             + rpCoef * phi[OY_LEN_1 * OX_LEN];
                if (fabs(density[OY_LEN_1 * OX_LEN]) < fabs(DBL_MIN_TRIM))
                    density[OY_LEN_1 * OX_LEN] = 0;
            }

            // point (0,1)
            if (CP01 == 1) {
                density[OY_LEN] = -1. / 3. * (prev_density[OY_LEN - 1]
                                              + prev_density[OY_LEN_1 + OY_LEN])
                                  - 1. / 9. * prev_density[OY_LEN_1 + OY_LEN - 1]
                                  + rpCoef * phi[OY_LEN];
                if (fabs(density[OY_LEN]) < fabs(DBL_MIN_TRIM))
                    density[OY_LEN] = 0;
            }

            // point (1,1)
            if (CP11 == 1) {
                density[OY_LEN_1 * OX_LEN + OY_LEN] = -1. / 3. * (prev_density[OY_LEN_1 * OX_LEN + OY_LEN - 1]
                                                                  + prev_density[OY_LEN_1 * (OX_LEN - 1) + OY_LEN])
                                                      - 1. / 9. * prev_density[OY_LEN_1 * (OX_LEN - 1) + OY_LEN - 1]
                                                      + rpCoef * phi[OY_LEN_1 * OX_LEN + OY_LEN];
                if (fabs(density[OY_LEN_1 * OX_LEN + OY_LEN]) < fabs(DBL_MIN_TRIM))
                    density[OY_LEN_1 * OX_LEN + OY_LEN] = 0;
            }

            rpCoef = 16. / (9. * HX * HY);
            for (int i = 1; i < OX_LEN; ++i) {
                for (int j = 1; j < OY_LEN; ++j) {
                    density[OY_LEN_1 * i + j] = -1. / 6. * (
                            prev_density[OY_LEN_1 * i + j - 1] + // left
                            prev_density[OY_LEN_1 * (i - 1) + j] + // upper
                            prev_density[OY_LEN_1 * i + j + 1] + // right
                            prev_density[OY_LEN_1 * (i + 1) + j] // bottom
                    ) - 1. / 36. * (
                            prev_density[OY_LEN_1 * (i + 1) + j + 1] + // bottom right
                            prev_density[OY_LEN_1 * (i + 1) + j - 1] + // bottom left
                            prev_density[OY_LEN_1 * (i - 1) + j + 1] + // upper right
                            prev_density[OY_LEN_1 * (i - 1) + j - 1]  // upper left
                    ) + rpCoef * phi[OY_LEN_1 * i + j];
                    if (fabs(density[OY_LEN_1 * i + j]) < fabs(DBL_MIN_TRIM)) density[OY_LEN_1 * i + j] = 0;
                }
            }
            ++ic;

            maxErr = FLT_MIN;
            for (int i = 0; i < OX_LEN_1; ++i)
                for (int j = 0; j < OY_LEN_1; ++j) {
                    double val = fabs(density[i * OY_LEN_1 + j] - prev_density[i * OY_LEN_1 + j]);
                    if (val > maxErr) maxErr = val;
                }

            //</editor-fold>

            //<editor-fold desc="Calculate residual">

            // !!!!!! КАЖЕТСЯ ЧТО НУЖНО ЗАДАВАТЬ rp_coef правильно для каждой из границ !!!!!!!

            rpCoef = HX * HY / 64.;

            // G1 -- (x_i, 0=C) -- bottom boundary
            for (int i = 1; i < OX_LEN; ++i) {
                if (G1[i] == 1) {
                    residual[OY_LEN_1 * i] = rpCoef * (
                            18. * density[OY_LEN_1 * i] +
                            6. * density[OY_LEN_1 * i + 1] +
                            3. * (
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
                            18. * density[OY_LEN_1 * OX_LEN + j] +
                            6. * density[OY_LEN_1 * (OX_LEN - 1) + j] +
                            3. * (
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
                            18. * density[OY_LEN_1 * i + OY_LEN] +
                            6. * density[OY_LEN_1 * i + OY_LEN - 1] +
                            3. * (
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
                            18. * density[j] +
                            6. * density[OY_LEN_1 + j] +
                            3. * (
                                    density[j + 1] +
                                    density[j - 1]
                            ) +
                            density[OY_LEN_1 + j + 1] +
                            density[OY_LEN_1 + j - 1]
                    ) - phi[j];
                }
            }

            // point (0,0)
            if (CP00 == 1) {
                residual[0] = rpCoef * (
                        9. * density[0] +
                        3. * (
                                density[1] +
                                density[OY_LEN_1]
                        ) +
                        density[OY_LEN_1 + 1]
                ) - phi[0];
            }

            // point (1,0)
            if (CP10 == 1) {
                residual[OY_LEN_1 * OX_LEN] = rpCoef * (
                        9. * density[OY_LEN_1 * OX_LEN] +
                        3. * (
                                density[OY_LEN_1 * (OX_LEN - 1)] +
                                density[OY_LEN_1 * OX_LEN + 1]
                        ) +
                        density[OY_LEN_1 * (OX_LEN - 1) + 1]
                ) - phi[OY_LEN_1 * OX_LEN];
            }

            // point (0,1)
            if (CP01 == 1) {
                residual[OY_LEN] = rpCoef * (
                        9. * density[OY_LEN] +
                        3. * (
                                density[OY_LEN - 1] +
                                density[OY_LEN_1 + OY_LEN]
                        ) +
                        density[OY_LEN_1 + OY_LEN - 1]
                ) - phi[OY_LEN];
            }

            // point (1,1)
            if (CP11 == 1) {
                residual[OY_LEN_1 * OX_LEN + OY_LEN] = rpCoef * (
                        9. * density[OY_LEN_1 * OX_LEN + OY_LEN] +
                        3. * (
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
                            36. * density[OY_LEN_1 * i + j] +
                            6. * (
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

            maxRes = FLT_MIN;
            for (int i = 0; i < OX_LEN_1; ++i) {
                for (int j = 0; j < OY_LEN_1; ++j) {
                    double val = fabs(residual[i * OY_LEN_1 + j]);
                    if (val > maxRes) maxRes = val;
                }
            }

            //</editor-fold>

            memcpy(prev_density, density, XY_LEN * sizeof(double));
        }

        printf("tl = %d IterCount = %d Max(Residual) = %le Sum(Rho) = %le Sum(absRho) = %le\n",
               tl, ic, maxRes, calc_array_sum(density, OX_LEN_1, OY_LEN_1, 0),
               calc_array_sum(density, OX_LEN_1, OY_LEN_1, 1));
        fflush(stdout);

        if (tl % 5 == 0) {
            print_data_to_files(phi, density, residual, tl);
            /*int fixed_x = (int) (get_center_x() / HX);
            int fixed_y = (int) (get_center_y() / HY);
            print_line_along_x("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                               U_VELOCITY, V_VELOCITY, density, fixed_y);
            print_line_along_y("rho", OX_LEN, OY_LEN, HX, HY, tl, A, C, get_center_x(), get_center_y(), TAU,
                               U_VELOCITY, V_VELOCITY, density, fixed_x);*/
        }
    }

    double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
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

double *calc_error_3(double hx, double hy, double tt, double *solution) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(solution[i * OY_LEN_1 + j]
                                         - analytical_solution_circle(tt, A + hx * i, C + hy * j));
    return res;
}

double *get_exact_solution_3(double hx, double hy, double t) {
    double *res = new double[XY_LEN];
    for (int i = 0; i < OX_LEN_1; i++)
        for (int j = 0; j < OY_LEN_1; j++)
            res[i * OY_LEN_1 + j] = fabs(analytical_solution_circle(t, A + hx * i, C + hy * j));
    return res;
}