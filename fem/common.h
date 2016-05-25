#ifndef FEM_CIRCLE_COMMON_H
#define FEM_CIRCLE_COMMON_H

#include "consts.h"

inline double get_center_x() {
    return A + CENTER_OFFSET_X;
}

inline double get_center_y() {
    return C + CENTER_OFFSET_Y;
}

inline void get_coordinates_on_curr(
        int ii, int jj,
        double &x1, double &y1,
        double &x2, double &y2,
        double &x3, double &y3,
        double &x4, double &y4
) {
    // проверить условия вылета
    if (ii > 0 && ii < OX_LEN && jj > 0 && jj < OY_LEN) {
        // p1 (x_{i-1/2}, y_{j-1/2})
        x1 = A + ii * HX - HX / 2.;
        y1 = C + jj * HY - HY / 2.;
        // p2 (x_{i+1/2}, y_{j-1/2})
        x2 = A + ii * HX + HX / 2.;
        y2 = C + jj * HY - HY / 2.;
        // p3 (x_{i+1/2}, y_{j+1/2})
        x3 = A + ii * HX + HX / 2.;
        y3 = C + jj * HY + HY / 2.;
        // p4 (x_{i-1/2}, y_{j+1/2})
        x4 = A + ii * HX - HX / 2.;
        y4 = C + jj * HY + HY / 2.;
        if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
            || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
            printf("1. Inner point, ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii == OX_LEN && jj == 0) { // point (1,0)  omega_{i-1,j}
        // p1 (x_{OX_LEN-1/2}, C)
        x1 = B - HX / 2.;
        y1 = C;
        // p2 (B, C)
        x2 = B;
        y2 = C;
        // p3 (B, y_{1/2})
        x3 = B;
        y3 = C + HY / 2.;
        // p4 (x_{OX_LEN-1/2}, y_{1/2})
        x4 = B - HX / 2.;
        y4 = C + HY / 2.;
        if (x1 <= A || x1 >= B || x2 <= A || x2 > B || x3 <= A || x3 > B || x4 <= A || x4 >= B
            || y1 < C || y1 >= D || y2 < C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
            printf("2. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii == 0 && jj == OY_LEN) { // point (0,1)  omega_{i,j-1}
        // p1 (A, y_{OY_LEN-1/2})
        x1 = A;
        y1 = D - HY / 2.;
        // p2 (x_{1/2}, y_{OY_LEN-1/2})
        x2 = A + HX / 2.;
        y2 = D - HY / 2.;
        // p3 (x_{1/2}, D)
        x3 = A + HX / 2.;
        y3 = D;
        // p4 (A, D)
        x4 = A;
        y4 = D;
        if (x1 < A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 < A || x4 >= B
            || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 > D || y4 <= C || y4 > D)
            printf("3. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii == 0 && jj == 0) { // point (0,0)  omega_{i,j}
        // p1 (A, C)
        x1 = A;
        y1 = C;
        // p2 (x_{1/2}, C)
        x2 = A + HX / 2.;
        y2 = C;
        // p3 (x_{1/2}, y_{1/2})
        x3 = A + HX / 2.;
        y3 = C + HY / 2.;
        // p4 (A, y_{1/2})
        x4 = A;
        y4 = C + HY / 2.;
        if (x1 < A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 < A || x4 >= B
            || y1 < C || y1 >= D || y2 < C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
            printf("4. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii == OX_LEN && jj == OY_LEN) { // point (1,1)  omega_{i-1,j-1}
        // p1 (x_{OX_LEN-1/2}, y_{OY_LEN-1/2})
        x1 = B - HX / 2.;
        y1 = D - HY / 2.;
        // p2 (B, y_{OY_LEN-1/2})
        x2 = B;
        y2 = D - HY / 2.;
        // p3 (B, D)
        x3 = B;
        y3 = D;
        // p4 (x_{OX_LEN-1/2}, D)
        x4 = B - HX / 2.;
        y4 = D;
        if (x1 <= A || x1 >= B || x2 <= A || x2 > B || x3 <= A || x3 > B || x4 <= A || x4 >= B
            || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 > D || y4 <= C || y4 > D)
            printf("5. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii > 0 && ii < OX_LEN && jj == 0) { // G1 -- bottom boundary
        // p1 (x_{i-1/2}, C)
        x1 = A + ii * HX - HX / 2.;
        y1 = C;
        // p2 (x_{i+1/2}, C)
        x2 = A + ii * HX + HX / 2.;
        y2 = C;
        // p3 (x_{i+1/2}, y_{1/2})
        x3 = A + ii * HX + HX / 2.;
        y3 = C + HY / 2.;
        // p4 (x_{i-1/2}, y_{1/2})
        x4 = A + ii * HX - HX / 2.;
        y4 = C + HY / 2.;
        if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
            || y1 < C || y1 >= D || y2 < C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
            printf("6. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii == OX_LEN && jj > 0 && jj < OY_LEN) { // G2 -- right boundary
        // p1 (x_{OX_LEN-1/2}, y_{j-1/2})
        x1 = B - HX / 2.;
        y1 = C + jj * HY - HY / 2.;
        // p2 (B, y_{j-1/2})
        x2 = B;
        y2 = C + jj * HY - HY / 2.;
        // p3 (B, y_{j+1/2})
        x3 = B;
        y3 = C + jj * HY + HY / 2.;
        // p4 (x_{OX_LEN-1/2}, y_{j+1/2})
        x4 = B - HX / 2.;
        y4 = C + jj * HY + HY / 2.;
        if (x1 <= A || x1 >= B || x2 <= A || x2 > B || x3 <= A || x3 > B || x4 <= A || x4 >= B
            || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
            printf("7. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (jj == OY_LEN && ii > 0 && ii < OX_LEN) { // G3 -- top boundary
        // p1 (x_{i-1/2}, y_{OY_LEN-1/2})
        x1 = A + ii * HX - HX / 2.;
        y1 = D - HY / 2.;
        // p2 (x_{i+1/2}, y_{OY_LEN-1/2})
        x2 = A + ii * HX + HX / 2.;
        y2 = D - HY / 2.;
        //p3 (x_{i+1/2}, D)
        x3 = A + ii * HX + HX / 2.;
        y3 = D;
        //p4 (x_{i-1/2}, D)
        x4 = A + ii * HX - HX / 2.;
        y4 = D;
        if (x1 <= A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 <= A || x4 >= B
            || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 > D || y4 <= C || y4 > D)
            printf("8. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else if (ii == 0 && jj > 0 && jj < OY_LEN) { // G4 -- left boundary
        // p1 (A, y_{j-1/2})
        x1 = A;
        y1 = C + jj * HY - HY / 2.;
        // p2 (x_{1/2}, y_{j-1/2})
        x2 = A + HX / 2.;
        y2 = C + jj * HY - HY / 2.;
        //p3 (x_{1/2}, y_{j+1/2})
        x3 = A + HX / 2.;
        y3 = C + jj * HY + HY / 2.;
        //p4 (A, y_{j+1/2})
        x4 = A;
        y4 = C + jj * HY + HY / 2.;
        if (x1 < A || x1 >= B || x2 <= A || x2 >= B || x3 <= A || x3 >= B || x4 < A || x4 >= B
            || y1 <= C || y1 >= D || y2 <= C || y2 >= D || y3 <= C || y3 >= D || y4 <= C || y4 >= D)
            printf("9. ERROR INDEX i=%d j=%d : x1=%.8le * y1=%.8le ** x2=%.8le * y2=%.8le ** x3=%.8le * y3=%.8le ** "
                           "x4=%.8le * y4%.8le\n ", ii, jj, x1, y1, x2, y2, x3, y3, x4, y4);
    }
    else {
        printf("ERROR! INDEX i=%d j=%d ", ii, jj);
    }
}

double *solve_3(double &tme);

double *calc_error_3(double hx, double hy, double tt, double *solution);

double *get_exact_solution_3(double hx, double hy, double t);

double *solve_4(double &tme);

double *calc_error_4(double hx, double hy, double tt, double *solution);

double *get_exact_solution_4(double hx, double hy, double t);

double *solve_5(double &tme);

double *calc_error_5(double hx, double hy, double tt, double *solution);

double *get_exact_solution_5(double hx, double hy, double t);

#endif //FEM_CIRCLE_COMMON_H
