#include <utils.h>
#include <cmath>
#include <common.h>
#include <dirent.h>
#include <fstream>
#include "gtest/gtest.h"
#include <algorithm>
#include <locale>

class FemFixture : public ::testing::Test {
protected:
    virtual void TearDown() {
        if (G1 != NULL)
            delete[] G1;
        if (G2 != NULL)
            delete[] G2;
        if (G3 != NULL)
            delete[] G3;
        if (G4 != NULL)
            delete[] G4;
    }

    virtual void SetUp() {
    }

public:
    FemFixture() : Test() {
    }
};

void print_params() {
    printf("\nOX_LENxOY_LEN = %dx%d\n", OX_LEN, OY_LEN);
    printf("(U, V) = (%le, %le)\n", U_VELOCITY, V_VELOCITY);
    printf("(HX, HY) = (%le, %le)\n", HX, HY);
    printf("TAU = %le\n", TAU);
    printf("TIME_STEP_CNT = %d\n", TIME_STEP_CNT);
    printf("INTEGR_TYPE = %d\n", INTEGR_TYPE);
    printf("IDEAL_SQ_SIZE = %dx%d\n", IDEAL_SQ_SIZE_X, IDEAL_SQ_SIZE_Y);
    printf("CENTER_OFFSET = %le, %le\n", CENTER_OFFSET_X, CENTER_OFFSET_Y);
}

void init_boundary_arrays_and_cp() {
    G1 = new int[OX_LEN_1];
    G2 = new int[OY_LEN_1];
    G3 = new int[OX_LEN_1];
    G4 = new int[OY_LEN_1];
    for (int i = 0; i < OX_LEN_1; ++i) {
        G1[i] = 0;
        G3[i] = 0;
    }
    for (int j = 0; j < OY_LEN_1; ++j) {
        G2[j] = 0;
        G4[j] = 0;
    }
    CP00 = 0;
    CP10 = 0;
    CP01 = 0;
    CP11 = 0;
}

// тестируем третий случай - движение по кругу
TEST_F(FemFixture, test3_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 2.5e-3;

            TIME_STEP_CNT = 2700;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            /*
            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 1;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 1;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
                    G4[j] = 1;
                }
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 0;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);
            printf("G1\n");
            print_vector(G1, OX_LEN_1);
            printf("G2\n");
            print_vector(G2, OY_LEN_1);
            printf("G3\n");
            print_vector(G3, OX_LEN_1);
            printf("G4\n");
            print_vector(G4, OY_LEN_1);

            double *density = solve_3(tme);
            double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_3(HX, HY, 0);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            //           double l1 = get_l1_norm_int_middle(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
             */

        }
    }
}

// тестируем третий случай - движение по кругу
TEST_F(FemFixture, test3_2) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 2; i < 3; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 128 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 128 * (iter + 1);

            CENTER_OFFSET_X = 0.5;
            CENTER_OFFSET_Y = 0.5;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 2.5e-3;

            TIME_STEP_CNT = 50;
            XY_LEN = OX_LEN_1 * OY_LEN_1;
/*

            double x0 = get_center_x();
            double y0 = get_center_y();

            double *exact0 = get_exact_solution_3(HX, HY, 0.);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);
            for (int j = 0; j < TIME_STEP_CNT; ++j) {
                double *exact0 = get_exact_solution_3(HX, HY, j * TAU);
                print_surface("exact", OX_LEN, OY_LEN, HX, HY, j, A, C, x0, y0, TAU, U_VELOCITY,
                              V_VELOCITY, exact0);
            }

            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);


            delete[] exact0;
            delete[] exactT;
            */
        }
    }
}

// тестируем третий случай - движение по кругу
// сходимость
TEST_F(FemFixture, test3_3) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 0; i < 1; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            d = 100.;

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            XY_LEN = OX_LEN_1 * OY_LEN_1;
            //TYPE_EXACT = 1;

            //init_boundary_arrays_and_cp();


            print_params();


            double *u = solve_3(tme);

            /*
            double *err = calc_error_3(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_3(HX, HY, 0);
            double *exactT = get_exact_solution_3(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);

            delete[] exact0;
            delete[] exactT;
            delete[] err;

            */
            delete[] u;
        }
    }
}

// тестируем третий случай - движение по кругу
// сходимость
TEST_F(FemFixture, test4_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 32 * (iter + 1);
            IDEAL_SQ_SIZE_Y = 32 * (iter + 1);

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;

            //TAU = 16. / pow(2., (i + 1));
            TAU = 1.e-3;

            //TIME_STEP_CNT = (int) pow(2., i);
            TIME_STEP_CNT = 2;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;
            for (int i = 1; i < OX_LEN; ++i) {
                G1[i] = 0;
                G3[i] = 1;
            }
            for (int j = 1; j < OY_LEN; ++j) {
                G2[j] = 1;
                G4[j] = 0;
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 1;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);
            printf("G1\n");
            print_vector(G1, OX_LEN_1);
            printf("G2\n");
            print_vector(G2, OY_LEN_1);
            printf("G3\n");
            print_vector(G3, OX_LEN_1);
            printf("G4\n");
            print_vector(G4, OY_LEN_1);

            double *density = solve_4(tme);
            double *err = calc_error_4(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_4(HX, HY, 0);
            double *exactT = get_exact_solution_4(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm(HX, HY, OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}

// тестируем солвер 5 движение из угла в угол
// подход с накоплением коэффициентов
// оказался медленным и размазывает решение
TEST_F(FemFixture, test5_1) {
    double tme = 0.;
    for (int iter = 0; iter < 1; ++iter) {

        double d = 0;
        for (int i = 3; i < 4; ++i) {
            switch (i) {
                case 0:
                    d = 50.;
                    break;
                case 1:
                    d = 100.;
                    break;
                case 2:
                    d = 200.;
                    break;
                case 3:
                    d = 400.;
                    break;
                case 4:
                    d = 800.;
                    break;
                case 5:
                    d = 1600.;
                    break;
                default:
                    return;
            }

            A = 0.;
            B = 1.;
            C = 0.;
            D = 1.;
            R_SQ = 0.099 * 0.099;
            INN_DENSITY = 1.;
            OUT_DENSITY = 0.;

            OX_LEN = (int) d;
            OY_LEN = (int) d;
            OX_LEN_1 = OX_LEN + 1;
            OY_LEN_1 = OY_LEN + 1;
            HX = (B - A) / OX_LEN;
            HY = (D - C) / OY_LEN;
            IDEAL_SQ_SIZE_X = 64;
            IDEAL_SQ_SIZE_Y = 64;

            CENTER_OFFSET_X = 0.3;
            CENTER_OFFSET_Y = 0.3;

            INTEGR_TYPE = 1;

            U_VELOCITY = 1.;
            V_VELOCITY = 1.;
            OMEGA = 1.;
            TAU = 1.2475e-3;

            TIME_STEP_CNT = 322;
            XY_LEN = OX_LEN_1 * OY_LEN_1;

            init_boundary_arrays_and_cp();

            int midIndexX = OX_LEN_1 / 2;
            int midIndexY = OY_LEN_1 / 2;

            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A < .5 && i > 0)
                    G1[i] = 0;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C < .5 && j > 0)
                    G2[j] = 0;
            }
            for (int i = 0; i < OX_LEN_1; ++i) {
                if (i * HX + A > .5 && i < OX_LEN_1 - 1)
                    G3[i] = 0;
            }
            for (int j = 0; j < OY_LEN_1; ++j) {
                if (j * HY + C > .5 && j > 0 && j < OX_LEN_1 - 1) {
                    G4[j] = 0;
                }
            }

            CP00 = 0;
            CP10 = 0;
            CP01 = 0;
            CP11 = 0;

            print_params();
            printf("rel = %le\n", HX / (-HY + 1.));
            printf("midIndexX = %d\n", midIndexX);
            printf("midIndexY = %d\n", midIndexY);

            double *density = solve_5(tme);
            double *err = calc_error_5(HX, HY, TAU * TIME_STEP_CNT, density);
            double *exact0 = get_exact_solution_5(HX, HY, 0);
            double *exactT = get_exact_solution_5(HX, HY, TAU * TIME_STEP_CNT);

            double x0 = get_center_x();
            double y0 = get_center_y();
            print_surface("rho", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, density);
            print_surface("err", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, err);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, 0, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exact0);
            print_surface("exact", OX_LEN, OY_LEN, HX, HY, TIME_STEP_CNT, A, C, x0, y0, TAU, U_VELOCITY,
                          V_VELOCITY, exactT);

            double l1 = get_l1_norm_vec(OX_LEN_1, OY_LEN_1, err);
            double l_inf = get_l_inf_norm(OX_LEN_1, OY_LEN_1, err);
            printf("l1 %le \n", l1);
            printf("l_inf %le\n", l_inf);
            delete[] density;
            delete[] exact0;
            delete[] exactT;
            delete[] err;
        }
    }
}