#include "consts.h"
#include <algorithm>

#ifdef WIN32
__int64 timerStart = 0;
#else
struct timeval timerStart;
#endif

double A = 0.;
double B = 0.;
double C = 0.;
double D = 0.;
int OX_LEN = 0;
int OY_LEN = 0;
int OX_LEN_1 = 0;
int OY_LEN_1 = 0;
int XY_LEN = 0;
double TAU = 0.;
int TIME_STEP_CNT = 0;
int JAK_ITER_CNT = 0;
double HX = 0.;
double HY = 0.;
double R_SQ = 0.;
double INN_DENSITY = 0.;
double OUT_DENSITY = 0.;
double U_VELOCITY = 0.;
double V_VELOCITY = 0.;
double DBL_MIN_TRIM = 1.e-16;
double RES_EPS = 1.e-14;
double EPS = 1.e-8;
int INTEGR_TYPE = 1;
int IDEAL_SQ_SIZE_X = 0;
int IDEAL_SQ_SIZE_Y = 0;
double CENTER_OFFSET_X = 0.;
double CENTER_OFFSET_Y = 0.;
double OMEGA = 0.;
int* G1 = NULL;
int* G2 = NULL;
int* G3 = NULL;
int* G4 = NULL;
int CP00 = 0;
int CP10 = 0;
int CP11 = 0;
int CP01 = 0;
