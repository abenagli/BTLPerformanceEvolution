#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include "TMath.h"

// parameters for alpha model (sum of 20 decay times)
// from Sasha: https://indico.cern.ch/event/957609/contributions/4025071/attachments/2113221/3554839/200930_mtd_ledovskoy.pdf
#define alpha_p0 39.17
#define alpha_p1 -3.243
#define alpha_p2 -0.1545
#define alpha_p3 0.01374

#define kB 8.617E-05  // Boltzmann constant

#define minInDay 1440.
#define minInMonth 43200.
#define minInYear 525600.

double funcAlpha(double* x, double* par);

double timeScale(double T_a, double T_r);

double myfunc_DCR(double* x, double* par);

double myfunc_amp(double* x, double* par);

#endif
