//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: get_chi.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 05-Mar-2019 16:54:46
//

// Include Files
#include <cmath>
#include "rt_nonfinite.h"
#include "get_chi.h"
#include "gammaincinv.h"
#include "gammaln.h"

// Function Definitions

//
// Arguments    : double p
//                double n
// Return Type  : double
//
double get_chi(double p, double n)
{
  double X;
  double a;
  double d0;
  double d1;
  a = n / 2.0;
  if ((0.0 <= a) && (!rtIsInf(a)) && (p >= 0.0) && (p <= 1.0)) {
    if ((p > 0.0) && (p < 1.0) && (a > 0.0)) {
      d0 = a;
      gammaln(&d0);
      d1 = a + 1.0;
      gammaln(&d1);
      X = eml_gammaincinv(p, a, std::log(a), d0, d1) * 2.0;
    } else if ((a == 0.0) || (p == 0.0)) {
      X = 0.0;
    } else if (p == 1.0) {
      X = rtInf;
    } else {
      X = rtNaN;
    }
  } else {
    X = rtNaN;
  }

  return X;
}

//
// File trailer for get_chi.cpp
//
// [EOF]
//
