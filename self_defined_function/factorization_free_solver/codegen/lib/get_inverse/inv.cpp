//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: inv.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 01-Mar-2019 11:25:12
//

// Include Files
#include "get_inverse.h"
#include "inv.h"
#include "xtrsm.h"
#include "xgetrf.h"

// Function Declarations
static void eml_ipiv2perm(const int ipiv[1000], int perm[1000]);

// Function Definitions

//
// Arguments    : const int ipiv[1000]
//                int perm[1000]
// Return Type  : void
//
static void eml_ipiv2perm(const int ipiv[1000], int perm[1000])
{
  int k;
  int pipk;
  for (k = 0; k < 1000; k++) {
    perm[k] = 1 + k;
  }

  for (k = 0; k < 999; k++) {
    if (ipiv[k] > 1 + k) {
      pipk = perm[ipiv[k] - 1];
      perm[ipiv[k] - 1] = perm[k];
      perm[k] = pipk;
    }
  }
}

//
// Arguments    : const double x[1000000]
//                double y[1000000]
// Return Type  : void
//
void invNxN(const double x[1000000], double y[1000000])
{
  int i0;
  static double b_x[1000000];
  int ipiv[1000];
  int p[1000];
  int k;
  int j;
  int i1;
  int i;
  for (i0 = 0; i0 < 1000000; i0++) {
    y[i0] = 0.0;
    b_x[i0] = x[i0];
  }

  xgetrf(b_x, ipiv);
  eml_ipiv2perm(ipiv, p);
  for (k = 0; k < 1000; k++) {
    i0 = p[k];
    y[k + 1000 * (p[k] - 1)] = 1.0;
    for (j = k + 1; j < 1001; j++) {
      if (y[(j + 1000 * (i0 - 1)) - 1] != 0.0) {
        i1 = j + 1;
        for (i = i1; i < 1001; i++) {
          y[(i + 1000 * (i0 - 1)) - 1] -= y[(j + 1000 * (i0 - 1)) - 1] * b_x[(i
            + 1000 * (j - 1)) - 1];
        }
      }
    }
  }

  xtrsm(b_x, y);
}

//
// File trailer for inv.cpp
//
// [EOF]
//
