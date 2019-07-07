//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xgetrf.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 01-Mar-2019 11:25:12
//

// Include Files
#include <cmath>
#include "get_inverse.h"
#include "xgetrf.h"
#include "xswap.h"

// Function Definitions

//
// Arguments    : double A[1000000]
//                int ipiv[1000]
// Return Type  : void
//
void xgetrf(double A[1000000], int ipiv[1000])
{
  int i2;
  int j;
  int b;
  int jj;
  int jp1j;
  int n;
  int jy;
  int ix;
  double smax;
  int jA;
  double s;
  int ijA;
  for (i2 = 0; i2 < 1000; i2++) {
    ipiv[i2] = 1 + i2;
  }

  for (j = 0; j < 999; j++) {
    b = j * 1001;
    jj = j * 1001;
    jp1j = b + 2;
    n = 1000 - j;
    jy = 1;
    ix = b;
    smax = std::abs(A[b]);
    for (jA = 2; jA <= n; jA++) {
      ix++;
      s = std::abs(A[ix]);
      if (s > smax) {
        jy = jA;
        smax = s;
      }
    }

    if (A[(jj + jy) - 1] != 0.0) {
      if (jy - 1 != 0) {
        ipiv[j] = j + jy;
        xswap(A, j + 1, j + jy);
      }

      i2 = (jj - j) + 1000;
      for (jy = jp1j; jy <= i2; jy++) {
        A[jy - 1] /= A[jj];
      }
    }

    n = 998 - j;
    jy = b + 1000;
    jA = jj;
    for (b = 0; b <= n; b++) {
      smax = A[jy];
      if (A[jy] != 0.0) {
        ix = jj + 1;
        i2 = jA + 1002;
        jp1j = (jA - j) + 2000;
        for (ijA = i2; ijA <= jp1j; ijA++) {
          A[ijA - 1] += A[ix] * -smax;
          ix++;
        }
      }

      jy += 1000;
      jA += 1000;
    }
  }
}

//
// File trailer for xgetrf.cpp
//
// [EOF]
//
