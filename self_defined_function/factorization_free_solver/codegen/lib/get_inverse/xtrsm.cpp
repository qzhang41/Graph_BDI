//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xtrsm.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 01-Mar-2019 11:25:12
//

// Include Files
#include "get_inverse.h"
#include "xtrsm.h"

// Function Definitions

//
// Arguments    : const double A[1000000]
//                double B[1000000]
// Return Type  : void
//
void xtrsm(const double A[1000000], double B[1000000])
{
  int j;
  int jBcol;
  int k;
  int kAcol;
  int i3;
  double d0;
  int i;
  int i4;
  for (j = 0; j < 1000; j++) {
    jBcol = 1000 * j;
    for (k = 999; k >= 0; k--) {
      kAcol = 1000 * k;
      i3 = k + jBcol;
      d0 = B[i3];
      if (d0 != 0.0) {
        B[i3] = d0 / A[k + kAcol];
        for (i = 0; i < k; i++) {
          i4 = i + jBcol;
          B[i4] -= B[i3] * A[i + kAcol];
        }
      }
    }
  }
}

//
// File trailer for xtrsm.cpp
//
// [EOF]
//
