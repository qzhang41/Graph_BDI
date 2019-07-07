//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xswap.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 01-Mar-2019 11:25:12
//

// Include Files
#include "get_inverse.h"
#include "xswap.h"

// Function Definitions

//
// Arguments    : double x[1000000]
//                int ix0
//                int iy0
// Return Type  : void
//
void xswap(double x[1000000], int ix0, int iy0)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 0; k < 1000; k++) {
    temp = x[ix];
    x[ix] = x[iy];
    x[iy] = temp;
    ix += 1000;
    iy += 1000;
  }
}

//
// File trailer for xswap.cpp
//
// [EOF]
//
