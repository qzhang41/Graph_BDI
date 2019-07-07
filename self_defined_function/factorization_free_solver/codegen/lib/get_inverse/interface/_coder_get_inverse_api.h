/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_get_inverse_api.h
 *
 * MATLAB Coder version            : 4.1
 * C/C++ source code generated on  : 01-Mar-2019 11:25:12
 */

#ifndef _CODER_GET_INVERSE_API_H
#define _CODER_GET_INVERSE_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_get_inverse_api.h"

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void get_inverse(real_T A[1000000], real_T B[1000000]);
extern void get_inverse_api(const mxArray * const prhs[1], int32_T nlhs, const
  mxArray *plhs[1]);
extern void get_inverse_atexit(void);
extern void get_inverse_initialize(void);
extern void get_inverse_terminate(void);
extern void get_inverse_xil_terminate(void);

#endif

/*
 * File trailer for _coder_get_inverse_api.h
 *
 * [EOF]
 */
