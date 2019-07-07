//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: main.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 01-Mar-2019 11:25:12
//

//***********************************************************************
// This automatically generated example C main file shows how to call
// entry-point functions that MATLAB Coder generated. You must customize
// this file for your application. Do not modify this file directly.
// Instead, make a copy of this file, modify it, and integrate it into
// your development environment.
//
// This file initializes entry-point function arguments to a default
// size and value before calling the entry-point functions. It does
// not store or use any values returned from the entry-point functions.
// If necessary, it does pre-allocate memory for returned values.
// You can use this file as a starting point for a main function that
// you can deploy in your application.
//
// After you copy the file, and before you deploy it, you must make the
// following changes:
// * For variable-size function arguments, change the example sizes to
// the sizes that your application requires.
// * Change the example values of function arguments to the values that
// your application requires.
// * If the entry-point functions return values, store these values or
// otherwise use them as required by your application.
//
//***********************************************************************
// Include Files
#include "get_inverse.h"
#include "main.h"
#include "get_inverse_terminate.h"
#include "get_inverse_initialize.h"

// Function Declarations
static void argInit_1000x1000_real_T(double result[1000000]);
static double argInit_real_T();
static void main_get_inverse();

// Function Definitions

//
// Arguments    : double result[1000000]
// Return Type  : void
//
static void argInit_1000x1000_real_T(double result[1000000])
{
  int idx0;
  int idx1;

  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < 1000; idx0++) {
    for (idx1 = 0; idx1 < 1000; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result[idx0 + 1000 * idx1] = argInit_real_T();
    }
  }
}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_get_inverse()
{
  static double dv0[1000000];
  static double B[1000000];

  // Initialize function 'get_inverse' input arguments.
  // Initialize function input argument 'A'.
  // Call the entry-point 'get_inverse'.
  argInit_1000x1000_real_T(dv0);
  get_inverse(dv0, B);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // Initialize the application.
  // You do not need to do this more than one time.
  get_inverse_initialize();

  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_get_inverse();

  // Terminate the application.
  // You do not need to do this more than one time.
  get_inverse_terminate();
  return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
