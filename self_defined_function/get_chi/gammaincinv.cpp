//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: gammaincinv.cpp
//
// MATLAB Coder version            : 4.1
// C/C++ source code generated on  : 05-Mar-2019 16:54:46
//

// Include Files
#include <cmath>
#include "rt_nonfinite.h"
#include <math.h>
#include "get_chi.h"
#include "gammaincinv.h"
#include "log.h"
#include "sqrt.h"

// Function Declarations
static double PHIinv(double p);
static double eml_PHIc(double z);
static double eml_gammainc(double x, double a, double la, double lgap1,
  boolean_T upper);
static double rt_powd_snf(double u0, double u1);

// Function Definitions

//
// Arguments    : double p
// Return Type  : double
//
static double PHIinv(double p)
{
  double z;
  double r;
  if (p <= 0.0) {
    if (p == 0.0) {
      z = rtMinusInf;
    } else {
      z = rtNaN;
    }
  } else if (p >= 1.0) {
    if (p == 1.0) {
      z = rtInf;
    } else {
      z = rtNaN;
    }
  } else if (std::abs(p - 0.5) <= 0.425) {
    r = 0.180625 - (p - 0.5) * (p - 0.5);
    z = (p - 0.5) * (((((((2509.0809287301227 * r + 33430.575583588128) * r +
                          67265.7709270087) * r + 45921.95393154987) * r +
                        13731.693765509461) * r + 1971.5909503065513) * r +
                      133.14166789178438) * r + 3.3871328727963665) /
      (((((((5226.4952788528544 * r + 28729.085735721943) * r +
            39307.895800092709) * r + 21213.794301586597) * r +
          5394.1960214247511) * r + 687.18700749205789) * r + 42.313330701600911)
       * r + 1.0);
  } else {
    if (p - 0.5 < 0.0) {
      r = std::sqrt(-std::log(p));
    } else {
      r = std::sqrt(-std::log(1.0 - p));
    }

    if (r <= 5.0) {
      r -= 1.6;
      z = (((((((0.00077454501427834139 * r + 0.022723844989269184) * r +
                0.24178072517745061) * r + 1.2704582524523684) * r +
              3.6478483247632045) * r + 5.769497221460691) * r +
            4.6303378461565456) * r + 1.4234371107496835) /
        (((((((1.0507500716444169E-9 * r + 0.00054759380849953455) * r +
              0.015198666563616457) * r + 0.14810397642748008) * r +
            0.6897673349851) * r + 1.6763848301838038) * r + 2.053191626637759) *
         r + 1.0);
    } else {
      r -= 5.0;
      z = (((((((2.0103343992922881E-7 * r + 2.7115555687434876E-5) * r +
                0.0012426609473880784) * r + 0.026532189526576124) * r +
              0.29656057182850487) * r + 1.7848265399172913) * r +
            5.4637849111641144) * r + 6.6579046435011033) /
        (((((((2.0442631033899397E-15 * r + 1.4215117583164459E-7) * r +
              1.8463183175100548E-5) * r + 0.00078686913114561329) * r +
            0.014875361290850615) * r + 0.13692988092273581) * r +
          0.599832206555888) * r + 1.0);
    }

    if (p - 0.5 < 0.0) {
      z = -z;
    }
  }

  return z;
}

//
// Arguments    : double z
// Return Type  : double
//
static double eml_PHIc(double z)
{
  double y;
  double x;
  double t;
  x = 0.70710678118654746 * z;
  t = 3.97886080735226 / (std::abs(x) + 3.97886080735226);
  y = 0.5 * ((((((((((((((((((((((0.0012710976495261409 * (t - 0.5) +
    0.00011931402283834095) * (t - 0.5) + -0.0039638509736051354) * (t - 0.5) +
    -0.00087077963531729586) * (t - 0.5) + 0.0077367252831352668) * (t - 0.5) +
    0.0038333512626488732) * (t - 0.5) + -0.012722381378212275) * (t - 0.5) +
    -0.013382364453346007) * (t - 0.5) + 0.016131532973325226) * (t - 0.5) +
    0.039097684558848406) * (t - 0.5) + 0.0024936720005350331) * (t - 0.5) +
                        -0.0838864557023002) * (t - 0.5) + -0.11946395996432542)
                      * (t - 0.5) + 0.016620792496936737) * (t - 0.5) +
                     0.35752427444953105) * (t - 0.5) + 0.80527640875291062) *
                   (t - 0.5) + 1.1890298290927332) * (t - 0.5) +
                  1.3704021768233816) * (t - 0.5) + 1.313146538310231) * (t -
    0.5) + 1.0792551515585667) * (t - 0.5) + 0.77436819911953858) * (t - 0.5) +
              0.49016508058531844) * (t - 0.5) + 0.27537474159737679) * t * std::
    exp(-x * x);
  if (x < 0.0) {
    y = 1.0 - y;
  }

  return y;
}

//
// Arguments    : double x
//                double a
//                double la
//                double lgap1
//                boolean_T upper
// Return Type  : double
//
static double eml_gammainc(double x, double a, double la, double lgap1,
  boolean_T upper)
{
  double rval;
  double LOGREALMAX;
  double asq;
  double stirlerr;
  static const double dv0[31] = { 0.0, 0.15342640972002736, 0.081061466795327261,
    0.054814121051917651, 0.0413406959554093, 0.033162873519936291,
    0.027677925684998338, 0.023746163656297496, 0.020790672103765093,
    0.018488450532673187, 0.016644691189821193, 0.015134973221917378,
    0.013876128823070748, 0.012810465242920227, 0.01189670994589177,
    0.011104559758206917, 0.010411265261972096, 0.0097994161261588039,
    0.0092554621827127329, 0.0087687001341393862, 0.00833056343336287,
    0.00793411456431402, 0.0075736754879518406, 0.007244554301320383,
    0.00694284010720953, 0.0066652470327076821, 0.0064089941880042071,
    0.0061717122630394576, 0.0059513701127588475, 0.0057462165130101155,
    0.0055547335519628011 };

  double a1;
  double afrac_tmp;
  double xD0;
  double vsq;
  double old;
  double logpax;
  double term;
  double fac;
  int exitg1;
  double n;
  int i0;
  int i;
  LOGREALMAX = 1.7976931348623157E+308;
  b_log(&LOGREALMAX);
  if (!(a > 0.0)) {
    if (a == 0.0) {
      if (x == x) {
        rval = 1.0 - (double)upper;
      } else {
        rval = rtNaN;
      }
    } else {
      rval = rtNaN;
    }
  } else if (!(x > 0.0)) {
    if (x == 0.0) {
      rval = upper;
    } else {
      rval = rtNaN;
    }
  } else if (rtIsInf(a)) {
    if (rtIsInf(x)) {
      rval = rtNaN;
    } else {
      rval = upper;
    }
  } else if (rtIsInf(x)) {
    rval = 1.0 - (double)upper;
  } else {
    if (a <= 15.0) {
      asq = 2.0 * a;
      if (asq == std::floor(asq)) {
        stirlerr = dv0[(int)(asq + 1.0) - 1];
      } else {
        stirlerr = ((lgap1 - (a + 0.5) * la) + a) - 0.91893853320467267;
      }
    } else {
      asq = a * a;
      stirlerr = (0.083333333333333329 + (-0.0027777777777777779 +
        (0.00079365079365079365 + (-0.00059523809523809529 +
        0.00084175084175084171 / asq) / asq) / asq) / asq) / a;
    }

    a1 = a - x;
    afrac_tmp = a + x;
    if (std::abs(a1) > 0.1 * afrac_tmp) {
      if (a < 2.2250738585072014E-308 * x) {
        xD0 = x;
      } else if ((x < 1.0) && (a > 1.7976931348623157E+308 * x)) {
        asq = x;
        b_log(&asq);
        xD0 = (a * la - a * asq) - a;
      } else {
        asq = a / x;
        b_log(&asq);
        xD0 = (a * asq + x) - a;
      }
    } else {
      asq = x / a;
      asq = (1.0 - asq) / (1.0 + asq);
      vsq = asq * asq;
      xD0 = a1 * asq;
      old = xD0;
      term = 2.0 * (a * asq);
      asq = 3.0;
      do {
        exitg1 = 0;
        term *= vsq;
        xD0 += term / asq;
        if (xD0 == old) {
          exitg1 = 1;
        } else {
          old = xD0;
          asq += 2.0;
        }
      } while (exitg1 == 0);
    }

    logpax = (-0.5 * (1.8378770664093453 + la) - stirlerr) - xD0;
    if (x > 1.0E+6) {
      stirlerr = x;
      b_sqrt(&stirlerr);
      asq = std::abs(a1 - 1.0) / stirlerr;
      vsq = asq * asq;
      if (asq < 15.0) {
        a1 = 6.2831853071795862;
        b_sqrt(&a1);
        a1 = eml_PHIc(asq) * a1 * std::exp(0.5 * vsq);
        xD0 = (a1 * ((vsq - 3.0) * asq) - (vsq - 4.0)) / 6.0;
        old = (a1 * ((vsq * vsq - 9.0) * vsq + 6.0) - ((vsq - 1.0) * vsq - 6.0) *
               asq) / 72.0;
        asq = (a1 * (((((5.0 * vsq + 45.0) * vsq - 81.0) * vsq - 315.0) * vsq +
                      270.0) * asq) - ((((5.0 * vsq + 40.0) * vsq - 111.0) * vsq
                 - 174.0) * vsq + 192.0)) / 6480.0;
      } else {
        a1 = (1.0 + (-1.0 + (3.0 - 15.0 / vsq) / vsq) / vsq) / asq;
        xD0 = (1.0 + (-4.0 + (25.0 - 210.0 / vsq) / vsq) / vsq) / vsq;
        old = (1.0 + (-11.0 + (130.0 - 1750.0 / vsq) / vsq) / vsq) / (vsq * asq);
        asq = (1.0 + (-26.0 + (546.0 - 11368.0 / vsq) / vsq) / vsq) / (vsq * vsq);
      }

      if (x < a - 1.0) {
        asq = a * (((a1 / stirlerr - xD0 / x) + old / (x * stirlerr)) - asq / (x
                    * x));
        if (logpax < LOGREALMAX) {
          rval = std::exp(logpax) * asq;
        } else {
          b_log(&asq);
          rval = std::exp(logpax + asq);
        }

        if (upper) {
          rval = 1.0 - rval;
        }
      } else {
        asq = a * (((a1 / stirlerr + xD0 / x) + old / (x * stirlerr)) + asq / (x
                    * x));
        if (logpax < LOGREALMAX) {
          rval = std::exp(logpax) * asq;
        } else {
          b_log(&asq);
          rval = std::exp(logpax + asq);
        }

        if (!upper) {
          rval = 1.0 - rval;
        }
      }
    } else if (upper && (x < 1.0) && (a < 1.0)) {
      fac = a * -x;
      xD0 = fac / (a + 1.0);
      n = 2.0;
      do {
        exitg1 = 0;
        fac = -x * fac / n;
        term = fac / (a + n);
        xD0 += term;
        if (std::abs(term) <= std::abs(xD0) * 2.2204460492503131E-16) {
          exitg1 = 1;
        } else {
          n++;
        }
      } while (exitg1 == 0);

      a1 = x;
      b_log(&a1);
      asq = a * a1 - lgap1;
      vsq = std::exp(asq);
      if (!(vsq == 1.0)) {
        if (vsq - 1.0 == -1.0) {
          asq = -1.0;
        } else {
          a1 = vsq;
          b_log(&a1);
          asq = (vsq - 1.0) * asq / a1;
        }
      }

      rval = -((xD0 + asq) + xD0 * asq);
      if (rval < 0.0) {
        rval = 0.0;
      }
    } else if ((x < a) || (x < 1.0)) {
      n = 1.0;
      if ((!(x < a)) && (a < 2.2250738585072014E-308) && (x >
           1.7976931348623157E+308 * a)) {
        rval = 1.0 - (double)upper;
      } else {
        if (x > 2.2204460492503131E-16 * a) {
          fac = x / a;
          n = 2.0;
          do {
            exitg1 = 0;
            fac = x * fac / (a + (n - 1.0));
            if (fac < 2.2204460492503131E-16) {
              exitg1 = 1;
            } else {
              n++;
            }
          } while (exitg1 == 0);
        }

        asq = 0.0;
        i0 = (int)((1.0 + (-1.0 - (n - 1.0))) / -1.0);
        for (i = 0; i < i0; i++) {
          asq = x * (1.0 + asq) / (a + ((n - 1.0) + -(double)i));
        }

        asq++;
        if (logpax < LOGREALMAX) {
          rval = std::exp(logpax) * asq;
        } else {
          b_log(&asq);
          rval = std::exp(logpax + asq);
        }

        if (rval > 1.0) {
          rval = 1.0;
        }

        if (upper) {
          rval = 1.0 - rval;
        }
      }
    } else {
      fac = 1.0;
      n = 1.0;
      do {
        exitg1 = 0;
        a1 = std::floor(afrac_tmp);
        if (n <= a1) {
          fac = (a - n) * fac / x;
          if (std::abs(fac) < 2.2204460492503131E-16) {
            exitg1 = 1;
          } else {
            n++;
          }
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      if (n <= a1) {
        asq = 1.0;
      } else {
        afrac_tmp = std::floor(a);
        vsq = a - afrac_tmp;
        if (vsq == 0.0) {
          asq = 1.0;
          n = afrac_tmp;
        } else if (vsq == 0.5) {
          a1 = 3.1415926535897931 * x;
          b_sqrt(&a1);
          afrac_tmp = 2.0 * x;
          b_sqrt(&afrac_tmp);
          asq = a1 * std::exp(x) * (2.0 * eml_PHIc(afrac_tmp));
          n = std::floor(a) + 1.0;
        } else {
          xD0 = 1.0;
          a1 = x;
          old = 0.0;
          stirlerr = 1.0;
          fac = 1.0 / x;
          n = 1.0;
          asq = fac;
          term = 0.0;
          while (std::abs(asq - term) > 2.2204460492503131E-16 * asq) {
            term = asq;
            asq = n - vsq;
            xD0 = (a1 + xD0 * asq) * fac;
            old = (stirlerr + old * asq) * fac;
            asq = n * fac;
            a1 = x * xD0 + asq * a1;
            stirlerr = x * old + asq * stirlerr;
            fac = 1.0 / a1;
            asq = stirlerr * fac;
            n++;
          }

          asq *= x;
          n = afrac_tmp + 1.0;
        }
      }

      i0 = (int)((1.0 + (-1.0 - (n - 1.0))) / -1.0);
      for (i = 0; i < i0; i++) {
        asq = 1.0 + (a - ((n - 1.0) + -(double)i)) * asq / x;
      }

      asq = asq * a / x;
      if (logpax < 709.782712893384) {
        rval = std::exp(logpax) * asq;
      } else {
        rval = std::exp(logpax + std::log(asq));
      }

      if (rval > 1.0) {
        rval = 1.0;
      }

      if (!upper) {
        rval = 1.0 - rval;
      }
    }
  }

  return rval;
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d2;
  double d3;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d2 = std::abs(u0);
    d3 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d2 == 1.0) {
        y = 1.0;
      } else if (d2 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d3 == 0.0) {
      y = 1.0;
    } else if (d3 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

//
// Arguments    : double p
//                double a
//                double la
//                double lga
//                double lgap1
// Return Type  : double
//
double eml_gammaincinv(double p, double a, double la, double lga, double lgap1)
{
  double rval;
  boolean_T upper;
  double nu;
  double pLower;
  double log1mpLower;
  double chi2;
  double p1;
  int i;
  boolean_T exitg1;
  double oldz;
  int sgn;
  boolean_T guard1 = false;
  upper = false;
  if (!(a > 0.0)) {
    if ((a == 0.0) && (0.0 <= p) && (p <= 1.0)) {
      rval = 0.0;
    }
  } else if ((0.0 < p) && (p < 1.0)) {
    if (a == rtInf) {
      if (p == 0.0) {
        rval = 0.0;
      } else {
        rval = rtInf;
      }
    } else if (lga == rtInf) {
      rval = a;
    } else {
      nu = 2.0 * a;
      if (p > 0.5) {
        p = 1.0 - p;
        upper = true;
      }

      if (upper) {
        pLower = 1.0 - p;
        if (1.0 - p == 1.0) {
          pLower = 0.9999999999999778;
        }

        log1mpLower = std::log(p);
      } else {
        pLower = p;
        if (1.0 - p != 1.0) {
          log1mpLower = std::log(1.0 - p) * (-p / ((1.0 - p) - 1.0));
        } else {
          log1mpLower = -p;
        }
      }

      if (nu < -1.24 * std::log(pLower)) {
        chi2 = rt_powd_snf(pLower * std::exp(lgap1 + a * 0.693147180559945), 1.0
                           / a);
        if (chi2 < 2.2250738585072014E-306) {
          chi2 = 2.2250738585072014E-306;
        }
      } else if (nu <= 0.32) {
        chi2 = 0.4;
        i = 0;
        exitg1 = false;
        while ((!exitg1) && (i < 200)) {
          oldz = chi2;
          p1 = 1.0 + chi2 * (4.67 + chi2);
          pLower = chi2 * (6.73 + chi2 * (6.66 + chi2));
          chi2 -= (1.0 - std::exp(((log1mpLower + lga) + 0.5 * chi2) + (a - 1.0)
                    * 0.693147180559945) * pLower / p1) / ((-0.5 + (4.67 + 2.0 *
            chi2) / p1) - (6.73 + chi2 * (13.32 + 3.0 * chi2)) / pLower);
          if (std::abs(oldz - chi2) < 0.01 * chi2) {
            exitg1 = true;
          } else {
            i++;
          }
        }
      } else {
        p1 = 0.222222 / nu;
        chi2 = (PHIinv(pLower) * std::sqrt(p1) + 1.0) - p1;
        chi2 *= nu * chi2 * chi2;
        if (chi2 > 2.2 * nu + 6.0) {
          chi2 = -2.0 * ((log1mpLower - (a - 1.0) * std::log(0.5 * chi2)) + lga);
        }
      }

      rval = 0.5 * chi2;
      pLower = rtInf;
      oldz = rtInf;
      if (p > 1.0021E-294) {
        nu = 2.2204460492503131E-14 * p;
      } else {
        nu = 2.2251089859537388E-308;
      }

      if (upper) {
        sgn = -1;
      } else {
        sgn = 1;
      }

      chi2 = 0.0;
      log1mpLower = 1.7976931348623157E+308;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 1000)) {
        p1 = (double)sgn * (eml_gammainc(rval, a, la, lgap1, upper) - p);
        if ((p1 * pLower < 0.0) && (std::abs(pLower) <= std::abs(p1))) {
          rval = 0.5 * rval + 0.5 * oldz;
          p1 = (double)sgn * (eml_gammainc(rval, a, la, lgap1, upper) - p);
        }

        if (p1 > 0.0) {
          log1mpLower = rval;
        } else {
          chi2 = rval;
        }

        if ((std::abs(p1) < nu) || (std::abs(rval - oldz) <
             2.2204460492503131E-16 * rval + 2.2250738585072014E-308)) {
          exitg1 = true;
        } else {
          oldz = rval;
          pLower = p1;
          guard1 = false;
          if (i < 500) {
            rval *= 1.0 - p1 / (rval * std::exp(((a - 1.0) * std::log(rval) -
              rval) - lga) + p1 * ((rval + 1.0) - a) / 2.0);
            if (rval <= chi2) {
              if (chi2 == 0.0) {
                if (std::abs((double)upper - p) < std::abs(eml_gammainc
                     (2.2250738585072014E-308, a, la, lgap1, upper) - p)) {
                  rval = 0.0;
                  exitg1 = true;
                } else {
                  chi2 = 2.2250738585072014E-308;
                  guard1 = true;
                }
              } else {
                guard1 = true;
              }
            } else {
              if (rval >= log1mpLower) {
                rval = 0.01 * chi2 + 0.99 * log1mpLower;
              }

              i++;
            }
          } else {
            if (1.0E+8 * chi2 < log1mpLower) {
              oldz = 1.0E+8 * chi2;
              pLower = (double)sgn * (eml_gammainc(oldz, a, la, lgap1, upper) -
                p);
              if (pLower > 0.0) {
                log1mpLower = oldz;
              } else {
                chi2 = oldz;
              }
            }

            rval = 0.5 * chi2 + 0.5 * log1mpLower;
            i++;
          }

          if (guard1) {
            rval = 0.99 * chi2 + 0.01 * log1mpLower;
            i++;
          }
        }
      }
    }
  } else if (p == 0.0) {
    rval = 0.0;
  } else {
    if (p == 1.0) {
      rval = rtInf;
    }
  }

  return rval;
}

//
// File trailer for gammaincinv.cpp
//
// [EOF]
//
