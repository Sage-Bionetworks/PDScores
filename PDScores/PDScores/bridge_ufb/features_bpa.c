/*
 * features_bpa.c
 *
 * Code generation for function 'features_bpa'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "extractaudiophon.h"
#include "fastdfa.h"
#include "features_bga.h"
#include "features_bpa.h"
#include "features_bta.h"
#include "features_bvav2.h"
#include "features_ufb.h"
#include "lomb.h"
#include "mfcc.h"
#include "swipep.h"
#include "vadsplitphon.h"
#include "bridge_ufb_emxutil.h"
#include "polyfit.h"
#include "log10.h"
#include "sortrows.h"
#include "log2.h"
#include "quantile.h"
#include "sqrt.h"
#include "sum.h"
#include "power.h"
#include "repmat.h"
#include "diff.h"
#include "combine_vector_elements.h"
#include "bridge_ufb_rtwutil.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Definitions */
void features_bpa(const emxArray_real_T *post, double ft[3])
{
  int vstride;
  emxArray_real_T *t;
  int idx;
  double dT;
  int ixstart;
  int n;
  int ix;
  boolean_T exitg3;
  emxArray_boolean_T *x;
  int ii_data[1];
  int ii;
  boolean_T exitg2;
  int istart_data[1];
  boolean_T exitg1;
  unsigned int iend_data[1];
  emxArray_real_T *posttrim;
  double b_x[3];
  double c_x[3];
  emxArray_real_T *d_x;
  emxArray_real_T *b_t;
  emxArray_real_T *dt;
  int k;
  int npages;
  double xlast;
  emxArray_real_T *b;
  emxArray_real_T *r11;
  double nscales;
  emxArray_real_T *ypts;
  emxArray_real_T *e_x;
  emxArray_real_T *b_dt;
  emxArray_real_T *c_dt;
  double coeffs[2];

  /*  Computes basic posture test features. */
  /*  Inputs: */
  /*   post - posture accelerometry vector: post(:,1) - time points, */
  /*          post(:,2:4) - X,Y,Z acceleration data */
  /*  */
  /*  (CC BY-SA 3.0) Max Little, 2014 */
  /*  Output feature vector */
  for (vstride = 0; vstride < 3; vstride++) {
    ft[vstride] = rtNaN;
  }

  /*  Ignore zero-length inputs */
  if (post->size[0] == 0) {
  } else {
    b_emxInit_real_T(&t, 1);

    /*  Calculate relative time */
    idx = post->size[0];
    dT = post->data[0];
    vstride = t->size[0];
    t->size[0] = idx;
    emxEnsureCapacity((emxArray__common *)t, vstride, (int)sizeof(double));
    for (vstride = 0; vstride < idx; vstride++) {
      t->data[vstride] = post->data[vstride] - dT;
    }

    /*  Ignore posture tests which do not contain enough data */
    ixstart = 1;
    n = t->size[0];
    dT = t->data[0];
    if (t->size[0] > 1) {
      if (rtIsNaN(t->data[0])) {
        ix = 2;
        exitg3 = false;
        while ((!exitg3) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(t->data[ix - 1])) {
            dT = t->data[ix - 1];
            exitg3 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < t->size[0]) {
        while (ixstart + 1 <= n) {
          if (t->data[ixstart] > dT) {
            dT = t->data[ixstart];
          }

          ixstart++;
        }
      }
    }

    if (dT < 19.0) {
    } else {
      emxInit_boolean_T(&x, 1);

      /*  Trim away start/end of test */
      vstride = x->size[0];
      x->size[0] = t->size[0];
      emxEnsureCapacity((emxArray__common *)x, vstride, (int)sizeof(boolean_T));
      idx = t->size[0];
      for (vstride = 0; vstride < idx; vstride++) {
        x->data[vstride] = (t->data[vstride] >= 3.0);
      }

      idx = 0;
      n = 1;
      ii = 1;
      exitg2 = false;
      while ((!exitg2) && (ii <= x->size[0])) {
        if (x->data[ii - 1]) {
          idx = 1;
          ii_data[0] = ii;
          exitg2 = true;
        } else {
          ii++;
        }
      }

      if (idx == 0) {
        n = 0;
      }

      vstride = 0;
      while (vstride <= n - 1) {
        istart_data[0] = ii_data[0];
        vstride = 1;
      }

      vstride = x->size[0];
      x->size[0] = t->size[0];
      emxEnsureCapacity((emxArray__common *)x, vstride, (int)sizeof(boolean_T));
      idx = t->size[0];
      for (vstride = 0; vstride < idx; vstride++) {
        x->data[vstride] = (t->data[vstride] <= 19.0);
      }

      idx = 0;
      n = 1;
      ii = x->size[0];
      exitg1 = false;
      while ((!exitg1) && (ii > 0)) {
        if (x->data[ii - 1]) {
          idx = 1;
          ii_data[0] = ii;
          exitg1 = true;
        } else {
          ii--;
        }
      }

      emxFree_boolean_T(&x);
      if (idx == 0) {
        n = 0;
      }

      vstride = 0;
      while (vstride <= n - 1) {
        iend_data[0] = (unsigned int)ii_data[0];
        vstride = 1;
      }

      if ((int)iend_data[0] < istart_data[0]) {
      } else {
        if (istart_data[0] > (int)iend_data[0]) {
          vstride = 1;
          ii = 1;
        } else {
          vstride = istart_data[0];
          ii = (int)iend_data[0] + 1;
        }

        emxInit_real_T(&posttrim, 2);
        ixstart = posttrim->size[0] * posttrim->size[1];
        posttrim->size[0] = ii - vstride;
        posttrim->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)posttrim, ixstart, (int)sizeof
                          (double));
        for (ixstart = 0; ixstart < 3; ixstart++) {
          idx = ii - vstride;
          for (ix = 0; ix < idx; ix++) {
            posttrim->data[ix + posttrim->size[0] * ixstart] = post->data
              [((vstride + ix) + post->size[0] * (1 + ixstart)) - 1];
          }
        }

        if (istart_data[0] > (int)iend_data[0]) {
          ixstart = 0;
          ix = 0;
        } else {
          ixstart = istart_data[0] - 1;
          ix = (int)iend_data[0];
        }

        if (0 == ii - vstride) {
          n = -1;
        } else if (ii - vstride > 3) {
          n = (ii - vstride) - 1;
        } else {
          n = 2;
        }

        dT = t->data[((ixstart + ix) - ixstart) - 1] - t->data[ixstart];

        /*  Orientation */
        combine_vector_elements(posttrim, b_x);
        idx = ii - vstride;

        /*  Orientation-corrected force signals */
        for (vstride = 0; vstride < 3; vstride++) {
          c_x[vstride] = b_x[vstride] / (double)idx;
        }

        emxInit_real_T(&d_x, 2);
        repmat(c_x, n + 1, d_x);
        vstride = posttrim->size[0] * posttrim->size[1];
        posttrim->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)posttrim, vstride, (int)sizeof
                          (double));
        idx = posttrim->size[0];
        ii = posttrim->size[1];
        idx *= ii;
        for (vstride = 0; vstride < idx; vstride++) {
          posttrim->data[vstride] -= d_x->data[vstride];
        }

        b_emxInit_real_T(&b_t, 1);

        /*  Scaled velocity signals */
        vstride = b_t->size[0];
        b_t->size[0] = ix - ixstart;
        emxEnsureCapacity((emxArray__common *)b_t, vstride, (int)sizeof(double));
        idx = ix - ixstart;
        for (vstride = 0; vstride < idx; vstride++) {
          b_t->data[vstride] = t->data[ixstart + vstride];
        }

        b_emxInit_real_T(&dt, 1);
        diff(b_t, dt);
        dt->data[dt->size[0]] = dt->data[dt->size[0] - 1];
        b_repmat(dt, d_x);
        vstride = d_x->size[0] * d_x->size[1];
        d_x->size[0] = posttrim->size[0];
        d_x->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)d_x, vstride, (int)sizeof(double));
        idx = posttrim->size[0] * posttrim->size[1];
        emxFree_real_T(&b_t);
        for (vstride = 0; vstride < idx; vstride++) {
          d_x->data[vstride] *= posttrim->data[vstride];
        }

        if (d_x->size[0] != 1) {
          idx = 0;
        } else {
          idx = 1;
        }

        n = d_x->size[idx];
        if ((!(d_x->size[0] == 0)) && (d_x->size[idx] > 1)) {
          vstride = 1;
          k = 1;
          while (k <= idx) {
            vstride *= d_x->size[0];
            k = 2;
          }

          npages = 1;
          k = idx + 2;
          while (k < 3) {
            npages *= 3;
            k = 3;
          }

          ix = 0;
          for (idx = 1; idx <= npages; idx++) {
            ixstart = ix;
            for (ii = 1; ii <= vstride; ii++) {
              ixstart++;
              ix = ixstart;
              xlast = d_x->data[ixstart - 1];
              for (k = 0; k <= n - 2; k++) {
                ix += vstride;
                xlast += d_x->data[ix - 1];
                d_x->data[ix - 1] = xlast;
              }
            }
          }
        }

        emxInit_real_T(&b, 2);
        emxInit_real_T(&r11, 2);

        /*  Average scaled power X,Y,Z */
        power(d_x, b);
        vstride = r11->size[0] * r11->size[1];
        r11->size[0] = b->size[0];
        r11->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)r11, vstride, (int)sizeof(double));
        idx = b->size[0] * b->size[1];
        for (vstride = 0; vstride < idx; vstride++) {
          r11->data[vstride] = 35.0 * b->data[vstride];
        }

        emxFree_real_T(&b);
        sum(r11, b_x);
        emxFree_real_T(&r11);
        for (vstride = 0; vstride < 3; vstride++) {
          b_x[vstride] /= dT;
        }

        dT = b_x[0];
        for (k = 0; k < 2; k++) {
          dT += b_x[k + 1];
        }

        /*  Force vector magnitude signal */
        power(posttrim, d_x);
        b_sum(d_x, t);
        b_sqrt(t);

        /*  Maximum force */
        xlast = quantile(t);

        /*  Detrended fluctuation analysis scaling exponent */
        /*  Performs fast detrended fluctuation analysis on a nonstationary input signal to */
        /*  obtain an estimate for the scaling exponent. */
        /*  */
        /*  Useage: */
        /*  [alpha, intervals, flucts] = fastdfa(x) */
        /*  [alpha, intervals, flucts] = fastdfa(x, intervals) */
        /*  Inputs */
        /*     x          - input signal: must be a row vector */
        /*  Optional inputs */
        /*     intervals  - List of sample interval widths at each scale */
        /*                  (If not specified, then a binary subdivision is constructed) */
        /*  */
        /*  Outputs: */
        /*     alpha      - Estimated scaling exponent */
        /*     intervals  - List of sample interval widths at each scale */
        /*     flucts     - List of fluctuation amplitudes at each scale */
        /*  */
        /*  (CC BY-SA 3.0) Max Little, 2006-2014. */
        /*  If you use this code for academic publication, please cite: */
        /*  M. Little, P. McSharry, I. Moroz, S. Roberts (2006), */
        /*  Nonlinear, biophysically-informed speech pathology detection */
        /*  in Proceedings of ICASSP 2006, IEEE Publishers: Toulouse, France. */
        nscales = floor(scalar_log2(t->size[0]));
        emxFree_real_T(&d_x);
        emxFree_real_T(&posttrim);
        if (rt_powd_snf(2.0, nscales - 1.0) > (double)t->size[0] / 2.5) {
          nscales--;
        }

        vstride = dt->size[0];
        dt->size[0] = (int)nscales;
        emxEnsureCapacity((emxArray__common *)dt, vstride, (int)sizeof(double));
        idx = (int)nscales;
        for (vstride = 0; vstride < idx; vstride++) {
          dt->data[vstride] = 0.0;
        }

        b_emxInit_real_T(&ypts, 1);
        vstride = ypts->size[0];
        ypts->size[0] = (int)nscales;
        emxEnsureCapacity((emxArray__common *)ypts, vstride, (int)sizeof(double));
        idx = (int)nscales;
        for (vstride = 0; vstride < idx; vstride++) {
          ypts->data[vstride] = 0.0;
        }

        b_emxInit_real_T(&e_x, 1);
        vstride = e_x->size[0];
        e_x->size[0] = t->size[0];
        emxEnsureCapacity((emxArray__common *)e_x, vstride, (int)sizeof(double));
        idx = t->size[0];
        for (vstride = 0; vstride < idx; vstride++) {
          e_x->data[vstride] = t->data[vstride];
        }

        emxInit_real_T(&b_dt, 2);
        fastdfa_core_nomex(&dt->data[0], &ypts->data[0], &e_x->data[0], (double)
                           t->size[0], 1.0);

        /*  Sort the intervals, and produce a log-log straight line fit */
        idx = dt->size[0];
        ii = ypts->size[0];
        vstride = b_dt->size[0] * b_dt->size[1];
        b_dt->size[0] = idx;
        b_dt->size[1] = 2;
        emxEnsureCapacity((emxArray__common *)b_dt, vstride, (int)sizeof(double));
        emxFree_real_T(&e_x);
        for (vstride = 0; vstride < idx; vstride++) {
          b_dt->data[vstride] = dt->data[vstride];
        }

        for (vstride = 0; vstride < ii; vstride++) {
          b_dt->data[vstride + b_dt->size[0]] = ypts->data[vstride];
        }

        emxInit_real_T(&c_dt, 2);
        sortrows(b_dt, c_dt);
        b_log10(dt);
        b_log10(ypts);
        polyfit(dt, ypts, coeffs);

        /*  Output posture test feature vector */
        ft[0] = xlast / 10.0;
        ft[1] = dT / 3.0 / 10000.0;
        ft[2] = coeffs[0];
        emxFree_real_T(&c_dt);
        emxFree_real_T(&b_dt);
        emxFree_real_T(&ypts);
        emxFree_real_T(&dt);
      }
    }

    emxFree_real_T(&t);
  }
}

/* End of code generation (features_bpa.c) */
