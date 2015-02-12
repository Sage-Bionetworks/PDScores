/*
 * fastdfa.c
 *
 * Code generation for function 'fastdfa'
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
#include "bridge_ufb_rtwutil.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Definitions */
void fastdfa(const emxArray_real_T *x, double *alpha, emxArray_real_T *intervals,
             emxArray_real_T *flucts)
{
  double nscales;
  emxArray_real_T *xpts;
  int i5;
  int loop_ub;
  emxArray_real_T *ypts;
  emxArray_real_T *b_x;
  emxArray_real_T *datapts;
  int ypts_idx_0;
  double coeffs[2];

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
  nscales = floor(scalar_log2(x->size[1]));
  if (rt_powd_snf(2.0, nscales - 1.0) > (double)x->size[1] / 2.5) {
    nscales--;
  }

  b_emxInit_real_T(&xpts, 1);
  i5 = xpts->size[0];
  xpts->size[0] = (int)nscales;
  emxEnsureCapacity((emxArray__common *)xpts, i5, (int)sizeof(double));
  loop_ub = (int)nscales;
  for (i5 = 0; i5 < loop_ub; i5++) {
    xpts->data[i5] = 0.0;
  }

  b_emxInit_real_T(&ypts, 1);
  i5 = ypts->size[0];
  ypts->size[0] = (int)nscales;
  emxEnsureCapacity((emxArray__common *)ypts, i5, (int)sizeof(double));
  loop_ub = (int)nscales;
  for (i5 = 0; i5 < loop_ub; i5++) {
    ypts->data[i5] = 0.0;
  }

  emxInit_real_T(&b_x, 2);
  i5 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = 1;
  b_x->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)b_x, i5, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i5 = 0; i5 < loop_ub; i5++) {
    b_x->data[i5] = x->data[i5];
  }

  emxInit_real_T(&datapts, 2);
  fastdfa_core_nomex(&xpts->data[0], &ypts->data[0], &b_x->data[0], 1.0, (double)
                     x->size[1]);

  /*  Sort the intervals, and produce a log-log straight line fit */
  loop_ub = xpts->size[0];
  ypts_idx_0 = ypts->size[0];
  i5 = datapts->size[0] * datapts->size[1];
  datapts->size[0] = loop_ub;
  datapts->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)datapts, i5, (int)sizeof(double));
  emxFree_real_T(&b_x);
  for (i5 = 0; i5 < loop_ub; i5++) {
    datapts->data[i5] = xpts->data[i5];
  }

  for (i5 = 0; i5 < ypts_idx_0; i5++) {
    datapts->data[i5 + datapts->size[0]] = ypts->data[i5];
  }

  b_sortrows(datapts);
  loop_ub = datapts->size[0];
  i5 = intervals->size[0];
  intervals->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)intervals, i5, (int)sizeof(double));
  for (i5 = 0; i5 < loop_ub; i5++) {
    intervals->data[i5] = datapts->data[i5];
  }

  loop_ub = datapts->size[0];
  i5 = flucts->size[0];
  flucts->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)flucts, i5, (int)sizeof(double));
  for (i5 = 0; i5 < loop_ub; i5++) {
    flucts->data[i5] = datapts->data[i5 + datapts->size[0]];
  }

  emxFree_real_T(&datapts);
  b_log10(xpts);
  b_log10(ypts);
  polyfit(xpts, ypts, coeffs);
  *alpha = coeffs[0];
  emxFree_real_T(&ypts);
  emxFree_real_T(&xpts);
}

/* End of code generation (fastdfa.c) */
