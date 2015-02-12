/*
 * abs.c
 *
 * Code generation for function 'abs'
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
#include "abs.h"
#include "bridge_ufb_emxutil.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Definitions */
void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int unnamed_idx_0;
  int k;
  unnamed_idx_0 = (unsigned int)x->size[0];
  k = y->size[0];
  y->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k < x->size[0]; k++) {
    y->data[k] = fabs(x->data[k]);
  }
}

void c_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int uv4[2];
  int i18;
  int k;
  for (i18 = 0; i18 < 2; i18++) {
    uv4[i18] = (unsigned int)x->size[i18];
  }

  i18 = y->size[0] * y->size[1];
  y->size[0] = (int)uv4[0];
  y->size[1] = (int)uv4[1];
  emxEnsureCapacity((emxArray__common *)y, i18, (int)sizeof(double));
  i18 = x->size[0] * x->size[1];
  for (k = 0; k < i18; k++) {
    y->data[k] = fabs(x->data[k]);
  }
}

void d_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int uv6[2];
  int k;
  for (k = 0; k < 2; k++) {
    uv6[k] = (unsigned int)x->size[k];
  }

  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv6[1];
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k < x->size[1]; k++) {
    y->data[k] = fabs(x->data[k]);
  }
}

/* End of code generation (abs.c) */
