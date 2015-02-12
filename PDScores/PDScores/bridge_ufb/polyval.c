/*
 * polyval.c
 *
 * Code generation for function 'polyval'
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
#include "polyval.h"
#include "bridge_ufb_emxutil.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Definitions */
void polyval(const double p[3], const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int uv5[2];
  int i25;
  int loop_ub;
  int k;
  for (i25 = 0; i25 < 2; i25++) {
    uv5[i25] = (unsigned int)x->size[i25];
  }

  i25 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int)uv5[1];
  emxEnsureCapacity((emxArray__common *)y, i25, (int)sizeof(double));
  if (!((int)uv5[1] == 0)) {
    i25 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)y, i25, (int)sizeof(double));
    i25 = y->size[0] * y->size[1];
    y->size[1] = (int)uv5[1];
    emxEnsureCapacity((emxArray__common *)y, i25, (int)sizeof(double));
    loop_ub = (int)uv5[1];
    for (i25 = 0; i25 < loop_ub; i25++) {
      y->data[i25] = p[0];
    }

    for (k = 0; k < 2; k++) {
      i25 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity((emxArray__common *)y, i25, (int)sizeof(double));
      loop_ub = x->size[0] * x->size[1];
      for (i25 = 0; i25 < loop_ub; i25++) {
        y->data[i25] = x->data[i25] * y->data[i25] + p[k + 1];
      }
    }
  }
}

/* End of code generation (polyval.c) */
