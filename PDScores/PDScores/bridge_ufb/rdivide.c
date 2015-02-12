/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
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
#include "rdivide.h"
#include "bridge_ufb_emxutil.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Definitions */
double b_rdivide(double x, double y)
{
  return x / y;
}

void c_rdivide(const emxArray_real_T *x, double y, emxArray_real_T *z)
{
  int i4;
  int loop_ub;
  i4 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)z, i4, (int)sizeof(double));
  loop_ub = x->size[0] * x->size[1];
  for (i4 = 0; i4 < loop_ub; i4++) {
    z->data[i4] = x->data[i4] / y;
  }
}

void d_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
               emxArray_real_T *z)
{
  int i12;
  int loop_ub;
  i12 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i12, (int)sizeof(double));
  loop_ub = x->size[0];
  for (i12 = 0; i12 < loop_ub; i12++) {
    z->data[i12] = x->data[i12] / y->data[i12];
  }
}

void e_rdivide(double x, const double y[2], double z[2])
{
  int i16;
  for (i16 = 0; i16 < 2; i16++) {
    z[i16] = x / y[i16];
  }
}

void f_rdivide(double x, const emxArray_real_T *y, emxArray_real_T *z)
{
  int i17;
  int loop_ub;
  i17 = z->size[0] * z->size[1];
  z->size[0] = 1;
  z->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)z, i17, (int)sizeof(double));
  loop_ub = y->size[0] * y->size[1];
  for (i17 = 0; i17 < loop_ub; i17++) {
    z->data[i17] = x / y->data[i17];
  }
}

void g_rdivide(const emxArray_real_T *y, emxArray_real_T *z)
{
  int i20;
  int loop_ub;
  i20 = z->size[0];
  z->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)z, i20, (int)sizeof(double));
  loop_ub = y->size[0];
  for (i20 = 0; i20 < loop_ub; i20++) {
    z->data[i20] = 1.0 / y->data[i20];
  }
}

void rdivide(const emxArray_real_T *x, double y, emxArray_real_T *z)
{
  int i3;
  int loop_ub;
  i3 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i3, (int)sizeof(double));
  loop_ub = x->size[0];
  for (i3 = 0; i3 < loop_ub; i3++) {
    z->data[i3] = x->data[i3] / y;
  }
}

/* End of code generation (rdivide.c) */
