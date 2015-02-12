/*
 * sortrows.c
 *
 * Code generation for function 'sortrows'
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
#include "sortrows.h"
#include "bridge_ufb_emxutil.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Declarations */
static boolean_T eml_sort_le(const emxArray_real_T *v, int irow1, int irow2);

/* Function Definitions */
static boolean_T eml_sort_le(const emxArray_real_T *v, int irow1, int irow2)
{
  boolean_T p;
  boolean_T b0;
  p = true;
  if ((v->data[irow1 - 1] == v->data[irow2 - 1]) || (rtIsNaN(v->data[irow1 - 1])
       && rtIsNaN(v->data[irow2 - 1]))) {
    b0 = true;
  } else {
    b0 = false;
  }

  if (!b0) {
    if ((v->data[irow1 - 1] <= v->data[irow2 - 1]) || rtIsNaN(v->data[irow2 - 1]))
    {
      p = true;
    } else {
      p = false;
    }
  }

  return p;
}

void b_sortrows(emxArray_real_T *y)
{
  emxArray_int32_T *idx;
  int n;
  int i2;
  int k;
  emxArray_int32_T *idx0;
  int i;
  int j;
  int pEnd;
  int p;
  int q;
  int qEnd;
  int kEnd;
  emxArray_real_T *ycol;
  emxInit_int32_T(&idx, 1);
  n = y->size[0];
  i2 = idx->size[0];
  idx->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)idx, i2, (int)sizeof(int));
  if (y->size[0] == 0) {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }
  } else {
    for (k = 1; k <= n; k++) {
      idx->data[k - 1] = k;
    }

    for (k = 1; k <= n - 1; k += 2) {
      if (eml_sort_le(y, k, k + 1)) {
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    emxInit_int32_T(&idx0, 1);
    i2 = idx0->size[0];
    idx0->size[0] = y->size[0];
    emxEnsureCapacity((emxArray__common *)idx0, i2, (int)sizeof(int));
    k = y->size[0];
    for (i2 = 0; i2 < k; i2++) {
      idx0->data[i2] = 1;
    }

    i = 2;
    while (i < n) {
      i2 = i << 1;
      j = 1;
      for (pEnd = 1 + i; pEnd < n + 1; pEnd = qEnd + i) {
        p = j;
        q = pEnd;
        qEnd = j + i2;
        if (qEnd > n + 1) {
          qEnd = n + 1;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if (eml_sort_le(y, idx->data[p - 1], idx->data[q - 1])) {
            idx0->data[k] = idx->data[p - 1];
            p++;
            if (p == pEnd) {
              while (q < qEnd) {
                k++;
                idx0->data[k] = idx->data[q - 1];
                q++;
              }
            }
          } else {
            idx0->data[k] = idx->data[q - 1];
            q++;
            if (q == qEnd) {
              while (p < pEnd) {
                k++;
                idx0->data[k] = idx->data[p - 1];
                p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k + 1 <= kEnd; k++) {
          idx->data[(j + k) - 1] = idx0->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&idx0);
  }

  b_emxInit_real_T(&ycol, 1);
  pEnd = y->size[0];
  k = y->size[0];
  i2 = ycol->size[0];
  ycol->size[0] = k;
  emxEnsureCapacity((emxArray__common *)ycol, i2, (int)sizeof(double));
  for (j = 0; j < 2; j++) {
    for (i = 0; i + 1 <= pEnd; i++) {
      ycol->data[i] = y->data[(idx->data[i] + y->size[0] * j) - 1];
    }

    for (i = 0; i + 1 <= pEnd; i++) {
      y->data[i + y->size[0] * j] = ycol->data[i];
    }
  }

  emxFree_real_T(&ycol);
  emxFree_int32_T(&idx);
}

void sortrows(const emxArray_real_T *y, emxArray_real_T *b_y)
{
  int i6;
  int loop_ub;
  i6 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = y->size[0];
  b_y->size[1] = 2;
  emxEnsureCapacity((emxArray__common *)b_y, i6, (int)sizeof(double));
  loop_ub = y->size[0] * y->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    b_y->data[i6] = y->data[i6];
  }

  b_sortrows(b_y);
}

/* End of code generation (sortrows.c) */
