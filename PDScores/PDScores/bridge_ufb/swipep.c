/*
 * swipep.c
 *
 * Code generation for function 'swipep'
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
#include "repmat.h"
#include "abs.h"
#include "interp1.h"
#include "rdivide.h"
#include "sqrt.h"
#include "flipud.h"
#include "iqr.h"
#include "log2.h"
#include "polyval.h"
#include "power.h"
#include "polyfit.h"
#include "round.h"
#include "primes.h"
#include "bridge_ufb_rtwutil.h"
#include "bridge_ufb_data.h"

/* Custom Source Code */
#include "fastdfa_core_nomex.h"
#include "buffer.h"
#include "signalprocessing.h"

/* Function Declarations */
static void pitchStrengthAllCandidates(const emxArray_real_T *f, const
  emxArray_real_T *L, const emxArray_real_T *pc, emxArray_real_T *S);
static void pitchStrengthOneCandidate(const emxArray_real_T *f, const
  emxArray_real_T *NL, double pc, emxArray_real_T *S);

/* Function Definitions */
static void pitchStrengthAllCandidates(const emxArray_real_T *f, const
  emxArray_real_T *L, const emxArray_real_T *pc, emxArray_real_T *S)
{
  int u1;
  int dim;
  int i;
  int absb;
  int n;
  emxArray_real_T *k;
  int j;
  emxArray_boolean_T *x;
  int vstride;
  double xlast;
  int npages;
  int b_k;
  int ii_data[1];
  int cdiff;
  boolean_T exitg1;
  int tmp_data[1];
  emxArray_real_T *c_k;
  emxArray_real_T *N;
  emxArray_real_T *b_n;
  emxArray_int32_T *r20;
  emxArray_int32_T *r21;
  emxArray_boolean_T *b_x;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  emxArray_real_T *b_f;
  emxArray_real_T *b_L;

  /*  Create pitch strength matrix */
  if ((0 == pc->size[0]) || (0 == pc->size[1])) {
    u1 = 0;
  } else {
    dim = pc->size[0];
    u1 = pc->size[1];
    if (dim >= u1) {
      u1 = dim;
    }
  }

  i = S->size[0] * S->size[1];
  S->size[0] = u1;
  emxEnsureCapacity((emxArray__common *)S, i, (int)sizeof(double));
  dim = L->size[1];
  i = S->size[0] * S->size[1];
  S->size[1] = dim;
  emxEnsureCapacity((emxArray__common *)S, i, (int)sizeof(double));
  absb = u1 * L->size[1];
  for (i = 0; i < absb; i++) {
    S->data[i] = 0.0;
  }

  /*  Define integration regions */
  if ((0 == pc->size[0]) || (0 == pc->size[1])) {
    n = 1;
  } else if (pc->size[0] > pc->size[1]) {
    n = pc->size[0] + 1;
  } else {
    n = pc->size[1] + 1;
  }

  emxInit_real_T(&k, 2);
  i = k->size[0] * k->size[1];
  k->size[0] = 1;
  k->size[1] = n;
  emxEnsureCapacity((emxArray__common *)k, i, (int)sizeof(double));
  for (i = 0; i < n; i++) {
    k->data[i] = 1.0;
  }

  j = 0;
  emxInit_boolean_T(&x, 1);
  while (j <= n - 2) {
    if (k->data[j] > f->size[0]) {
      i = 0;
      vstride = 0;
    } else {
      i = (int)k->data[j] - 1;
      vstride = f->size[0];
    }

    xlast = pc->data[j] / 4.0;
    npages = x->size[0];
    x->size[0] = vstride - i;
    emxEnsureCapacity((emxArray__common *)x, npages, (int)sizeof(boolean_T));
    absb = vstride - i;
    for (vstride = 0; vstride < absb; vstride++) {
      x->data[vstride] = (f->data[i + vstride] > xlast);
    }

    if (1 <= x->size[0]) {
      b_k = 1;
    } else {
      b_k = 0;
    }

    dim = 0;
    cdiff = 1;
    exitg1 = false;
    while ((!exitg1) && (cdiff <= x->size[0])) {
      if (x->data[cdiff - 1]) {
        dim = 1;
        ii_data[0] = cdiff;
        exitg1 = true;
      } else {
        cdiff++;
      }
    }

    if (b_k == 1) {
      if (dim == 0) {
        b_k = 0;
      }
    } else if (1 > dim) {
      b_k = 0;
    } else {
      b_k = 1;
    }

    i = 0;
    while (i <= b_k - 1) {
      tmp_data[0] = ii_data[0];
      i = 1;
    }

    k->data[j + 1] = (k->data[j] - 1.0) + (double)tmp_data[0];
    j++;
  }

  emxFree_boolean_T(&x);
  if (2 > k->size[1]) {
    i = 0;
    vstride = 0;
  } else {
    i = 1;
    vstride = k->size[1];
  }

  emxInit_real_T(&c_k, 2);
  npages = c_k->size[0] * c_k->size[1];
  c_k->size[0] = 1;
  c_k->size[1] = vstride - i;
  emxEnsureCapacity((emxArray__common *)c_k, npages, (int)sizeof(double));
  absb = vstride - i;
  for (vstride = 0; vstride < absb; vstride++) {
    c_k->data[c_k->size[0] * vstride] = k->data[i + vstride];
  }

  i = k->size[0] * k->size[1];
  k->size[0] = 1;
  k->size[1] = c_k->size[1];
  emxEnsureCapacity((emxArray__common *)k, i, (int)sizeof(double));
  absb = c_k->size[1];
  for (i = 0; i < absb; i++) {
    k->data[k->size[0] * i] = c_k->data[c_k->size[0] * i];
  }

  emxFree_real_T(&c_k);
  emxInit_real_T(&N, 2);

  /*  Create loudness normalization matrix */
  i = N->size[0] * N->size[1];
  N->size[0] = L->size[0];
  N->size[1] = L->size[1];
  emxEnsureCapacity((emxArray__common *)N, i, (int)sizeof(double));
  absb = L->size[0] * L->size[1];
  for (i = 0; i < absb; i++) {
    N->data[i] = L->data[i] * L->data[i];
  }

  flipud(N);
  dim = 1;
  if (N->size[0] != 1) {
    dim = 0;
  }

  absb = N->size[dim];
  if ((!((N->size[0] == 0) || (N->size[1] == 0))) && (N->size[dim] > 1)) {
    vstride = 1;
    b_k = 1;
    while (b_k <= dim) {
      vstride *= N->size[0];
      b_k = 2;
    }

    npages = 1;
    b_k = dim + 2;
    while (b_k < 3) {
      npages *= N->size[1];
      b_k = 3;
    }

    dim = 0;
    for (i = 1; i <= npages; i++) {
      cdiff = dim;
      for (j = 1; j <= vstride; j++) {
        cdiff++;
        dim = cdiff;
        xlast = N->data[cdiff - 1];
        for (b_k = 0; b_k <= absb - 2; b_k++) {
          dim += vstride;
          xlast += N->data[dim - 1];
          N->data[dim - 1] = xlast;
        }
      }
    }
  }

  flipud(N);
  c_sqrt(N);
  if ((0 == pc->size[0]) || (0 == pc->size[1])) {
    u1 = 0;
  } else {
    dim = pc->size[0];
    u1 = pc->size[1];
    if (dim >= u1) {
      u1 = dim;
    }
  }

  j = 0;
  emxInit_real_T(&b_n, 2);
  b_emxInit_int32_T(&r20, 2);
  emxInit_int32_T(&r21, 1);
  b_emxInit_boolean_T(&b_x, 2);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&b_y, 2);
  b_emxInit_real_T(&b_f, 1);
  emxInit_real_T(&b_L, 2);
  while (j <= u1 - 1) {
    /*  Normalize loudness */
    absb = N->size[1];
    b_k = (int)k->data[j];
    i = b_n->size[0] * b_n->size[1];
    b_n->size[0] = 1;
    b_n->size[1] = absb;
    emxEnsureCapacity((emxArray__common *)b_n, i, (int)sizeof(double));
    for (i = 0; i < absb; i++) {
      b_n->data[b_n->size[0] * i] = N->data[(b_k + N->size[0] * i) - 1];
    }

    i = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1;
    b_x->size[1] = b_n->size[1];
    emxEnsureCapacity((emxArray__common *)b_x, i, (int)sizeof(boolean_T));
    absb = b_n->size[0] * b_n->size[1];
    for (i = 0; i < absb; i++) {
      b_x->data[i] = (b_n->data[i] == 0.0);
    }

    n = b_x->size[1];
    b_k = 0;
    for (i = 1; i <= n; i++) {
      if (b_x->data[i - 1]) {
        b_k++;
      }
    }

    i = r20->size[0] * r20->size[1];
    r20->size[0] = 1;
    r20->size[1] = b_k;
    emxEnsureCapacity((emxArray__common *)r20, i, (int)sizeof(int));
    dim = 0;
    for (i = 1; i <= n; i++) {
      if (b_x->data[i - 1]) {
        r20->data[dim] = i;
        dim++;
      }
    }

    absb = r20->size[0] * r20->size[1];
    for (i = 0; i < absb; i++) {
      b_n->data[r20->data[i] - 1] = rtInf;
    }

    /*  to make zero-loudness equal zero after normalization */
    /*      NL = L(k(j):end,:) ./ repmat( n, size(L,1)-k(j)+1, 1); */
    if (b_n->size[1] < 1) {
      n = -1;
      vstride = 0;
    } else {
      dim = (int)floor(((double)b_n->size[1] - 1.0) + 0.5);
      vstride = dim + 1;
      cdiff = (dim - b_n->size[1]) + 1;
      absb = b_n->size[1];
      if (fabs(cdiff) < 4.4408920985006262E-16 * (double)absb) {
        dim++;
        vstride = b_n->size[1];
      } else if (cdiff > 0) {
        vstride = dim;
      } else {
        dim++;
      }

      n = dim - 1;
    }

    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
    if (n + 1 > 0) {
      y->data[0] = 1.0;
      if (n + 1 > 1) {
        y->data[n] = vstride;
        dim = (n + (n < 0)) >> 1;
        for (b_k = 1; b_k < dim; b_k++) {
          y->data[b_k] = 1.0 + (double)b_k;
          y->data[n - b_k] = vstride - b_k;
        }

        if (dim << 1 == n) {
          y->data[dim] = (1.0 + (double)vstride) / 2.0;
        } else {
          y->data[dim] = 1.0 + (double)dim;
          y->data[dim + 1] = vstride - dim;
        }
      }
    }

    /*      NL = L(k(j):end,:) ./ n(rowIdx(:,ones(size(L,1)-k(j)+1,1)), colIdx(:,1)); */
    if (k->data[j] > L->size[0]) {
      i = 0;
      vstride = 0;
    } else {
      i = (int)k->data[j] - 1;
      vstride = L->size[0];
    }

    npages = r21->size[0];
    r21->size[0] = (int)(((double)L->size[0] - k->data[j]) + 1.0);
    emxEnsureCapacity((emxArray__common *)r21, npages, (int)sizeof(int));
    absb = (int)(((double)L->size[0] - k->data[j]) + 1.0);
    for (npages = 0; npages < absb; npages++) {
      r21->data[npages] = 0;
    }

    dim = r21->size[0];
    npages = b_y->size[0] * b_y->size[1];
    b_y->size[0] = dim;
    b_y->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, npages, (int)sizeof(double));
    absb = y->size[1];
    for (npages = 0; npages < absb; npages++) {
      for (cdiff = 0; cdiff < dim; cdiff++) {
        b_y->data[cdiff + b_y->size[0] * npages] = b_n->data[b_n->size[0] *
          ((int)y->data[y->size[0] * npages] - 1)];
      }
    }

    /*  Compute pitch strength */
    if (k->data[j] > f->size[0]) {
      npages = 0;
      cdiff = 0;
    } else {
      npages = (int)k->data[j] - 1;
      cdiff = f->size[0];
    }

    dim = b_f->size[0];
    b_f->size[0] = cdiff - npages;
    emxEnsureCapacity((emxArray__common *)b_f, dim, (int)sizeof(double));
    absb = cdiff - npages;
    for (cdiff = 0; cdiff < absb; cdiff++) {
      b_f->data[cdiff] = f->data[npages + cdiff];
    }

    absb = L->size[1];
    npages = b_L->size[0] * b_L->size[1];
    b_L->size[0] = vstride - i;
    b_L->size[1] = absb;
    emxEnsureCapacity((emxArray__common *)b_L, npages, (int)sizeof(double));
    for (npages = 0; npages < absb; npages++) {
      dim = vstride - i;
      for (cdiff = 0; cdiff < dim; cdiff++) {
        b_L->data[cdiff + b_L->size[0] * npages] = L->data[(i + cdiff) + L->
          size[0] * npages] / b_y->data[cdiff + b_y->size[0] * npages];
      }
    }

    pitchStrengthOneCandidate(b_f, b_L, pc->data[j], b_n);
    absb = b_n->size[1];
    for (i = 0; i < absb; i++) {
      S->data[j + S->size[0] * i] = b_n->data[b_n->size[0] * i];
    }

    j++;
  }

  emxFree_real_T(&b_L);
  emxFree_real_T(&b_f);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  emxFree_boolean_T(&b_x);
  emxFree_int32_T(&r21);
  emxFree_int32_T(&r20);
  emxFree_real_T(&b_n);
  emxFree_real_T(&N);
  emxFree_real_T(&k);
}

static void pitchStrengthOneCandidate(const emxArray_real_T *f, const
  emxArray_real_T *NL, double pc, emxArray_real_T *S)
{
  double n;
  int i19;
  emxArray_real_T *k;
  unsigned int idx;
  int ar;
  emxArray_real_T *q;
  int numPrimes;
  double i;
  emxArray_real_T *a;
  emxArray_boolean_T *p;
  emxArray_real_T *x;
  emxArray_real_T *b;
  emxArray_real_T *r22;
  emxArray_int32_T *r23;
  emxArray_int32_T *r24;
  emxArray_real_T *b_q;
  int b_k;
  emxArray_boolean_T *c_k;
  double absxk;
  double t;
  emxArray_real_T *d_k;
  emxArray_real_T *b_a;
  int br;
  int b_n;
  int ic;
  int ib;
  int ia;
  n = trunc(f->data[f->size[0] - 1] / pc - 0.75);

  /*  Number of harmonics */
  if (n == 0.0) {
    i19 = S->size[0] * S->size[1];
    S->size[0] = 1;
    S->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)S, i19, (int)sizeof(double));
    S->data[0] = rtNaN;
  } else {
    b_emxInit_real_T(&k, 1);
    idx = (unsigned int)f->size[0];
    i19 = k->size[0];
    k->size[0] = (int)idx;
    emxEnsureCapacity((emxArray__common *)k, i19, (int)sizeof(double));
    ar = (int)idx;
    for (i19 = 0; i19 < ar; i19++) {
      k->data[i19] = 0.0;
    }

    b_emxInit_real_T(&q, 1);

    /*  Kernel */
    /*  Normalize frequency w.r.t. candidate */
    rdivide(f, pc, q);

    /*  Create kernel */
    /*  for i = [ 1 primes(n) ] */
    /*  Coder error: "FOR loop index expressions of unknown size are only */
    /*  supported if they are of the form A:B or A:B:C." */
    /* for i = primetable(primetable <= n) */
    idx = 1U;
    numPrimes = primetable->size[1];
    i = primetable->data[0];
    b_emxInit_real_T(&a, 1);
    emxInit_boolean_T(&p, 1);
    b_emxInit_real_T(&x, 1);
    b_emxInit_real_T(&b, 1);
    b_emxInit_real_T(&r22, 1);
    emxInit_int32_T(&r23, 1);
    emxInit_int32_T(&r24, 1);
    b_emxInit_real_T(&b_q, 1);
    while (((int)idx <= numPrimes) && (i <= n)) {
      i19 = b_q->size[0];
      b_q->size[0] = q->size[0];
      emxEnsureCapacity((emxArray__common *)b_q, i19, (int)sizeof(double));
      ar = q->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        b_q->data[i19] = q->data[i19] - i;
      }

      b_abs(b_q, a);

      /*  Peak's weigth */
      i19 = p->size[0];
      p->size[0] = a->size[0];
      emxEnsureCapacity((emxArray__common *)p, i19, (int)sizeof(boolean_T));
      ar = a->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        p->data[i19] = (a->data[i19] < 0.25);
      }

      eml_li_find(p, r23);
      eml_li_find(p, r24);
      i19 = b->size[0];
      b->size[0] = r24->size[0];
      emxEnsureCapacity((emxArray__common *)b, i19, (int)sizeof(double));
      ar = r24->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        b->data[i19] = q->data[r24->data[i19] - 1];
      }

      i19 = b->size[0];
      emxEnsureCapacity((emxArray__common *)b, i19, (int)sizeof(double));
      ar = b->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        b->data[i19] *= 6.2831853071795862;
      }

      i19 = r22->size[0];
      r22->size[0] = b->size[0];
      emxEnsureCapacity((emxArray__common *)r22, i19, (int)sizeof(double));
      ar = b->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        r22->data[i19] = b->data[i19];
      }

      for (b_k = 0; b_k < b->size[0]; b_k++) {
        r22->data[b_k] = cos(r22->data[b_k]);
      }

      ar = r22->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        k->data[r23->data[i19] - 1] = r22->data[i19];
      }

      /*  Valleys' weights */
      i19 = p->size[0];
      p->size[0] = a->size[0];
      emxEnsureCapacity((emxArray__common *)p, i19, (int)sizeof(boolean_T));
      ar = a->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        p->data[i19] = ((0.25 < a->data[i19]) && (a->data[i19] < 0.75));
      }

      eml_li_find(p, r23);
      i19 = r22->size[0];
      r22->size[0] = r23->size[0];
      emxEnsureCapacity((emxArray__common *)r22, i19, (int)sizeof(double));
      ar = r23->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        r22->data[i19] = k->data[r23->data[i19] - 1];
      }

      eml_li_find(p, r23);
      i19 = b->size[0];
      b->size[0] = r23->size[0];
      emxEnsureCapacity((emxArray__common *)b, i19, (int)sizeof(double));
      ar = r23->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        b->data[i19] = q->data[r23->data[i19] - 1];
      }

      i19 = b->size[0];
      emxEnsureCapacity((emxArray__common *)b, i19, (int)sizeof(double));
      ar = b->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        b->data[i19] *= 6.2831853071795862;
      }

      i19 = x->size[0];
      x->size[0] = b->size[0];
      emxEnsureCapacity((emxArray__common *)x, i19, (int)sizeof(double));
      ar = b->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        x->data[i19] = b->data[i19];
      }

      for (b_k = 0; b_k < b->size[0]; b_k++) {
        x->data[b_k] = cos(x->data[b_k]);
      }

      eml_li_find(p, r23);
      rdivide(x, 2.0, a);
      ar = r22->size[0];
      for (i19 = 0; i19 < ar; i19++) {
        k->data[r23->data[i19] - 1] = r22->data[i19] + a->data[i19];
      }

      idx++;
      i = primetable->data[(int)idx - 1];
    }

    emxFree_real_T(&b_q);
    emxFree_int32_T(&r24);
    emxFree_real_T(&b);
    emxFree_boolean_T(&p);
    emxFree_real_T(&a);
    emxFree_real_T(&q);

    /*  Apply envelope */
    g_rdivide(f, r22);
    b_sqrt(r22);
    i19 = k->size[0];
    emxEnsureCapacity((emxArray__common *)k, i19, (int)sizeof(double));
    ar = k->size[0];
    for (i19 = 0; i19 < ar; i19++) {
      k->data[i19] *= r22->data[i19];
    }

    emxFree_real_T(&r22);
    emxInit_boolean_T(&c_k, 1);

    /*  K+-normalize kernel */
    i19 = c_k->size[0];
    c_k->size[0] = k->size[0];
    emxEnsureCapacity((emxArray__common *)c_k, i19, (int)sizeof(boolean_T));
    ar = k->size[0];
    for (i19 = 0; i19 < ar; i19++) {
      c_k->data[i19] = (k->data[i19] > 0.0);
    }

    eml_li_find(c_k, r23);
    i19 = x->size[0];
    x->size[0] = r23->size[0];
    emxEnsureCapacity((emxArray__common *)x, i19, (int)sizeof(double));
    ar = r23->size[0];
    emxFree_boolean_T(&c_k);
    for (i19 = 0; i19 < ar; i19++) {
      x->data[i19] = k->data[r23->data[i19] - 1];
    }

    emxFree_int32_T(&r23);
    n = 0.0;
    if (x->size[0] < 1) {
    } else if (x->size[0] == 1) {
      n = fabs(x->data[0]);
    } else {
      i = 2.2250738585072014E-308;
      for (b_k = 1; b_k <= x->size[0]; b_k++) {
        absxk = fabs(x->data[b_k - 1]);
        if (absxk > i) {
          t = i / absxk;
          n = 1.0 + n * t * t;
          i = absxk;
        } else {
          t = absxk / i;
          n += t * t;
        }
      }

      n = i * sqrt(n);
    }

    emxFree_real_T(&x);
    b_emxInit_real_T(&d_k, 1);
    i19 = d_k->size[0];
    d_k->size[0] = k->size[0];
    emxEnsureCapacity((emxArray__common *)d_k, i19, (int)sizeof(double));
    ar = k->size[0];
    for (i19 = 0; i19 < ar; i19++) {
      d_k->data[i19] = k->data[i19];
    }

    emxInit_real_T(&b_a, 2);
    rdivide(d_k, n, k);

    /*  Compute pitch strength */
    i19 = b_a->size[0] * b_a->size[1];
    b_a->size[0] = 1;
    b_a->size[1] = k->size[0];
    emxEnsureCapacity((emxArray__common *)b_a, i19, (int)sizeof(double));
    ar = k->size[0];
    emxFree_real_T(&d_k);
    for (i19 = 0; i19 < ar; i19++) {
      b_a->data[b_a->size[0] * i19] = k->data[i19];
    }

    emxFree_real_T(&k);
    if ((b_a->size[1] == 1) || (NL->size[0] == 1)) {
      i19 = S->size[0] * S->size[1];
      S->size[0] = 1;
      S->size[1] = NL->size[1];
      emxEnsureCapacity((emxArray__common *)S, i19, (int)sizeof(double));
      ar = NL->size[1];
      for (i19 = 0; i19 < ar; i19++) {
        S->data[S->size[0] * i19] = 0.0;
        numPrimes = b_a->size[1];
        for (br = 0; br < numPrimes; br++) {
          S->data[S->size[0] * i19] += b_a->data[b_a->size[0] * br] * NL->
            data[br + NL->size[0] * i19];
        }
      }
    } else {
      b_k = b_a->size[1];
      idx = (unsigned int)NL->size[1];
      b_n = NL->size[1] - 1;
      i19 = S->size[0] * S->size[1];
      S->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)S, i19, (int)sizeof(double));
      i19 = S->size[0] * S->size[1];
      S->size[1] = (int)idx;
      emxEnsureCapacity((emxArray__common *)S, i19, (int)sizeof(double));
      ar = (int)idx;
      for (i19 = 0; i19 < ar; i19++) {
        S->data[i19] = 0.0;
      }

      if (NL->size[1] == 0) {
      } else {
        for (numPrimes = 0; numPrimes <= b_n; numPrimes++) {
          for (ic = numPrimes; ic + 1 <= numPrimes + 1; ic++) {
            S->data[ic] = 0.0;
          }
        }

        br = 0;
        for (numPrimes = 0; numPrimes <= b_n; numPrimes++) {
          ar = 0;
          i19 = br + b_k;
          for (ib = br; ib + 1 <= i19; ib++) {
            if (NL->data[ib] != 0.0) {
              ia = ar;
              for (ic = numPrimes; ic + 1 <= numPrimes + 1; ic++) {
                ia++;
                S->data[ic] += NL->data[ib] * b_a->data[ia - 1];
              }
            }

            ar++;
          }

          br += b_k;
        }
      }
    }

    emxFree_real_T(&b_a);
  }
}

void b_swipep(const emxArray_real_T *x, double fs, emxArray_real_T *p,
              emxArray_real_T *t, emxArray_real_T *s)
{
  double y;
  int n;
  double o;
  double apnd;
  double ndbl;
  double cdiff;
  emxArray_real_T *b_y;
  int i15;
  int nm1d2;
  int b_cdiff;
  double kd;
  int b_ndbl;
  int tmp_size[2];
  double tmp_data[1229];
  double y_data[160];
  unsigned char y_size[2];
  static const unsigned char uv2[2] = { 1U, 160U };

  int log2pc_size[1];
  double log2pc_data[160];
  emxArray_real_T *xzp;
  emxArray_real_T b_log2pc_data;
  int pc_size_idx_0;
  double pc_data[160];
  emxArray_real_T *S;
  double b_x[2];
  double logWs[2];
  double absb;
  emxArray_real_T *ws;
  emxArray_real_T *pO;
  int ixstart;
  boolean_T exitg9;
  emxArray_real_T *c_y;
  emxArray_real_T *b_xzp;
  emxArray_real_T *fERBs;
  int i;
  emxArray_real_T *w;
  emxArray_real_T *X;
  emxArray_real_T *f;
  emxArray_real_T *ti;
  emxArray_real_T *L;
  emxArray_real_T *Si;
  emxArray_real_T *mu;
  emxArray_real_T *j;
  emxArray_real_T *k;
  emxArray_real_T *c_xzp;
  emxArray_real_T *r17;
  emxArray_real_T *r18;
  emxArray_boolean_T *c_x;
  emxArray_boolean_T *d_x;
  emxArray_int32_T *ii;
  emxArray_real_T *r19;
  emxArray_real_T *log2pc;
  emxArray_real_T *b_Si;
  emxArray_real_T *c_Si;
  emxArray_real_T *d_Si;
  emxArray_real_T *pc;
  emxArray_real_T *b_fERBs;
  emxArray_real_T *b_L;
  emxArray_real_T *b_S;
  boolean_T exitg8;
  boolean_T guard5 = false;
  int b_tmp_data[160];
  boolean_T exitg7;
  boolean_T guard4 = false;
  boolean_T exitg6;
  boolean_T guard3 = false;
  boolean_T exitg5;
  boolean_T guard2 = false;
  double c_log2pc_data[160];
  int b_log2pc_size[1];
  emxArray_real_T d_log2pc_data;
  boolean_T exitg4;
  boolean_T guard1 = false;
  int b_apnd;
  int ii_data[1];
  boolean_T exitg3;
  int found_data[1];
  int b_j;
  emxArray_real_T *b_pO;
  emxArray_real_T *d_xzp;
  boolean_T exitg2;
  double b_pc_data[161];
  int pc_size[1];
  emxArray_real_T c_pc_data;
  int tc_size[1];
  double tc_data[161];
  emxArray_real_T b_tc_data;
  double S_data[161];
  int S_size[1];
  emxArray_real_T b_S_data;
  double c[3];
  boolean_T exitg1;

  /*  SWIPEP Pitch estimation using SWIPE'. */
  /*     P = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,STHR) estimates the pitch  */
  /*     of the vector signal X every DT seconds. The sampling frequency of */
  /*     the signal is Fs (in Hertz). The spectrum is computed using a Hann */
  /*     window with an overlap WOVERLAP between 0 and 1. The spectrum is */
  /*     sampled uniformly in the ERB scale with a step size of DERBS ERBs. The */
  /*     pitch is searched within the range [PMIN PMAX] (in Hertz) with samples */
  /*     distributed every DLOG2P units on a base-2 logarithmic scale of Hertz.  */
  /*     The pitch is fine-tuned using parabolic interpolation with a resolution */
  /*     of 1 cent. Pitch estimates with a strength lower than STHR are treated */
  /*     as undefined. */
  /*      */
  /*     [P,T,S] = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,WOVERLAP,STHR)  */
  /*     returns the times T at which the pitch was estimated and the pitch  */
  /*     strength S of every pitch estimate. */
  /*  */
  /*     P = SWIPEP(X,Fs) estimates the pitch using the default settings PMIN = */
  /*     30 Hz, PMAX = 5000 Hz, DT = 0.001 s, DLOG2P = 1/48 (48 steps per  */
  /*     octave), DERBS = 0.1 ERBs, WOVERLAP = 0.5, and STHR = -Inf. */
  /*  */
  /*     P = SWIPEP(X,Fs,...,[],...) uses the default setting for the parameter */
  /*     replaced with the placeholder []. */
  /*  */
  /*     REMARKS: (1) For better results, make DLOG2P and DERBS as small as  */
  /*     possible and WOVERLAP as large as possible. However, take into account */
  /*     that the computational complexity of the algorithm is inversely  */
  /*     proportional to DLOG2P, DERBS and 1-WOVERLAP, and that the  default  */
  /*     values have been found empirically to produce good results. Consider  */
  /*     also that the computational complexity is directly proportional to the */
  /*     number of octaves in the pitch search range, and therefore , it is  */
  /*     recommendable to restrict the search range to the expected range of */
  /*     pitch, if any. (2) This code implements SWIPE', which uses only the */
  /*     first and prime harmonics of the signal. To convert it into SWIPE, */
  /*     which uses all the harmonics of the signal, replace the word */
  /*     PRIMES with a colon (it is located almost at the end of the code). */
  /*     However, this may not be recommendable since SWIPE' is reported to  */
  /*     produce on average better results than SWIPE (Camacho and Harris, */
  /*     2008). */
  /*  */
  /*     EXAMPLE: Estimate the pitch of the signal X every 10 ms within the */
  /*     range 75-500 Hz using the default resolution (i.e., 48 steps per */
  /*     octave), sampling the spectrum every 1/20th of ERB, using a window  */
  /*     overlap factor of 50%, and discarding samples with pitch strength  */
  /*     lower than 0.2. Plot the pitch trace. */
  /*        [x,Fs] = wavread(filename); */
  /*        [p,t,s] = swipep(x,Fs,[75 500],0.01,[],1/20,0.5,0.2); */
  /*        plot(1000*t,p) */
  /*        xlabel('Time (ms)') */
  /*        ylabel('Pitch (Hz)') */
  /*  */
  /*     REFERENCES: Camacho, A., Harris, J.G, (2008) "A sawtooth waveform  */
  /*     inspired pitch estimator for speech and music," J. Acoust. Soc. Am. */
  /*     124, 1638-1652. */
  /*  */
  /*     MAINTENANCE HISTORY: */
  /*     - Added line 153 to avoid division by zero in line 154 if loudness */
  /*       equals zero (06/23/2010). */
  y = (double)x->size[0] / fs;
  if (rtIsNaN(y)) {
    n = 0;
    o = rtNaN;
    apnd = y;
  } else if (y < 0.0) {
    n = -1;
    o = 0.0;
    apnd = y;
  } else if (rtIsInf(y)) {
    n = 0;
    o = rtNaN;
    apnd = y;
  } else {
    o = 0.0;
    ndbl = floor(y / 0.02 + 0.5);
    apnd = ndbl * 0.02;
    cdiff = apnd - y;
    if (fabs(cdiff) < 4.4408920985006262E-16 * y) {
      ndbl++;
      apnd = y;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * 0.02;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  emxInit_real_T(&b_y, 2);
  i15 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i15, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = o;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (b_cdiff = 1; b_cdiff < nm1d2; b_cdiff++) {
        kd = (double)b_cdiff * 0.02;
        b_y->data[b_cdiff] = o + kd;
        b_y->data[n - b_cdiff] = apnd - kd;
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (o + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * 0.02;
        b_y->data[nm1d2] = o + kd;
        b_y->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  i15 = t->size[0];
  t->size[0] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)t, i15, (int)sizeof(double));
  b_ndbl = b_y->size[1];
  for (i15 = 0; i15 < b_ndbl; i15++) {
    t->data[i15] = b_y->data[b_y->size[0] * i15];
  }

  /*  Times */
  if (primetable->size[1] == 0) {
    primes(tmp_data, tmp_size);
    i15 = primetable->size[0] * primetable->size[1];
    primetable->size[0] = 1;
    primetable->size[1] = 1 + tmp_size[1];
    emxEnsureCapacity((emxArray__common *)primetable, i15, (int)sizeof(double));
    primetable->data[0] = 1.0;
    b_ndbl = tmp_size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      primetable->data[primetable->size[0] * (i15 + 1)] = tmp_data[tmp_size[0] *
        i15];
    }
  }

  /*  Define pitch candidates */
  for (i15 = 0; i15 < 2; i15++) {
    y_size[i15] = uv2[i15];
  }

  y_data[0] = 5.6438561897747244;
  y_data[159] = 8.9563561897747235;
  for (b_cdiff = 0; b_cdiff < 78; b_cdiff++) {
    kd = ((double)b_cdiff + 1.0) * 0.020833333333333332;
    y_data[b_cdiff + 1] = 5.6438561897747244 + kd;
    y_data[158 - b_cdiff] = 8.9563561897747235 - kd;
  }

  y_data[79] = 7.2896895231080574;
  y_data[80] = 7.31052285644139;
  log2pc_size[0] = y_size[1];
  b_ndbl = y_size[1];
  for (i15 = 0; i15 < b_ndbl; i15++) {
    log2pc_data[i15] = y_data[y_size[0] * i15];
  }

  b_emxInit_real_T(&xzp, 1);
  b_log2pc_data.data = (double *)&log2pc_data;
  b_log2pc_data.size = (int *)&log2pc_size;
  b_log2pc_data.allocatedSize = 160;
  b_log2pc_data.numDimensions = 1;
  b_log2pc_data.canFreeData = false;
  e_power(&b_log2pc_data, xzp);
  pc_size_idx_0 = xzp->size[0];
  b_ndbl = xzp->size[0];
  for (i15 = 0; i15 < b_ndbl; i15++) {
    pc_data[i15] = xzp->data[i15];
  }

  emxInit_real_T(&S, 2);
  i15 = S->size[0] * S->size[1];
  S->size[0] = pc_size_idx_0;
  emxEnsureCapacity((emxArray__common *)S, i15, (int)sizeof(double));
  nm1d2 = t->size[0];
  i15 = S->size[0] * S->size[1];
  S->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)S, i15, (int)sizeof(double));
  b_ndbl = pc_size_idx_0 * t->size[0];
  for (i15 = 0; i15 < b_ndbl; i15++) {
    S->data[i15] = 0.0;
  }

  /*  Pitch strength matrix */
  /*  Determine P2-WSs */
  kd = 8.0 * fs;
  for (i15 = 0; i15 < 2; i15++) {
    b_x[i15] = kd / (50.0 + 450.0 * (double)i15);
  }

  b_log2(b_x, logWs);
  b_round(logWs);
  if (rtIsNaN(logWs[0]) || rtIsNaN(logWs[1])) {
    n = 0;
    o = rtNaN;
    apnd = logWs[1];
  } else if (logWs[0] < logWs[1]) {
    n = -1;
    o = logWs[0];
    apnd = logWs[1];
  } else if (rtIsInf(logWs[0]) || rtIsInf(logWs[1])) {
    n = 0;
    o = rtNaN;
    apnd = logWs[1];
  } else {
    o = logWs[0];
    ndbl = floor(-(logWs[1] - logWs[0]) + 0.5);
    apnd = logWs[0] + -ndbl;
    cdiff = logWs[1] - apnd;
    kd = fabs(logWs[0]);
    absb = fabs(logWs[1]);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
      ndbl++;
      apnd = logWs[1];
    } else if (cdiff > 0.0) {
      apnd = logWs[0] + -(ndbl - 1.0);
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  i15 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i15, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = o;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (b_cdiff = 1; b_cdiff < nm1d2; b_cdiff++) {
        b_y->data[b_cdiff] = o + -(double)b_cdiff;
        b_y->data[n - b_cdiff] = apnd - (-(double)b_cdiff);
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (o + apnd) / 2.0;
      } else {
        b_y->data[nm1d2] = o + -(double)nm1d2;
        b_y->data[nm1d2 + 1] = apnd - (-(double)nm1d2);
      }
    }
  }

  emxInit_real_T(&ws, 2);
  emxInit_real_T(&pO, 2);
  f_power(b_y, ws);

  /*  P2-WSs */
  f_rdivide(8.0 * fs, ws, pO);

  /*  Optimal pitches for P2-WSs */
  /*  Determine window sizes used by each pitch candidate */
  y = scalar_log2(b_rdivide(8.0 * fs, ws->data[0]));
  b_ndbl = log2pc_size[0];
  for (i15 = 0; i15 < b_ndbl; i15++) {
    log2pc_data[i15] = (1.0 + log2pc_data[i15]) - y;
  }

  /*  Create ERB-scale uniformly-spaced frequencies (in Hertz) */
  ixstart = 1;
  kd = pc_data[0];
  if (pc_size_idx_0 > 1) {
    if (rtIsNaN(pc_data[0])) {
      nm1d2 = 2;
      exitg9 = false;
      while ((!exitg9) && (nm1d2 <= pc_size_idx_0)) {
        ixstart = nm1d2;
        if (!rtIsNaN(pc_data[nm1d2 - 1])) {
          kd = pc_data[nm1d2 - 1];
          exitg9 = true;
        } else {
          nm1d2++;
        }
      }
    }

    if (ixstart < pc_size_idx_0) {
      while (ixstart + 1 <= pc_size_idx_0) {
        if (pc_data[ixstart] < kd) {
          kd = pc_data[ixstart];
        }

        ixstart++;
      }
    }
  }

  y = 6.44 * (scalar_log2(229.0 + kd / 4.0) - 7.84);
  o = 6.44 * (scalar_log2(229.0 + fs / 2.0) - 7.84);
  if (rtIsNaN(y) || rtIsNaN(o)) {
    n = 0;
    y = rtNaN;
    apnd = o;
  } else if (o < y) {
    n = -1;
    apnd = o;
  } else if (rtIsInf(y) || rtIsInf(o)) {
    n = 0;
    y = rtNaN;
    apnd = o;
  } else {
    ndbl = floor((o - y) / 0.1 + 0.5);
    apnd = y + ndbl * 0.1;
    cdiff = apnd - o;
    kd = fabs(y);
    absb = fabs(o);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
      ndbl++;
      apnd = o;
    } else if (cdiff > 0.0) {
      apnd = y + (ndbl - 1.0) * 0.1;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  i15 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i15, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = y;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (b_cdiff = 1; b_cdiff < nm1d2; b_cdiff++) {
        kd = (double)b_cdiff * 0.1;
        b_y->data[b_cdiff] = y + kd;
        b_y->data[n - b_cdiff] = apnd - kd;
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (y + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * 0.1;
        b_y->data[nm1d2] = y + kd;
        b_y->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  b_emxInit_real_T(&c_y, 1);
  i15 = c_y->size[0];
  c_y->size[0] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)c_y, i15, (int)sizeof(double));
  b_ndbl = b_y->size[1];
  for (i15 = 0; i15 < b_ndbl; i15++) {
    c_y->data[i15] = b_y->data[b_y->size[0] * i15];
  }

  b_emxInit_real_T(&b_xzp, 1);
  rdivide(c_y, 6.44, xzp);
  i15 = b_xzp->size[0];
  b_xzp->size[0] = xzp->size[0];
  emxEnsureCapacity((emxArray__common *)b_xzp, i15, (int)sizeof(double));
  b_ndbl = xzp->size[0];
  emxFree_real_T(&c_y);
  for (i15 = 0; i15 < b_ndbl; i15++) {
    b_xzp->data[i15] = xzp->data[i15] + 7.84;
  }

  b_emxInit_real_T(&fERBs, 1);
  e_power(b_xzp, fERBs);
  i15 = fERBs->size[0];
  emxEnsureCapacity((emxArray__common *)fERBs, i15, (int)sizeof(double));
  b_ndbl = fERBs->size[0];
  emxFree_real_T(&b_xzp);
  for (i15 = 0; i15 < b_ndbl; i15++) {
    fERBs->data[i15] -= 229.0;
  }

  i = 0;
  b_emxInit_real_T(&w, 1);
  emxInit_real_T(&X, 2);
  b_emxInit_real_T(&f, 1);
  b_emxInit_real_T(&ti, 1);
  emxInit_real_T(&L, 2);
  emxInit_real_T(&Si, 2);
  emxInit_real_T(&mu, 2);
  emxInit_real_T(&j, 2);
  emxInit_real_T(&k, 2);
  b_emxInit_real_T(&c_xzp, 1);
  emxInit_real_T(&r17, 2);
  emxInit_real_T(&r18, 2);
  emxInit_boolean_T(&c_x, 1);
  b_emxInit_boolean_T(&d_x, 2);
  emxInit_int32_T(&ii, 1);
  emxInit_real_T(&r19, 2);
  emxInit_real_T(&log2pc, 2);
  emxInit_real_T(&b_Si, 2);
  emxInit_real_T(&c_Si, 2);
  emxInit_real_T(&d_Si, 2);
  emxInit_real_T(&pc, 2);
  b_emxInit_real_T(&b_fERBs, 1);
  emxInit_real_T(&b_L, 2);
  emxInit_real_T(&b_S, 2);
  while (i <= ws->size[1] - 1) {
    kd = fmax(1.0, rt_roundd_snf(4.0 * fs / pO->data[i]));

    /*  Hop size */
    /*  Zero pad signal */
    y = ws->data[i] / 2.0;
    o = ws->data[i] / 2.0;
    i15 = xzp->size[0];
    xzp->size[0] = ((int)y + x->size[0]) + (int)(kd + o);
    emxEnsureCapacity((emxArray__common *)xzp, i15, (int)sizeof(double));
    b_ndbl = (int)y;
    for (i15 = 0; i15 < b_ndbl; i15++) {
      xzp->data[i15] = 0.0;
    }

    b_ndbl = x->size[0];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      xzp->data[i15 + (int)y] = x->data[i15];
    }

    b_ndbl = (int)(kd + o);
    for (i15 = 0; i15 < b_ndbl; i15++) {
      xzp->data[(i15 + (int)y) + x->size[0]] = 0.0;
    }

    /*  Compute spectrum */
    i15 = w->size[0];
    w->size[0] = (int)ws->data[i];
    emxEnsureCapacity((emxArray__common *)w, i15, (int)sizeof(double));
    b_ndbl = (int)ws->data[i];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      w->data[i15] = 0.0;
    }

    hanning(&w->data[0], ws->data[i]);
    o = fmax(0.0, rt_roundd_snf(ws->data[i] - kd));

    /*  Window overlap */
    kd = ws->data[i] - o;
    y = ws->data[i] / 2.0;
    kd = (((double)xzp->size[0] + kd) - 1.0) / kd;
    i15 = X->size[0] * X->size[1];
    X->size[0] = (int)(y + 1.0);
    X->size[1] = (int)kd;
    emxEnsureCapacity((emxArray__common *)X, i15, (int)sizeof(double));
    b_ndbl = (int)(y + 1.0) * (int)kd;
    for (i15 = 0; i15 < b_ndbl; i15++) {
      X->data[i15] = 0.0;
    }

    i15 = f->size[0];
    f->size[0] = (int)(y + 1.0);
    emxEnsureCapacity((emxArray__common *)f, i15, (int)sizeof(double));
    b_ndbl = (int)(y + 1.0);
    for (i15 = 0; i15 < b_ndbl; i15++) {
      f->data[i15] = 0.0;
    }

    i15 = ti->size[0];
    ti->size[0] = (int)kd;
    emxEnsureCapacity((emxArray__common *)ti, i15, (int)sizeof(double));
    b_ndbl = (int)kd;
    for (i15 = 0; i15 < b_ndbl; i15++) {
      ti->data[i15] = 0.0;
    }

    /*  use argument order of newer spectrogram function, for which we */
    /*  have documentation */
    i15 = c_xzp->size[0];
    c_xzp->size[0] = xzp->size[0];
    emxEnsureCapacity((emxArray__common *)c_xzp, i15, (int)sizeof(double));
    b_ndbl = xzp->size[0];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      c_xzp->data[i15] = xzp->data[i15];
    }

    spectrogram(&X->data[0], &f->data[0], &ti->data[0], &c_xzp->data[0], (double)
                xzp->size[0], &w->data[0], o, ws->data[i], fs);

    /*  Select candidates that use this window size */
    if (ws->size[1] == 1) {
      i15 = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = pc_size_idx_0;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      for (i15 = 0; i15 < pc_size_idx_0; i15++) {
        j->data[j->size[0] * i15] = pc_data[i15];
      }

      i15 = k->size[0] * k->size[1];
      k->size[0] = 0;
      k->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)k, i15, (int)sizeof(double));
    } else if (1 + i == ws->size[1]) {
      i15 = c_x->size[0];
      c_x->size[0] = log2pc_size[0];
      emxEnsureCapacity((emxArray__common *)c_x, i15, (int)sizeof(boolean_T));
      b_ndbl = log2pc_size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        c_x->data[i15] = (log2pc_data[i15] - (1.0 + (double)i) > -1.0);
      }

      b_ndbl = 0;
      i15 = ii->size[0];
      ii->size[0] = 160;
      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      ixstart = 1;
      exitg8 = false;
      while ((!exitg8) && (ixstart <= 160)) {
        guard5 = false;
        if (c_x->data[ixstart - 1]) {
          b_ndbl++;
          ii->data[b_ndbl - 1] = ixstart;
          if (b_ndbl >= 160) {
            exitg8 = true;
          } else {
            guard5 = true;
          }
        } else {
          guard5 = true;
        }

        if (guard5) {
          ixstart++;
        }
      }

      i15 = ii->size[0];
      if (1 > b_ndbl) {
        ii->size[0] = 0;
      } else {
        ii->size[0] = b_ndbl;
      }

      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      i15 = xzp->size[0];
      xzp->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)xzp, i15, (int)sizeof(double));
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        xzp->data[i15] = ii->data[i15];
      }

      ixstart = ii->size[0];
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        b_tmp_data[i15] = ii->data[i15];
      }

      nm1d2 = xzp->size[0];
      i15 = j->size[0] * j->size[1];
      j->size[0] = nm1d2;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      i15 = j->size[0] * j->size[1];
      j->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        j->data[i15] = xzp->data[i15];
      }

      i15 = d_x->size[0] * d_x->size[1];
      d_x->size[0] = ixstart;
      emxEnsureCapacity((emxArray__common *)d_x, i15, (int)sizeof(boolean_T));
      i15 = d_x->size[0] * d_x->size[1];
      d_x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)d_x, i15, (int)sizeof(boolean_T));
      for (i15 = 0; i15 < ixstart; i15++) {
        d_x->data[i15] = (log2pc_data[b_tmp_data[i15] - 1] - (1.0 + (double)i) <
                          0.0);
      }

      nm1d2 = d_x->size[0];
      b_ndbl = 0;
      i15 = ii->size[0];
      ii->size[0] = d_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      ixstart = 1;
      exitg7 = false;
      while ((!exitg7) && (ixstart <= nm1d2)) {
        guard4 = false;
        if (d_x->data[ixstart - 1]) {
          b_ndbl++;
          ii->data[b_ndbl - 1] = ixstart;
          if (b_ndbl >= nm1d2) {
            exitg7 = true;
          } else {
            guard4 = true;
          }
        } else {
          guard4 = true;
        }

        if (guard4) {
          ixstart++;
        }
      }

      if (d_x->size[0] == 1) {
        if (b_ndbl == 0) {
          i15 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
        }
      } else {
        i15 = ii->size[0];
        if (1 > b_ndbl) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = b_ndbl;
        }

        emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      }

      i15 = xzp->size[0];
      xzp->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)xzp, i15, (int)sizeof(double));
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        xzp->data[i15] = ii->data[i15];
      }

      nm1d2 = xzp->size[0];
      i15 = k->size[0] * k->size[1];
      k->size[0] = nm1d2;
      emxEnsureCapacity((emxArray__common *)k, i15, (int)sizeof(double));
      i15 = k->size[0] * k->size[1];
      k->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)k, i15, (int)sizeof(double));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        k->data[i15] = xzp->data[i15];
      }
    } else if (1 + i == 1) {
      i15 = c_x->size[0];
      c_x->size[0] = log2pc_size[0];
      emxEnsureCapacity((emxArray__common *)c_x, i15, (int)sizeof(boolean_T));
      b_ndbl = log2pc_size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        c_x->data[i15] = (log2pc_data[i15] - 1.0 < 1.0);
      }

      b_ndbl = 0;
      i15 = ii->size[0];
      ii->size[0] = 160;
      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      ixstart = 1;
      exitg6 = false;
      while ((!exitg6) && (ixstart <= 160)) {
        guard3 = false;
        if (c_x->data[ixstart - 1]) {
          b_ndbl++;
          ii->data[b_ndbl - 1] = ixstart;
          if (b_ndbl >= 160) {
            exitg6 = true;
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }

        if (guard3) {
          ixstart++;
        }
      }

      i15 = ii->size[0];
      if (1 > b_ndbl) {
        ii->size[0] = 0;
      } else {
        ii->size[0] = b_ndbl;
      }

      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      i15 = xzp->size[0];
      xzp->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)xzp, i15, (int)sizeof(double));
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        xzp->data[i15] = ii->data[i15];
      }

      ixstart = ii->size[0];
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        b_tmp_data[i15] = ii->data[i15];
      }

      nm1d2 = xzp->size[0];
      i15 = j->size[0] * j->size[1];
      j->size[0] = nm1d2;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      i15 = j->size[0] * j->size[1];
      j->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        j->data[i15] = xzp->data[i15];
      }

      i15 = d_x->size[0] * d_x->size[1];
      d_x->size[0] = ixstart;
      emxEnsureCapacity((emxArray__common *)d_x, i15, (int)sizeof(boolean_T));
      i15 = d_x->size[0] * d_x->size[1];
      d_x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)d_x, i15, (int)sizeof(boolean_T));
      for (i15 = 0; i15 < ixstart; i15++) {
        d_x->data[i15] = (log2pc_data[b_tmp_data[i15] - 1] - 1.0 > 0.0);
      }

      nm1d2 = d_x->size[0];
      b_ndbl = 0;
      i15 = ii->size[0];
      ii->size[0] = d_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      ixstart = 1;
      exitg5 = false;
      while ((!exitg5) && (ixstart <= nm1d2)) {
        guard2 = false;
        if (d_x->data[ixstart - 1]) {
          b_ndbl++;
          ii->data[b_ndbl - 1] = ixstart;
          if (b_ndbl >= nm1d2) {
            exitg5 = true;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }

        if (guard2) {
          ixstart++;
        }
      }

      if (d_x->size[0] == 1) {
        if (b_ndbl == 0) {
          i15 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
        }
      } else {
        i15 = ii->size[0];
        if (1 > b_ndbl) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = b_ndbl;
        }

        emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      }

      i15 = xzp->size[0];
      xzp->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)xzp, i15, (int)sizeof(double));
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        xzp->data[i15] = ii->data[i15];
      }

      nm1d2 = xzp->size[0];
      i15 = k->size[0] * k->size[1];
      k->size[0] = nm1d2;
      emxEnsureCapacity((emxArray__common *)k, i15, (int)sizeof(double));
      i15 = k->size[0] * k->size[1];
      k->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)k, i15, (int)sizeof(double));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        k->data[i15] = xzp->data[i15];
      }
    } else {
      b_log2pc_size[0] = log2pc_size[0];
      b_ndbl = log2pc_size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        c_log2pc_data[i15] = log2pc_data[i15] - (1.0 + (double)i);
      }

      d_log2pc_data.data = (double *)&c_log2pc_data;
      d_log2pc_data.size = (int *)&b_log2pc_size;
      d_log2pc_data.allocatedSize = 160;
      d_log2pc_data.numDimensions = 1;
      d_log2pc_data.canFreeData = false;
      b_abs(&d_log2pc_data, xzp);
      i15 = c_x->size[0];
      c_x->size[0] = xzp->size[0];
      emxEnsureCapacity((emxArray__common *)c_x, i15, (int)sizeof(boolean_T));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        c_x->data[i15] = (xzp->data[i15] < 1.0);
      }

      nm1d2 = c_x->size[0];
      b_ndbl = 0;
      i15 = ii->size[0];
      ii->size[0] = c_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      ixstart = 1;
      exitg4 = false;
      while ((!exitg4) && (ixstart <= nm1d2)) {
        guard1 = false;
        if (c_x->data[ixstart - 1]) {
          b_ndbl++;
          ii->data[b_ndbl - 1] = ixstart;
          if (b_ndbl >= nm1d2) {
            exitg4 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          ixstart++;
        }
      }

      if (c_x->size[0] == 1) {
        if (b_ndbl == 0) {
          i15 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
        }
      } else {
        i15 = ii->size[0];
        if (1 > b_ndbl) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = b_ndbl;
        }

        emxEnsureCapacity((emxArray__common *)ii, i15, (int)sizeof(int));
      }

      i15 = xzp->size[0];
      xzp->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)xzp, i15, (int)sizeof(double));
      b_ndbl = ii->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        xzp->data[i15] = ii->data[i15];
      }

      nm1d2 = xzp->size[0];
      i15 = j->size[0] * j->size[1];
      j->size[0] = nm1d2;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      i15 = j->size[0] * j->size[1];
      j->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)j, i15, (int)sizeof(double));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        j->data[i15] = xzp->data[i15];
      }

      nm1d2 = j->size[0];
      if (nm1d2 >= 1) {
      } else {
        nm1d2 = 1;
      }

      if (0 == j->size[0]) {
        n = 0;
      } else {
        n = nm1d2;
      }

      if (n < 1) {
        n = -1;
        b_apnd = 0;
      } else {
        b_ndbl = (int)floor(((double)n - 1.0) + 0.5);
        b_apnd = b_ndbl + 1;
        b_cdiff = (b_ndbl - n) + 1;
        if (fabs(b_cdiff) < 4.4408920985006262E-16 * (double)n) {
          b_ndbl++;
          b_apnd = n;
        } else if (b_cdiff > 0) {
          b_apnd = b_ndbl;
        } else {
          b_ndbl++;
        }

        n = b_ndbl - 1;
      }

      i15 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = n + 1;
      emxEnsureCapacity((emxArray__common *)b_y, i15, (int)sizeof(double));
      if (n + 1 > 0) {
        b_y->data[0] = 1.0;
        if (n + 1 > 1) {
          b_y->data[n] = b_apnd;
          nm1d2 = (n + (n < 0)) >> 1;
          for (b_cdiff = 1; b_cdiff < nm1d2; b_cdiff++) {
            b_y->data[b_cdiff] = 1.0 + (double)b_cdiff;
            b_y->data[n - b_cdiff] = b_apnd - b_cdiff;
          }

          if (nm1d2 << 1 == n) {
            b_y->data[nm1d2] = (1.0 + (double)b_apnd) / 2.0;
          } else {
            b_y->data[nm1d2] = 1.0 + (double)nm1d2;
            b_y->data[nm1d2 + 1] = b_apnd - nm1d2;
          }
        }
      }

      i15 = k->size[0] * k->size[1];
      k->size[0] = 1;
      k->size[1] = b_y->size[1];
      emxEnsureCapacity((emxArray__common *)k, i15, (int)sizeof(double));
      b_ndbl = b_y->size[0] * b_y->size[1];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        k->data[i15] = b_y->data[i15];
      }
    }

    /*  Compute loudness at ERBs uniformly-spaced frequencies */
    y = pc_data[(int)j->data[0] - 1] / 4.0;
    i15 = c_x->size[0];
    c_x->size[0] = fERBs->size[0];
    emxEnsureCapacity((emxArray__common *)c_x, i15, (int)sizeof(boolean_T));
    b_ndbl = fERBs->size[0];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      c_x->data[i15] = (fERBs->data[i15] > y);
    }

    if (1 <= c_x->size[0]) {
      b_cdiff = 1;
    } else {
      b_cdiff = 0;
    }

    b_ndbl = 0;
    ixstart = 1;
    exitg3 = false;
    while ((!exitg3) && (ixstart <= c_x->size[0])) {
      if (c_x->data[ixstart - 1]) {
        b_ndbl = 1;
        ii_data[0] = ixstart;
        exitg3 = true;
      } else {
        ixstart++;
      }
    }

    if (b_cdiff == 1) {
      if (b_ndbl == 0) {
        b_cdiff = 0;
      }
    } else if (1 > b_ndbl) {
      b_cdiff = 0;
    } else {
      b_cdiff = 1;
    }

    i15 = 0;
    while (i15 <= b_cdiff - 1) {
      found_data[0] = ii_data[0];
      i15 = 1;
    }

    if (found_data[0] > fERBs->size[0]) {
      i15 = 0;
      b_cdiff = 0;
    } else {
      i15 = found_data[0] - 1;
      b_cdiff = fERBs->size[0];
    }

    nm1d2 = b_fERBs->size[0];
    b_fERBs->size[0] = b_cdiff - i15;
    emxEnsureCapacity((emxArray__common *)b_fERBs, nm1d2, (int)sizeof(double));
    b_ndbl = b_cdiff - i15;
    for (b_cdiff = 0; b_cdiff < b_ndbl; b_cdiff++) {
      b_fERBs->data[b_cdiff] = fERBs->data[i15 + b_cdiff];
    }

    i15 = fERBs->size[0];
    fERBs->size[0] = b_fERBs->size[0];
    emxEnsureCapacity((emxArray__common *)fERBs, i15, (int)sizeof(double));
    b_ndbl = b_fERBs->size[0];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      fERBs->data[i15] = b_fERBs->data[i15];
    }

    c_abs(X, r19);
    interp1(f, r19, fERBs, X);
    for (i15 = 0; i15 < 2; i15++) {
      logWs[i15] = X->size[i15];
    }

    i15 = L->size[0] * L->size[1];
    L->size[0] = (int)logWs[0];
    L->size[1] = (int)logWs[1];
    emxEnsureCapacity((emxArray__common *)L, i15, (int)sizeof(double));
    i15 = (int)logWs[0] * (int)logWs[1];
    for (b_cdiff = 0; b_cdiff + 1 <= i15; b_cdiff++) {
      L->data[b_cdiff] = fmax(0.0, X->data[b_cdiff]);
    }

    c_sqrt(L);

    /*  Compute pitch strength */
    i15 = pc->size[0] * pc->size[1];
    pc->size[0] = j->size[0];
    pc->size[1] = j->size[1];
    emxEnsureCapacity((emxArray__common *)pc, i15, (int)sizeof(double));
    b_ndbl = j->size[0] * j->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      pc->data[i15] = pc_data[(int)j->data[i15] - 1];
    }

    pitchStrengthAllCandidates(fERBs, L, pc, r19);
    i15 = Si->size[0] * Si->size[1];
    Si->size[0] = r19->size[0];
    Si->size[1] = r19->size[1];
    emxEnsureCapacity((emxArray__common *)Si, i15, (int)sizeof(double));
    b_ndbl = r19->size[0] * r19->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      Si->data[i15] = r19->data[i15];
    }

    /*  Interpolate pitch strength at desired times */
    if (Si->size[1] > 1) {
      i15 = d_Si->size[0] * d_Si->size[1];
      d_Si->size[0] = Si->size[1];
      d_Si->size[1] = Si->size[0];
      emxEnsureCapacity((emxArray__common *)d_Si, i15, (int)sizeof(double));
      b_ndbl = Si->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        nm1d2 = Si->size[1];
        for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
          d_Si->data[b_cdiff + d_Si->size[0] * i15] = Si->data[i15 + Si->size[0]
            * b_cdiff];
        }
      }

      b_interp1(ti, d_Si, t, r19);
      i15 = c_Si->size[0] * c_Si->size[1];
      c_Si->size[0] = Si->size[1];
      c_Si->size[1] = Si->size[0];
      emxEnsureCapacity((emxArray__common *)c_Si, i15, (int)sizeof(double));
      b_ndbl = Si->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        nm1d2 = Si->size[1];
        for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
          c_Si->data[b_cdiff + c_Si->size[0] * i15] = Si->data[i15 + Si->size[0]
            * b_cdiff];
        }
      }

      b_interp1(ti, c_Si, t, X);
      i15 = b_Si->size[0] * b_Si->size[1];
      b_Si->size[0] = Si->size[1];
      b_Si->size[1] = Si->size[0];
      emxEnsureCapacity((emxArray__common *)b_Si, i15, (int)sizeof(double));
      b_ndbl = Si->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        nm1d2 = Si->size[1];
        for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
          b_Si->data[b_cdiff + b_Si->size[0] * i15] = Si->data[i15 + Si->size[0]
            * b_cdiff];
        }
      }

      b_interp1(ti, b_Si, t, L);
      i15 = b_L->size[0] * b_L->size[1];
      b_L->size[0] = L->size[1];
      b_L->size[1] = L->size[0];
      emxEnsureCapacity((emxArray__common *)b_L, i15, (int)sizeof(double));
      b_ndbl = L->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        nm1d2 = L->size[1];
        for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
          b_L->data[b_cdiff + b_L->size[0] * i15] = L->data[i15 + L->size[0] *
            b_cdiff];
        }
      }

      nm1d2 = r19->size[1];
      ixstart = X->size[0];
      i15 = Si->size[0] * Si->size[1];
      Si->size[0] = nm1d2;
      Si->size[1] = ixstart;
      emxEnsureCapacity((emxArray__common *)Si, i15, (int)sizeof(double));
      for (i15 = 0; i15 < ixstart; i15++) {
        for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
          Si->data[b_cdiff + Si->size[0] * i15] = b_L->data[b_cdiff + nm1d2 *
            i15];
        }
      }
    } else {
      if ((0 == Si->size[0]) || (0 == Si->size[1])) {
        ixstart = 0;
      } else {
        nm1d2 = Si->size[0];
        ixstart = Si->size[1];
        if (nm1d2 >= ixstart) {
          ixstart = nm1d2;
        }
      }

      d_repmat(ixstart, t->size[0], r19);
      i15 = Si->size[0] * Si->size[1];
      Si->size[0] = r19->size[0];
      Si->size[1] = r19->size[1];
      emxEnsureCapacity((emxArray__common *)Si, i15, (int)sizeof(double));
      b_ndbl = r19->size[0] * r19->size[1];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        Si->data[i15] = r19->data[i15];
      }
    }

    /*  Add pitch strength to combination */
    for (i15 = 0; i15 < 2; i15++) {
      logWs[i15] = j->size[i15];
    }

    i15 = mu->size[0] * mu->size[1];
    mu->size[0] = (int)logWs[0];
    emxEnsureCapacity((emxArray__common *)mu, i15, (int)sizeof(double));
    i15 = mu->size[0] * mu->size[1];
    mu->size[1] = (int)logWs[1];
    emxEnsureCapacity((emxArray__common *)mu, i15, (int)sizeof(double));
    b_ndbl = (int)logWs[0] * (int)logWs[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      mu->data[i15] = 1.0;
    }

    i15 = log2pc->size[0] * log2pc->size[1];
    log2pc->size[0] = k->size[0];
    log2pc->size[1] = k->size[1];
    emxEnsureCapacity((emxArray__common *)log2pc, i15, (int)sizeof(double));
    b_ndbl = k->size[0] * k->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      log2pc->data[i15] = log2pc_data[(int)j->data[(int)k->data[i15] - 1] - 1] -
        (1.0 + (double)i);
    }

    c_abs(log2pc, r19);
    i15 = r17->size[0] * r17->size[1];
    r17->size[0] = r19->size[0];
    r17->size[1] = r19->size[1];
    emxEnsureCapacity((emxArray__common *)r17, i15, (int)sizeof(double));
    b_ndbl = r19->size[0] * r19->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      r17->data[i15] = r19->data[i15];
    }

    b_ndbl = r17->size[0] * r17->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      mu->data[(int)k->data[i15] - 1] = 1.0 - r17->data[i15];
    }

    e_repmat(mu, Si->size[1], r19);
    i15 = r18->size[0] * r18->size[1];
    r18->size[0] = r19->size[0];
    r18->size[1] = r19->size[1];
    emxEnsureCapacity((emxArray__common *)r18, i15, (int)sizeof(double));
    b_ndbl = r19->size[0] * r19->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      r18->data[i15] = r19->data[i15];
    }

    nm1d2 = j->size[0] * j->size[1];
    ixstart = S->size[1];
    i15 = b_S->size[0] * b_S->size[1];
    b_S->size[0] = nm1d2;
    b_S->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)b_S, i15, (int)sizeof(double));
    for (i15 = 0; i15 < ixstart; i15++) {
      for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
        b_S->data[b_cdiff + b_S->size[0] * i15] = S->data[((int)j->data[b_cdiff]
          + S->size[0] * i15) - 1] + r18->data[b_cdiff + r18->size[0] * i15] *
          Si->data[b_cdiff + Si->size[0] * i15];
      }
    }

    b_ndbl = b_S->size[1];
    for (i15 = 0; i15 < b_ndbl; i15++) {
      nm1d2 = b_S->size[0];
      for (b_cdiff = 0; b_cdiff < nm1d2; b_cdiff++) {
        S->data[((int)j->data[b_cdiff] + S->size[0] * i15) - 1] = b_S->
          data[b_cdiff + b_S->size[0] * i15];
      }
    }

    i++;
  }

  emxFree_real_T(&b_S);
  emxFree_real_T(&b_L);
  emxFree_real_T(&b_fERBs);
  emxFree_real_T(&pc);
  emxFree_real_T(&d_Si);
  emxFree_real_T(&c_Si);
  emxFree_real_T(&b_Si);
  emxFree_real_T(&log2pc);
  emxFree_real_T(&r19);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&d_x);
  emxFree_boolean_T(&c_x);
  emxFree_real_T(&r18);
  emxFree_real_T(&r17);
  emxFree_real_T(&c_xzp);
  emxFree_real_T(&k);
  emxFree_real_T(&j);
  emxFree_real_T(&mu);
  emxFree_real_T(&Si);
  emxFree_real_T(&L);
  emxFree_real_T(&ti);
  emxFree_real_T(&f);
  emxFree_real_T(&X);
  emxFree_real_T(&w);
  emxFree_real_T(&fERBs);

  /*  Fine tune pitch using parabolic interpolation */
  f_repmat(S->size[1], p);
  f_repmat(S->size[1], s);
  b_j = 0;
  emxInit_real_T(&b_pO, 2);
  b_emxInit_real_T(&d_xzp, 1);
  while (b_j <= S->size[1] - 1) {
    ixstart = 1;
    n = S->size[0];
    kd = S->data[S->size[0] * b_j];
    i = 1;
    i15 = S->size[0];
    if (i15 > 1) {
      if (rtIsNaN(kd)) {
        nm1d2 = 2;
        exitg2 = false;
        while ((!exitg2) && (nm1d2 <= n)) {
          ixstart = nm1d2;
          if (!rtIsNaN(S->data[(nm1d2 + S->size[0] * b_j) - 1])) {
            kd = S->data[(nm1d2 + S->size[0] * b_j) - 1];
            i = nm1d2;
            exitg2 = true;
          } else {
            nm1d2++;
          }
        }
      }

      i15 = S->size[0];
      if (ixstart < i15) {
        while (ixstart + 1 <= n) {
          if (S->data[ixstart + S->size[0] * b_j] > kd) {
            kd = S->data[ixstart + S->size[0] * b_j];
            i = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    s->data[b_j] = kd;
    if ((i == 1) || (i == pc_size_idx_0)) {
      p->data[b_j] = pc_data[i - 1];
    } else {
      if (i + 1 < i - 1) {
        n = -1;
        o = (double)i - 1.0;
        b_apnd = i + 1;
      } else {
        o = (double)i - 1.0;
        b_ndbl = 2;
        b_apnd = i + 1;
        b_cdiff = (b_apnd - i) - 1;
        nm1d2 = i - 1;
        ixstart = i + 1;
        if (nm1d2 >= ixstart) {
          ixstart = nm1d2;
        }

        if (fabs(b_cdiff) < 4.4408920985006262E-16 * (double)ixstart) {
          b_ndbl = 3;
          b_apnd = i + 1;
        } else if (b_cdiff > 0) {
          b_apnd = i;
        } else {
          b_ndbl = 3;
        }

        n = b_ndbl - 1;
      }

      i15 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = n + 1;
      emxEnsureCapacity((emxArray__common *)b_y, i15, (int)sizeof(double));
      if (n + 1 > 0) {
        b_y->data[0] = o;
        if (n + 1 > 1) {
          b_y->data[n] = b_apnd;
          nm1d2 = (n + (n < 0)) >> 1;
          if (nm1d2 << 1 == n) {
            b_y->data[nm1d2] = (o + (double)b_apnd) / 2.0;
          } else {
            b_y->data[nm1d2] = o + (double)nm1d2;
            b_y->data[nm1d2 + 1] = b_apnd - nm1d2;
          }
        }
      }

      pc_size[0] = b_y->size[1];
      b_ndbl = b_y->size[1];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        b_pc_data[i15] = pc_data[(int)b_y->data[b_y->size[0] * i15] - 1];
      }

      c_pc_data.data = (double *)&b_pc_data;
      c_pc_data.size = (int *)&pc_size;
      c_pc_data.allocatedSize = 161;
      c_pc_data.numDimensions = 1;
      c_pc_data.canFreeData = false;
      g_rdivide(&c_pc_data, xzp);
      tc_size[0] = xzp->size[0];
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        tc_data[i15] = xzp->data[i15];
      }

      b_tc_data.data = (double *)&tc_data;
      b_tc_data.size = (int *)&tc_size;
      b_tc_data.allocatedSize = 161;
      b_tc_data.numDimensions = 1;
      b_tc_data.canFreeData = false;
      rdivide(&b_tc_data, tc_data[1], xzp);
      i15 = d_xzp->size[0];
      d_xzp->size[0] = xzp->size[0];
      emxEnsureCapacity((emxArray__common *)d_xzp, i15, (int)sizeof(double));
      b_ndbl = xzp->size[0];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        d_xzp->data[i15] = (xzp->data[i15] - 1.0) * 2.0 * 3.1415926535897931;
      }

      S_size[0] = b_y->size[1];
      b_ndbl = b_y->size[1];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        S_data[i15] = S->data[((int)b_y->data[b_y->size[0] * i15] + S->size[0] *
          b_j) - 1];
      }

      b_S_data.data = (double *)&S_data;
      b_S_data.size = (int *)&S_size;
      b_S_data.allocatedSize = 161;
      b_S_data.numDimensions = 1;
      b_S_data.canFreeData = false;
      b_polyfit(d_xzp, &b_S_data, c);
      y = scalar_log2(pc_data[(int)b_y->data[0] - 1]);
      o = scalar_log2(pc_data[(int)b_y->data[2] - 1]);
      if (rtIsNaN(y) || rtIsNaN(o)) {
        n = 0;
        y = rtNaN;
        apnd = o;
      } else if (o < y) {
        n = -1;
        apnd = o;
      } else if (rtIsInf(y) || rtIsInf(o)) {
        n = 0;
        y = rtNaN;
        apnd = o;
      } else {
        ndbl = floor((o - y) / 0.00083333333333333328 + 0.5);
        apnd = y + ndbl * 0.00083333333333333328;
        cdiff = apnd - o;
        kd = fabs(y);
        absb = fabs(o);
        if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
          ndbl++;
          apnd = o;
        } else if (cdiff > 0.0) {
          apnd = y + (ndbl - 1.0) * 0.00083333333333333328;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl - 1;
        } else {
          n = -1;
        }
      }

      i15 = ws->size[0] * ws->size[1];
      ws->size[0] = 1;
      ws->size[1] = n + 1;
      emxEnsureCapacity((emxArray__common *)ws, i15, (int)sizeof(double));
      if (n + 1 > 0) {
        ws->data[0] = y;
        if (n + 1 > 1) {
          ws->data[n] = apnd;
          nm1d2 = (n + (n < 0)) >> 1;
          for (b_cdiff = 1; b_cdiff < nm1d2; b_cdiff++) {
            kd = (double)b_cdiff * 0.00083333333333333328;
            ws->data[b_cdiff] = y + kd;
            ws->data[n - b_cdiff] = apnd - kd;
          }

          if (nm1d2 << 1 == n) {
            ws->data[nm1d2] = (y + apnd) / 2.0;
          } else {
            kd = (double)nm1d2 * 0.00083333333333333328;
            ws->data[nm1d2] = y + kd;
            ws->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }

      f_power(ws, pO);
      f_rdivide(1.0, pO, ws);
      c_rdivide(ws, tc_data[1], pO);
      i15 = b_pO->size[0] * b_pO->size[1];
      b_pO->size[0] = 1;
      b_pO->size[1] = pO->size[1];
      emxEnsureCapacity((emxArray__common *)b_pO, i15, (int)sizeof(double));
      b_ndbl = pO->size[0] * pO->size[1];
      for (i15 = 0; i15 < b_ndbl; i15++) {
        b_pO->data[i15] = (pO->data[i15] - 1.0) * 2.0 * 3.1415926535897931;
      }

      polyval(c, b_pO, ws);
      ixstart = 1;
      n = ws->size[1];
      kd = ws->data[0];
      i = 1;
      if (ws->size[1] > 1) {
        if (rtIsNaN(ws->data[0])) {
          nm1d2 = 2;
          exitg1 = false;
          while ((!exitg1) && (nm1d2 <= n)) {
            ixstart = nm1d2;
            if (!rtIsNaN(ws->data[nm1d2 - 1])) {
              kd = ws->data[nm1d2 - 1];
              i = nm1d2;
              exitg1 = true;
            } else {
              nm1d2++;
            }
          }
        }

        if (ixstart < ws->size[1]) {
          while (ixstart + 1 <= n) {
            if (ws->data[ixstart] > kd) {
              kd = ws->data[ixstart];
              i = ixstart + 1;
            }

            ixstart++;
          }
        }
      }

      s->data[b_j] = kd;
      kd = pc_data[(int)b_y->data[0] - 1];
      p->data[b_j] = rt_powd_snf(2.0, scalar_log2(kd) + ((double)i - 1.0) / 12.0
        / 100.0);
    }

    b_j++;
  }

  emxFree_real_T(&d_xzp);
  emxFree_real_T(&b_pO);
  emxFree_real_T(&b_y);
  emxFree_real_T(&xzp);
  emxFree_real_T(&pO);
  emxFree_real_T(&ws);
  emxFree_real_T(&S);
}

void swipep(const emxArray_real_T *x, double fs, const double plim[2], double dt,
            double dlog2p, double dERBs, double woverlap, double sTHR,
            emxArray_real_T *p, emxArray_real_T *t, emxArray_real_T *s)
{
  double y;
  int n;
  double o;
  double apnd;
  double ndbl;
  double cdiff;
  emxArray_real_T *b_y;
  int i28;
  int nm1d2;
  int k;
  double kd;
  int ixstart;
  int tmp_size[2];
  double tmp_data[1229];
  double absb;
  emxArray_real_T *log2pc;
  emxArray_real_T *pc;
  emxArray_real_T *S;
  double dv4[2];
  double logWs[2];
  emxArray_real_T *ws;
  emxArray_real_T *pO;
  boolean_T exitg9;
  emxArray_real_T *c_y;
  emxArray_real_T *w;
  emxArray_real_T *b_w;
  emxArray_real_T *fERBs;
  int i;
  emxArray_real_T *X;
  emxArray_real_T *f;
  emxArray_real_T *ti;
  emxArray_real_T *L;
  emxArray_real_T *Si;
  emxArray_real_T *tc;
  emxArray_real_T *j;
  emxArray_real_T *b_k;
  emxArray_real_T *xzp;
  emxArray_boolean_T *b_x;
  emxArray_boolean_T *c_x;
  emxArray_int32_T *ii;
  emxArray_real_T *b_log2pc;
  emxArray_real_T *b_Si;
  emxArray_real_T *b_pc;
  emxArray_real_T *c_log2pc;
  emxArray_real_T *b_fERBs;
  emxArray_real_T *b_S;
  int idx;
  boolean_T exitg8;
  boolean_T guard5 = false;
  boolean_T exitg7;
  boolean_T guard4 = false;
  boolean_T exitg6;
  boolean_T guard3 = false;
  boolean_T exitg5;
  boolean_T guard2 = false;
  boolean_T exitg4;
  boolean_T guard1 = false;
  int b_apnd;
  int ii_data[1];
  boolean_T exitg3;
  int found_data[1];
  emxArray_real_T *I;
  emxArray_real_T *b_pO;
  emxArray_real_T *c_w;
  emxArray_real_T *c_S;
  emxArray_real_T *c_pc;
  boolean_T exitg2;
  unsigned int b_absb;
  unsigned int u0;
  double c[3];
  boolean_T exitg1;

  /*  SWIPEP Pitch estimation using SWIPE'. */
  /*     P = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,STHR) estimates the pitch  */
  /*     of the vector signal X every DT seconds. The sampling frequency of */
  /*     the signal is Fs (in Hertz). The spectrum is computed using a Hann */
  /*     window with an overlap WOVERLAP between 0 and 1. The spectrum is */
  /*     sampled uniformly in the ERB scale with a step size of DERBS ERBs. The */
  /*     pitch is searched within the range [PMIN PMAX] (in Hertz) with samples */
  /*     distributed every DLOG2P units on a base-2 logarithmic scale of Hertz.  */
  /*     The pitch is fine-tuned using parabolic interpolation with a resolution */
  /*     of 1 cent. Pitch estimates with a strength lower than STHR are treated */
  /*     as undefined. */
  /*      */
  /*     [P,T,S] = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,WOVERLAP,STHR)  */
  /*     returns the times T at which the pitch was estimated and the pitch  */
  /*     strength S of every pitch estimate. */
  /*  */
  /*     P = SWIPEP(X,Fs) estimates the pitch using the default settings PMIN = */
  /*     30 Hz, PMAX = 5000 Hz, DT = 0.001 s, DLOG2P = 1/48 (48 steps per  */
  /*     octave), DERBS = 0.1 ERBs, WOVERLAP = 0.5, and STHR = -Inf. */
  /*  */
  /*     P = SWIPEP(X,Fs,...,[],...) uses the default setting for the parameter */
  /*     replaced with the placeholder []. */
  /*  */
  /*     REMARKS: (1) For better results, make DLOG2P and DERBS as small as  */
  /*     possible and WOVERLAP as large as possible. However, take into account */
  /*     that the computational complexity of the algorithm is inversely  */
  /*     proportional to DLOG2P, DERBS and 1-WOVERLAP, and that the  default  */
  /*     values have been found empirically to produce good results. Consider  */
  /*     also that the computational complexity is directly proportional to the */
  /*     number of octaves in the pitch search range, and therefore , it is  */
  /*     recommendable to restrict the search range to the expected range of */
  /*     pitch, if any. (2) This code implements SWIPE', which uses only the */
  /*     first and prime harmonics of the signal. To convert it into SWIPE, */
  /*     which uses all the harmonics of the signal, replace the word */
  /*     PRIMES with a colon (it is located almost at the end of the code). */
  /*     However, this may not be recommendable since SWIPE' is reported to  */
  /*     produce on average better results than SWIPE (Camacho and Harris, */
  /*     2008). */
  /*  */
  /*     EXAMPLE: Estimate the pitch of the signal X every 10 ms within the */
  /*     range 75-500 Hz using the default resolution (i.e., 48 steps per */
  /*     octave), sampling the spectrum every 1/20th of ERB, using a window  */
  /*     overlap factor of 50%, and discarding samples with pitch strength  */
  /*     lower than 0.2. Plot the pitch trace. */
  /*        [x,Fs] = wavread(filename); */
  /*        [p,t,s] = swipep(x,Fs,[75 500],0.01,[],1/20,0.5,0.2); */
  /*        plot(1000*t,p) */
  /*        xlabel('Time (ms)') */
  /*        ylabel('Pitch (Hz)') */
  /*  */
  /*     REFERENCES: Camacho, A., Harris, J.G, (2008) "A sawtooth waveform  */
  /*     inspired pitch estimator for speech and music," J. Acoust. Soc. Am. */
  /*     124, 1638-1652. */
  /*  */
  /*     MAINTENANCE HISTORY: */
  /*     - Added line 153 to avoid division by zero in line 154 if loudness */
  /*       equals zero (06/23/2010). */
  if (rtIsNaN(dt)) {
    dt = 0.001;
  }

  if (rtIsNaN(dlog2p)) {
    dlog2p = 0.020833333333333332;
  }

  if (rtIsNaN(dERBs)) {
    dERBs = 0.1;
  }

  if (rtIsNaN(woverlap)) {
    woverlap = 0.5;
  }

  if (rtIsNaN(sTHR)) {
    sTHR = rtMinusInf;
  }

  y = (double)x->size[0] / fs;
  if (rtIsNaN(dt) || rtIsNaN(y)) {
    n = 0;
    o = rtNaN;
    apnd = y;
  } else if ((dt == 0.0) || ((0.0 < y) && (dt < 0.0)) || ((y < 0.0) && (dt > 0.0)))
  {
    n = -1;
    o = 0.0;
    apnd = y;
  } else if (rtIsInf(y)) {
    n = 0;
    o = rtNaN;
    apnd = y;
  } else if (rtIsInf(dt)) {
    n = 0;
    o = 0.0;
    apnd = y;
  } else {
    o = 0.0;
    ndbl = floor(y / dt + 0.5);
    apnd = ndbl * dt;
    if (dt > 0.0) {
      cdiff = apnd - y;
    } else {
      cdiff = y - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(y)) {
      ndbl++;
      apnd = y;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * dt;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  emxInit_real_T(&b_y, 2);
  i28 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i28, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = o;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (k = 1; k < nm1d2; k++) {
        kd = (double)k * dt;
        b_y->data[k] = o + kd;
        b_y->data[n - k] = apnd - kd;
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (o + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * dt;
        b_y->data[nm1d2] = o + kd;
        b_y->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  i28 = t->size[0];
  t->size[0] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)t, i28, (int)sizeof(double));
  ixstart = b_y->size[1];
  for (i28 = 0; i28 < ixstart; i28++) {
    t->data[i28] = b_y->data[b_y->size[0] * i28];
  }

  /*  Times */
  if (primetable->size[1] == 0) {
    primes(tmp_data, tmp_size);
    i28 = primetable->size[0] * primetable->size[1];
    primetable->size[0] = 1;
    primetable->size[1] = 1 + tmp_size[1];
    emxEnsureCapacity((emxArray__common *)primetable, i28, (int)sizeof(double));
    primetable->data[0] = 1.0;
    ixstart = tmp_size[1];
    for (i28 = 0; i28 < ixstart; i28++) {
      primetable->data[primetable->size[0] * (i28 + 1)] = tmp_data[tmp_size[0] *
        i28];
    }
  }

  /*  Define pitch candidates */
  y = scalar_log2(plim[0]);
  o = scalar_log2(plim[1]);
  if (rtIsNaN(y) || rtIsNaN(dlog2p) || rtIsNaN(o)) {
    n = 0;
    y = rtNaN;
    apnd = o;
  } else if ((dlog2p == 0.0) || ((y < o) && (dlog2p < 0.0)) || ((o < y) &&
              (dlog2p > 0.0))) {
    n = -1;
    apnd = o;
  } else if (rtIsInf(y) || rtIsInf(o)) {
    n = 0;
    y = rtNaN;
    apnd = o;
  } else if (rtIsInf(dlog2p)) {
    n = 0;
    apnd = o;
  } else {
    ndbl = floor((o - y) / dlog2p + 0.5);
    apnd = y + ndbl * dlog2p;
    if (dlog2p > 0.0) {
      cdiff = apnd - o;
    } else {
      cdiff = o - apnd;
    }

    kd = fabs(y);
    absb = fabs(o);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
      ndbl++;
      apnd = o;
    } else if (cdiff > 0.0) {
      apnd = y + (ndbl - 1.0) * dlog2p;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  i28 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i28, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = y;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (k = 1; k < nm1d2; k++) {
        kd = (double)k * dlog2p;
        b_y->data[k] = y + kd;
        b_y->data[n - k] = apnd - kd;
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (y + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * dlog2p;
        b_y->data[nm1d2] = y + kd;
        b_y->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  b_emxInit_real_T(&log2pc, 1);
  i28 = log2pc->size[0];
  log2pc->size[0] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)log2pc, i28, (int)sizeof(double));
  ixstart = b_y->size[1];
  for (i28 = 0; i28 < ixstart; i28++) {
    log2pc->data[i28] = b_y->data[b_y->size[0] * i28];
  }

  b_emxInit_real_T(&pc, 1);
  emxInit_real_T(&S, 2);
  e_power(log2pc, pc);
  nm1d2 = pc->size[0];
  i28 = S->size[0] * S->size[1];
  S->size[0] = nm1d2;
  emxEnsureCapacity((emxArray__common *)S, i28, (int)sizeof(double));
  nm1d2 = t->size[0];
  i28 = S->size[0] * S->size[1];
  S->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)S, i28, (int)sizeof(double));
  ixstart = pc->size[0] * t->size[0];
  for (i28 = 0; i28 < ixstart; i28++) {
    S->data[i28] = 0.0;
  }

  /*  Pitch strength matrix */
  /*  Determine P2-WSs */
  e_rdivide(8.0 * fs, plim, dv4);
  b_log2(dv4, logWs);
  b_round(logWs);
  if (rtIsNaN(logWs[0]) || rtIsNaN(logWs[1])) {
    n = 0;
    o = rtNaN;
    apnd = logWs[1];
  } else if (logWs[0] < logWs[1]) {
    n = -1;
    o = logWs[0];
    apnd = logWs[1];
  } else if (rtIsInf(logWs[0]) || rtIsInf(logWs[1])) {
    n = 0;
    o = rtNaN;
    apnd = logWs[1];
  } else {
    o = logWs[0];
    ndbl = floor(-(logWs[1] - logWs[0]) + 0.5);
    apnd = logWs[0] + -ndbl;
    cdiff = logWs[1] - apnd;
    kd = fabs(logWs[0]);
    absb = fabs(logWs[1]);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
      ndbl++;
      apnd = logWs[1];
    } else if (cdiff > 0.0) {
      apnd = logWs[0] + -(ndbl - 1.0);
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  i28 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i28, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = o;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (k = 1; k < nm1d2; k++) {
        b_y->data[k] = o + -(double)k;
        b_y->data[n - k] = apnd - (-(double)k);
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (o + apnd) / 2.0;
      } else {
        b_y->data[nm1d2] = o + -(double)nm1d2;
        b_y->data[nm1d2 + 1] = apnd - (-(double)nm1d2);
      }
    }
  }

  emxInit_real_T(&ws, 2);
  emxInit_real_T(&pO, 2);
  f_power(b_y, ws);

  /*  P2-WSs */
  f_rdivide(8.0 * fs, ws, pO);

  /*  Optimal pitches for P2-WSs */
  /*  Determine window sizes used by each pitch candidate */
  y = scalar_log2(b_rdivide(8.0 * fs, ws->data[0]));
  i28 = log2pc->size[0];
  emxEnsureCapacity((emxArray__common *)log2pc, i28, (int)sizeof(double));
  ixstart = log2pc->size[0];
  for (i28 = 0; i28 < ixstart; i28++) {
    log2pc->data[i28] = (1.0 + log2pc->data[i28]) - y;
  }

  /*  Create ERB-scale uniformly-spaced frequencies (in Hertz) */
  ixstart = 1;
  n = pc->size[0];
  kd = pc->data[0];
  if (pc->size[0] > 1) {
    if (rtIsNaN(pc->data[0])) {
      nm1d2 = 2;
      exitg9 = false;
      while ((!exitg9) && (nm1d2 <= n)) {
        ixstart = nm1d2;
        if (!rtIsNaN(pc->data[nm1d2 - 1])) {
          kd = pc->data[nm1d2 - 1];
          exitg9 = true;
        } else {
          nm1d2++;
        }
      }
    }

    if (ixstart < pc->size[0]) {
      while (ixstart + 1 <= n) {
        if (pc->data[ixstart] < kd) {
          kd = pc->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  y = 6.44 * (scalar_log2(229.0 + kd / 4.0) - 7.84);
  o = 6.44 * (scalar_log2(229.0 + fs / 2.0) - 7.84);
  if (rtIsNaN(y) || rtIsNaN(dERBs) || rtIsNaN(o)) {
    n = 0;
    y = rtNaN;
    apnd = o;
  } else if ((dERBs == 0.0) || ((y < o) && (dERBs < 0.0)) || ((o < y) && (dERBs >
    0.0))) {
    n = -1;
    apnd = o;
  } else if (rtIsInf(y) || rtIsInf(o)) {
    n = 0;
    y = rtNaN;
    apnd = o;
  } else if (rtIsInf(dERBs)) {
    n = 0;
    apnd = o;
  } else {
    ndbl = floor((o - y) / dERBs + 0.5);
    apnd = y + ndbl * dERBs;
    if (dERBs > 0.0) {
      cdiff = apnd - o;
    } else {
      cdiff = o - apnd;
    }

    kd = fabs(y);
    absb = fabs(o);
    if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
      ndbl++;
      apnd = o;
    } else if (cdiff > 0.0) {
      apnd = y + (ndbl - 1.0) * dERBs;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl - 1;
    } else {
      n = -1;
    }
  }

  i28 = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  b_y->size[1] = n + 1;
  emxEnsureCapacity((emxArray__common *)b_y, i28, (int)sizeof(double));
  if (n + 1 > 0) {
    b_y->data[0] = y;
    if (n + 1 > 1) {
      b_y->data[n] = apnd;
      nm1d2 = (n + (n < 0)) >> 1;
      for (k = 1; k < nm1d2; k++) {
        kd = (double)k * dERBs;
        b_y->data[k] = y + kd;
        b_y->data[n - k] = apnd - kd;
      }

      if (nm1d2 << 1 == n) {
        b_y->data[nm1d2] = (y + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * dERBs;
        b_y->data[nm1d2] = y + kd;
        b_y->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  b_emxInit_real_T(&c_y, 1);
  i28 = c_y->size[0];
  c_y->size[0] = b_y->size[1];
  emxEnsureCapacity((emxArray__common *)c_y, i28, (int)sizeof(double));
  ixstart = b_y->size[1];
  for (i28 = 0; i28 < ixstart; i28++) {
    c_y->data[i28] = b_y->data[b_y->size[0] * i28];
  }

  b_emxInit_real_T(&w, 1);
  b_emxInit_real_T(&b_w, 1);
  rdivide(c_y, 6.44, w);
  i28 = b_w->size[0];
  b_w->size[0] = w->size[0];
  emxEnsureCapacity((emxArray__common *)b_w, i28, (int)sizeof(double));
  ixstart = w->size[0];
  emxFree_real_T(&c_y);
  for (i28 = 0; i28 < ixstart; i28++) {
    b_w->data[i28] = w->data[i28] + 7.84;
  }

  b_emxInit_real_T(&fERBs, 1);
  e_power(b_w, fERBs);
  i28 = fERBs->size[0];
  emxEnsureCapacity((emxArray__common *)fERBs, i28, (int)sizeof(double));
  ixstart = fERBs->size[0];
  emxFree_real_T(&b_w);
  for (i28 = 0; i28 < ixstart; i28++) {
    fERBs->data[i28] -= 229.0;
  }

  i = 0;
  emxInit_real_T(&X, 2);
  b_emxInit_real_T(&f, 1);
  b_emxInit_real_T(&ti, 1);
  emxInit_real_T(&L, 2);
  emxInit_real_T(&Si, 2);
  b_emxInit_real_T(&tc, 1);
  emxInit_real_T(&j, 2);
  emxInit_real_T(&b_k, 2);
  b_emxInit_real_T(&xzp, 1);
  emxInit_boolean_T(&b_x, 1);
  b_emxInit_boolean_T(&c_x, 2);
  emxInit_int32_T(&ii, 1);
  emxInit_real_T(&b_log2pc, 2);
  emxInit_real_T(&b_Si, 2);
  emxInit_real_T(&b_pc, 2);
  b_emxInit_real_T(&c_log2pc, 1);
  b_emxInit_real_T(&b_fERBs, 1);
  emxInit_real_T(&b_S, 2);
  while (i <= ws->size[1] - 1) {
    kd = fmax(1.0, rt_roundd_snf(8.0 * (1.0 - woverlap) * fs / pO->data[i]));

    /*  Hop size */
    /*  Zero pad signal */
    y = ws->data[i] / 2.0;
    o = ws->data[i] / 2.0;
    i28 = tc->size[0];
    tc->size[0] = ((int)y + x->size[0]) + (int)(kd + o);
    emxEnsureCapacity((emxArray__common *)tc, i28, (int)sizeof(double));
    ixstart = (int)y;
    for (i28 = 0; i28 < ixstart; i28++) {
      tc->data[i28] = 0.0;
    }

    ixstart = x->size[0];
    for (i28 = 0; i28 < ixstart; i28++) {
      tc->data[i28 + (int)y] = x->data[i28];
    }

    ixstart = (int)(kd + o);
    for (i28 = 0; i28 < ixstart; i28++) {
      tc->data[(i28 + (int)y) + x->size[0]] = 0.0;
    }

    /*  Compute spectrum */
    i28 = w->size[0];
    w->size[0] = (int)ws->data[i];
    emxEnsureCapacity((emxArray__common *)w, i28, (int)sizeof(double));
    ixstart = (int)ws->data[i];
    for (i28 = 0; i28 < ixstart; i28++) {
      w->data[i28] = 0.0;
    }

    hanning(&w->data[0], ws->data[i]);
    o = fmax(0.0, rt_roundd_snf(ws->data[i] - kd));

    /*  Window overlap */
    kd = ws->data[i] - o;
    y = ws->data[i] / 2.0;
    kd = (((double)tc->size[0] + kd) - 1.0) / kd;
    i28 = X->size[0] * X->size[1];
    X->size[0] = (int)(y + 1.0);
    X->size[1] = (int)kd;
    emxEnsureCapacity((emxArray__common *)X, i28, (int)sizeof(double));
    ixstart = (int)(y + 1.0) * (int)kd;
    for (i28 = 0; i28 < ixstart; i28++) {
      X->data[i28] = 0.0;
    }

    i28 = f->size[0];
    f->size[0] = (int)(y + 1.0);
    emxEnsureCapacity((emxArray__common *)f, i28, (int)sizeof(double));
    ixstart = (int)(y + 1.0);
    for (i28 = 0; i28 < ixstart; i28++) {
      f->data[i28] = 0.0;
    }

    i28 = ti->size[0];
    ti->size[0] = (int)kd;
    emxEnsureCapacity((emxArray__common *)ti, i28, (int)sizeof(double));
    ixstart = (int)kd;
    for (i28 = 0; i28 < ixstart; i28++) {
      ti->data[i28] = 0.0;
    }

    /*  use argument order of newer spectrogram function, for which we */
    /*  have documentation */
    i28 = xzp->size[0];
    xzp->size[0] = tc->size[0];
    emxEnsureCapacity((emxArray__common *)xzp, i28, (int)sizeof(double));
    ixstart = tc->size[0];
    for (i28 = 0; i28 < ixstart; i28++) {
      xzp->data[i28] = tc->data[i28];
    }

    spectrogram(&X->data[0], &f->data[0], &ti->data[0], &xzp->data[0], (double)
                tc->size[0], &w->data[0], o, ws->data[i], fs);

    /*  Select candidates that use this window size */
    if (ws->size[1] == 1) {
      i28 = j->size[0] * j->size[1];
      j->size[0] = 1;
      j->size[1] = pc->size[0];
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      ixstart = pc->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        j->data[j->size[0] * i28] = pc->data[i28];
      }

      i28 = b_k->size[0] * b_k->size[1];
      b_k->size[0] = 0;
      b_k->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)b_k, i28, (int)sizeof(double));
    } else if (1 + i == ws->size[1]) {
      i28 = b_x->size[0];
      b_x->size[0] = log2pc->size[0];
      emxEnsureCapacity((emxArray__common *)b_x, i28, (int)sizeof(boolean_T));
      ixstart = log2pc->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_x->data[i28] = (log2pc->data[i28] - (1.0 + (double)i) > -1.0);
      }

      nm1d2 = b_x->size[0];
      idx = 0;
      i28 = ii->size[0];
      ii->size[0] = b_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      ixstart = 1;
      exitg8 = false;
      while ((!exitg8) && (ixstart <= nm1d2)) {
        guard5 = false;
        if (b_x->data[ixstart - 1]) {
          idx++;
          ii->data[idx - 1] = ixstart;
          if (idx >= nm1d2) {
            exitg8 = true;
          } else {
            guard5 = true;
          }
        } else {
          guard5 = true;
        }

        if (guard5) {
          ixstart++;
        }
      }

      if (b_x->size[0] == 1) {
        if (idx == 0) {
          i28 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
        }
      } else {
        i28 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      }

      i28 = w->size[0];
      w->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)w, i28, (int)sizeof(double));
      ixstart = ii->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        w->data[i28] = ii->data[i28];
      }

      idx = w->size[0];
      i28 = j->size[0] * j->size[1];
      j->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      i28 = j->size[0] * j->size[1];
      j->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        j->data[i28] = w->data[i28];
      }

      idx = w->size[0];
      i28 = c_x->size[0] * c_x->size[1];
      c_x->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)c_x, i28, (int)sizeof(boolean_T));
      i28 = c_x->size[0] * c_x->size[1];
      c_x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)c_x, i28, (int)sizeof(boolean_T));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        c_x->data[i28] = (log2pc->data[(int)w->data[i28] - 1] - (1.0 + (double)i)
                          < 0.0);
      }

      nm1d2 = c_x->size[0];
      idx = 0;
      i28 = ii->size[0];
      ii->size[0] = c_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      ixstart = 1;
      exitg7 = false;
      while ((!exitg7) && (ixstart <= nm1d2)) {
        guard4 = false;
        if (c_x->data[ixstart - 1]) {
          idx++;
          ii->data[idx - 1] = ixstart;
          if (idx >= nm1d2) {
            exitg7 = true;
          } else {
            guard4 = true;
          }
        } else {
          guard4 = true;
        }

        if (guard4) {
          ixstart++;
        }
      }

      if (c_x->size[0] == 1) {
        if (idx == 0) {
          i28 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
        }
      } else {
        i28 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      }

      i28 = w->size[0];
      w->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)w, i28, (int)sizeof(double));
      ixstart = ii->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        w->data[i28] = ii->data[i28];
      }

      idx = w->size[0];
      i28 = b_k->size[0] * b_k->size[1];
      b_k->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)b_k, i28, (int)sizeof(double));
      i28 = b_k->size[0] * b_k->size[1];
      b_k->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_k, i28, (int)sizeof(double));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_k->data[i28] = w->data[i28];
      }
    } else if (1 + i == 1) {
      i28 = b_x->size[0];
      b_x->size[0] = log2pc->size[0];
      emxEnsureCapacity((emxArray__common *)b_x, i28, (int)sizeof(boolean_T));
      ixstart = log2pc->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_x->data[i28] = (log2pc->data[i28] - 1.0 < 1.0);
      }

      nm1d2 = b_x->size[0];
      idx = 0;
      i28 = ii->size[0];
      ii->size[0] = b_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      ixstart = 1;
      exitg6 = false;
      while ((!exitg6) && (ixstart <= nm1d2)) {
        guard3 = false;
        if (b_x->data[ixstart - 1]) {
          idx++;
          ii->data[idx - 1] = ixstart;
          if (idx >= nm1d2) {
            exitg6 = true;
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }

        if (guard3) {
          ixstart++;
        }
      }

      if (b_x->size[0] == 1) {
        if (idx == 0) {
          i28 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
        }
      } else {
        i28 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      }

      i28 = w->size[0];
      w->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)w, i28, (int)sizeof(double));
      ixstart = ii->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        w->data[i28] = ii->data[i28];
      }

      idx = w->size[0];
      i28 = j->size[0] * j->size[1];
      j->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      i28 = j->size[0] * j->size[1];
      j->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        j->data[i28] = w->data[i28];
      }

      idx = w->size[0];
      i28 = c_x->size[0] * c_x->size[1];
      c_x->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)c_x, i28, (int)sizeof(boolean_T));
      i28 = c_x->size[0] * c_x->size[1];
      c_x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)c_x, i28, (int)sizeof(boolean_T));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        c_x->data[i28] = (log2pc->data[(int)w->data[i28] - 1] - 1.0 > 0.0);
      }

      nm1d2 = c_x->size[0];
      idx = 0;
      i28 = ii->size[0];
      ii->size[0] = c_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      ixstart = 1;
      exitg5 = false;
      while ((!exitg5) && (ixstart <= nm1d2)) {
        guard2 = false;
        if (c_x->data[ixstart - 1]) {
          idx++;
          ii->data[idx - 1] = ixstart;
          if (idx >= nm1d2) {
            exitg5 = true;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }

        if (guard2) {
          ixstart++;
        }
      }

      if (c_x->size[0] == 1) {
        if (idx == 0) {
          i28 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
        }
      } else {
        i28 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      }

      i28 = w->size[0];
      w->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)w, i28, (int)sizeof(double));
      ixstart = ii->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        w->data[i28] = ii->data[i28];
      }

      idx = w->size[0];
      i28 = b_k->size[0] * b_k->size[1];
      b_k->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)b_k, i28, (int)sizeof(double));
      i28 = b_k->size[0] * b_k->size[1];
      b_k->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_k, i28, (int)sizeof(double));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_k->data[i28] = w->data[i28];
      }
    } else {
      i28 = c_log2pc->size[0];
      c_log2pc->size[0] = log2pc->size[0];
      emxEnsureCapacity((emxArray__common *)c_log2pc, i28, (int)sizeof(double));
      ixstart = log2pc->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        c_log2pc->data[i28] = log2pc->data[i28] - (1.0 + (double)i);
      }

      b_abs(c_log2pc, w);
      i28 = b_x->size[0];
      b_x->size[0] = w->size[0];
      emxEnsureCapacity((emxArray__common *)b_x, i28, (int)sizeof(boolean_T));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_x->data[i28] = (w->data[i28] < 1.0);
      }

      nm1d2 = b_x->size[0];
      idx = 0;
      i28 = ii->size[0];
      ii->size[0] = b_x->size[0];
      emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      ixstart = 1;
      exitg4 = false;
      while ((!exitg4) && (ixstart <= nm1d2)) {
        guard1 = false;
        if (b_x->data[ixstart - 1]) {
          idx++;
          ii->data[idx - 1] = ixstart;
          if (idx >= nm1d2) {
            exitg4 = true;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          ixstart++;
        }
      }

      if (b_x->size[0] == 1) {
        if (idx == 0) {
          i28 = ii->size[0];
          ii->size[0] = 0;
          emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
        }
      } else {
        i28 = ii->size[0];
        if (1 > idx) {
          ii->size[0] = 0;
        } else {
          ii->size[0] = idx;
        }

        emxEnsureCapacity((emxArray__common *)ii, i28, (int)sizeof(int));
      }

      i28 = w->size[0];
      w->size[0] = ii->size[0];
      emxEnsureCapacity((emxArray__common *)w, i28, (int)sizeof(double));
      ixstart = ii->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        w->data[i28] = ii->data[i28];
      }

      idx = w->size[0];
      i28 = j->size[0] * j->size[1];
      j->size[0] = idx;
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      i28 = j->size[0] * j->size[1];
      j->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)j, i28, (int)sizeof(double));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        j->data[i28] = w->data[i28];
      }

      idx = w->size[0];
      if (0 == idx) {
        idx = 0;
      } else {
        idx = w->size[0];
        if (idx > 1) {
          idx = w->size[0];
        } else {
          idx = 1;
        }
      }

      if (idx < 1) {
        n = -1;
        b_apnd = 0;
      } else {
        nm1d2 = (int)floor(((double)idx - 1.0) + 0.5);
        b_apnd = nm1d2 + 1;
        ixstart = (nm1d2 - idx) + 1;
        if (fabs(ixstart) < 4.4408920985006262E-16 * (double)idx) {
          nm1d2++;
          b_apnd = idx;
        } else if (ixstart > 0) {
          b_apnd = nm1d2;
        } else {
          nm1d2++;
        }

        n = nm1d2 - 1;
      }

      i28 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = n + 1;
      emxEnsureCapacity((emxArray__common *)b_y, i28, (int)sizeof(double));
      if (n + 1 > 0) {
        b_y->data[0] = 1.0;
        if (n + 1 > 1) {
          b_y->data[n] = b_apnd;
          nm1d2 = (n + (n < 0)) >> 1;
          for (k = 1; k < nm1d2; k++) {
            b_y->data[k] = 1.0 + (double)k;
            b_y->data[n - k] = b_apnd - k;
          }

          if (nm1d2 << 1 == n) {
            b_y->data[nm1d2] = (1.0 + (double)b_apnd) / 2.0;
          } else {
            b_y->data[nm1d2] = 1.0 + (double)nm1d2;
            b_y->data[nm1d2 + 1] = b_apnd - nm1d2;
          }
        }
      }

      i28 = b_k->size[0] * b_k->size[1];
      b_k->size[0] = 1;
      b_k->size[1] = b_y->size[1];
      emxEnsureCapacity((emxArray__common *)b_k, i28, (int)sizeof(double));
      ixstart = b_y->size[0] * b_y->size[1];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_k->data[i28] = b_y->data[i28];
      }
    }

    /*  Compute loudness at ERBs uniformly-spaced frequencies */
    y = pc->data[(int)j->data[0] - 1] / 4.0;
    i28 = b_x->size[0];
    b_x->size[0] = fERBs->size[0];
    emxEnsureCapacity((emxArray__common *)b_x, i28, (int)sizeof(boolean_T));
    ixstart = fERBs->size[0];
    for (i28 = 0; i28 < ixstart; i28++) {
      b_x->data[i28] = (fERBs->data[i28] > y);
    }

    if (1 <= b_x->size[0]) {
      k = 1;
    } else {
      k = 0;
    }

    idx = 0;
    ixstart = 1;
    exitg3 = false;
    while ((!exitg3) && (ixstart <= b_x->size[0])) {
      if (b_x->data[ixstart - 1]) {
        idx = 1;
        ii_data[0] = ixstart;
        exitg3 = true;
      } else {
        ixstart++;
      }
    }

    if (k == 1) {
      if (idx == 0) {
        k = 0;
      }
    } else if (1 > idx) {
      k = 0;
    } else {
      k = 1;
    }

    i28 = 0;
    while (i28 <= k - 1) {
      found_data[0] = ii_data[0];
      i28 = 1;
    }

    if (found_data[0] > fERBs->size[0]) {
      i28 = 0;
      idx = 0;
    } else {
      i28 = found_data[0] - 1;
      idx = fERBs->size[0];
    }

    nm1d2 = b_fERBs->size[0];
    b_fERBs->size[0] = idx - i28;
    emxEnsureCapacity((emxArray__common *)b_fERBs, nm1d2, (int)sizeof(double));
    ixstart = idx - i28;
    for (idx = 0; idx < ixstart; idx++) {
      b_fERBs->data[idx] = fERBs->data[i28 + idx];
    }

    i28 = fERBs->size[0];
    fERBs->size[0] = b_fERBs->size[0];
    emxEnsureCapacity((emxArray__common *)fERBs, i28, (int)sizeof(double));
    ixstart = b_fERBs->size[0];
    for (i28 = 0; i28 < ixstart; i28++) {
      fERBs->data[i28] = b_fERBs->data[i28];
    }

    c_abs(X, L);
    interp1(f, L, fERBs, X);
    for (i28 = 0; i28 < 2; i28++) {
      tmp_size[i28] = X->size[i28];
    }

    i28 = L->size[0] * L->size[1];
    L->size[0] = tmp_size[0];
    L->size[1] = tmp_size[1];
    emxEnsureCapacity((emxArray__common *)L, i28, (int)sizeof(double));
    i28 = tmp_size[0] * tmp_size[1];
    for (k = 0; k + 1 <= i28; k++) {
      L->data[k] = fmax(0.0, X->data[k]);
    }

    c_sqrt(L);

    /*  Compute pitch strength */
    i28 = b_pc->size[0] * b_pc->size[1];
    b_pc->size[0] = j->size[0];
    b_pc->size[1] = j->size[1];
    emxEnsureCapacity((emxArray__common *)b_pc, i28, (int)sizeof(double));
    ixstart = j->size[0] * j->size[1];
    for (i28 = 0; i28 < ixstart; i28++) {
      b_pc->data[i28] = pc->data[(int)j->data[i28] - 1];
    }

    pitchStrengthAllCandidates(fERBs, L, b_pc, Si);

    /*  Interpolate pitch strength at desired times */
    if (Si->size[1] > 1) {
      i28 = b_Si->size[0] * b_Si->size[1];
      b_Si->size[0] = Si->size[1];
      b_Si->size[1] = Si->size[0];
      emxEnsureCapacity((emxArray__common *)b_Si, i28, (int)sizeof(double));
      ixstart = Si->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        nm1d2 = Si->size[1];
        for (idx = 0; idx < nm1d2; idx++) {
          b_Si->data[idx + b_Si->size[0] * i28] = Si->data[i28 + Si->size[0] *
            idx];
        }
      }

      b_interp1(ti, b_Si, t, L);
      i28 = Si->size[0] * Si->size[1];
      Si->size[0] = L->size[1];
      Si->size[1] = L->size[0];
      emxEnsureCapacity((emxArray__common *)Si, i28, (int)sizeof(double));
      ixstart = L->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        nm1d2 = L->size[1];
        for (idx = 0; idx < nm1d2; idx++) {
          Si->data[idx + Si->size[0] * i28] = L->data[i28 + L->size[0] * idx];
        }
      }
    } else {
      if ((0 == Si->size[0]) || (0 == Si->size[1])) {
        nm1d2 = 0;
      } else {
        nm1d2 = Si->size[0];
        if (nm1d2 >= 1) {
        } else {
          nm1d2 = 1;
        }
      }

      d_repmat(nm1d2, t->size[0], Si);
    }

    /*  Add pitch strength to combination */
    for (i28 = 0; i28 < 2; i28++) {
      tmp_size[i28] = j->size[i28];
    }

    i28 = X->size[0] * X->size[1];
    X->size[0] = tmp_size[0];
    emxEnsureCapacity((emxArray__common *)X, i28, (int)sizeof(double));
    i28 = X->size[0] * X->size[1];
    X->size[1] = tmp_size[1];
    emxEnsureCapacity((emxArray__common *)X, i28, (int)sizeof(double));
    ixstart = tmp_size[0] * tmp_size[1];
    for (i28 = 0; i28 < ixstart; i28++) {
      X->data[i28] = 1.0;
    }

    i28 = b_log2pc->size[0] * b_log2pc->size[1];
    b_log2pc->size[0] = b_k->size[0];
    b_log2pc->size[1] = b_k->size[1];
    emxEnsureCapacity((emxArray__common *)b_log2pc, i28, (int)sizeof(double));
    ixstart = b_k->size[0] * b_k->size[1];
    for (i28 = 0; i28 < ixstart; i28++) {
      b_log2pc->data[i28] = log2pc->data[(int)j->data[(int)b_k->data[i28] - 1] -
        1] - (1.0 + (double)i);
    }

    c_abs(b_log2pc, L);
    ixstart = L->size[0] * L->size[1];
    for (i28 = 0; i28 < ixstart; i28++) {
      X->data[(int)b_k->data[i28] - 1] = 1.0 - L->data[i28];
    }

    e_repmat(X, Si->size[1], L);
    nm1d2 = j->size[0] * j->size[1];
    ixstart = S->size[1];
    i28 = b_S->size[0] * b_S->size[1];
    b_S->size[0] = nm1d2;
    b_S->size[1] = ixstart;
    emxEnsureCapacity((emxArray__common *)b_S, i28, (int)sizeof(double));
    for (i28 = 0; i28 < ixstart; i28++) {
      for (idx = 0; idx < nm1d2; idx++) {
        b_S->data[idx + b_S->size[0] * i28] = S->data[((int)j->data[idx] +
          S->size[0] * i28) - 1] + L->data[idx + L->size[0] * i28] * Si->
          data[idx + Si->size[0] * i28];
      }
    }

    ixstart = b_S->size[1];
    for (i28 = 0; i28 < ixstart; i28++) {
      nm1d2 = b_S->size[0];
      for (idx = 0; idx < nm1d2; idx++) {
        S->data[((int)j->data[idx] + S->size[0] * i28) - 1] = b_S->data[idx +
          b_S->size[0] * i28];
      }
    }

    i++;
  }

  emxFree_real_T(&b_S);
  emxFree_real_T(&b_fERBs);
  emxFree_real_T(&c_log2pc);
  emxFree_real_T(&b_pc);
  emxFree_real_T(&b_Si);
  emxFree_real_T(&b_log2pc);
  emxFree_int32_T(&ii);
  emxFree_boolean_T(&c_x);
  emxFree_boolean_T(&b_x);
  emxFree_real_T(&xzp);
  emxFree_real_T(&b_k);
  emxFree_real_T(&j);
  emxFree_real_T(&Si);
  emxFree_real_T(&L);
  emxFree_real_T(&ti);
  emxFree_real_T(&f);
  emxFree_real_T(&X);
  emxFree_real_T(&fERBs);
  emxFree_real_T(&log2pc);

  /*  Fine tune pitch using parabolic interpolation */
  f_repmat(S->size[1], p);
  f_repmat(S->size[1], s);
  b_apnd = 0;
  emxInit_real_T(&I, 2);
  emxInit_real_T(&b_pO, 2);
  b_emxInit_real_T(&c_w, 1);
  b_emxInit_real_T(&c_S, 1);
  b_emxInit_real_T(&c_pc, 1);
  while (b_apnd <= S->size[1] - 1) {
    ixstart = 1;
    n = S->size[0];
    kd = S->data[S->size[0] * b_apnd];
    idx = 1;
    i28 = S->size[0];
    if (i28 > 1) {
      if (rtIsNaN(kd)) {
        nm1d2 = 2;
        exitg2 = false;
        while ((!exitg2) && (nm1d2 <= n)) {
          ixstart = nm1d2;
          if (!rtIsNaN(S->data[(nm1d2 + S->size[0] * b_apnd) - 1])) {
            kd = S->data[(nm1d2 + S->size[0] * b_apnd) - 1];
            idx = nm1d2;
            exitg2 = true;
          } else {
            nm1d2++;
          }
        }
      }

      i28 = S->size[0];
      if (ixstart < i28) {
        while (ixstart + 1 <= n) {
          if (S->data[ixstart + S->size[0] * b_apnd] > kd) {
            kd = S->data[ixstart + S->size[0] * b_apnd];
            idx = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    s->data[b_apnd] = kd;
    if (s->data[b_apnd] < sTHR) {
    } else if ((idx == 1) || (idx == pc->size[0])) {
      p->data[b_apnd] = pc->data[idx - 1];
    } else {
      if ((double)idx + 1.0 < (double)idx - 1.0) {
        n = -1;
        o = (double)idx - 1.0;
        apnd = (double)idx + 1.0;
      } else {
        o = (double)idx - 1.0;
        ndbl = floor((((double)idx + 1.0) - ((double)idx - 1.0)) + 0.5);
        apnd = ((double)idx - 1.0) + ndbl;
        cdiff = apnd - ((double)idx + 1.0);
        nm1d2 = (int)fabs((double)idx - 1.0);
        b_absb = idx + 1U;
        u0 = (unsigned int)nm1d2;
        if (u0 >= b_absb) {
          b_absb = u0;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * (double)b_absb) {
          ndbl++;
          apnd = (double)idx + 1.0;
        } else if (cdiff > 0.0) {
          apnd = ((double)idx - 1.0) + (ndbl - 1.0);
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl - 1;
        } else {
          n = -1;
        }
      }

      i28 = I->size[0] * I->size[1];
      I->size[0] = 1;
      I->size[1] = n + 1;
      emxEnsureCapacity((emxArray__common *)I, i28, (int)sizeof(double));
      if (n + 1 > 0) {
        I->data[0] = o;
        if (n + 1 > 1) {
          I->data[n] = apnd;
          nm1d2 = (n + (n < 0)) >> 1;
          for (k = 1; k < nm1d2; k++) {
            I->data[k] = o + (double)k;
            I->data[n - k] = apnd - (double)k;
          }

          if (nm1d2 << 1 == n) {
            I->data[nm1d2] = (o + apnd) / 2.0;
          } else {
            I->data[nm1d2] = o + (double)nm1d2;
            I->data[nm1d2 + 1] = apnd - (double)nm1d2;
          }
        }
      }

      i28 = c_pc->size[0];
      c_pc->size[0] = I->size[1];
      emxEnsureCapacity((emxArray__common *)c_pc, i28, (int)sizeof(double));
      ixstart = I->size[1];
      for (i28 = 0; i28 < ixstart; i28++) {
        c_pc->data[i28] = pc->data[(int)I->data[I->size[0] * i28] - 1];
      }

      g_rdivide(c_pc, tc);
      rdivide(tc, tc->data[1], w);
      i28 = c_w->size[0];
      c_w->size[0] = w->size[0];
      emxEnsureCapacity((emxArray__common *)c_w, i28, (int)sizeof(double));
      ixstart = w->size[0];
      for (i28 = 0; i28 < ixstart; i28++) {
        c_w->data[i28] = (w->data[i28] - 1.0) * 2.0 * 3.1415926535897931;
      }

      i28 = c_S->size[0];
      c_S->size[0] = I->size[1];
      emxEnsureCapacity((emxArray__common *)c_S, i28, (int)sizeof(double));
      ixstart = I->size[1];
      for (i28 = 0; i28 < ixstart; i28++) {
        c_S->data[i28] = S->data[((int)I->data[I->size[0] * i28] + S->size[0] *
          b_apnd) - 1];
      }

      b_polyfit(c_w, c_S, c);
      y = scalar_log2(pc->data[(int)I->data[0] - 1]);
      o = scalar_log2(pc->data[(int)I->data[2] - 1]);
      if (rtIsNaN(y) || rtIsNaN(o)) {
        n = 0;
        y = rtNaN;
        apnd = o;
      } else if (o < y) {
        n = -1;
        apnd = o;
      } else if (rtIsInf(y) || rtIsInf(o)) {
        n = 0;
        y = rtNaN;
        apnd = o;
      } else {
        ndbl = floor((o - y) / 0.00083333333333333328 + 0.5);
        apnd = y + ndbl * 0.00083333333333333328;
        cdiff = apnd - o;
        kd = fabs(y);
        absb = fabs(o);
        if (fabs(cdiff) < 4.4408920985006262E-16 * fmax(kd, absb)) {
          ndbl++;
          apnd = o;
        } else if (cdiff > 0.0) {
          apnd = y + (ndbl - 1.0) * 0.00083333333333333328;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl - 1;
        } else {
          n = -1;
        }
      }

      i28 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = n + 1;
      emxEnsureCapacity((emxArray__common *)b_y, i28, (int)sizeof(double));
      if (n + 1 > 0) {
        b_y->data[0] = y;
        if (n + 1 > 1) {
          b_y->data[n] = apnd;
          nm1d2 = (n + (n < 0)) >> 1;
          for (k = 1; k < nm1d2; k++) {
            kd = (double)k * 0.00083333333333333328;
            b_y->data[k] = y + kd;
            b_y->data[n - k] = apnd - kd;
          }

          if (nm1d2 << 1 == n) {
            b_y->data[nm1d2] = (y + apnd) / 2.0;
          } else {
            kd = (double)nm1d2 * 0.00083333333333333328;
            b_y->data[nm1d2] = y + kd;
            b_y->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }

      f_power(b_y, pO);
      f_rdivide(1.0, pO, ws);
      c_rdivide(ws, tc->data[1], pO);
      i28 = b_pO->size[0] * b_pO->size[1];
      b_pO->size[0] = 1;
      b_pO->size[1] = pO->size[1];
      emxEnsureCapacity((emxArray__common *)b_pO, i28, (int)sizeof(double));
      ixstart = pO->size[0] * pO->size[1];
      for (i28 = 0; i28 < ixstart; i28++) {
        b_pO->data[i28] = (pO->data[i28] - 1.0) * 2.0 * 3.1415926535897931;
      }

      polyval(c, b_pO, ws);
      ixstart = 1;
      n = ws->size[1];
      kd = ws->data[0];
      idx = 1;
      if (ws->size[1] > 1) {
        if (rtIsNaN(ws->data[0])) {
          nm1d2 = 2;
          exitg1 = false;
          while ((!exitg1) && (nm1d2 <= n)) {
            ixstart = nm1d2;
            if (!rtIsNaN(ws->data[nm1d2 - 1])) {
              kd = ws->data[nm1d2 - 1];
              idx = nm1d2;
              exitg1 = true;
            } else {
              nm1d2++;
            }
          }
        }

        if (ixstart < ws->size[1]) {
          while (ixstart + 1 <= n) {
            if (ws->data[ixstart] > kd) {
              kd = ws->data[ixstart];
              idx = ixstart + 1;
            }

            ixstart++;
          }
        }
      }

      s->data[b_apnd] = kd;
      p->data[b_apnd] = rt_powd_snf(2.0, scalar_log2(pc->data[(int)I->data[0] -
        1]) + ((double)idx - 1.0) / 12.0 / 100.0);
    }

    b_apnd++;
  }

  emxFree_real_T(&c_pc);
  emxFree_real_T(&c_S);
  emxFree_real_T(&c_w);
  emxFree_real_T(&b_pO);
  emxFree_real_T(&b_y);
  emxFree_real_T(&tc);
  emxFree_real_T(&I);
  emxFree_real_T(&w);
  emxFree_real_T(&pO);
  emxFree_real_T(&ws);
  emxFree_real_T(&S);
  emxFree_real_T(&pc);
}

/* End of code generation (swipep.c) */
