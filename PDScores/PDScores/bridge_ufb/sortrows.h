/*
 * sortrows.h
 *
 * Code generation for function 'sortrows'
 *
 */

#ifndef __SORTROWS_H__
#define __SORTROWS_H__

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "bridge_ufb_types.h"

/* Function Declarations */
extern void b_sortrows(emxArray_real_T *y);
extern void sortrows(const emxArray_real_T *y, emxArray_real_T *b_y);

#endif

/* End of code generation (sortrows.h) */
