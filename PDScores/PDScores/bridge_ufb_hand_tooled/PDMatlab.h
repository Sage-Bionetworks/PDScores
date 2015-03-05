//
//  PDMatlab.h
//  PDScores
//
//  Created by Erin Mounts on 3/3/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "PDArray.h"

extern PDRealArray *zeros(size_t rows, size_t columns);
extern PDRealArray *NaN(size_t rows, size_t columns);
extern PDRealArray *sortrows(PDRealArray *table, size_t column);
extern PDRealArray *polyfit(PDRealArray *x, PDRealArray *y, int order);
extern PDRealArray *repmat(PDRealArray *x, size_t rowsreps, size_t colsreps);
extern PDRealArray *quantile(PDRealArray *x, double p);