//
//  bridge_ufb_hand_tooled.h
//  PDScores
//
//  Created by Erin Mounts on 3/3/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "PDMatlab.h"

void fastdfa(const PDRealArray *x, double *alpha, PDRealArray **intervals, PDRealArray **flucts);
void features_bpa(PDRealArray *post, PDRealArray **ft);
void features_bga(PDRealArray *gait, PDRealArray **ft);
void features_bta(PDRealArray *tap, PDRealArray **ft);
void lomb(PDRealArray *t, PDRealArray *h, double ofac, double hifac, PDRealArray **f, PDRealArray **P, PDRealArray **prob);
