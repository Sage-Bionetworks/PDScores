//
//  bridge_ufb_hand_tooled.h
//  PDScores
//
//  Created by Erin Mounts on 3/3/15.
//  (CC BY-SA 3.0) Sage Bionetworks, 2015.
//

#import "PDMatlab.h"

void fastdfa(const PDRealArray *x, double *alpha, PDRealArray **intervals, PDRealArray **flucts);

void features_bpa(PDRealArray *post, PDRealArray **ft);
void features_bga(PDRealArray *gait, PDRealArray **ft);
void features_bta(PDRealArray *tap, PDRealArray **ft);
void features_bvav2(PDRealArray *audio, double srate, PDRealArray **ft);

double features_ufb(PDRealArray *ftvec, PDRealArray *wvec, PDIntArray *ilog, PDRealArray *ftmin, PDRealArray *ftmax, double fbmin, double fbmax);

void lomb(PDRealArray *t, PDRealArray *h, double ofac, double hifac, PDRealArray **f, PDRealArray **P, PDRealArray **prob);
void mfcc(PDRealArray *samples, PDRealArray **cepstra);
void swipep(PDRealArray *x, double fs, double plim[2], double dt, double dlog2p, double dERBs, double woverlap, double sTHR,
            PDRealArray **p, PDRealArray **t, PDRealArray **s);