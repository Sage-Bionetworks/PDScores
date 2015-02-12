//
//  signalprocessing.h
//  PDScores
//
//  Created by Erin Mounts on 2/9/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#ifndef __PDScores__signalprocessing__
#define __PDScores__signalprocessing__

void hanning(double *outBuf, unsigned long windowSize);
void spectrogram(double *outFourierTransform, double *outFrequencies, double *outTimes, double *inSignal, unsigned long signalSize, double *window, unsigned long overlap, unsigned long windowSize, double samplingRate);

#endif /* defined(__PDScores__signalprocessing__) */
