//
//  signalprocessing.c
//  PDScores
//
//  Created by Erin Mounts on 2/9/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#include "signalprocessing.h"
#include <Accelerate/Accelerate.h>

void hanning(double *outBuf, unsigned long windowSize)
{
    vDSP_hann_windowD(outBuf, windowSize, vDSP_HANN_DENORM);
}

void spectrogram(double *outFourierTransform, double *outFrequencies, double *outTimes, double *inSignal, unsigned long signalSize, double *window, unsigned long overlap, unsigned long windowSize, double samplingRate)
{
    unsigned long fftSize = windowSize;					// sample size
    unsigned long fftSizeOver2 = fftSize/2;
    unsigned long log2n = log2f(fftSize);			// bins
    
    double *in_real = (double *) malloc(fftSize * sizeof(double));
    DOUBLE_COMPLEX_SPLIT split_data;
    split_data.realp = (double *) malloc(fftSizeOver2 * sizeof(double));
    split_data.imagp = (double *) malloc(fftSizeOver2 * sizeof(double));
    
    FFTSetupD fftSetup = vDSP_create_fftsetupD(log2n, FFT_RADIX2);
    
    unsigned long framestep = windowSize - overlap;
    unsigned long frames = (signalSize + framestep - 1) / framestep;
    
    double frameDuration = (double)framestep / samplingRate;
    double frameTime = frameDuration / 2.0;
    
    double frequencyStep = samplingRate / fftSize;
    double frequencyForBin = 0; // DC
    for (unsigned long bin = 0; bin <= fftSizeOver2; ++bin) {
        outFrequencies[bin] = frequencyForBin;
        frequencyForBin += frequencyStep;
    }
    
    for (unsigned long frame = 0; frame < frames; ++frame) {
        double *frameSignal = inSignal + (frame * framestep);
        
        //multiply by window
        vDSP_vmulD(frameSignal, 1, window, 1, in_real, 1, fftSize);
        
        //convert to split complex format with evens in real and odds in imag
        vDSP_ctozD((DOUBLE_COMPLEX *) in_real, 2, &split_data, 1, fftSizeOver2);
        
        //calc fft
        vDSP_fft_zripD(fftSetup, &split_data, 1, log2n, FFT_FORWARD);
        
        // convert back to interleaved complex format
        vDSP_ztocD(&split_data, 1, (DOUBLE_COMPLEX *) in_real, 2, fftSizeOver2);
        
        // TODO: Divide all coefficients by 2 due to scaling?
        // https://developer.apple.com/library/ios/documentation/Performance/Conceptual/vDSP_Programming_Guide/UsingFourierTransforms/UsingFourierTransforms.html#//apple_ref/doc/uid/TP40005147-CH202-15952
        unsigned long colSize = fftSizeOver2 + 1;
        double *outColStart = outFourierTransform + frame * colSize;
        cblas_dcopy((int)fftSizeOver2, in_real, 2, outColStart, 1);
        outColStart[fftSizeOver2] = in_real[1]; // Nyquist value packed in complex part alongside DC value in real part
        
        outTimes[frame] = frameTime;
        frameTime += frameDuration;
    }
    
}