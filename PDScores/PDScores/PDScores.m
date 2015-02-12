//
//  PDScores.m
//  PDScores
//
//  Created by Erin Mounts on 2/9/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "PDScores.h"
#import "bridge_ufb_types.h"
#import "features_bga.h"
#import "features_bpa.h"
#import "features_bta.h"
#import "features_bvav2.h"
#import "features_ufb.h"
#import "bridge_ufb_initialize.h"
@import AudioToolbox;
@import CoreMotion;

const double fbmin = -2.6912;
const double fbmax = 0.1125;

@implementation PDScores

+ (void)initialize
{
    bridge_ufb_initialize();
}

+ (double)scoreFromGaitTest:(NSArray *)gaitData
{
    // format the gait data for the MATLAB code
    emxArray_real_T gait;
    int sizeOfGait[2] = {(int)[gaitData count], 4};
    gait.size = sizeOfGait;
    gait.data = malloc(sizeof(double) * sizeOfGait[0] * sizeOfGait[1]);
    double *pT = gait.data;
    double *pX = pT + sizeOfGait[0];
    double *pY = pX + sizeOfGait[0];
    double *pZ = pY + sizeOfGait[0];
    for (int i = 0; i < sizeOfGait[0]; ++i) {
        NSDictionary *accel = [gaitData objectAtIndex:i];
        *pT++ = [[accel objectForKey:@"timestamp"] doubleValue];
        *pX++ = [[accel objectForKey:@"x"] doubleValue];
        *pY++ = [[accel objectForKey:@"y"] doubleValue];
        *pZ++ = [[accel objectForKey:@"z"] doubleValue];
    }
    
    // get the features array
    double ft[7];
    features_bga(&gait, ft);
    for (int i = 0; i < 7; ++i) {
        if (isnan(ft[i]))
            return NAN;
    }

    static emxArray_real_T ftvec;
    static int sizeOfFtvec[2] = {1, 7};
    static emxArray_real_T wvec;
    static int sizeOfWvec[1] = {7};
    static double wvecData[7] = {
        0.1744,
        -0.0529,
        0.4017,
        0.4050,
        0.2125,
        -0.0343,
        0.2773
    };
    static emxArray_real_T ilog;
    static int sizeOfIlog[2] = {1, 1};
    static double ilogData[1] = {2};
    static emxArray_real_T ftmin;
    static int sizeOfFtmin[1] = {7};
    static double ftminData[7] = {0.0144, -3.8249, 0.1433, 0.0475, 0.0156, 0, 10.7954};
    static emxArray_real_T ftmax;
    static int sizeOfFtmax[1] = {7};
    static double ftmaxData[7] = {1.7500, 2.1518, 1.8985, 2.9557, 3.9239, 2.1887, 208.9349};
    
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        ftvec.size = sizeOfFtvec;
        wvec.size = sizeOfWvec;
        wvec.data = wvecData;
        ilog.size = sizeOfIlog;
        ilog.data = ilogData;
        ftmin.size = sizeOfFtmin;
        ftmin.data = ftminData;
        ftmax.size = sizeOfFtmax;
        ftmax.data = ftmaxData;
    });
    
    ftvec.data = ft;
    double score = features_ufb(&ftvec, &wvec, &ilog, &ftmin, &ftmax, fbmin, fbmax);
    
    return score;
}

+ (double)scoreFromPhonationTest:(NSURL *)phonationAudioFile
{
    OSStatus err = noErr;
    
    ExtAudioFileRef audiofile;
    err = ExtAudioFileOpenURL((__bridge CFURLRef)phonationAudioFile, &audiofile);
    if (err || !audiofile) {
        NSLog(@"Failed to open audio file %@, OSStatus error: %d", phonationAudioFile.absoluteString, (int)err);
        return NAN;
    }
    
    // get some info about the file's format.
    AudioStreamBasicDescription fileFormat;
    UInt32 size = sizeof(fileFormat);
    err = ExtAudioFileGetProperty(audiofile, kExtAudioFileProperty_FileDataFormat, &size, &fileFormat);
    if (err) {
        NSLog(@"OSStatus error getting audio file property FileDataFormat: %d", (int)err);
        return NAN;
    }

    // tell the ExtAudioFile API what format we want samples back in
    AudioStreamBasicDescription clientFormat;
    bzero(&clientFormat, sizeof(clientFormat));
    clientFormat.mChannelsPerFrame = fileFormat.mChannelsPerFrame;
    clientFormat.mBytesPerFrame = 4;
    clientFormat.mBytesPerPacket = clientFormat.mBytesPerFrame;
    clientFormat.mFramesPerPacket = 1;
    clientFormat.mBitsPerChannel = 32;
    clientFormat.mFormatID = kAudioFormatLinearPCM;
    clientFormat.mSampleRate = fileFormat.mSampleRate;
    clientFormat.mFormatFlags = kLinearPCMFormatFlagIsFloat | kAudioFormatFlagIsNonInterleaved;
    err = ExtAudioFileSetProperty(audiofile, kExtAudioFileProperty_ClientDataFormat, sizeof(clientFormat), &clientFormat);
    if (err) {
        NSLog(@"OSStatus error setting audio file property ClientDataFormat: %d", (int)err);
        return NAN;
    }
    
    // find out how many frames we need to read
    SInt64 numFrames = 0;
    size = sizeof(numFrames);
    err = ExtAudioFileGetProperty(audiofile, kExtAudioFileProperty_FileLengthFrames, &size, &numFrames);
    if (err) {
        NSLog(@"OSStatus error getting audio file property FileLengthFrames: %d", (int)err);
        return NAN;
    }
    
    // create the buffers for reading in data
    AudioBufferList *bufferList = malloc(sizeof(AudioBufferList) + sizeof(AudioBuffer) * (clientFormat.mChannelsPerFrame - 1));
    bufferList->mNumberBuffers = clientFormat.mChannelsPerFrame;
    for (int ii=0; ii < bufferList->mNumberBuffers; ++ii) {
        bufferList->mBuffers[ii].mDataByteSize = (UInt32)(sizeof(float) * numFrames);
        bufferList->mBuffers[ii].mNumberChannels = 1;
        bufferList->mBuffers[ii].mData = malloc(bufferList->mBuffers[ii].mDataByteSize);
    }
    
    // read in the data
    UInt32 rFrames = (UInt32)numFrames;
    err = ExtAudioFileRead(audiofile, &rFrames, bufferList);
    if (err) {
        NSLog(@"OSStatus error %d reading audio file %@", (int)err, phonationAudioFile.absoluteString);
        // release the AudioToolbox buffers
        for (int ii=0; ii < bufferList->mNumberBuffers; ++ii) {
            free(bufferList->mBuffers[ii].mData);
        }
        free(bufferList);
        return NAN;
    }
    
    // close the file
    err = ExtAudioFileDispose(audiofile);
    if (err) {
        NSLog(@"OSStatus error %d closing audio file %@", (int)err, phonationAudioFile.absoluteString);
    }
    
    // convert to emxArray_real_T for use with MATLAB converted code
    emxArray_real_T audio;
    int sizeOfAudio[1] = {(int)numFrames};
    audio.size = sizeOfAudio;
    audio.numDimensions = 1;
    audio.canFreeData = true;
    audio.allocatedSize = (int)(sizeof(double) * numFrames);
    audio.data = malloc(audio.allocatedSize);
    float *pSrc = bufferList->mBuffers[0].mData;
    float *pSrcEnd = pSrc + numFrames;
    double *pDst = audio.data;
    while (pSrc < pSrcEnd) {
        *pDst++ = (double)*pSrc++;
    }
    
    // release the AudioToolbox buffers
    for (int ii=0; ii < bufferList->mNumberBuffers; ++ii) {
        free(bufferList->mBuffers[ii].mData);
    }
    free(bufferList);
    bufferList = NULL;
    
    // get the features array
    double ft[13];
    double srate = fileFormat.mSampleRate;
    features_bvav2(&audio, srate, ft);
    for (int i = 0; i < 13; ++i) {
        if (isnan(ft[i]))
            return NAN;
    }
    
    // clean up
    free(audio.data);
    
    // generate a score
    static emxArray_real_T ftvec;
    static int sizeOfFtvec[2] = {1, 13};
    static emxArray_real_T wvec;
    static int sizeOfWvec[1] = {13};
    static double wvecData[13] = {
        0.2862,
        0.2612,
        0.0989,
        -0.0045,
        0.2607,
        0.0151,
        -0.1823,
        0.0008,
        0.0663,
        0.0353,
        0.0258,
        0.0086,
        -0.1278
    };
    static emxArray_real_T ilog;
    static int sizeOfIlog[2] = {1, 7};
    static double ilogData[7] = {2, 3, 4, 10, 11, 12, 13};
    static emxArray_real_T ftmin;
    static int sizeOfFtmin[1] = {13};
    static double ftminData[13] = {75.8734,  0.9091, 0.7501, -1.2153, 0.0284, -1.1089, -1.7412, -10.3110, -17.7297, 0.7543, 0.6865, 0.8511, 0.7226};
    static emxArray_real_T ftmax;
    static int sizeOfFtmax[1] = {13};
    static double ftmaxData[13] = {466.9711, 3.5756, 2.2598, 0.1496, 0.9342, 35.4322, 16.9814, 11.8234, 11.5353, 2.2795, 2.0381, 2.0521, 1.9541};

    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        ftvec.size = sizeOfFtvec;
        wvec.size = sizeOfWvec;
        wvec.data = wvecData;
        ilog.size = sizeOfIlog;
        ilog.data = ilogData;
        ftmin.size = sizeOfFtmin;
        ftmin.data = ftminData;
        ftmax.size = sizeOfFtmax;
        ftmax.data = ftmaxData;
    });
    
    ftvec.data = ft;
    double score = features_ufb(&ftvec, &wvec, &ilog, &ftmin, &ftmax, fbmin, fbmax);
    
    return score;
}

+ (double)scoreFromPostureTest:(NSArray *)postureData
{
    // format the posture data for the MATLAB code
    emxArray_real_T posture;
    int sizeOfPosture[2] = {(int)[postureData count], 4};
    posture.size = sizeOfPosture;
    posture.data = malloc(sizeof(double) * sizeOfPosture[0] * sizeOfPosture[1]);
    double *pT = posture.data;
    double *pX = pT + sizeOfPosture[0];
    double *pY = pX + sizeOfPosture[0];
    double *pZ = pY + sizeOfPosture[0];
    for (int i = 0; i < sizeOfPosture[0]; ++i) {
        NSDictionary *accel = [postureData objectAtIndex:i];
        *pT++ = [[accel objectForKey:@"timestamp"] doubleValue];
        *pX++ = [[accel objectForKey:@"x"] doubleValue];
        *pY++ = [[accel objectForKey:@"y"] doubleValue];
        *pZ++ = [[accel objectForKey:@"z"] doubleValue];
    }
    
    // get the features array
    double ft[3];
    features_bpa(&posture, ft);
    for (int i = 0; i < 3; ++i) {
        if (isnan(ft[i]))
            return NAN;
    }

    static emxArray_real_T ftvec;
    static int sizeOfFtvec[2] = {1, 3};
    static emxArray_real_T wvec;
    static int sizeOfWvec[1] = {3};
    static double wvecData[3] = {
        0.3222,
        0.1695,
        0.1633
    };
    static emxArray_real_T ilog;
    static int sizeOfIlog[2] = {1, 2};
    static double ilogData[2] = {1, 2};
    static emxArray_real_T ftmin;
    static int sizeOfFtmin[1] = {3};
    static double ftminData[3] = {-1.9092, -3.9434, 0.7115};
    static emxArray_real_T ftmax;
    static int sizeOfFtmax[1] = {3};
    static double ftmaxData[3] = {0.1923, 2.1383, 1.5896};
    
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        ftvec.size = sizeOfFtvec;
        wvec.size = sizeOfWvec;
        wvec.data = wvecData;
        ilog.size = sizeOfIlog;
        ilog.data = ilogData;
        ftmin.size = sizeOfFtmin;
        ftmin.data = ftminData;
        ftmax.size = sizeOfFtmax;
        ftmax.data = ftmaxData;
    });
    
    ftvec.data = ft;
    double score = features_ufb(&ftvec, &wvec, &ilog, &ftmin, &ftmax, fbmin, fbmax);
    
    return score;
}

+ (double)scoreFromTappingTest:(NSArray *)tappingData
{
    // format the tapping data for the MATLAB code
    emxArray_real_T tapping;
    int sizeOfTapping[2] = {(int)[tappingData count], 3};
    tapping.size = sizeOfTapping;
    tapping.data = malloc(sizeof(double) * sizeOfTapping[0] * sizeOfTapping[1]);
    double *pT = tapping.data;
    double *pX = pT + sizeOfTapping[0];
    double *pY = pX + sizeOfTapping[0];
    for (int i = 0; i < sizeOfTapping[0]; ++i) {
        NSDictionary *tap = [tappingData objectAtIndex:i];
        NSString *tapCoord = [tap objectForKey:@"TapCoordinate"];
        NSCharacterSet *curlyBraces = [NSCharacterSet characterSetWithCharactersInString:@"{}"];
        tapCoord = [tapCoord stringByTrimmingCharactersInSet:curlyBraces];
        NSArray *tapCoords = [tapCoord componentsSeparatedByString:@", "];
        
        *pT++ = [[tap objectForKey:@"TapTimeStamp"] doubleValue];
        *pX++ = [[tapCoords objectAtIndex:0] doubleValue];
        *pY++ = [[tapCoords objectAtIndex:1] doubleValue];
    }
    
    // get the features array
    double ft[3];
    features_bta(&tapping, ft);
    for (int i = 0; i < 3; ++i) {
        if (isnan(ft[i]))
            return NAN;
    }
    
    static emxArray_real_T ftvec;
    static int sizeOfFtvec[2] = {1, 2};
    static emxArray_real_T wvec;
    static int sizeOfWvec[1] = {2};
    static double wvecData[2] = {
        0.0980,
        0.2554
    };
    static emxArray_real_T ilog;
    static int sizeOfIlog[2] = {1, 2};
    static double ilogData[2] = {1, 2};
    static emxArray_real_T ftmin;
    static int sizeOfFtmin[1] = {2};
    static double ftminData[2] = {-1.3010, -2.0000};
    static emxArray_real_T ftmax;
    static int sizeOfFtmax[1] = {2};
    static double ftmaxData[2] = {0.6781, 0.7158};
    
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        ftvec.size = sizeOfFtvec;
        wvec.size = sizeOfWvec;
        wvec.data = wvecData;
        ilog.size = sizeOfIlog;
        ilog.data = ilogData;
        ftmin.size = sizeOfFtmin;
        ftmin.data = ftminData;
        ftmax.size = sizeOfFtmax;
        ftmax.data = ftmaxData;
    });
    
    ftvec.data = ft;
    double score = features_ufb(&ftvec, &wvec, &ilog, &ftmin, &ftmax, fbmin, fbmax);
    
    return score;
}

@end
