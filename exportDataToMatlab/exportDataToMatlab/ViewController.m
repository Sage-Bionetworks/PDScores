//
//  ViewController.m
//  exportDataToMatlab
//
//  Created by Erin Mounts on 2/13/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "ViewController.h"
@import AudioToolbox;

@interface ViewController()

@property (weak) IBOutlet NSTextField *phonationTextField;
@property (weak) IBOutlet NSTextField *gaitTextField;
@property (weak) IBOutlet NSTextField *postureTextField;
@property (weak) IBOutlet NSTextField *tappingTextField;

@property (weak) IBOutlet NSProgressIndicator *phonationProgressIndicator;
@property (weak) IBOutlet NSProgressIndicator *gaitProgressIndicator;
@property (weak) IBOutlet NSProgressIndicator *postureProgressIndicator;
@property (weak) IBOutlet NSProgressIndicator *tappingProgressIndicator;

@end

@implementation ViewController

void dispatchToMain(dispatch_block_t block) {
    dispatch_async(dispatch_get_main_queue(), block);
}

- (NSURL *)outputURLForFileName:(NSString *)fileName
{
    NSFileManager *fm = [NSFileManager defaultManager];
    
    static NSURL *dirURL = nil;
    
    static dispatch_once_t once;
    dispatch_once(&once, ^{
        dirURL = [fm URLForDirectory:NSDesktopDirectory inDomain:NSUserDomainMask appropriateForURL:nil create:NO error:nil];
        dirURL = [dirURL URLByAppendingPathComponent:@"For Matlab" isDirectory:YES];
        [fm createDirectoryAtURL:dirURL withIntermediateDirectories:YES attributes:nil error:nil];
    });
    
    NSString *outName = [NSString stringWithFormat:@"%@.csv", fileName];
    NSURL *outURL = [dirURL URLByAppendingPathComponent:outName];
    
    return outURL;
}

- (NSURL *)outputURLForFileURL:(NSURL *)fileURL
{
    NSString *bareName = [[fileURL lastPathComponent] stringByDeletingPathExtension];
    
    return [self outputURLForFileName:bareName];
}

- (void)exportPhonationData
{
    NSString *successString = @"Phonation data exported \u2713";
    NSString *failureString = @"Problem exporting phonation data ¯\\_(ツ)_/¯";
    
    NSString *doneString = successString;
    
    dispatchToMain(^{
        [_phonationProgressIndicator startAnimation:nil];
    });
    
    NSURL *audioFileURL = [[NSBundle mainBundle] URLForResource:@"phonation" withExtension:@"m4a"];
    NSURL *outURL = [self outputURLForFileURL:audioFileURL];
    if (!outURL) {
        NSLog(@"Failed to generate output file URL");
        doneString = failureString;
        goto done;
    }

    OSStatus err = noErr;
    
    ExtAudioFileRef audiofile;
    err = ExtAudioFileOpenURL((__bridge CFURLRef)audioFileURL, &audiofile);
    if (err || !audiofile) {
        NSLog(@"Failed to open audio file %@, OSStatus error: %d", audioFileURL.absoluteString, (int)err);
        doneString = failureString;
        goto done;
    }
    
    // get some info about the file's format.
    AudioStreamBasicDescription fileFormat;
    UInt32 size = sizeof(fileFormat);
    err = ExtAudioFileGetProperty(audiofile, kExtAudioFileProperty_FileDataFormat, &size, &fileFormat);
    if (err) {
        NSLog(@"OSStatus error getting audio file property FileDataFormat: %d", (int)err);
        doneString = failureString;
        goto done;
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
        doneString = failureString;
        goto done;
    }
    
    // find out how many frames we need to read
    SInt64 numFrames = 0;
    size = sizeof(numFrames);
    err = ExtAudioFileGetProperty(audiofile, kExtAudioFileProperty_FileLengthFrames, &size, &numFrames);
    if (err) {
        NSLog(@"OSStatus error getting audio file property FileLengthFrames: %d", (int)err);
        doneString = failureString;
        goto done;
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
        NSLog(@"OSStatus error %d reading audio file %@", (int)err, audioFileURL.absoluteString);
        doneString = failureString;
        goto done;
    }
    
    // close the file
    err = ExtAudioFileDispose(audiofile);
    if (err) {
        NSLog(@"OSStatus error %d closing audio file %@", (int)err, audioFileURL.absoluteString);
    }
    
    // write the buffer to the output file as a comma-delimited array of floats
    float *pSrc = bufferList->mBuffers[0].mData;
    float *pSrcEnd = pSrc + numFrames;
    // open file
    const char *fname = [outURL  fileSystemRepresentation];
    FILE *f = fopen(fname, "w");
    char outBuf[1100];
    sprintf(outBuf, "%.1080g", *pSrc++);
    // write to file
    fwrite(outBuf, sizeof(char), strlen(outBuf), f);
    while (pSrc < pSrcEnd) {
        sprintf(outBuf, "\n%.1080g", *pSrc++);
        // write to file
        fwrite(outBuf, sizeof(char), strlen(outBuf), f);
    }
    // close file
    fclose(f);
    
done:
    // release the AudioToolbox buffers
    if (bufferList) {
        for (int ii=0; ii < bufferList->mNumberBuffers; ++ii) {
            free(bufferList->mBuffers[ii].mData);
        }
        free(bufferList);
        bufferList = NULL;
    }
    
    dispatchToMain(^{
        [_phonationProgressIndicator stopAnimation:nil];
        _phonationProgressIndicator.hidden = YES;
        _phonationTextField.stringValue = doneString;
    });
}

- (void)exportGaitPostureData
{
    NSString *successString = @"Gait data exported \u2713";
    NSString *failureString = @"Problem exporting gait data ¯\\_(ツ)_/¯";
    
    NSString *doneString = successString;
    
    NSURL *gpFileURL = [[NSBundle mainBundle] URLForResource:@"gait-posture" withExtension:@"json"];
    NSData *jsonData = [NSData dataWithContentsOfURL:gpFileURL];
    NSDictionary *json = [NSJSONSerialization JSONObjectWithData:jsonData options:0 error:NULL];
    NSArray *keys = [json allKeys];
    
    dispatchToMain(^{
        [_gaitProgressIndicator startAnimation:nil];
    });
    
    NSPredicate *gaitPred = [NSPredicate predicateWithFormat:@"SELF contains %@ OR SELF contains %@", @"accel_walking.outbound-", @"accel_walking.return-"];
    NSArray *gaitKeys = [keys filteredArrayUsingPredicate:gaitPred];
    NSDictionary *gaitData = [json objectForKey:[gaitKeys firstObject]];
    NSArray *gaitItems = [gaitData objectForKey:@"items"];
    NSURL *outGaitURL = [self outputURLForFileName:@"gait"];
    
    if (!gaitItems) {
        NSLog(@"Error reading gait items from gait-posture.json");
        doneString = failureString;
        goto doneGait;
    }
    
    if (!gaitItems.count) {
        NSLog(@"Empty array reading gait items from gait-posture.json");
        doneString = failureString;
        goto doneGait;
    }
    
    if (!outGaitURL) {
        NSLog(@"Failed to generate output file URL");
        doneString = failureString;
        goto doneGait;
    }
    
    // open output file
    const char *fnamegait = [outGaitURL  fileSystemRepresentation];
    FILE *fg = fopen(fnamegait, "w");
    // export one sample per line as timestamp, x, y, z
    for (NSDictionary *item in gaitItems) {
        NSString *timestamp = [item objectForKey:@"timestamp"];
        NSString *x = [item objectForKey:@"x"];
        NSString *y = [item objectForKey:@"y"];
        NSString *z = [item objectForKey:@"z"];
        
        NSString *sampleString = [NSString stringWithFormat:@"%@, %@, %@, %@\n", timestamp, x, y, z];
        const char *outBuf = [sampleString cStringUsingEncoding:NSUTF8StringEncoding];
        fwrite(outBuf, sizeof(char), strlen(outBuf), fg);
    }
    // close file
    fclose(fg);

doneGait:
    dispatchToMain(^{
        [_gaitProgressIndicator stopAnimation:nil];
        _gaitProgressIndicator.hidden = YES;
        _gaitTextField.stringValue = doneString;
    });
    
    successString = @"Posture data exported \u2713";
    failureString = @"Problem exporting posture data ¯\\_(ツ)_/¯";
    
    doneString = successString;
    
    dispatchToMain(^{
        [_postureProgressIndicator startAnimation:nil];
    });
    
    NSPredicate *postPred = [NSPredicate predicateWithFormat:@"SELF contains %@", @"accel_walking.rest-"];
    NSArray *postKeys = [keys filteredArrayUsingPredicate:postPred];
    NSDictionary *postData = [json objectForKey:[postKeys firstObject]];
    NSArray *postItems = [postData objectForKey:@"items"];
    NSURL *outPostURL = [self outputURLForFileName:@"posture"];
    
    if (!postItems) {
        NSLog(@"Error reading posture items from gait-posture.json");
        doneString = failureString;
        goto donePosture;
    }
    
    if (!postItems.count) {
        NSLog(@"Empty array reading posture items from gait-posture.json");
        doneString = failureString;
        goto donePosture;
    }
    
    if (!outPostURL) {
        NSLog(@"Failed to generate output file URL");
        doneString = failureString;
        goto donePosture;
    }
    
    // open output file
    const char *fnamepost = [outPostURL  fileSystemRepresentation];
    FILE *fp = fopen(fnamepost, "w");
    // export one sample per line as timestamp, x, y, z
    for (NSDictionary *item in postItems) {
        NSString *timestamp = [item objectForKey:@"timestamp"];
        NSString *x = [item objectForKey:@"x"];
        NSString *y = [item objectForKey:@"y"];
        NSString *z = [item objectForKey:@"z"];
        
        NSString *sampleString = [NSString stringWithFormat:@"%@, %@, %@, %@\n", timestamp, x, y, z];
        const char *outBuf = [sampleString cStringUsingEncoding:NSUTF8StringEncoding];
        fwrite(outBuf, sizeof(char), strlen(outBuf), fp);
    }
    // close file
    fclose(fp);
    
donePosture:
    dispatchToMain(^{
        [_postureProgressIndicator stopAnimation:nil];
        _postureProgressIndicator.hidden = YES;
        _postureTextField.stringValue = doneString;
    });
}

- (void)exportTappingData
{
    NSString *successString = @"Tapping data exported \u2713";
    NSString *failureString = @"Problem exporting tapping data ¯\\_(ツ)_/¯";
    
    NSString *doneString = successString;
    
    dispatchToMain(^{
        [_tappingProgressIndicator startAnimation:nil];
    });

    NSURL *tappingFileURL = [[NSBundle mainBundle] URLForResource:@"tapping" withExtension:@"json"];
    NSData *jsonData = [NSData dataWithContentsOfURL:tappingFileURL];
    NSDictionary *json = [NSJSONSerialization JSONObjectWithData:jsonData options:0 error:NULL];
    NSDictionary *noName = [json objectForKey:@"NoName.json"];
    NSArray *tapSamples = [noName objectForKey:@"TappingSamples"];
    
    NSURL *outURL = [self outputURLForFileURL:tappingFileURL];
    if (!outURL) {
        NSLog(@"Failed to generate output file URL");
        doneString = failureString;
        goto done;
    }

    // open output file
    const char *fname = [outURL  fileSystemRepresentation];
    FILE *f = fopen(fname, "w");
    // export one sample per line as timestamp, x, y
    for (NSDictionary *tap in tapSamples) {
        NSString *tapCoord = [tap objectForKey:@"TapCoordinate"];
        NSCharacterSet *curlyBraces = [NSCharacterSet characterSetWithCharactersInString:@"{}"];
        tapCoord = [tapCoord stringByTrimmingCharactersInSet:curlyBraces];
        NSNumber *tapTime = [tap objectForKey:@"TapTimeStamp"];
        NSString *sampleString = [NSString stringWithFormat:@"%@, %@\n",tapTime, tapCoord];
        const char *outBuf = [sampleString cStringUsingEncoding:NSUTF8StringEncoding];
        fwrite(outBuf, sizeof(char), strlen(outBuf), f);
    }
    // close file
    fclose(f);
    
done:
    dispatchToMain(^{
        [_tappingProgressIndicator stopAnimation:nil];
        _tappingProgressIndicator.hidden = YES;
        _tappingTextField.stringValue = doneString;
    });
}

- (void)viewDidLoad {
    [super viewDidLoad];

    // Do any additional setup after loading the view.
    dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
        [self exportPhonationData];
        [self exportGaitPostureData];
        [self exportTappingData];
    });
}

- (void)setRepresentedObject:(id)representedObject {
    [super setRepresentedObject:representedObject];

    // Update the view, if already loaded.
}

@end
