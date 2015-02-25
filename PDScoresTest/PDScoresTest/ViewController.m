//
//  ViewController.m
//  PDScoresTest
//
//  Created by Erin Mounts on 2/11/15.
//  Copyright (c) 2015 Sage Bionetworks. All rights reserved.
//

#import "ViewController.h"
#import "PDScores.h"

@interface ViewController ()
@property (weak, nonatomic) IBOutlet UILabel *phonationScoreLabel;
@property (weak, nonatomic) IBOutlet UILabel *gaitScoreLabel;
@property (weak, nonatomic) IBOutlet UILabel *postureScoreLabel;
@property (weak, nonatomic) IBOutlet UILabel *tappingScoreLabel;
@property (weak, nonatomic) IBOutlet UILabel *combinedScoreLabel;

@property (nonatomic) double phonationScore;
@property (nonatomic) double gaitScore;
@property (nonatomic) double postureScore;
@property (nonatomic) double tappingScore;

@end

@implementation ViewController

- (void)scorePhonation
{
    NSURL *audioFileURL = [[NSBundle mainBundle] URLForResource:@"phonation" withExtension:@"m4a"];
    _phonationScore = [PDScores scoreFromPhonationTest:audioFileURL];
    dispatch_async(dispatch_get_main_queue(), ^{
        if (isnan(_phonationScore)) {
            _phonationScoreLabel.text = @"Unable to compute a score";
        } else {
            _phonationScoreLabel.text = [NSString stringWithFormat:@"%f", _phonationScore];
        }
    });
}

- (void)scoreGaitAndPosture
{
    NSURL *gpFileURL = [[NSBundle mainBundle] URLForResource:@"gait-posture" withExtension:@"json"];
    NSData *jsonData = [NSData dataWithContentsOfURL:gpFileURL];
    NSDictionary *json = [NSJSONSerialization JSONObjectWithData:jsonData options:0 error:NULL];
    NSArray *keys = [json allKeys];
    NSPredicate *gaitPred = [NSPredicate predicateWithFormat:@"SELF contains %@ OR SELF contains %@", @"accel_walking.outbound-", @"accel_walking.return-"];
    NSPredicate *postPred = [NSPredicate predicateWithFormat:@"SELF contains %@", @"accel_walking.rest-"];
    NSArray *gaitKeys = [keys filteredArrayUsingPredicate:gaitPred];
    NSArray *postKeys = [keys filteredArrayUsingPredicate:postPred];
    NSDictionary *gaitData = [json objectForKey:[gaitKeys firstObject]];
    NSDictionary *postData = [json objectForKey:[postKeys firstObject]];
    NSArray *gaitItems = [gaitData objectForKey:@"items"];
    NSArray *postItems = [postData objectForKey:@"items"];
    _gaitScore = [PDScores scoreFromGaitTest:gaitItems];
    _postureScore = [PDScores scoreFromPostureTest:postItems];
    dispatch_async(dispatch_get_main_queue(), ^{
        if (isnan(_gaitScore)) {
            _gaitScoreLabel.text = @"Unable to compute a score";
        } else {
            _gaitScoreLabel.text = [NSString stringWithFormat:@"%f", _gaitScore];
        }
        if (isnan(_postureScore)) {
            _postureScoreLabel.text = @"Unable to compute a score";
        } else {
            _postureScoreLabel.text = [NSString stringWithFormat:@"%f", _postureScore];
        }
    });
}

- (void)scoreTapping
{
    NSURL *tappingFileURL = [[NSBundle mainBundle] URLForResource:@"tapping" withExtension:@"json"];
    NSData *jsonData = [NSData dataWithContentsOfURL:tappingFileURL];
    NSDictionary *json = [NSJSONSerialization JSONObjectWithData:jsonData options:0 error:NULL];
    NSDictionary *noName = [json objectForKey:@"NoName.json"];
    NSArray *tapSamples = [noName objectForKey:@"TappingSamples"];
    _tappingScore = [PDScores scoreFromTappingTest:tapSamples];
    dispatch_async(dispatch_get_main_queue(), ^{
        if (isnan(_tappingScore)) {
            _tappingScoreLabel.text = @"Unable to compute a score";
        } else {
            _tappingScoreLabel.text = [NSString stringWithFormat:@"%f", _tappingScore];
        }
    });
}

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    
    dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0), ^{
        [self scorePhonation];
        [self scoreGaitAndPosture];
        [self scoreTapping];
        double combinedScore = [PDScores combinedScoreFromPhonationScore:_phonationScore gaitScore:_gaitScore postureScore:_postureScore tappingScore:_tappingScore];
        dispatch_async(dispatch_get_main_queue(), ^{
            if (isnan(combinedScore)) {
                _combinedScoreLabel.text = @"Unable to compute a score";
            } else {
                _combinedScoreLabel.text = [NSString stringWithFormat:@"%f", combinedScore];
            }
        });
    });
}

- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
