#import "bridge_ufb_hand_tooled.h"

PDRealArray *pitchStrengthAllCandidates(PDRealArray *f, PDRealArray *L, PDRealArray *pc);
PDRealArray *pitchStrengthOneCandidate(PDRealArray *f, PDRealArray *NL, double pc);

static PDRealArray *primetable;

static double gPrimes_to_10000[1230] = { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
    107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
    191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
    271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359,
    367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449,
    457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557,
    563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
    647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
    751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853,
    857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953,
    967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049,
    1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129,
    1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231,
    1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319,
    1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439,
    1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523,
    1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609,
    1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709,
    1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
    1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913,
    1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017,
    2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113,
    2129, 2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237,
    2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333,
    2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411,
    2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539,
    2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657,
    2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719,
    2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
    2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 2909, 2917, 2927,
    2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041,
    3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169,
    3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271,
    3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371,
    3373, 3389, 3391, 3407, 3413, 3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491,
    3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581,
    3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671, 3673, 3677,
    3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793,
    3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
    3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007,
    4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111,
    4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229,
    4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337,
    4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451,
    4457, 4463, 4481, 4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561,
    4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663,
    4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 4759, 4783, 4787, 4789,
    4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919,
    4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
    5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107,
    5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233,
    5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381,
    5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471,
    5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569,
    5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683,
    5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801,
    5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879,
    5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037,
    6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
    6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257,
    6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343,
    6353, 6359, 6361, 6367, 6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469,
    6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581,
    6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701, 6703,
    6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827,
    6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947,
    6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027,
    7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177,
    7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
    7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451,
    7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541,
    7547, 7549, 7559, 7561, 7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639,
    7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741,
    7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853, 7867, 7873, 7877,
    7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009,
    8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117,
    8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, 8221, 8231, 8233, 8237,
    8243, 8263, 8269, 8273, 8287, 8291, 8293, 8297, 8311, 8317, 8329, 8353, 8363,
    8369, 8377, 8387, 8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
    8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, 8599, 8609, 8623,
    8627, 8629, 8641, 8647, 8663, 8669, 8677, 8681, 8689, 8693, 8699, 8707, 8713,
    8719, 8731, 8737, 8741, 8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821,
    8831, 8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, 8933, 8941,
    8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, 9013, 9029, 9041, 9043, 9049,
    9059, 9067, 9091, 9103, 9109, 9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181,
    9187, 9199, 9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, 9293,
    9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, 9391, 9397, 9403, 9413,
    9419, 9421, 9431, 9433, 9437, 9439, 9461, 9463, 9467, 9473, 9479, 9491, 9497,
    9511, 9521, 9533, 9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
    9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733, 9739, 9743, 9749,
    9767, 9769, 9781, 9787, 9791, 9803, 9811, 9817, 9829, 9833, 9839, 9851, 9857,
    9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973
};

//function erbs = hz2erbs(hz)
//erbs = 6.44 * ( log2( 229 + hz ) - 7.84 );
double hz2erbs(double hz)
{
    return 6.44 * (log2(229.0 + hz) - 7.84);
}

//function hz = erbs2hz(erbs)
//hz = ( 2 .^ ( erbs./6.44 + 7.84) ) - 229;
double erbs2hz(double erbs)
{
    return pow(2.0, erbs / 6.44 + 7.84) - 229.0;
}


//function [p,t,s] = swipep(x,fs,plim,dt,dlog2p,dERBs,woverlap,sTHR)
void swipep(PDRealArray *x, double fs, double plim[2], double dt, double dlog2p, double dERBs, double woverlap, double sTHR,
            PDRealArray **p, PDRealArray **t, PDRealArray **s)
{
    //% SWIPEP Pitch estimation using SWIPE'.
    //%    P = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,STHR) estimates the pitch
    //%    of the vector signal X every DT seconds. The sampling frequency of
    //%    the signal is Fs (in Hertz). The spectrum is computed using a Hann
    //%    window with an overlap WOVERLAP between 0 and 1. The spectrum is
    //%    sampled uniformly in the ERB scale with a step size of DERBS ERBs. The
    //%    pitch is searched within the range [PMIN PMAX] (in Hertz) with samples
    //%    distributed every DLOG2P units on a base-2 logarithmic scale of Hertz.
    //%    The pitch is fine-tuned using parabolic interpolation with a resolution
    //%    of 1 cent. Pitch estimates with a strength lower than STHR are treated
    //%    as undefined.
    //%
    //%    [P,T,S] = SWIPEP(X,Fs,[PMIN PMAX],DT,DLOG2P,DERBS,WOVERLAP,STHR)
    //%    returns the times T at which the pitch was estimated and the pitch
    //%    strength S of every pitch estimate.
    //%
    //%    P = SWIPEP(X,Fs) estimates the pitch using the default settings PMIN =
    //%    30 Hz, PMAX = 5000 Hz, DT = 0.001 s, DLOG2P = 1/48 (48 steps per
    //%    octave), DERBS = 0.1 ERBs, WOVERLAP = 0.5, and STHR = -Inf.
    //%
    //%    P = SWIPEP(X,Fs,...,[],...) uses the default setting for the parameter
    //%    replaced with the placeholder [].
    //%
    //%    REMARKS: (1) For better results, make DLOG2P and DERBS as small as
    //%    possible and WOVERLAP as large as possible. However, take into account
    //%    that the computational complexity of the algorithm is inversely
    //%    proportional to DLOG2P, DERBS and 1-WOVERLAP, and that the  default
    //%    values have been found empirically to produce good results. Consider
    //%    also that the computational complexity is directly proportional to the
    //%    number of octaves in the pitch search range, and therefore , it is
    //%    recommendable to restrict the search range to the expected range of
    //%    pitch, if any. (2) This code implements SWIPE', which uses only the
    //%    first and prime harmonics of the signal. To convert it into SWIPE,
    //%    which uses all the harmonics of the signal, replace the word
    //%    PRIMES with a colon (it is located almost at the end of the code).
    //%    However, this may not be recommendable since SWIPE' is reported to
    //%    produce on average better results than SWIPE (Camacho and Harris,
    //%    2008).
    //%
    //%    EXAMPLE: Estimate the pitch of the signal X every 10 ms within the
    //%    range 75-500 Hz using the default resolution (i.e., 48 steps per
    //%    octave), sampling the spectrum every 1/20th of ERB, using a window
    //%    overlap factor of 50%, and discarding samples with pitch strength
    //%    lower than 0.2. Plot the pitch trace.
    //%       [x,Fs] = wavread(filename);
    //%       [p,t,s] = swipep(x,Fs,[75 500],0.01,[],1/20,0.5,0.2);
    //%       plot(1000*t,p)
    //%       xlabel('Time (ms)')
    //%       ylabel('Pitch (Hz)')
    //%
    //%    REFERENCES: Camacho, A., Harris, J.G, (2008) "A sawtooth waveform
    //%    inspired pitch estimator for speech and music," J. Acoust. Soc. Am.
    //%    124, 1638-1652.
    //%
    //%    MAINTENANCE HISTORY:
    //%    - Added line 153 to avoid division by zero in line 154 if loudness
    //%      equals zero (06/23/2010).
    //if ~ exist( 'plim', 'var' ) || isempty(plim), plim = [30 5000]; end
    //if ~ exist( 'dt', 'var' ) || isempty(dt), dt = 0.001; end
    //if ~ exist( 'dlog2p', 'var' ) || isempty(dlog2p), dlog2p = 1/48; end
    //if ~ exist( 'dERBs', 'var' ) || isempty(dERBs), dERBs = 0.1; end
    //if ~ exist( 'woverlap', 'var' ) || isempty(woverlap)
    //    woverlap = 0.5;
    //elseif woverlap>1 || woverlap<0
    //    error('Window overlap must be between 0 and 1.')
    //end
    //if ~ exist( 'sTHR', 'var' ) || isempty(sTHR), sTHR = -Inf; end
    // --ignoring for conversion, just make sure callers include all params and use the above defaults as appropriate
    
    //t = [ 0: dt: length(x)/fs ]'; % Times
    double lengthx_over_fs = (double)MAX(x.rows, x.cols) / fs;
    *t = [PDRealArray rowVectorWithStart:0 step:dt cap:lengthx_over_fs];

    //global primetable;
    //if (isempty(primetable))
    //    primetable = [1 primes(10000)];
    //end
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        primetable = [PDRealArray new];
        primetable.rows = 1;
        primetable.cols = 1230;
        primetable.data = gPrimes_to_10000;
    });

    //% Define pitch candidates
    //log2pc = [ log2(plim(1)): dlog2p: log2(plim(2)) ]';
    PDRealArray *log2pc = [[PDRealArray rowVectorWithStart:log2(plim[0]) step:dlog2p cap:log2(plim[1])] transpose];
    //pc = 2 .^ log2pc;
    PDRealArray *pc = [log2pc exp2];
    //S = zeros( length(pc), length(t) ); % Pitch strength matrix
    size_t length_pc = MAX(pc.rows, pc.cols);
    size_t length_t = MAX((*t).rows, (*t).cols);
    PDRealArray *S = zeros(length_pc, length_t);
    //% Determine P2-WSs
    //logWs = round( log2( 8*fs ./ plim ) );
    PDRealArray *plimArray = [PDRealArray new];
    plimArray.rows = 1;
    plimArray.cols = 2;
    plimArray.data = plim;
    PDRealArray *logWs = [plimArray applyReal:^double(const double element) {
        return round(log2(8.0 * fs / element));
    }];
    //ws = 2.^[ logWs(1): -1: logWs(2) ]; % P2-WSs
    PDRealArray *ws = [[PDRealArray rowVectorWithStart:logWs.data[0] step:-1 cap:logWs.data[1]] exp2];
    //pO = 8 * fs ./ ws; % Optimal pitches for P2-WSs
    PDRealArray *p0 = [ws applyReal:^double(const double element) {
        return 8.0 * fs / element;
    }];
    //% Determine window sizes used by each pitch candidate
    //d = 1 + log2pc - log2( 8*fs./ws(1) );
    PDRealArray *d = [log2pc applyReal:^double(const double element) {
        return 1.0 + element - log2(8.0 * fs / ws.data[0]);
    }];
    
    //% Create ERB-scale uniformly-spaced frequencies (in Hertz)
    //fERBs = erbs2hz([ hz2erbs(min(pc)/4): dERBs: hz2erbs(fs/2) ]');
    PDRealArray *hzerbsderbs = [[PDRealArray rowVectorWithStart:hz2erbs([pc min].data[0] / 4.0) step:dERBs cap:hz2erbs(fs / 2.0)] transpose];
    PDRealArray *fERBs = [hzerbsderbs applyReal:^double(const double element) {
        return erbs2hz(element);
    }];
    //for i = 1 : length(ws)
    size_t length_ws = MAX(ws.rows, ws.cols);
    for (size_t i = 0; i < length_ws; ++i) {
        @autoreleasepool {
            //    dn = max( 1, round( 8*(1-woverlap) * fs / pO(i) ) ); % Hop size
            double dn = MAX(1.0, round(8.0 * (1.0 - woverlap) * fs / p0.data[i]));
            //    % Zero pad signal
            //    xzp = [ zeros( ws(i)/2, 1 ); x(:); zeros( dn + ws(i)/2, 1 ) ];
            double ws_i = ws.data[i];
            PDRealArray *f;
            PDRealArray *ti;
            PDComplexArray *X;
            @autoreleasepool {
                PDRealArray *xzp = [PDRealArray new];
                [xzp concatenateColumnVectors:@[ zeros(ws_i / 2.0, 1), x, zeros(dn + ws_i / 2.0, 1)]];
                //    % Compute spectrum
                //    w = hanning( ws(i) ); % Hann window
                PDRealArray *w = hanning(ws_i);
                //    o = max( 0, round( ws(i) - dn ) ); % Window overlap
                size_t o = MAX(0, round(ws_i - dn));
                //    [ X, f, ti ] = specgram( xzp, ws(i), fs, w, o );
                X = specgram(xzp, ws_i, fs, w, o, &f, &ti);
            }
            //    % Select candidates that use this window size
            //    if length(ws) == 1
            //        j=[(pc)]'; k = [];
            //    elseif i == length(ws)
            //        j=find(d-i>-1); k=find(d(j)-i<0);
            //    elseif i==1
            //        j=find(d-i<1); k=find(d(j)-i>0);
            //    else
            //        j=find(abs(d-i)<1); k=1:length(j);
            //    end
            PDIntArray *j;
            PDIntArray *k;
            if (length_ws == 1) {
                // The Matlab code for this case makes no sense to me--the output would be a real array, not an
                // integer or logical one, so couldn't be used to index into or select rows from other arrays. And
                // even if there's some magic in the (apparently superfluous) [(..)] construction, it won't affect that
                // pc will have values in the range 50-500 while find(d...) would have values in the 1-length(d) range.
                // Be that as it may... I don't think it matters as we never get in here anyway.
                j = [[pc applyInt:^size_t(const double element) {
                    return element;
                }] transpose];
                k = [PDIntArray new];
            } else if (i + 1 == length_ws) {
                j = [[d applyInt:^size_t(const double element) {
                    return element - (i + 1) > -1;
                }] find];
                k = [[[d elementsWithIndices:j] applyInt:^size_t(const double element) {
                    return element - (i + 1) < 0;
                }] find];
            } else if (i + 1 == 1) {
                j = [[d applyInt:^size_t(const double element) {
                    return element - (i + 1) < 1;
                }] find];
                k = [[[d elementsWithIndices:j] applyInt:^size_t(const double element) {
                    return element - (i + 1) > 0;
                }] find];
            } else {
                j = [[d applyInt:^size_t(const double element) {
                    return fabs(element - (i + 1)) < 1;
                }] find];
                size_t length_j = MAX(j.rows, j.cols);
                k = [PDIntArray rowVectorFrom:1 to:length_j];
            }
            //    % Compute loudness at ERBs uniformly-spaced frequencies
            //    fERBs = fERBs( find( fERBs > pc(j(1))/4, 1, 'first' ) : end );
            size_t foundfERBs = [[fERBs applyInt:^size_t(const double element) {
                return element > pc.data[j.data[0] - 1] / 4.0;
            }] findFirst:1].data[0];
            PDIntArray *ifERBs = [PDIntArray rowVectorFrom:foundfERBs to:fERBs.rows];
            fERBs = [fERBs elementsWithIndices:ifERBs];
            //    L = sqrt( max( 0, interp1( f, abs(X), fERBs, 'spline', 0) ) );
            PDRealArray *L;
            @autoreleasepool {
                PDRealArray *Linterp = interp1(f, [X abs], fERBs, PDInterp1MethodSpline, 0.0);
                L = [[Linterp applyReal:^double(const double element) {
                    return isnan(element) || element < 0.0 ? 0.0 : element;
                }] sqrt];
                X = nil;
            }
            //    % Compute pitch strength
            //    Si = pitchStrengthAllCandidates( fERBs, L, pc(j) );
            PDRealArray *Si;
            @autoreleasepool {
                Si = pitchStrengthAllCandidates(fERBs, L, [pc elementsWithIndices:j]);
                L = nil;
            }
            //    % Interpolate pitch strength at desired times
            //    if size(Si,2) > 1
            if (Si.cols > 1) {
                //        warning off MATLAB:interp1:NaNinY
                //        Si = interp1( ti, Si', t, 'linear', NaN )';
                Si = [interp1(ti, [Si transpose], *t, PDInterp1MethodLinear, NAN) transpose];
                //        warning on MATLAB:interp1:NaNinY
                //    else
            } else {
                
                //        Si = repmat( NaN, length(Si), length(t) );
                size_t length_Si = MAX(Si.rows, Si.cols);
                Si = NaN(length_Si, length_t);
                //    end
            }
            //    % Add pitch strength to combination
            //    lambda = d( j(k) ) - i;
            PDRealArray *lambda = [[d elementsWithIndices:[j elementsWithIndices:k]] subtract:i + 1];
            //    mu = ones( size(j) );
            PDRealArray *mu = ones(j.rows, j.cols);
            //    mu(k) = 1 - abs( lambda );
            PDRealArray *one_minus_abs_lambda = [[lambda abs] applyReal:^double(const double element) {
                return 1.0 - element;
            }];
            [mu setElementsWithIndices:k fromArray:one_minus_abs_lambda];
            //    S(j,:) = S(j,:) + repmat(mu,1,size(Si,2)) .* Si;
            PDRealArray *reppedmu_times_si = [repmat(mu, 1, Si.cols) applyReal:^double(const double element, const double otherArrayElement) {
                return element * otherArrayElement;
            } withRealArray:Si];
            PDIntArray *S_cols = [PDIntArray rowVectorFrom:1 to:S.cols];
            PDRealArray *S_j = [S subarrayWithRowIndices:j columnIndices:S_cols];
            PDRealArray *S_j_plus_reppedmu = [S_j applyReal:^double(const double element, const double otherArrayElement) {
                return element + otherArrayElement;
            } withRealArray:reppedmu_times_si];
            [S setElementsWithRowIndices:j columnIndices:S_cols fromArray:S_j_plus_reppedmu];
            //        size_t j_size = j.rows * j.cols;
            //        for (size_t jIdx = 0; jIdx < j_size; ++jIdx) {
            //            PDIntArray *SColumns = [PDIntArray rowVectorFrom:1 to:S.cols];
            //            PDRealArray *reppedmu = [repmat(mu, 1, Si.cols) applyReal:^double(const double element, const double otherArrayElement) {
            //                return element * otherArrayElement;
            //            } withRealArray:Si];
            //            PDRealArray *S_j = [S subarrayWithRowIndices:j columnIndices:SColumns];
            //            PDRealArray *S_plus_repmat = [S_j applyReal:^double(const double element, const double otherArrayElement) {
            //                return element + otherArrayElement;
            //            } withRealArray:reppedmu];
            //            [S setElementsWithRowIndices:j columnIndices:SColumns fromArray:S_plus_repmat];
            //        }
            //end
        }
    }

    //% Fine tune pitch using parabolic interpolation
    //p = repmat( NaN, size(S,2), 1 );
    *p = NaN(S.cols, 1);
    //s = repmat( NaN, size(S,2), 1 );
    *s = NaN(S.cols, 1);
    //for j = 1 : size(S,2)
    for (size_t j = 0; j < S.cols; ++j) {
        @autoreleasepool {
            //    [ s(j), i ] = max( S(:,j), [], 1 );
            NSRange Srows = NSMakeRange(0, S.rows);
            NSRange Sjcol = NSMakeRange(j, 1);
            PDIntArray *i;
            PDRealArray *S_j = [S subarrayWithRows:Srows columns:Sjcol];
            (*s).data[j] = [S_j maxAndIndices:&i].data[0];
            //    if s(j) < sTHR, continue, end
            if ((*s).data[j] < sTHR) {
                continue;
            }
            //    if i == 1 || i == length(pc)
            if (i.data[0] == 1 || i.data[0] == length_pc) {
                //        p(j) = pc(i);
                (*p).data[j] = pc.data[i.data[0] - 1];
                //    else
            } else {
                //        I = i-1 : i+1;
                PDIntArray *I = [PDIntArray rowVectorFrom:i.data[0]-1 to:i.data[0]+1];
                // Turns out this can be done without all the element-by-element inversions (i.e faster)
                // and also without the subtractions (i.e. more accurately). Unfortunately then it differs from Matlab's output
                // from the original code by up to 5 parts in 10^-9. "Improved" code left in, commented out, for reference...
                //        tc = 1 ./ pc(I);
                PDRealArray *tc = [[pc elementsWithIndices:I] oneOverX];
//                PDRealArray *tc = [pc elementsWithIndices:I];
                //        ntc = ( tc/tc(2) - 1 ) * 2*pi;
                PDRealArray *ntc = [[[tc divide:tc.data[1]] subtract:1.0] multiply:2.0 * M_PI];
//                PDRealArray *ntc = [[tc under:tc.data[1]] multiply:2.0 * M_PI];
                //        c = polyfit( ntc, S(I,j), 2 );
                PDRealArray *c = polyfit(ntc, [S_j elementsWithIndices:I], 2);
                // OK, this whole section appears to be finding the local maximum of a parabola by evaluating it over a hundred
                // points per semitone and taking the maximum result. Much simpler, faster, and more precise and accurate would be
                // to simply evaluate the polynomial (c(0)x^2 + c(1)x + c(2)) at the point where the first derivative is 0, i.e.
                // -c(1)/(2 * c(0)) since we know by construction we've got a local maximum and not a minimum. In other words,
                // the s(j) we're looking for is in fact polyval(c, -c(1)/(2 * c(0))) and no need to go searching for it.
                // However, in the interest of maintaining as close equivalence with the original as possible (since the scoring
                // weights etc. are hand-tuned to that calculation) and the results differ ever so slightly, we'll stick with this.
                //        ftc = 1 ./ 2.^[ log2(pc(I(1))): 1/12/100: log2(pc(I(3))) ];
                PDRealArray *ftc = [[[PDRealArray rowVectorWithStart:log2(pc.data[I.data[0] - 1]) step:1./12./100. cap:log2(pc.data[I.data[2] - 1])] exp2] oneOverX];
//                PDRealArray *ftc = [[PDRealArray rowVectorWithStart:log2(pc.data[I.data[0] - 1]) step:1./12./100. cap:log2(pc.data[I.data[2] - 1])] exp2];
                //        nftc = ( ftc/tc(2) - 1 ) * 2*pi;
                PDRealArray *nftc = [[[ftc divide:tc.data[1]] subtract:1.] multiply:2. * M_PI];
//                PDRealArray *nftc = [[ftc under:tc.data[1]] multiply:2. * M_PI];
                //        [s(j) k] = max( polyval( c, nftc ) );
                PDIntArray *k;
                (*s).data[j] = [polyval(c, nftc) maxAndIndices:&k].data[0];
                //        p(j) = 2 ^ ( log2(pc(I(1))) + (k-1)/12/100 );
                (*p).data[j] = exp2( log2(pc.data[I.data[0] - 1]) + (k.data[0] - 1.)/12./100.);
                //    end
            }
        }
    //end
    }
}

    //function S = pitchStrengthAllCandidates( f, L, pc )
PDRealArray *pitchStrengthAllCandidates(PDRealArray *f, PDRealArray *L, PDRealArray *pc)
{
    //% Create pitch strength matrix
    //S = zeros( length(pc), size(L,2) );
    size_t length_pc = MAX(pc.rows, pc.cols);
    PDRealArray *S = zeros(length_pc, L.cols);
    @autoreleasepool {
        //% Define integration regions
        //k = ones( 1, length(pc)+1 );
        PDRealArray *k = ones(1, length_pc + 1);
        //for j = 1 : length(k)-1
        //    k(j+1) = k(j) - 1 + find( f(k(j):end) > pc(j)/4, 1, 'first' );
        //end
        size_t length_k = MAX(k.rows, k.cols);
        NSRange fcols = NSMakeRange(0, 1);
        for (size_t j = 0; j < length_k - 1; ++j) {
            NSRange frows = NSMakeRange(k.data[j] - 1, f.rows - k.data[j] + 1);
            double *p = f.data + frows.location;
            double *pEnd = p + frows.length;
            int idx = 1;
            double pcj_4 = pc.data[j] / 4.0;
            while (p < pEnd) {
                double fx = *p++;
                if (fx > pcj_4) {
                    break;
                }
                ++idx;
            }
            k.data[j + 1] = k.data[j] - 1.0 + idx;
//            PDRealArray *f_k_j_end = [f subarrayWithRows:frows columns:fcols];
//            k.data[j + 1] = k.data[j] - 1.0 + [[f_k_j_end applyInt:^size_t(const double element) {
//                return element > pc.data[j] / 4.0;
//            }] findFirst:1].data[0];
        }
        //k = k(2:end);
        k = [k subarrayWithRows:NSMakeRange(0, 1) columns:NSMakeRange(1, k.cols - 1)];
        //% Create loudness normalization matrix
        //N = sqrt( flipud( cumsum( flipud(L.*L) ) ) );
//        PDRealArray *Nslow = [[[[[L  square] flipud] cumsum] flipud] sqrt]; // uses way too much peak memory
        PDRealArray *N;
        @autoreleasepool {
            N = [L copy];
            // square-flipud-cumsum-flipud-sqrt in place, in one pass
            for (size_t col = 0; col < N.cols; ++col) {
                double *pCol = N.data + col * N.rows;
                double *pBottom = pCol + N.rows - 1;
                double prevfsumfsq = *pBottom * *pBottom;
                *pBottom = sqrt(prevfsumfsq); // square-sqrt (this one not affected by flipud-cumsum-flipud)
                while (--pBottom >= pCol) {
                    double value = *pBottom;
                    value *= value; // square
                    value += prevfsumfsq; // flipud-cumsum-flipud
                    prevfsumfsq = value;
                    *pBottom = sqrt(value); // sqrt
                }
            }
            
            // only keep around the rows of N we're actually gonna end up using
            PDRealArray *kmax = [k max];
            double maxRowIdx = kmax.data[0];
            N = [N subarrayWithRows:NSMakeRange(0, maxRowIdx) columns:NSMakeRange(0, N.cols)];
            
            //for j = 1 : length(pc)
            // do this here because it's more convenient to do it all at once before we use each row
            //    % Normalize loudness
            //    n(n==0) = Inf; % to make zero-loudness equal zero after normalization
            N = [N applyReal:^double(const double element) {
                double value = element;
                if (value == 0.0) {
                    value = INFINITY;
                }
                return value;
            }];
       }
        
        // precompute NL arrays for each value of k[j]--there's only a few
        // on second thought, they're huge so even a few is too many at once
//        NSMutableDictionary *NLbyK = [NSMutableDictionary dictionary];
//        for (size_t j = 0; j < length_pc; ++j) {
//            PDRealArray *NLforkj = [NLbyK objectForKey:@(k.data[j])];
//            if (!NLforkj) {
//                NLforkj = [L divideRows:NSMakeRange(k.data[j] - 1, L.rows - k.data[j] + 1) byRow:k.data[j] - 1 ofRealArray:N];
//                [NLbyK setObject:NLforkj forKey:@(k.data[j])];
//            }
//        }
        PDRealArray *NL = nil;
        double this_k_j = -1; // force a compute of NL first time through
        for (size_t j = 0; j < length_pc; ++j) {
            //    % Normalize loudness
            //    n = N(k(j),:);
//            PDRealArray *n = [N subarrayWithRows:NSMakeRange(k.data[j] - 1, 1) columns:NSMakeRange(0, N.cols)];
//            //    n(n==0) = Inf; % to make zero-loudness equal zero after normalization
//            [n applyReal:^double(const double element) {
//                double value = element;
//                if (value == 0.0) {
//                    value = INFINITY;
//                }
//                return value;
//            }];
            //
            //%     NL = L(k(j):end,:) ./ repmat( n, size(L,1)-k(j)+1, 1);
            //
            //    rowIdx = (1:size(n,1))';
            //        PDIntArray *rowIdx = [[PDIntArray rowVectorFrom:1 to:n.rows] transpose];
            //    colIdx = (1:size(n,2))';
            //        PDIntArray *colIdx = [[PDIntArray rowVectorFrom:1 to:n.cols] transpose];
            //%     NL = L(k(j):end,:) ./ n(rowIdx(:,ones(size(L,1)-k(j)+1,1)), colIdx(:,1));
            //    NL = L(k(j):end,:) ./ n(rowIdx(:,ones(size(L,1)-k(j)+1,1)), colIdx);
//            PDRealArray *L_k_j_end = [L subarrayWithRows:NSMakeRange(k.data[j] - 1, L.rows - k.data[j] + 1) columns:NSMakeRange(0, L.cols)];
            //        PDIntArray *onesLRows_k = ones(1, L.rows - k.data[j] + 1);
            //        PDRealArray *nSubset = [n subarrayWithRowIndices:onesLRows_k columnIndices:colIdx];
            ////        PDRealArray *nSubset = [n subarrayWithRows:NSMakeRange(0, L.rows - k.data[j] + 1) columns:NSMakeRange(0, n.cols)];
            // only recompute NL when it changes
            @autoreleasepool {
                if (k.data[j] != this_k_j) {
                    this_k_j = k.data[j];
                    NL = [L divideRows:NSMakeRange(k.data[j] - 1, L.rows - k.data[j] + 1) byRow:k.data[j] - 1 ofRealArray:N];
                }
//            PDRealArray *NL = [L divideRows:NSMakeRange(k.data[j] - 1, L.rows - k.data[j] + 1) byRow:k.data[j] - 1 ofRealArray:N];
//            PDRealArray *nSubset = repmat(n, L.rows - k.data[j] + 1, 1);
//            PDRealArray *NL = [L_k_j_end divideElementByElement:nSubset];
//            PDRealArray *NL = [L_k_j_end applyReal:^double(const double element, const double otherArrayElement) {
//                return element / otherArrayElement;
//            } withRealArray:nSubset];
                //
                //    % Compute pitch strength
                //    S(j,:) = pitchStrengthOneCandidate( f(k(j):end), NL, pc(j) );
                NSRange jRow = NSMakeRange(j, 1);
                NSRange Scols = NSMakeRange(0, S.cols);
                PDRealArray *psoc = pitchStrengthOneCandidate([f subarrayWithRows:NSMakeRange(k.data[j] - 1, f.rows - k.data[j] + 1) columns:fcols], NL, pc.data[j]);
                [S setSubarrayRows:jRow columns:Scols fromArray:psoc];
            }
            //end
        }
    }

    return S;
}
    //
    //function S = pitchStrengthOneCandidate( f, NL, pc )
PDRealArray *pitchStrengthOneCandidate(PDRealArray *f, PDRealArray *NL, double pc)
{
    PDRealArray *S;
    @autoreleasepool {
        //global primetable;
        //
        //n = fix( f(end)/pc - 0.75 ); % Number of harmonics
        double unfixed_n = f.data[f.rows * f.cols - 1] / pc - 0.75;
        double n = unfixed_n < 0.0 ? ceil(unfixed_n) : floor(unfixed_n);
        
        //if n==0, S=NaN; return, end
        if (n == 0.0) {
            return NaN(1, NL.cols);
        }
        //k = zeros( size(f) ); % Kernel
        PDRealArray *k = zeros(f.rows, f.cols);
        //% Normalize frequency w.r.t. candidate
        //q = f / pc;
        PDRealArray *q = [f divide:pc];
        //% Create kernel
        //% for i = [ 1 primes(n) ]
        //for i = primetable(primetable <= n)
        size_t idx = 0;
        size_t numPrimes = primetable.rows * primetable.cols;
        double i;
        while (idx < numPrimes && (i = primetable.data[idx]) <= n) {
            //    a = abs( q - i );
            PDRealArray *a = [[q subtract:i] abs];
            //    % Peak's weigth
            //    p = a < .25;
            PDIntArray *p = [[a applyInt:^size_t(const double element) {
                return element < .25;
            }] find];
            //    k(p) = cos( 2*pi * q(p) );
            [k setElementsWithIndices:p fromArray:[[[q elementsWithIndices:p] multiply:2.0] cospi]];
            //    % Valleys' weights
            //    v = .25 < a & a < .75;
            PDIntArray *v = [[a applyInt:^size_t(const double element) {
                return .25 < element && element < .75;
            }] find];
            //    k(v) = k(v) + cos( 2*pi * q(v) ) / 2;
            PDRealArray *kv_plus_cos = [[k elementsWithIndices:v] applyReal:^double(const double element, const double otherArrayElement) {
                return element + cos(2 * M_PI * otherArrayElement) / 2.0;
            } withRealArray:[q elementsWithIndices:v]];
            [k setElementsWithIndices:v fromArray:kv_plus_cos];
            //end
            ++idx;
        }
        //% Apply envelope
        //k = k .* sqrt( 1./f  );
        k = [k applyReal:^double(const double element, const double otherArrayElement) {
            return element * sqrt(1. / otherArrayElement);
        } withRealArray:f];
        //% K+-normalize kernel
        //k = k / norm( k(k>0) );
        PDIntArray *pos_ks = [[k applyInt:^size_t(const double element) {
            return element > 0;
        }] find];
        k = [k divide:[[k elementsWithIndices:pos_ks] norm]];
        //% Compute pitch strength
        //S = k' * NL;
        S = [[k transpose] matmult:NL];
    }
    //
    //function erbs = hz2erbs(hz)
    //erbs = 6.44 * ( log2( 229 + hz ) - 7.84 );
    //
    //function hz = erbs2hz(erbs)
    //hz = ( 2 .^ ( erbs./6.44 + 7.84) ) - 229;
    
    return S;
}
