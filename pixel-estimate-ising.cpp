/* Enlarges source image trying to estimate missing pixels - Ising Model
   (c) 2016, Eduardo Valle, blog.eduardovalle.com/
*/
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include <gd.h>

using namespace std;

/* Main function */
int main(int argc, char *argv[]) {

    // Check command-line and open source image
    const int nParamsMin = 2;
    const int nParamsMax = 9;
    if ( argc <= nParamsMin || argc>(nParamsMax+1) ) {
        cerr << "Usage: " << argv[0]
             << " <in.png> <out.png> [<dom>] [<asp>] [H up] [V up] [gweight] [mindelta] [maxiter]"
             << endl
             << endl
             << "   in.png    Input image in PNG format (mandatory)" << endl
             << "   out.png   Output image in PNG format (mandatory)" << endl
             << "   dom       Pixel connectedness domain: 4 or 8 (def: 4)" << endl
             << "   asp       Pixel aspect ratio width/height (def: 1.0)" << endl
             << "   H up      How much to upscale the image horizontally (def: 3)" << endl
             << "   V up      How much to upscale the image vertically (def: 3)" << endl
             << "   gweight   Relative importance of respecting data (def: 1.0 = neutral)" << endl
             << "   mindelta  Minimum total energy delta for convergence (def: 0.001)" << endl
             << "   maxiter   Maximum epochs for convergence (def: 100)" << endl
             << endl
             << "   Optional parameters must be entered in the right sequence" << endl
             << endl;
        return 1;
    }
    FILE *infile = fopen(argv[1], "r");
    if (infile == NULL) {
        cerr << "Error opening:" << argv[1] << endl;
        return 1;
    }
    gdImagePtr simage = gdImageCreateFromPng(infile);
    fclose(infile);
    if (simage == NULL) {
        cerr << "Error reading:" << argv[1] << endl;
        return 1;
    }
    const int w = simage->sx;
    const int h = simage->sy;
    const float maxGray = 255.0;

    // Are pixel and sample neighborhoods 4-connected or 8-connected ?
    const int connectedness = (argc <= 3) ? 4 : atoi(argv[3]);
    if (connectedness!=4 && connectedness!=8) {
        cerr << "Connectedness domain must be 4 or 8: " << argv[3] << endl;
        return 1;
    }
    const bool  four = (connectedness == 4);

    // What is the relative importance of neighbors
    const float aspect     = (argc <= 4) ? 1.0 : atof(argv[4]);
    if (aspect <= 0.0f) {
        cerr << "Aspect ratio must be >= 0.0: " << argv[4] << endl;
        return 1;
    }
    const float weightV = 1.0f;
    const float weightH = 1.0f/aspect;
    const float weightD = 1.0f/sqrt(1.0f + aspect*aspect);
    const float weightN = 1.0f /
        (four ? (2*weightH + 2*weightV) : (4*weightD + 2*weightV + 2*weightH));

    // Horizontal, vertical, and surface ampliation
    const int upX = (argc <= 5) ? 3 : atoi(argv[5]);
    const int upY = (argc <= 6) ? 3 : atoi(argv[6]);
    if (upX < 1) {
        cerr << "Horizontal upscaling must be 1 or greater: " << argv[5] << endl;
        return 1;
    }
    if (upY < 1) {
        cerr << "Vertical upscaling must be 1 or greater: " << argv[6] << endl;
        return 1;
    }
    const int up2D = upX*upY;

    // Relative importance of attachment to graylevels
    const float grayWeight  = (argc <= 7) ? 1.0 : atof(argv[7]);
    if (grayWeight <= 0.0) {
        cerr << "Grayscale attachment weight must be > 0.0: " << argv[7] << endl;
        return 1;
    }

    // Convergence parameters
    const float minDeltaDiff = (argc <= 8) ? 0.001f : atof(argv[8]);
    if (minDeltaDiff < 0.0f) {
        cerr << "Minimum delta for convergence must be >= 0.0: " << argv[8] << endl;
        return 1;
    }
    const int maxEpochs = (argc <= 9) ? 100 : atoi(argv[9]);
    if (maxEpochs <= 0) {
        cerr << "Maximum epochs for convergence must be >0: " << argv[9] << endl;
        return 1;
    }


    // First pass : reads pixels matrix and count the number of needed variables
    const int ws = w+2;
    const int hs = h+2;
    float samples[hs][ws];
    for (int y=0, ys=1; y<h; y++, ys++) {
        for (int x=0, xs=1; x<w; x++, xs++) {
            const int rgb = gdImageGetTrueColorPixel(simage, x, y);
            const float r = rgb >> 16 & 0xff;
            const float g = rgb >> 8  & 0xff;
            const float b = rgb       & 0xff;
            // Luminance coefficient for each color
            const float rl = 0.21f;
            const float gl = 0.72f;
            const float bl = 0.07f;
            float l = (r*rl + g*gl + b*bl) / maxGray;
            // Quantizes value into one of the values representable by the number of samples ?
            // l <= value quantized between 0 and 1.0 in 1/up2D levels ?
            // l = round(l*up2D) / up2D;
            // Destination pixels
            samples[ys][xs] = l;
        }
    }
    // ...correct borders (tries to create uniform border) to ease processing
    for (int ys=0; ys<hs; ys++) {
        samples[ys][0   ] = samples[   1][   1];
        samples[ys][ws-1] = samples[hs-2][ws-2];
    }
    for (int xs=0; xs<ws; xs++) {
        samples[   0][xs] = samples[   1][   1];
        samples[hs-1][xs] = samples[hs-2][ws-2];
    }

    gdImageDestroy(simage);

    // Assembles optimization matrix...

    // The optimization we'll attempt tries to minimize the
    // difference between neighboring subsamples, while constraining all
    // subsamples inside a pixel to average to the given pixel.
    //
    // The total model energy will be
    //      gw * Sum_pi |pi - {Avg of subsamples in pi}| + Sum_{si,sj neighbors}|si-sj|
    // where gw weights the relative importance of the terms
    //
    // In the model the samples may assume one of a small number of discrete values
    // in our case 0 or 1


    default_random_engine            generator;
    uniform_real_distribution<float> uniform(0.0, 1.0);

    // Initializes model with random samples in the right proportion
    cerr << "Initializing subsamples..." << endl;
    int  counts[hs][ws];
    char subsamps[hs][ws][upY][upX];
    for (int ys=0; ys<hs; ys++) {
        for (int xs=0; xs<ws; xs++) {
            const int g = round(samples[ys][xs]*up2D);
            for (int ss=0; ss<up2D; ss++) {
                subsamps[ys][xs][0][ss] = (ss < g) ? 1 : 0;
            }
            if (g>0 && g<up2D) {
                shuffle(&subsamps[ys][xs][0][0], &subsamps[ys][xs][0][up2D], generator);
            }
            counts[ys][xs] = g;
        }
    }

    // DEBUG
    // Outputs intermediary mask
    cerr << "Outputting initial mask..." << endl;
    gdImagePtr timage = gdImageCreateTrueColor(ws*upX, hs*upY);
    if (timage == NULL) {
        fprintf(stderr, "Not enough memory.\n");
        return 1;
    }
    int tgray[256];
    for (int l=0; l<256; l++) {
        tgray[l] = gdImageColorAllocate(timage, l, l, l);
    }
    for (int yu=0; yu<(hs*upY); yu++) {
        for (int xu=0; xu<(ws*upX); xu++) {
            int l = 255 * subsamps[yu/upY][xu/upX][yu%upY][xu%upX];
            gdImageSetPixel(timage, xu, yu, tgray[l]);
        }
    }
    // Save results to output file
    FILE *tfile = fopen("debug.png", "w");
    if (tfile == NULL) {
        cerr << "Error opening: debug.png" << endl;
        return 1;
    }
    gdImagePng(timage, tfile);
    fclose(tfile);
    gdImageDestroy(timage);


    // Decides on a visitation order for the subsamples
    cerr << "Initializing route..." << endl;
    vector<tuple<int, int, int, int> > route;
    for (int ys=1; ys<hs-1; ys++) {
        for (int xs=1; xs<ws-1; xs++) {
            for (int v=0; v<upY; v++) {
                for (int u=0; u<upX; u++) {
                    route.push_back(make_tuple(xs, ys, u, v));
                }
            }
        }
    }
    shuffle(route.begin(), route.end(), generator);

    // Finds neighbor n, m of subsample at [y][x][v][u]
    // i.e., subsample at [y][x][v+n][u+m] (n, m = {-1, 0, 1}) compensating for pixel frotiers
    const auto neighbor = [&](int x, int y, int u, int v, int m, int n) -> char {
        assert(u>=0 && u<upX && v>=0 && v<upY && m>=-1 && m<=1 && n>=-1 && n<=1);
        int nu = u+m, nv = v+n;
        // Safe cases : neighbor is in same pixel
        if (nu>=0 && nu<upX && nv>=0 && nv<upY) {
            return subsamps[y][x][nv][nu];
        }
        // Problematic cases : neighbor is accross pixels
        int nx = x,
            ny = y;
        if (nu < 0) {
            if (nv <    0) { nx--; ny--; nu=upX-1; nv=upY-1; } else
            if (nv >= upY) { nx--; ny++; nu=upX-1; nv=    0; } else
                           { nx--;       nu=upX-1; }
        }
        else if (nu >= upX) {
            if (nv <    0) { nx++; ny--; nu=    0; nv=upY-1; } else
            if (nv >= upY) { nx++; ny++; nu=    0; nv=    0; } else
                           { nx++;       nu=    0; }
        }
        else {
            if (nv <    0) { ny--; nv=upY-1; } else
            if (nv >= upY) { ny++; nv=    0; } else
                           { assert(false); }
        }
        assert(nx>=0 && nx<ws && ny>=0 && ny<hs);
        return subsamps[ny][nx][nv][nu];
    };

    // Cost functions
    const auto attachment = [&](float sample, int oldcount, int newcount) -> float {
        const float oldap = ((float) oldcount) / up2D;
        const float newap = ((float) newcount) / up2D;
        return pow(sample-newap, 2.0)-pow(sample-oldap, 2.0);
    };

    const auto neighborhood = [&](char nLU, char nL, char nLD, char nU, char nD,
        char nRU, char nR, char nRD, char oldsub, char newsub) -> float {
        return (( four ? 0.0f : (abs(nLU-newsub)-abs(nLU-oldsub))*weightD ) +
                (               (abs(nL -newsub)-abs(nL -oldsub))*weightH ) +
                ( four ? 0.0f : (abs(nLD-newsub)-abs(nLD-oldsub))*weightD ) +
                (               (abs(nU -newsub)-abs(nU -oldsub))*weightV ) +
                (               (abs(nD -newsub)-abs(nD -oldsub))*weightV ) +
                ( four ? 0.0f : (abs(nRU-newsub)-abs(nRU-oldsub))*weightD ) +
                (               (abs(nR -newsub)-abs(nR -oldsub))*weightH ) +
                ( four ? 0.0f : (abs(nRD-newsub)-abs(nRD-oldsub))*weightD ))*weightN/grayWeight;
    };

    // Iterates, trying to lower energy
    cerr << "Iterating model..." << endl;
    int epoch = 0;
    float totalDeltaEnergy;
    do {
        totalDeltaEnergy = 0.0f;
        int xs, ys, u, v;
        for (auto r : route) {
            tie(xs, ys, u, v) = r;
            const char sub  = subsamps[ys][xs][v][u];
            const char flip = (sub == 0) ? 1 : 0;
            const char nLU  = four ? 0 : neighbor(xs, ys, u, v, -1, -1);
            const char nL   =            neighbor(xs, ys, u, v, -1,  0);
            const char nLD  = four ? 0 : neighbor(xs, ys, u, v, -1,  1);
            const char nU   =            neighbor(xs, ys, u, v,  0, -1);
            const char nD   =            neighbor(xs, ys, u, v,  0,  1);
            const char nRU  = four ? 0 : neighbor(xs, ys, u, v,  1, -1);
            const char nR   =            neighbor(xs, ys, u, v,  1,  0);
            const char nRD  = four ? 0 : neighbor(xs, ys, u, v,  1,  1);
            const float deltaE =
                attachment(samples[ys][xs], counts[ys][xs], counts[ys][xs]+(flip-sub)) +
                neighborhood(nLU, nL, nLD, nU, nD, nRU, nR, nRD, sub, flip);
            if (deltaE < 0.0f) {
                subsamps[ys][xs][v][u] = flip;
                counts[ys][xs]        += flip-sub;
                totalDeltaEnergy      += deltaE;
            }
        }
        epoch++;
        cerr << "Epoch " << epoch << " - delta " << totalDeltaEnergy << endl;
    } while (epoch<maxEpochs && abs(totalDeltaEnergy)>minDeltaDiff);

    // Prepares destination image
    cerr << "Preparing output..." << endl;
    const int wu = (int) (w * upX);
    const int hu = (int) (h * upY);
    gdImagePtr dimage = gdImageCreateTrueColor(wu, hu);
    if (dimage == NULL) {
        fprintf(stderr, "Not enough memory.\n");
        return 1;
    }
    int gray[256];
    for (int l=0; l<256; l++) {
        gray[l] = gdImageColorAllocate(dimage, l, l, l);
    }
    for (int yu=0; yu<hu; yu++) {
        for (int xu=0; xu<wu; xu++) {
            int l = 255 * subsamps[(yu/upY)+1][(xu/upX)+1][yu%upY][xu%upX];
            gdImageSetPixel(dimage, xu, yu, gray[l]);
        }
    }
    // Save results to output file
    FILE *outfile = fopen(argv[2], "w");
    if (outfile == NULL) {
        cerr << "Error opening: " << argv[2] << endl;
        return 1;
    }
    gdImagePng(dimage, outfile);
    fclose(outfile);
    gdImageDestroy(dimage);
    return 0;
}
