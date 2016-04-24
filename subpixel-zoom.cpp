/* Enlarges source image with bicubic interpolation using subpixel information
   (c) 2016, Eduardo Valle, blog.eduardovalle.com/
*/
#include <iostream>

#include <gd.h>

using namespace std;


/* A possible bicubic interpolation consists in performing a convolution
   with the following kernel :
                            [ 0  2  0  0][s0]
                            [-1  0  1  0][s1]
    p(t) = 0.5 [1 t t^2 t^3][ 2 -5  4 -1][s2]
                            [-1  3 -3  1][s3]
    Where s0..s3 are the samples, and t is the range from 0.0 to 1.0 */
static inline float bicubic(float s0, float s1, float s2, float s3, float t) {
    float r0 = 0.5f * (         2.0f*s1               );
    float r1 = 0.5f * (   -s0           +      s2     );
    float r2 = 0.5f * (2.0*s0  -5.0f*s1 + 4.0f*s2  -s3);
    float r3 = 0.5f * (   -s0 + 3.0f*s1  -3.0f*s2 + s3);
    return r3*t*t*t + r2*t*t + r1*t + r0;
}


/* sYX is a matrix with the samples (row-major)
   xt and yt are the interpolation fractions in the two directions ranging from
   0.0 to 1.0 */
static inline float bicubic2D(
    float s00, float s01, float s02, float s03,
    float s10, float s11, float s12, float s13,
    float s20, float s21, float s22, float s23,
    float s30, float s31, float s32, float s33, float xt, float yt) {
    // The bicubic convolution consists in passing the bicubic kernel in x
    // and then in y (or vice-versa, really)
    float r0 = bicubic(s00, s01, s02, s03, xt);
    float r1 = bicubic(s10, s11, s12, s13, xt);
    float r2 = bicubic(s20, s21, s22, s23, xt);
    float r3 = bicubic(s30, s31, s32, s33, xt);
    return bicubic(r0, r1, r2, r3, yt);
}


/* Main function */
int main(int argc, char *argv[]) {

    // Check command-line and open source image
    const int nParamsMin = 2;
    const int nParamsMax = 4;
    if ( argc <= nParamsMin || argc>(nParamsMax+1) ) {
        cerr << "Usage: " << argv[0]
             << " <in.png> <out.png> [<upsampling> (default: 3.0)] [samples_out.png]" << endl
             << endl
             << "   Assumption: subpixels ordered as R-G-B, horizontally." << endl;
        return 1;
    }
    FILE *infile = fopen(argv[1], "r");
    if (infile == NULL) {
        cerr << "Error opening: " <<  argv[1] << endl;
        return 1;
    }
    gdImagePtr simage = gdImageCreateFromPng(infile);
    fclose(infile);
    if (simage == NULL) {
        cerr << "Error reading: " << argv[1] << endl;
        return 1;
    }
    const int w = simage->sx;
    const int h = simage->sy;
    if (w <= 4 || h <= 4) {
        cerr << "Source image must be at least 4-pixels wide and height" << endl;
        return 1;
    }
    const float upsampling = (argc <= 3) ? 3.0 : atof(argv[3]);
    if (upsampling < 3.0) {
        cerr << "Upsampling must be >= 3.0: " << argv[3] << endl;
        return 1;
    }
    const bool outputSamples = !(argc <= 4);

    // First pass : turns subpixels into pixels, widening 300%
    const int ws = w*3 + 3; // +3 pixels to accomodate 4x4 convolution
    const int hs = h   + 3;
    float samp[hs][ws];
    for (int y=0, ys=1; y<h; y++, ys++) {
        for (int x=0, xs=1; x<w; x++, xs+=3) {
            // Reads RGB value
            int rgb = gdImageGetTrueColorPixel(simage, x, y);
            float r = rgb >> 16 & 0xff;
            float g = rgb >> 8  & 0xff;
            float b = rgb       & 0xff;
            // Luminance coefficient for each color
            const float rl = 1.0f; // In original colorspace: 0.21f
            const float gl = 1.0f; // In original colorspace: 0.72f
            const float bl = 1.0f; // In original colorspace: 0.07f
            // Destination pixels
            samp[ys][xs]   = r*rl;
            samp[ys][xs+1] = g*gl;
            samp[ys][xs+2] = b*bl;
        }
    }
    // ...correct borders (copy pixels) to ease convolution
    for (int ys=0; ys<hs; ys++) {
        samp[ys][0]    = samp[ys][1];
        samp[ys][ws-2] = samp[ys][ws-3];
        samp[ys][ws-1] = samp[ys][ws-3];
    }
    for (int xs=0; xs<ws; xs++) {
        samp[0][xs]    = samp[1][xs];
        samp[hs-2][xs] = samp[hs-3][xs];
        samp[hs-1][xs] = samp[hs-3][xs];
    }

    gdImageDestroy(simage);

    if (outputSamples) {
        // Outputs temporary samples image
        gdImagePtr timage = gdImageCreateTrueColor(ws, hs);
        if (timage == NULL) {
            cerr << "Not enough memory" << endl;
            return 1;
        }
        int tgray[256];
        for (int l=0; l<256; l++) {
            tgray[l] = gdImageColorAllocate(timage, l, l, l);
        }
        for (int ys=0; ys<hs; ys++) {
            for (int xs=0; xs<ws; xs++) {
                gdImageSetPixel(timage, xs, ys, tgray[(int) samp[ys][xs]]);
            }
        }
        FILE *tempfile = fopen(argv[4], "w");
        if (tempfile == NULL) {
        cerr << "Error opening: " << argv[4] << endl;
            return 1;
        }
        gdImagePng(timage, tempfile);
        fclose(tempfile);
        gdImageDestroy(timage);
    }

    // Prepares destination image
    const int wu = (int) (w * upsampling);
    const int hu = (int) (h * upsampling);
    gdImagePtr dimage = gdImageCreateTrueColor(wu, hu);
    if (dimage == NULL) {
        cerr << "Not enough memory" << endl;
        return 1;
    }
    int gray[256];
    for (int l=0; l<256; l++) {
        gray[l] = gdImageColorAllocate(dimage, l, l, l);
    }

    // Second pass : resample image
    for (int yu=0; yu<hu; yu++) {
        for (int xu=0; xu<wu; xu++) {
            // Gets corresponding pixel in source samples, integer and fraction
            float xs  = ((float) xu) / wu * (ws-3) + 1.0;
            int   xsp = (int) xs;
            float xst = xs-xsp;
            float ys =  ((float) yu) / hu * (hs-3) + 1.0;
            int   ysp = (int) ys;
            float yst = ys-ysp;
            // Gets the interpolated value
            float pixel = bicubic2D(
                samp[ysp-1][xsp-1], samp[ysp-1][xsp  ], samp[ysp-1][xsp+1], samp[ysp-1][xsp+2],
                samp[ysp  ][xsp-1], samp[ysp  ][xsp  ], samp[ysp  ][xsp+1], samp[ysp  ][xsp+2],
                samp[ysp+1][xsp-1], samp[ysp+1][xsp  ], samp[ysp+1][xsp+1], samp[ysp+1][xsp+2],
                samp[ysp+2][xsp-1], samp[ysp+2][xsp  ], samp[ysp+2][xsp+1], samp[ysp+2][xsp+2],
                xst, yst);
            int l = (pixel < 0.0) ? 0 : ( (pixel > 255.0) ? 255 : ((int) pixel) );
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


