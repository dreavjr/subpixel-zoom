/* Enlarges source image trying to estimate missing pixels - Quadratic Model
   (c) 2016, Eduardo Valle, blog.eduardovalle.com/
*/
#include <cassert>
#include <iostream>
#include <map>

// #define SOLVER_SPARSE
#include <Eigen/Core>
#ifdef SOLVER_SPARSE
    #include <Eigen/Sparse>
    typedef Eigen::Triplet<double> Triplet;
#else
    #include <Eigen/Dense>
#endif

#include <gd.h>


using namespace std;

/* Main function */
int main(int argc, char *argv[]) {

    // Check command-line and open source image
    const int nParamsMin = 2;
    const int nParamsMax = 8;
    if ( argc <= nParamsMin || argc>(nParamsMax+1) ) {
        cerr << "Usage: " << argv[0]
             << " <in.png> <out.png> [<dom>] [<asp>] [H up] [V up] [gweight] [toler]" << endl
             << endl
             << "   in.png    Input image in PNG format (mandatory)" << endl
             << "   out.png   Output image in PNG format (mandatory)" << endl
             << "   dom       Pixel connectedness domain: 4 or 8 (def: 4)" << endl
             << "   asp       Pixel aspect ratio width/height (def: 1.0)" << endl
             << "   H up      How much to upscale the image horizontally (def: 3)" << endl
             << "   V up      How much to upscale the image vertically (def: 3)" << endl
             << "   gweight   Relative importance of respecting data (def: 1.0 = neutral)" << endl
             << "   toler     Graylevel difference for  equal pixels (def: 0, max: 255)" << endl
             << endl
             << "   Optional parameters must be entered in the right sequence" << endl
             << endl
             #ifdef SOLVER_SPARSE
                 << "   Solver: Eigen::SparseQR" << endl
             #else
                 << "   Solver: Eigen::HouseholderQR" << endl
             #endif
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
    const float weightV    = 1.0;
    const float weightH    = pow(1.0/aspect, 2.0);
    const float weightD    = 1.0/(1.0 + aspect*aspect);

    // Horizontal, vertical, and surface ampliation
    const int   upX        = (argc <= 5) ? 3 : atoi(argv[5]);
    const int   upY        = (argc <= 6) ? 3 : atoi(argv[6]);
    if (upX < 1) {
        cerr << "Horizontal upscaling must be 1 or greater: " << argv[5] << endl;
        return 1;
    }
    if (upY < 1) {
        cerr << "Vertical upscaling must be 1 or greater: " << argv[6] << endl;
        return 1;
    }
    const int   up2D       = upX*upY;

    // Relative importance of attachment to graylevels
    const float grayWeight = up2D * ((argc <= 7) ? 1.0 : atof(argv[7]));
    if (grayWeight <= 0.0) {
        cerr << "Grayscale attachment weight must be > 0.0: " << argv[7] << endl;
        return 1;
    }

    // Tolerance for pixels be considered equal
    const float tolerance = (argc <= 8) ? 1.0 : atof(argv[8]);
    if (tolerance < 0.0) {
        cerr << "Pixel difference tolerance must be >= 0.0: " << argv[8] << endl;
        return 1;
    }

    // First pass : reads pixels matrix and count the number of needed variables
    const int ws = w+2;
    const int hs = h+2;
    float samples[hs][ws],
          lp = 0.0f;
    int   hotEstimate = 0;
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
            const float l = r*rl + g*gl + b*bl;
            // Destination pixels
            samples[ys][xs] = l;
            hotEstimate += (abs(l-lp)>tolerance) ? 1 : 0;
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

    cerr << "Hot pixels estimate: " << hotEstimate << endl;

    // Second pass : assembles list of the hot pixels :
    // those which are not equal to their neighbors
    map<pair<int, int>, int> hotMap;
    int hotCount = 0;
    for (int ys=1; ys<hs-1; ys++) {
        for (int xs=1; xs<ws-1; xs++) {
            if ((abs(samples[ys][xs] - samples[ys+1][xs+1]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys+1][xs  ]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys+1][xs-1]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys  ][xs+1]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys  ][xs-1]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys-1][xs+1]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys-1][xs  ]) > tolerance) ||
                (abs(samples[ys][xs] - samples[ys-1][xs-1]) > tolerance)) {
                hotMap.insert(make_pair(make_pair(xs, ys), hotCount));
                hotCount++;
            }
        }
    }

    cerr << "Hot pixels count: " << hotCount << endl;

    // Assembles optimization matrix...

    // The optimization we'll attempt tries to minimize the quadradic
    // difference between neighboring subsamples, while constraining all
    // subsamples inside a pixel to average to the given pixel. I.e.:
    // If subsamples si,1 to si,n belong to pi, we are trying to minimize
    // Sum_{for all neighboring subsample pairs s_1 s_2 } (s_1 - s_2)^2
    // constrained to : For all pixels p_i :
    // Avg_{for all subsamples s_i,j of pixel pi} si,j = pi

    // The optimization above leads to an overdetermined linear system:
    // 1) One equation per pixel for the averaging constraint
    // 2) One equation per subsample for the quadratic minimization that
    //    attains its mininimum when each subsample is the average of its
    //    neighbors

    const int vars    =  hotCount*up2D;
    const int eqns    =  hotCount*up2D                 + hotCount;
    const int nonzero = (hotCount*up2D)*(four ? 4 : 8) + (hotCount*up2D);
    cerr << "Optimizing : " << vars
             << " variables and " << eqns
             << " equations. " << nonzero
             << " non-zero entries / " << ( ((float)vars)*eqns )
             << " total (" << ( ((float) nonzero)/vars/eqns*100.0 )
             << "%)"  << endl;

    #ifdef SOLVER_SPARSE
        cerr << "Using sparse matrices" << endl;
        std::vector<Triplet> triplets(nonzero);
        Eigen::SparseMatrix<float> optimization(eqns, vars);
    #else
        cerr << "Using dense matrices" << endl;
        Eigen::MatrixXf optimization(eqns, vars);
        optimization.setZero();
    #endif
    Eigen::VectorXf targets(eqns);


    // ...variable index for subsample at position (u,v) of the pixel (x,y)
    auto varindex = [&](int i, int u, int v) -> int  {
        assert(i>=0 && i<hotCount && u>=0 && u<upX && v>=0 && v<upY);
        return i*up2D + v*upX+u;
    };
    auto findpix = [&](int x, int y) -> int {
        assert(x>=0 && x<ws && y>=0 && y<hs);
        try {
            return hotMap.at(make_pair(x, y));
        }
        catch(out_of_range e) {
            return -1;
        }
    };

    // ...assemble the graylevel constraints
    cerr << "Assembling graylevel constraints... " << endl;
    for (auto pixel : hotMap) {
        const int x = pixel.first.first;
        const int y = pixel.first.second;
        const int i = pixel.second;
        // The constraint is for the average of the subsamples to
        // be the graylevel of the pixel
        for (int v=0; v<upY; v++) {
            for (int u=0; u<upX; u++) {
                const int var   = varindex(i, u, v);
                const float val = (1.0f / up2D) * grayWeight;
                #ifdef SOLVER_SPARSE
                    triplets.push_back(Triplet(i, var, val));
                #else
                    optimization(i, var) = val;
                #endif
            }
        }
        targets(i) = samples[y][x] * grayWeight;
    }
    int eqn = hotCount;

    // ...assemble the quadratic minimizations
    cerr << "Assembling quadratic minimizers... " << endl;
    for (auto pixel : hotMap) {
        const int x = pixel.first.first;
        const int y = pixel.first.second;
        const int i = pixel.second;

        // Finds neighboring pixels
        const int iLU = four ? -1 : findpix(x-1, y-1);
        const int iL  =             findpix(x-1, y  );
        const int iLD = four ? -1 : findpix(x-1, y+1);
        const int iU  =             findpix(x  , y-1);
        const int iD  =             findpix(x  , y+1);
        const int iRU = four ? -1 : findpix(x+1, y-1);
        const int iR  =             findpix(x+1, y  );
        const int iRD = four ? -1 : findpix(x+1, y+1);

        // Finds variable corresponding to subsample at u+m,v+n
        // n, m = {-1, 0, 1}
        auto neighbor = [&](int u, int v, int m, int n) -> int {
            assert(u>=0 && u<upX && v>=0 && v<upY && m>=-1 && m<=1 && n>=-1 && n<=1);
            int nu = u+m, nv = v+n;
            // Safe cases : neighbor is in same pixel
            if (nu>=0 && nu<upX && nv>=0 && nv<upY) {
                return varindex(i, nu, nv);
            }
            // Problematic cases : neighbor is accross pixels
            int ni = -1;
            if (nu < 0) {
                if (nv <    0) { ni = iLU; nu = upX-1; nv = upY-1; } else
                if (nv >= upY) { ni = iLD; nu = upX-1; nv =     0; } else
                               { ni = iL ; nu = upX-1; }
            }
            else if (nu >= upX) {
                if (nv <    0) { ni = iRU; nu =     0; nv = upY-1; } else
                if (nv >= upY) { ni = iRD; nu =     0; nv =     0; } else
                               { ni = iR ; nu =     0; }
            }
            else {
                if (nv <    0) { ni = iU ; nv = upY-1; } else
                if (nv >= upY) { ni = iD ; nv =     0; } else
                               { assert(false); }
            }
            if (ni == -1) {
                return -2;
            }
            const int si = varindex(ni, nu, nv);
            assert(si != -1);
            return si;
        };

        // Create constraint equations for all subsamples
        for (int v=0; v<upY; v++) {
            for (int u=0; u<upX; u++) {
                // Finds neighboring subsamples
                const int sLU = four ? -3 : neighbor(u, v, -1, -1);
                const int sL  =             neighbor(u, v, -1,  0);
                const int sLD = four ? -3 : neighbor(u, v, -1,  1);
                const int sU  =             neighbor(u, v,  0, -1);
                const int sD  =             neighbor(u, v,  0,  1);
                const int sRU = four ? -3 : neighbor(u, v,  1, -1);
                const int sR  =             neighbor(u, v,  1,  0);
                const int sRD = four ? -3 : neighbor(u, v,  1,  1);
                // Counts neighbors
                const float norm =
                    ( (sLU >= 0) ? weightD : 0.0f ) +
                    ( (sL  >= 0) ? weightH : 0.0f ) +
                    ( (sLD >= 0) ? weightD : 0.0f ) +
                    ( (sU  >= 0) ? weightV : 0.0f ) +
                    ( (sD  >= 0) ? weightV : 0.0f ) +
                    ( (sRU >= 0) ? weightD : 0.0f ) +
                    ( (sR  >= 0) ? weightH : 0.0f ) +
                    ( (sRD >= 0) ? weightD : 0.0f );
                // The constraint is for each subsample to be the average
                // of the neighboring subsamples
                #ifdef SOLVER_SPARSE
                    triplets.push_back( Triplet(eqn, varindex(i, u, v), -1.0f) );
                    if (sLU >= 0) { triplets.push_back( Triplet(eqn, sLU, weightD / norm) ); }
                    if (sL  >= 0) { triplets.push_back( Triplet(eqn, sL , weightH / norm) ); }
                    if (sLD >= 0) { triplets.push_back( Triplet(eqn, sLD, weightD / norm) ); }
                    if (sU  >= 0) { triplets.push_back( Triplet(eqn, sU , weightV / norm) ); }
                    if (sD  >= 0) { triplets.push_back( Triplet(eqn, sD , weightV / norm) ); }
                    if (sRU >= 0) { triplets.push_back( Triplet(eqn, sRU, weightD / norm) ); }
                    if (sR  >= 0) { triplets.push_back( Triplet(eqn, sR , weightH / norm) ); }
                    if (sRD >= 0) { triplets.push_back( Triplet(eqn, sRD, weightD / norm) ); }
                #else
                    optimization(eqn, varindex(i, u, v)) = -1.0f;
                    if (sLU >= 0) { optimization(eqn, sLU) = weightD / norm; }
                    if (sL  >= 0) { optimization(eqn, sL ) = weightH / norm; }
                    if (sLD >= 0) { optimization(eqn, sLD) = weightD / norm; }
                    if (sU  >= 0) { optimization(eqn, sU ) = weightV / norm; }
                    if (sD  >= 0) { optimization(eqn, sD ) = weightV / norm; }
                    if (sRU >= 0) { optimization(eqn, sRU) = weightD / norm; }
                    if (sR  >= 0) { optimization(eqn, sR ) = weightH / norm; }
                    if (sRD >= 0) { optimization(eqn, sRD) = weightD / norm; }
                #endif
                targets(eqn) = 0.0f;
                eqn++;
            }
        }
    }

    // ...showtime ! --- solve the humungous system
    cerr << "Showtime !" << endl;
    Eigen::VectorXf subsamples(vars);
    #ifdef SOLVER_SPARSE
        optimization.setFromTriplets(triplets.begin(), triplets.end());
        optimization.makeCompressed(); // Requirement of SparseQR
        Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int> > solver;
        solver.compute(optimization);
        if (solver.info() != Eigen::Success) {
            cerr << "Matrix decomposition failed" << endl;
            return 1;
        }
        subsamples = solver.solve(targets);
        if(solver.info() != Eigen::Success) {
            cerr << "Solver failed" << endl;
            return 1;
        }
    #else
        Eigen::HouseholderQR<Eigen::MatrixXf> solver;
        solver.compute(optimization);
        subsamples = solver.solve(targets);
    #endif
    // DEBUG :
    // Eigen::JacobiSVD<Eigen::MatrixXf>
    //    solver(optimization, Eigen::ComputeThinU | Eigen::ComputeThinV); // +precise, but ++slower
    // cout << optimization << endl;
    // cout << targets << endl;


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

    // DEBUG :
    // int gray2[256];
    // for (int l=0; l<256; l++) {
    //     gray2[l] = gdImageColorAllocate(dimage, l, (int) (0.8f*l), l);
    // }

    // Second pass : resample image
    for (int y=0, yu=0; y<h; y++, yu+=upY) {
        for (int x=0, xu=0; x<w; x++, xu+=upX) {
            int i = findpix(x+1, y+1);
            for (int v=0; v<upY; v++) {
                for (int u=0; u<upX; u++) {
                    int l;
                    if (i >= 0) {
                        // Hot pixels were estimated by the optimizer
                        l = subsamples( varindex(i, u, v) );
                    }
                    else {
                        // Cold pixels come from flat areas and are copied from the samples
                        l = samples[y+1][x+1];
                    }
                    l = (l < 0) ? 0 : ( (l>255) ? 255 : l );
                    gdImageSetPixel(dimage, xu+u, yu+v, gray[l]);
                    // gdImageSetPixel(dimage, xu+u, yu+v, i>=0 ? gray[l] : gray2[l]); // DEBUG
                }
            }
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
