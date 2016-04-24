# subpixel-zoom

This is a collection of tools for enlarging images of small fonts rendered on screen. The intended target application was OCR for screen grabbers. 

Currently, the only practical tool is subpixel-zoom, which extracts the subpixels and then performs a bicubic interpolation. All others are (very!) experimental, and currently not useful, but they are fun to tinker with. I am releasing them in the wild in the hope they might be helpful as a starting point for anyone interested in OCR, binarization, font rendering, or document restoration.

The projects also illustrate:
* How to use libgd to read, manipulate and write PNG files (all programs);
* How to use Eigen to solve large linear systems, using both dense and sparse matrices (pixel-estimate-quadratic);
* How to run a simple Ising model with iterated conditional models using just C++ (pixel-estimated-ising*).

# Dependencies

All three programs depends on libgd: http://libgd.github.io/ 

If you want to compile pixel-estimate-quadratic, you additionally need Eigen 3: http://eigen.tuxfamily.org/

# Building

Once you have the dependencies installed, compilation is a breeze:

```
g++ src/subpixel-zoom.cpp -o subpixel-zoom -std=c++11 -lgd -Wall -O3
g++ src/pixel-estimate-quadratic.cpp -o pixel-estimate-quadratic -std=c++11 -I/usr/local/include/eigen3 -lgd -Wall -O3 -ffast-math -fassociative-math -DSOLVER_SPARSE
g++ src/pixel-estimate-ising.cpp -o pixel-estimate-ising -std=c++11 -lgd -Wall -O3 -ffast-math -fassociative-math
g++ src/pixel-estimate-ising-smoothed.cpp -o pixel-estimate-ising-smoothed -std=c++11 -lgd -Wall -O3 -ffast-math -fassociative-math
```

Remember to change /usr/local/include/eigen3 for the actual path where you installed Eigen's header files.

I tested it on OS X/HomeBrew + GCC, but it should work on other unices (including Linux) without difficulty.

# Using

Running the programs without options show a short help screen. Most options are self-explicative. In doubt, you can use the default values.

# Credits

I based the initial version of subpixel-zoom from this blog post: http://ruletheweb.co.uk/blog/2014/02/subpixel-aware-image-scaling/  While the original is concerned with binarization, I wanted to preserve grayscales.

In my first attempt to understand Bicubic Interpolation in practice, this blog post helped a lot: http://blog.demofox.org/2015/08/15/resizing-images-with-bicubic-interpolation/ 

I based the final implementation on the description of Wikipedia for a convolution based algorithm: https://en.wikipedia.org/w/index.php?title=Bicubic_interpolation&oldid=715871752#Bicubic_convolution_algorithm

# Licence and Limitations

I am licensing this under a permissive Simplified FreeBSD license. I attempted this mostly for fun. No guarantees whatsoever. Should you decide to use the software, for any purpose, you are assuming all risks, including legal risks, and releasing me from all liabilities. Please check LICENSE.md for the legal text

If you need to credit me, you can point to this repo at github.com/dreavjr/subpixel-zoom, or to my blog at blog.eduardovalle.com/.
