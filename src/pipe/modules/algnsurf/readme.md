# algnsurf: align a still image to a reference image

this module warps an image so that it matches a reference image. it is
meant for burst photographs and exposure brackets taken by hand, for instance
in front of `hdrmerge`. it needs no manual input. the module detects features
in both images, matches them, estimates an affine matrix or a homography, and
resamples the image.

the feature detector and descriptor are a port of
[gpusurf](https://asrl.utias.utoronto.ca/code/gpusurf/) by paul furgale and chi
hay tong (SURF, speeded up robust features). the cuda kernels are rewritten as
vulkan compute shaders. the licence of gpusurf is in `LICENCE.gpusurf`.

the features do not depend on the exposure of the image. the module converts
the luminance to a logarithmic scale before it detects features, so a change
of exposure only adds a constant. the descriptors are also normalised. this
means you can align the frames of an exposure bracket.

the images for the feature detection are separate from the image that the
module warps. connect the same image to `input` and `alignsrc` to align it
directly. you can also detect features in a version of the image that is easy
to match, and warp a different version, for instance a raw image before
demosaic. the three images can have different sizes. the matrix is relative to
the image size.

if the module finds no valid transformation, it copies the input image to the
output. if `alignsrc` or `aligndst` is not connected, it copies the input image
as well.

## parameters

* `model` the transformation model: affine (6 degrees of freedom) or homography (8 degrees of freedom). use the homography if the camera rotated, use the affine matrix if you want the more stable result
* `maxfeat` the number of features per image. the module keeps the features with the strongest response
* `maxdim` the module reduces the image by an integer factor until its longer side is smaller than this value. the detection runs on the reduced image. a larger value is more accurate and slower
* `thresh` the minimum detector response in units of 1e-4. a larger value removes weak features in flat areas
* `ratio` a match is valid if the best distance is smaller than the second best distance times this ratio. a smaller value gives fewer and more reliable matches
* `rthresh` the maximum distance in pixels of the reduced image between the transformed reference feature and the target feature for the match to count as an inlier

## connectors

* `input` the image pixels that the module warps. the output has the size of this image
* `output` the warped input image, aligned to the reference image. the area where the input image has no data is black
* `alignsrc` the image of which the module detects features. the module aligns this image. it is usually a version of `input`, for instance with a different tone curve, or the same image
* `aligndst` the reference image of which the module detects features. it stays static

## technical

for each of the two images the module creates the following compute nodes.

* `iirow` and `iicol` reduce the image, compute the luminance on a log scale, and build an integral image. the integral image has unsigned 32 bit entries that wrap around. a difference of entries is exact as long as the true sum is smaller than 2^32. this is the reason why the box sums stay exact for large images
* `hessian` evaluates the box filter approximation of the hessian for 4 octaves and 4 intervals. it uses the filter shapes of gpusurf. the box edges are at real valued positions, and the code interpolates the integral image bilinearly, as the texture unit does in gpusurf
* `nms` finds the maxima in 3x3x3 neighbourhoods of the scale space and refines position and scale with a quadratic fit
* `select` ranks the candidates by response and keeps `maxfeat` of them. the rank does not depend on the order of the threads
* `describe` computes the orientation with haar wavelets in a sliding window of 60 degrees. then it computes the 64 dimensional descriptor from 4x4 sub regions with the sums of dx, |dx|, dy and |dy|

the remaining nodes work on both images.

* `match` finds the nearest and second nearest descriptor. it keeps a match if it passes the ratio test, has the same sign of the laplacian, and is the best match in both directions
* `solve` runs in one work group. each thread tests random minimal samples (RANSAC). the best hypothesis is refitted to its inliers with linear least squares, then with Gauss-Newton iterations that minimise the reprojection error. the matrix is converted to the pixel coordinates of the full resolution images. a result needs at least 8 inliers
* `align` samples the input image with bilinear filtering at the position the matrix gives for each output pixel

the matrix maps the pixel corner coordinates of the `aligndst` image to the pixel
corner coordinates of the `alignsrc` image.

## limitations

the module estimates one global transformation. it does not model parallax or
lens distortion. an object that moves in the scene breaks the model only if
most of the matches are on that object.

the module supports a feature buffer of at most 1024 features. the detector
keeps at most 16384 candidates before it selects the strongest features.

## tests

the script `generate_test_images.py` renders pairs of images with a known
homography. the script `check_alignment.py` runs the module on the pairs and
prints the error of the aligned image. it also prints the error of the image
that the true matrix gives. the module passes if both errors are about equal.
run the check from the root of the repository, where `bin/vkdt-cli` exists:

```
python3 src/pipe/modules/algnsurf/check_alignment.py
```
