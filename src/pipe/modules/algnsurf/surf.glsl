// Common definitions for the SURF feature pipeline.
// This code is a port of gpusurf by Paul Furgale and Chi Hay Tong.
// The licence of gpusurf is in the file LICENCE.gpusurf.
#ifndef IMGALIGN_SURF_GLSL
#define IMGALIGN_SURF_GLSL

#define SURF_OCTAVES   4
#define SURF_INTERVALS 4

// Maximum number of keypoint candidates before the selection of the strongest.
#define CAND_MAX    16384
// One keypoint record: x, y, size, response, angle, octave, laplacian sign, unused.
#define KP_STRIDE   8
#define DESC_DIM    64
// One feature record: the keypoint record, followed by the descriptor.
#define FEAT_STRIDE (KP_STRIDE + DESC_DIM)

// Position of the values in a keypoint record.
#define KP_X      0
#define KP_Y      1
#define KP_SIZE   2
#define KP_RESP   3
#define KP_ANGLE  4
#define KP_OCTAVE 5
#define KP_LAP    6

// Parameters of the box filters. These are the gpusurf default values.
const float L1 = 3.0/1.5;
const float L2 = 5.0/1.5;
const float L3 = 3.0/1.5;
const float L4 = 1.0/1.5;
const float EDGE_SCALE = 0.81;

// Module parameters. The layout must match shader_params_t in main.c.
layout(std140, set = 0, binding = 1) uniform params_t
{
  float thresh;   // detector threshold for the determinant of the hessian
  float ratio;    // maximum distance ratio of best and second best match
  float rthresh;  // maximum reprojection error of an inlier in pixels
  int   maxfeat;  // maximum number of features per image
  int   model;    // transformation model: 0 affine, 1 homography
} params;

// The host computes the layout of the octaves and sends it as push constants.
// A shader can define SURF_NO_PUSH to declare its own push constants.
#ifndef SURF_NO_PUSH
layout(push_constant, std140) uniform push_t
{
  ivec4 size[SURF_OCTAVES];   // x size, y size, border, step of each octave
  vec4  scale[SURF_OCTAVES];  // filter scale of each interval of each octave
  ivec4 offset;               // start of each octave in the hessian buffer
  ivec4 dim;                  // integral image width, height, reduction factor
} push;
#endif

#endif
