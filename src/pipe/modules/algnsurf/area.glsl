// Box sums on the integral image.
// The integral image has one extra row and one extra column of zero values.
// Each entry is the sum of all quantised pixels above and left of the entry.
// The sums use unsigned integers and wrap around at 2^32.
// A difference of entries is exact when the true result is smaller than 2^32.
#ifndef IMGALIGN_AREA_GLSL
#define IMGALIGN_AREA_GLSL

#include "surf.glsl"

// Pixel values in the integral image have this maximum.
const float II_QUANT = 4095.0;

layout(set = 1, binding = 0, std430) readonly buffer buf_ii { uint ii[]; };

struct corner_t
{
  uint  s00, s10, s01, s11; // the four integral image entries around the corner
  float fx, fy;             // position inside the pixel
};

corner_t ii_corner(float x, float y)
{
  const int w = push.dim.x;
  const int h = push.dim.y;
  x = clamp(x, 0.0, float(w));
  y = clamp(y, 0.0, float(h));
  const int ix = min(int(floor(x)), w-1);
  const int iy = min(int(floor(y)), h-1);
  const int stride = w + 1;
  const int i0 = iy * stride + ix;
  corner_t c;
  c.s00 = ii[i0];
  c.s10 = ii[i0 + 1];
  c.s01 = ii[i0 + stride];
  c.s11 = ii[i0 + stride + 1];
  c.fx = x - float(ix);
  c.fy = y - float(iy);
  return c;
}

float ii_pixel(corner_t c)
{
  return float(int(c.s11 - c.s10 - c.s01 + c.s00));
}

// Return the sum of the pixels in a box with real valued corners.
// The sum is exact for a piecewise constant image.
// The code splits the bilinear interpolation of the four corners into groups.
// Each group is a difference of integral image entries.
// This keeps all float values small.
float box_area(float x0, float y0, float x1, float y1)
{
  const corner_t a = ii_corner(x0, y0);
  const corner_t b = ii_corner(x1, y0);
  const corner_t c = ii_corner(x0, y1);
  const corner_t d = ii_corner(x1, y1);
  float r = float(int(a.s00 - b.s00 - c.s00 + d.s00));
  r += a.fx * float(int((a.s10 - a.s00) - (c.s10 - c.s00)));
  r -= b.fx * float(int((b.s10 - b.s00) - (d.s10 - d.s00)));
  r += a.fy * float(int((a.s01 - a.s00) - (b.s01 - b.s00)));
  r -= c.fy * float(int((c.s01 - c.s00) - (d.s01 - d.s00)));
  r += a.fx*a.fy*ii_pixel(a) - b.fx*b.fy*ii_pixel(b)
     - c.fx*c.fy*ii_pixel(c) + d.fx*d.fy*ii_pixel(d);
  return r / II_QUANT;
}

#endif
