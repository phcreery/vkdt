#!/usr/bin/env python3
"""
Generate test image pairs for the algnsurf module.
Each pair has a reference image and a target image.
The target shows the same scene with a known homography.
The homography maps pixel corner coordinates of the reference to the target.
The script writes PFM files and a JSON file with the matrices.
The first argument is the output directory. The default is /tmp/algnsurf_test.
"""

import json
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = "/tmp/algnsurf_test"
WIDTH, HEIGHT = 1024, 768
MARGIN = 200  # extra scene around the reference image


def smooth_noise(rng, w, h, cell):
    """Return value noise with the given cell size and cosine interpolation."""
    gw, gh = w // cell + 3, h // cell + 3
    grid = rng.standard_normal((gh, gw))
    ys, xs = np.mgrid[0:h, 0:w].astype(np.float64)
    gx, gy = xs / cell, ys / cell
    ix, iy = gx.astype(int), gy.astype(int)
    fx, fy = gx - ix, gy - iy
    fx, fy = 0.5 - 0.5 * np.cos(np.pi * fx), 0.5 - 0.5 * np.cos(np.pi * fy)
    a = grid[iy, ix] * (1 - fx) + grid[iy, ix + 1] * fx
    b = grid[iy + 1, ix] * (1 - fx) + grid[iy + 1, ix + 1] * fx
    return a * (1 - fy) + b * fy


def make_scene(w, h, seed):
    """Return a linear RGB scene with structure at many scales."""
    rng = np.random.default_rng(seed)
    lum = np.zeros((h, w))
    for cell, amp in [(160, 1.0), (80, 0.8), (40, 0.6), (20, 0.4), (10, 0.3), (5, 0.2)]:
        lum += amp * smooth_noise(rng, w, h, cell)
    ys, xs = np.mgrid[0:h, 0:w]
    for _ in range(160):  # rectangles and discs with sharp edges
        cx, cy = rng.uniform(0, w), rng.uniform(0, h)
        size = rng.uniform(6, 70)
        delta = rng.uniform(-1.5, 1.5)
        if rng.random() < 0.5:
            mask = (abs(xs - cx) < size) & (abs(ys - cy) < size * rng.uniform(0.4, 1.5))
        else:
            mask = (xs - cx) ** 2 + (ys - cy) ** 2 < size ** 2
        lum[mask] += delta
    lum = (lum - lum.min()) / (lum.max() - lum.min())
    lum = 0.02 + 0.9 * lum ** 1.5
    rgb = np.stack([lum * (1.0 + 0.15 * smooth_noise(rng, w, h, 200)) for _ in range(3)], axis=-1)
    rgb[..., 1] *= 1.0 + 0.1 * smooth_noise(rng, w, h, 150)
    return np.clip(rgb, 0.0, None)


def sample(img, x, y):
    """Bilinear sampling at array index coordinates."""
    h, w = img.shape[:2]
    x = np.clip(x, 0, w - 1.001)
    y = np.clip(y, 0, h - 1.001)
    ix, iy = x.astype(int), y.astype(int)
    fx, fy = (x - ix)[..., None], (y - iy)[..., None]
    a = img[iy, ix] * (1 - fx) + img[iy, ix + 1] * fx
    b = img[iy + 1, ix] * (1 - fx) + img[iy + 1, ix + 1] * fx
    return a * (1 - fy) + b * fy


def homography(angle_deg, scale, shift, persp):
    """Build a homography in pixel corner coordinates around the image centre."""
    a = np.radians(angle_deg)
    c = np.array([WIDTH / 2.0, HEIGHT / 2.0])
    lin = np.array([[scale * np.cos(a), -scale * np.sin(a), 0.0],
                    [scale * np.sin(a),  scale * np.cos(a), 0.0],
                    [persp[0], persp[1], 1.0]])
    t_in = np.array([[1, 0, -c[0]], [0, 1, -c[1]], [0, 0, 1.0]])
    t_out = np.array([[1, 0, c[0] + shift[0]], [0, 1, c[1] + shift[1]], [0, 0, 1.0]])
    h = t_out @ lin @ t_in
    return h / h[2, 2]


def render(scene, h, gain, noise, rng):
    """Render the scene seen through the homography h (reference to target)."""
    ys, xs = np.mgrid[0:HEIGHT, 0:WIDTH].astype(np.float64)
    q = np.stack([xs + 0.5, ys + 0.5, np.ones_like(xs)], axis=-1)
    p = q @ np.linalg.inv(h).T
    px, py = p[..., 0] / p[..., 2] - 0.5, p[..., 1] / p[..., 2] - 0.5
    img = sample(scene, px + MARGIN, py + MARGIN) * gain
    if noise > 0.0:
        img = img + rng.standard_normal(img.shape) * noise * np.sqrt(np.maximum(img, 1e-4))
    return np.clip(img, 0.0, None)


def save_pfm(path, img):
    with open(path, "wb") as f:
        f.write(b"PF\n%d %d\n-1.0\n" % (img.shape[1], img.shape[0]))
        f.write(np.flipud(img).astype("<f4").tobytes())


def load_pfm(path):
    with open(path, "rb") as f:
        assert f.readline().strip() == b"PF"
        w, h = map(int, f.readline().split())
        f.readline()
        data = np.frombuffer(f.read(w * h * 12), dtype="<f4").reshape(h, w, 3)
    return np.flipud(data)


# name: (angle in degrees, scale, shift in pixels, perspective terms, gain of the target, noise)
PAIRS = {
    "burst": (0.8, 1.000, (11.3, -7.6), (0.0, 0.0), 1.0, 0.01),
    "hdr":   (2.5, 1.020, (24.5, -17.2), (2.0e-5, -1.5e-5), 0.25, 0.01),
    "large": (7.0, 0.960, (-45.0, 30.0), (4.0e-5, 3.0e-5), 0.5, 0.01),
    "rotated": (30.0, 1.100, (-60.0, 40.0), (0.0, 0.0), 0.7, 0.01),
}


def main(out=DEFAULT_OUT):
    os.makedirs(out, exist_ok=True)
    scene = make_scene(WIDTH + 2 * MARGIN, HEIGHT + 2 * MARGIN, 7)
    truth = {}
    for name, (angle, scale, shift, persp, gain, noise) in PAIRS.items():
        rng = np.random.default_rng(sum(map(ord, name)))
        h = homography(angle, scale, shift, persp)
        save_pfm(os.path.join(out, name + "_ref.pfm"), render(scene, np.eye(3), 1.0, noise, rng))
        save_pfm(os.path.join(out, name + "_tgt.pfm"), render(scene, h, gain, noise, rng))
        truth[name] = {"h": h.tolist(), "gain": gain}
        print("wrote", name)
    with open(os.path.join(out, "truth.json"), "w") as f:
        json.dump(truth, f, indent=2)


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:2]))
