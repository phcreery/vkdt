#!/usr/bin/env python3
"""
Run the algnsurf module on the generated test pairs and measure the result.
Run this script where bin/vkdt-cli works, from the root of the repository.
The script compares three images with the reference image:
the unaligned target, the aligned target, and the target warped with the true matrix.
The error is the mean absolute difference to the reference, relative to its mean.
The comparison uses only the pixels where the target has data.
The exposure gain of the target is removed before the comparison.
"""

import argparse
import json
import os
import subprocess
import sys

import numpy as np

import generate_test_images as gen

CFG = """frames:1
fps:0
module:i-pfm:01:0:390
module:i-pfm:02:0:590
module:algnsurf:01:450:390.6
module:o-pfm:main:650:390.6
connect:i-pfm:01:output:algnsurf:01:aligndst
connect:i-pfm:02:output:algnsurf:01:alignsrc
connect:i-pfm:02:output:algnsurf:01:input
connect:algnsurf:01:output:o-pfm:main:input
param:i-pfm:01:filename:{ref}
param:i-pfm:01:startid:0
param:i-pfm:01:noise a:0
param:i-pfm:01:noise b:0
param:i-pfm:02:filename:{tgt}
param:i-pfm:02:startid:0
param:i-pfm:02:noise a:0
param:i-pfm:02:noise b:0
{params}
param:o-pfm:main:filename:output
"""

MARGIN = 8  # ignore pixels close to the border of the target image


def run_pair(out, name, params):
    cfg = os.path.join(out, name + ".cfg")
    with open(cfg, "w") as f:
        f.write(CFG.format(ref=os.path.join(out, name + "_ref.pfm"),
                           tgt=os.path.join(out, name + "_tgt.pfm"),
                           params="\n".join("param:algnsurf:01:%s" % p for p in params)))
    res = os.path.join(out, name + "_aligned")
    cmd = ["bin/vkdt-cli", "-g", cfg, "--last-frame-only", "--format", "o-pfm", "--filename", res]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(res + ".pfm"):
        print(r.stdout, r.stderr)
        raise RuntimeError("vkdt-cli failed for " + name)
    return gen.load_pfm(res + ".pfm")


def rel_error(img, ref, mask):
    return float(np.mean(np.abs(img[mask] - ref[mask])) / np.mean(ref[mask]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/tmp/algnsurf_test")
    ap.add_argument("--pairs", nargs="*", default=list(gen.PAIRS))
    ap.add_argument("--param", action="append", default=[], help="algnsurf parameter, name:value")
    ap.add_argument("--no-generate", action="store_true")
    args = ap.parse_args()
    if not args.no_generate or not os.path.exists(os.path.join(args.out, "truth.json")):
        gen.main(args.out)
    with open(os.path.join(args.out, "truth.json")) as f:
        truth = json.load(f)

    failed = 0
    for name in args.pairs:
        ref = gen.load_pfm(os.path.join(args.out, name + "_ref.pfm"))
        tgt = gen.load_pfm(os.path.join(args.out, name + "_tgt.pfm"))
        gain = truth[name]["gain"]
        h = np.array(truth[name]["h"])
        aligned = run_pair(args.out, name, args.param)
        # the target warped with the true matrix: sample the target at h * p
        ys, xs = np.mgrid[0:gen.HEIGHT, 0:gen.WIDTH].astype(np.float64)
        p = np.stack([xs + 0.5, ys + 0.5, np.ones_like(xs)], axis=-1) @ h.T
        ideal = gen.sample(tgt, p[..., 0] / p[..., 2] - 0.5, p[..., 1] / p[..., 2] - 0.5)
        # only compare pixels where the target has data
        qx, qy = p[..., 0] / p[..., 2], p[..., 1] / p[..., 2]
        mask = (qx > MARGIN) & (qx < gen.WIDTH - MARGIN) & (qy > MARGIN) & (qy < gen.HEIGHT - MARGIN)
        e_none  = rel_error(tgt / gain, ref, mask)
        e_align = rel_error(aligned / gain, ref, mask)
        e_ideal = rel_error(ideal / gain, ref, mask)
        ok = e_align < 0.5 * e_none and e_align < e_ideal * 1.5 + 0.005
        failed += 0 if ok else 1
        print("%-6s unaligned %.4f  aligned %.4f  true matrix %.4f  %s" %
              (name, e_none, e_align, e_ideal, "ok" if ok else "FAIL"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
