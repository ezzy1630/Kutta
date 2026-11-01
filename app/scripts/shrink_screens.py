"""Shrink captured screens for the repository: downscale, then a 256-colour
palette from libimagequant with Floyd-Steinberg dithering.

Pillow's own median-cut palette shifts contour colours visibly, so it is not
used. Needs `pip install pillow imagequant` (in a virtual environment).

Usage:
  python3 -I app/scripts/shrink_screens.py [--width 1500] FILE.png... [--out DIR]

Without --out, files are rewritten in place.
"""

import argparse
from pathlib import Path

import imagequant
from PIL import Image


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", type=Path)
    ap.add_argument("--width", type=int, default=1500)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    if args.out:
        args.out.mkdir(parents=True, exist_ok=True)
    before = after = 0
    for f in args.files:
        before += f.stat().st_size
        img = Image.open(f).convert("RGBA")
        if img.width > args.width:
            img = img.resize((args.width, round(img.height * args.width / img.width)), Image.LANCZOS)
        out = (args.out / f.name) if args.out else f
        imagequant.quantize_pil_image(img, dithering_level=1.0, max_colors=256).save(out, optimize=True)
        after += out.stat().st_size
    print(f"{len(args.files)} files: {before / 1e6:.1f} MB -> {after / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
