# ERCOFTAC T3 flat-plate transition (Coupland, Rolls-Royce ASL)

Zero-pressure-gradient flat-plate transitional boundary layers T3A, T3B and
T3A-, measured by J. Coupland at the Rolls-Royce Applied Science Laboratory
for the ERCOFTAC Special Interest Group on Transition (Savill, 1993).

- Page: ERCOFTAC Classic Collection, case 020 "Flat Plate Transitional
  Boundary Layers", http://cfd.mace.manchester.ac.uk/ercoftac/doku.php?id=cases:case020
- Files: `http://cfd.mace.manchester.ac.uk/ercoftac/lib/exe/fetch.php?media=cdata:case020:<name>`
  for each name below. Downloaded 2026-10-31.
- Licence: the wiki footer states "Except where otherwise noted, content on
  this wiki is licensed under ... CC Attribution-Noncommercial-Share Alike 4.0
  International". No other terms are given for the data files. Kutta is
  intended for commercial use, so the files are **not redistributed** here
  (gitignored, local V&V only). Reports quote only the handful of values a
  comparison needs.

| file | content | sha256 |
| --- | --- | --- |
| raw/readme.txt | file conventions | 0a0dbd3a869cf0e9667191b347a438240833a496c88047dcda7da7cfa8809ed9 |
| raw/t3ay.dat | T3A summary: x, Re_x, U0, Cf, H, Re_theta, delta, freestream Tu | 61b15b27511fb940409bab20b67c4776c170f2a1679e4a369a6536ab9043e8cc |
| raw/t3by.dat | T3B summary | 95e04aca4510ead26c80ede79182cbd8e97bda9f6f30dd8ca6ce8e1b1c21a394 |
| raw/t3amy.dat | T3A- summary | 65aca3ecfc9349abe8b727d1ae15c38d62289e12c0691e78d9629998fcbd90a9 |
| raw/t3as.dat | T3A single-wire profiles | 225095350fbb04cfb9c87782a3e8dacbf3ff76eed02d6a2661c5b0515420ccac |
| raw/t3av.dat | T3A U-V crossed wire | 749b02dc7677616b8a3ae970291d1f8eb4af4edafe339c5bc2e59e5a9eca837d |
| raw/t3aw.dat | T3A U-W crossed wire | 1920f8ce4fa17a04daae8d3c64df92a8673a993c51b7ef801d0cd6edb1bfdcda |
| raw/t3bs.dat | T3B single-wire profiles | 1bbc9c2fac2800877a15662c964d4e15aa50a47b8e38d1a55b74256d3ee7401d |
| raw/t3bv.dat | T3B U-V crossed wire | 8cc74b6411851ea6db738a2a1bd499b9abad4bd2b82deeb324adb9334c6e68c9 |
| raw/t3bw.dat | T3B U-W crossed wire | 08664b84e16d2d0418b0bf9eb23448853e59eeacc5a9f4d1478cb1c810dc7bf2 |
| raw/t3ams.dat | T3A- single-wire profiles | 87e6e29f0be773c39e43e54a4f8de5ff65e9e2807ac68b4968497eb9bb817404 |
| raw/t3amv.dat | T3A- U-V crossed wire | c75b49461feafbd079cc1d50087692d005d6aff1905b97334b70a0149eefc316 |
| raw/t3amw.dat | T3A- U-W crossed wire | ea2c2bbb5341c27e443f96d936b7c4473d64e7656ba0f5307b7a922148c187bd |

The summary (`*y.dat`) Cf is the experimenters' wall-shear estimate: a
Blasius-type fit in laminar profiles and the log law (kappa 0.41, B 5.2) or
Preston tube in turbulent ones (case page, "Experimental Methods"). Their
`TUFS` column is the measured freestream turbulence intensity in percent at
each station, used to set the inlet decay.

To re-fetch, run from this directory:

```sh
mkdir -p raw && cd raw
for f in readme.txt t3ay.dat t3by.dat t3amy.dat t3as.dat t3av.dat t3aw.dat \
         t3bs.dat t3bv.dat t3bw.dat t3ams.dat t3amv.dat t3amw.dat; do
  curl -sL -o "$f" "http://cfd.mace.manchester.ac.uk/ercoftac/lib/exe/fetch.php?media=cdata:case020:$f"
done
shasum -a 256 *
```
