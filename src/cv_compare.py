"""ADIL kiyas: bizim yontem (uydurma YOK) vs B-V polinomu (CAPRAZ DOGRULANMIS).

Onceki kiyasta B-V polinomu ayni veriye uydurulup ayni veride test ediliyordu
(in-sample) -> haksiz avantaj. Burada 5-kat capraz dogrulama kullanilir.
"""
import os
import sys
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import response as resp
import spectra as spec
from photometry import photon_count_integral as PCI
from validate_johnson13 import load_johnson, sed_ratio, MATCH_ARCSEC

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    grid = resp.GRID_NM
    ri, rv = resp.system_response(grid), resp.v_band_response(grid)
    lib = spec.PicklesLibrary()
    vfn, _, _ = lib.match(20.0, "V")
    vl, vf = lib.spectrum(vfn)
    fveg = np.interp(grid, vl, vf, left=0.0, right=0.0)
    vratio = PCI(vl, vf, ri, grid) / PCI(vl, vf, rv, grid)

    jr = load_johnson()
    vc = next(r for r in jr if r["hr"].strip() == "7001")["cols"]
    ours = list(csv.DictReader(open(os.path.join(HERE, "out", "master_catalog.csv"))))
    ora = np.array([float(r["ra_deg"]) for r in ours])
    ode = np.array([float(r["dec_deg"]) for r in ours])

    for label, clean_only in (("TUM ESLESENLER", False), ("TEMIZ (degisken+cift elendi)", True)):
        D = []
        for j in jr:
            d = np.hypot((ora - j["ra"]) * np.cos(np.radians(j["de"])), ode - j["de"]) * 3600.0
            i = int(np.argmin(d))
            if d[i] > MATCH_ARCSEC:
                continue
            o = ours[i]
            if clean_only:
                if str(o.get("varflag", "")).strip() not in ("", "--"):
                    continue
                if str(o.get("multflag", "")).strip() not in ("", "--"):
                    continue
            if o["bv"] == "":
                continue
            fm = fveg * sed_ratio(j["cols"], vc, grid)
            Ii, Iv = PCI(grid, fm, ri, grid), PCI(grid, fm, rv, grid)
            if Ii <= 0 or Iv <= 0:
                continue
            D.append((float(o["bv"]), -2.5 * np.log10((Ii / Iv) / vratio), float(o["delta_v"])))

        bv, dt, do = (np.array(x) for x in zip(*D))
        print("=" * 62)
        print("  %s — %d yildiz" % (label, len(bv)))
        print("=" * 62)
        print("  BIZIM (sentetik, uydurma YOK)  std = %.4f  ort = %+.4f"
              % ((do - dt).std(), (do - dt).mean()))
        print()
        print("  B-V polinomu:")
        print("  %-12s %12s %12s" % ("derece", "5-kat CV", "(in-sample)"))
        rng = np.random.default_rng(0)
        folds = np.array_split(rng.permutation(len(bv)), 5)
        for deg in (1, 2, 4):
            err = []
            for k in range(5):
                te = folds[k]
                tr = np.concatenate([folds[m] for m in range(5) if m != k])
                c = np.polyfit(bv[tr], dt[tr], deg)
                err.append(np.polyval(c, bv[te]) - dt[te])
            e = np.concatenate(err)
            ins = np.polyval(np.polyfit(bv, dt, deg), bv) - dt
            print("  %-12d %12.4f %12.4f" % (deg, e.std(), ins.std()))
        print()


if __name__ == "__main__":
    main()
