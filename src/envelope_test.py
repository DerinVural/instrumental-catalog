"""N8 — Renk terimi ZARF ve EGIM testi.

Soru: katalogdaki uc `delta_v` degerleri (min ~ -3.8) FIZIKSEL mi?

Yontem: bandimizin etkin dalgaboyu V ile I_C arasinda. Dolayisiyla her tip icin
renk terimi C, ayni spektrumdan olculen (V - I_C) rengiyle SINIRLI olmali:

    oran = C / (V - I_C)

Oran maviden kirmiziya MONOTON yukselmeli (band etkin dalgaboyu kirmizi
yildizlarda I_C'ye yaklasir) ama 1.0'i asmamali — asarsa zincirde sorun var.

Renkler Pickles Tablo 2'den transkribe EDILMEZ; ayni UVILIB spektrumlarindan
Bessell (1990) Cousins R_C/I_C egrileriyle HESAPLANIR -> kutuphane surumu
tutarsizligi (N6) bu teste sizmaz.
"""
import os
import sys
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import response as resp
import spectra as spec
from photometry import photon_count_integral as PCI

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Bessell (1990) Cousins R_C ve I_C bant tepkileri (yaklasik, 10 nm adim)
RC = [(550, .00), (560, .23), (570, .74), (580, .91), (590, .98), (600, 1.0),
      (610, .98), (620, .96), (630, .93), (640, .90), (650, .86), (660, .83),
      (670, .79), (680, .76), (690, .72), (700, .68), (710, .64), (720, .60),
      (730, .56), (740, .52), (750, .48), (760, .44), (770, .40), (780, .36),
      (790, .32), (800, .28), (810, .24), (820, .20), (830, .16), (840, .12),
      (850, .09), (860, .06), (870, .03), (880, .01), (890, .00)]

IC = [(700, .00), (710, .02), (720, .10), (730, .29), (740, .57), (750, .79),
      (760, .92), (770, .98), (780, 1.0), (790, 1.0), (800, .99), (810, .98),
      (820, .96), (830, .94), (840, .91), (850, .87), (860, .83), (870, .79),
      (880, .74), (890, .69), (900, .63), (910, .57), (920, .50), (930, .43),
      (940, .36), (950, .29), (960, .22), (970, .15), (980, .09), (990, .04),
      (1000, .01), (1010, .00)]


def band(pts, grid):
    lam = np.array([p[0] for p in pts], float)
    r = np.array([p[1] for p in pts], float)
    return np.clip(np.interp(grid, lam, r, left=0.0, right=0.0), 0.0, None)


def main():
    grid = resp.GRID_NM
    r_inst = resp.system_response(grid)
    r_v = resp.v_band_response(grid)
    r_rc, r_ic = band(RC, grid), band(IC, grid)

    lib = spec.PicklesLibrary()
    vfn, _, _ = lib.match(20.0, "V")
    vl, vf = lib.spectrum(vfn)

    def ratios(lam, flx):
        return (PCI(lam, flx, r_inst, grid), PCI(lam, flx, r_v, grid),
                PCI(lam, flx, r_rc, grid), PCI(lam, flx, r_ic, grid))

    Iiv, Ivv, Ircv, Iicv = ratios(vl, vf)
    ref_inst = Iiv / Ivv          # Vega: C  normalizasyonu
    ref_rc = Ivv / Ircv           # Vega: V-Rc normalizasyonu
    ref_ic = Ivv / Iicv           # Vega: V-Ic normalizasyonu

    rows = []
    for fn, sp, idx, lum in lib.entries:
        lam, flx = lib.spectrum(fn)
        Ii, Iv, Irc, Iic = ratios(lam, flx)
        if min(Ii, Iv, Irc, Iic) <= 0:
            continue
        C = -2.5 * np.log10((Ii / Iv) / ref_inst)
        v_rc = -2.5 * np.log10((Iv / Irc) / ref_rc)
        v_ic = -2.5 * np.log10((Iv / Iic) / ref_ic)
        rows.append(dict(sp=sp, idx=idx, C=C, v_rc=v_rc, v_ic=v_ic,
                         ratio=(-C / v_ic if abs(v_ic) > 0.05 else np.nan)))   # isaret: incelemeyle ayni konvansiyon (pozitif)

    rows.sort(key=lambda r: r["v_ic"])          # maviden kirmiziya

    print("=" * 74)
    print("  N8 — ZARF ve EGIM TESTI (%d Pickles tipi)" % len(rows))
    print("=" * 74)
    print("  %-10s %8s %8s %8s %8s" % ("tip", "V-Rc", "V-Ic", "C", "C/(V-Ic)"))
    show = [r for r in rows if abs(r["v_ic"]) > 0.05]
    for r in show[::max(1, len(show) // 14)]:
        print("  %-10s %8.3f %8.3f %8.3f %8.2f"
              % (r["sp"], r["v_rc"], r["v_ic"], r["C"], r["ratio"]))

    ok = True
    valid = [r for r in rows if np.isfinite(r["ratio"])]

    # 1) ZARF: C ile (V-Ic) arasinda kalmali; oran 1.0'i asmamali
    over = [r for r in valid if r["ratio"] > 0.75]
    print()
    print("  [%s] zarf: oran > 0.75 olan tip yok  (bulunan: %d)"
          % ("GECTI" if not over else "KALDI", len(over)))
    if over:
        ok = False
        for r in over[:5]:
            print("        %-10s V-Ic=%.3f C=%.3f oran=%.2f" % (r["sp"], r["v_ic"], r["C"], r["ratio"]))

    # 2) EGIM: mavi uc dusuk, kirmizi uc yuksek (monoton egilim)
    # V-Ic kirmizilastikca BUYUR (standart konvansiyon). Mavi yildiz: V-Ic <= 0.
    blue = [r["ratio"] for r in valid if r["v_ic"] < 0.5]
    red = [r["ratio"] for r in valid if r["v_ic"] > 2.0]
    if blue and red:
        mb, mr = float(np.median(blue)), float(np.median(red))
        cond = mr > mb
        ok &= cond
        print("  [%s] egim: kirmizi uc orani > mavi uc  (mavi %.2f -> kirmizi %.2f)"
              % ("GECTI" if cond else "KALDI", mb, mr))

    # 3) KATALOG UC DEGERI bu zarfla uyumlu mu?
    mpath = os.path.join(HERE, "out", "master_catalog.csv")
    if os.path.exists(mpath):
        cat = list(csv.DictReader(open(mpath)))
        worst = min(cat, key=lambda r: float(r["delta_v"]))
        mt = worst["matched_sptype"]
        ent = next((r for r in rows if r["sp"] == mt), None)
        print()
        print("  --- katalogun EN UC yildizi ---")
        print("  HIP %s  SpType=%s  ->  Pickles %s" % (worst["hip"], worst["sptype"][:14], mt))
        print("  delta_v = %.3f" % float(worst["delta_v"]))
        if ent:
            print("  o tipin: V-Ic = %.3f  C = %.3f  oran = %.2f"
                  % (ent["v_ic"], ent["C"], ent["ratio"]))
            consistent = abs(float(worst["delta_v"]) - ent["C"]) < 0.01
            ok &= consistent
            print("  [%s] katalog degeri, tipin C'siyle TUTARLI (fark %.4f)"
                  % ("GECTI" if consistent else "KALDI",
                     abs(float(worst["delta_v"]) - ent["C"])))
        else:
            print("  (tip zarf tablosunda bulunamadi)")

    print()
    print("  SONUC: %s" % ("ZARF TESTI GECTI" if ok else "TETIK — zincirde inceleme gerekli"))
    print("=" * 74)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
