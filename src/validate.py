"""Spec §8 dogrulamalari: Vega yuvarlak-gidis, renk terimi, bilinen yildizlar, S0."""
import os
import sys
import csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import response as resp
import spectra as spec
from photometry import Photometry

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MASTER = os.path.join(HERE, "out", "master_catalog.csv")

# Bilinen yildizlar (HIP): beklenen renk davranisi
KNOWN = {
    27989: ("Betelgeuse", "M2Iab", "cok kirmizi -> belirgin NEGATIF delta"),
    24436: ("Rigel", "B8Ia", "mavi -> delta ~0 veya pozitif"),
    32349: ("Sirius", "A1V", "A tipi -> delta ~0'a yakin"),
    69673: ("Arcturus", "K1.5III", "turuncu -> negatif delta"),
    91262: ("Vega", "A0V", "referans -> delta ~0"),
}

ok_all = True


def check(name, cond, detail):
    global ok_all
    print("  [%s] %s — %s" % ("GECTI" if cond else "KALDI", name, detail))
    if not cond:
        ok_all = False


print("=" * 68)
print("  DOGRULAMA")
print("=" * 68)

# ── Test 1: Vega yuvarlak-gidis ──
grid = resp.GRID_NM
r_inst = resp.system_response(grid)
r_v = resp.v_band_response(grid)
area, f_len, d_m = resp.aperture_area_cm2()
lib = spec.PicklesLibrary()
vfn, vsp, _ = lib.match(20.0, "V")
vlam, vflux = lib.spectrum(vfn)
phot = Photometry(grid, r_inst, r_v, vlam, vflux, area)

V_TEST = 3.14159
mi = phot.m_inst(vlam, vflux, V_TEST)
check("Test1 Vega yuvarlak-gidis", abs(mi - V_TEST) < 1e-3,
      "A0V spektrumu V=%.5f -> m_inst=%.5f (fark %.2e)" % (V_TEST, mi, abs(mi - V_TEST)))

# ── master yukle ──
rows = list(csv.DictReader(open(MASTER)))
dv = np.array([float(r["delta_v"]) for r in rows])
bv = np.array([float(r["bv"]) for r in rows if r["bv"] != ""])
dvb = np.array([float(r["delta_v"]) for r in rows if r["bv"] != ""])
corr = float(np.corrcoef(bv, dvb)[0, 1])

# ── Test 2: renk terimi yonu ──
check("Test2 renk terimi yonu", corr < -0.5,
      "corr(B-V, m_inst-V) = %.3f (guclu negatif beklenir)" % corr)

# monotonluk: B-V binlerinde ortalama delta azalmali
bins = [(-0.4, 0.0), (0.0, 0.4), (0.4, 0.8), (0.8, 1.2), (1.2, 2.5)]
means = []
for lo, hi in bins:
    sel = dvb[(bv >= lo) & (bv < hi)]
    means.append(sel.mean() if len(sel) else np.nan)
mono = all(means[i] > means[i + 1] for i in range(len(means) - 1)
           if np.isfinite(means[i]) and np.isfinite(means[i + 1]))
check("Test2b monotonluk", mono,
      "B-V bin ortalamalari: " + ", ".join("%.3f" % m for m in means))

# ── Test 3: bilinen yildizlar ──
by_hip = {int(r["hip"]): r for r in rows}
print("\n  Bilinen yildizlar:")
print("  %-12s %-8s %6s %6s %7s %7s  %s" % ("ad", "sptype", "V", "B-V", "m_inst", "delta", "beklenti"))
for hip, (name, sp_exp, note) in KNOWN.items():
    r = by_hip.get(hip)
    if not r:
        print("  %-12s (katalogda yok)" % name)
        continue
    print("  %-12s %-8s %6.2f %6.2f %7.3f %7.3f  %s"
          % (name, r["sptype"][:8], float(r["vmag"]),
             float(r["bv"]) if r["bv"] else float("nan"),
             float(r["m_inst"]), float(r["delta_v"]), note))

bet = by_hip.get(27989)
rig = by_hip.get(24436)
if bet and rig:
    check("Test3 Betelgeuse < Rigel", float(bet["delta_v"]) < float(rig["delta_v"]) - 0.3,
          "delta: Betelgeuse %.3f < Rigel %.3f" % (float(bet["delta_v"]), float(rig["delta_v"])))

# ── Test 4: S0 mertebesi ──
s0_new = phot.s0_electrons_per_s()
# YI_SIL mevcut kaba formulu
phi_spectral, band_width, tau_optics, qe_avg = 8.3e7, 320.0, 0.80, 0.50
area_m2 = area * 1e-4
s0_old = phi_spectral * band_width * area_m2 * tau_optics * qe_avg
ratio = s0_new / s0_old
check("Test4 S0 mertebesi", 0.1 < ratio < 10.0,
      "yeni %.4g / eski %.4g = %.2fx" % (s0_new, s0_old, ratio))

# ── Test 5: dagilim saglik ──
# Ust sinir: hicbir yildiz V'den belirgin SONUK olmamali (silisyum NIR'de genis bant).
# Alt sinir: gec-M devleri (M7/M8III) akisinin cogunu NIR'de yayar -> -4 kadire kadar
# NEGATIF delta FIZIKSELDIR. Kriter: uc degerler yalnizca gec-M tiplerinde olmali.
extreme = [r for r in rows if float(r["delta_v"]) < -2.0]
# gec-M esigi: tayf indeksi >= 66  (M0=60, M6=66) — 'M10III' gibi tipleri de kapsar
late_m = []
for r in extreme:
    p = spec.parse_sptype(r["matched_sptype"])
    if p is not None and p[0] >= 66.0:
        late_m.append(r)
check("Test5 delta araligi", dv.max() < 0.2 and dv.min() > -5.0,
      "delta_v min=%.3f max=%.3f ort=%.3f medyan=%.3f"
      % (dv.min(), dv.max(), dv.mean(), float(np.median(dv))))
check("Test5b uc degerler gec-M", len(extreme) == len(late_m),
      "delta<-2: %d yildiz, hepsi gec-M mi? %s (%.2f%% of catalog)"
      % (len(extreme), len(extreme) == len(late_m), 100.0 * len(extreme) / len(rows)))

print("\n" + "=" * 68)
print("  SONUC: %s" % ("TUM TESTLER GECTI" if ok_all else "BAZI TESTLER KALDI"))
print("=" * 68)
sys.exit(0 if ok_all else 1)
