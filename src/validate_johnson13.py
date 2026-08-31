"""BAGIMSIZ DOGRULUK DOGRULAMASI — Johnson & Mitchell 13-renk fotometrisi (VizieR II/84).

Lu & Wu (2019) makalesinin dogrulama yontemi:
  1. Her yildiz icin OLCULEN 13-renk fotometrisinden SED yeniden kurulur
     (13 bant, 337-1104 nm; sistem A0V=0 renklerine normalize).
  2. Bu OLCULEN SED, sistem tepkisiyle konvolue edilir -> "gercek" dM.
  3. Bizim Pickles+SpType tabanli dM'imizle karsilastirilir.
  4. Ayni gercek degere karsi B-V regresyonu da olculur (makale Tablo 1 karsiligi).

Cikti: hata ortalamasi/std — makalenin 0.058 std degeriyle kiyaslanabilir.
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

# Johnson-Mitchell 13-renk sistemi ortalama dalgaboylari [nm]
BANDS = {
    "33-52": 337.0, "35-52": 353.0, "37-52": 375.0, "40-52": 402.0, "45-52": 452.0,
    "52-58": 586.0, "52-63": 636.0, "52-72": 719.0, "52-80": 802.0,
    "52-86": 858.0, "52-99": 990.0, "52-110": 1104.0,
}
LAM_52 = 526.0
MATCH_ARCSEC = 10.0


def load_johnson():
    rows = []
    with open(os.path.join(HERE, "data", "johnson13.csv")) as f:
        for r in csv.DictReader(f):
            try:
                ra, de = float(r["_RA"]), float(r["_DE"])
            except (ValueError, KeyError):
                continue
            cols = {}
            ok = True
            for k in BANDS:
                v = r.get(k, "").strip()
                if v in ("", "--"):
                    ok = False
                    break
                cols[k] = float(v)
            if ok:
                rows.append(dict(hr=r.get("HR", ""), spt=r.get("SpT", ""),
                                 ra=ra, de=de, cols=cols))
    return rows


def sed_ratio(cols, vega_cols, grid_nm):
    """Olculen renklerden Vega'ya GORELI aki orani r(lambda) uretir (log-lineer interp)."""
    lam, logr = [LAM_52], [0.0]         # 52 bandi referans -> oran 1
    for k, lm in BANDS.items():
        dc = cols[k] - vega_cols[k]
        # '33-52' = m33-m52  -> ustel -0.4*dc ; '52-58' = m52-m58 -> ustel +0.4*dc
        e = -0.4 * dc if k.endswith("-52") else +0.4 * dc
        lam.append(lm)
        logr.append(e)
    lam = np.array(lam); logr = np.array(logr)
    o = np.argsort(lam); lam, logr = lam[o], logr[o]
    return 10.0 ** np.interp(grid_nm, lam, logr, left=logr[0], right=logr[-1])


def main():
    grid = resp.GRID_NM
    r_inst = resp.system_response(grid)
    r_v = resp.v_band_response(grid)
    lib = spec.PicklesLibrary()
    vfn, _, _ = lib.match(20.0, "V")
    vlam, vflux = lib.spectrum(vfn)
    f_vega = np.interp(grid, vlam, vflux, left=0.0, right=0.0)

    # Vega referans oranlari (bizim yontemin paydasi)
    Iiv = PCI(vlam, vflux, r_inst, grid)
    Ivv = PCI(vlam, vflux, r_v, grid)
    vega_ratio = Iiv / Ivv

    jr = load_johnson()
    vega = next((r for r in jr if r["hr"].strip() == "7001"), None)
    if vega is None:
        raise RuntimeError("Vega (HR 7001) 13-renk verisinde yok")
    vega_cols = vega["cols"]
    print("Johnson13: %d yildiz (tam renk seti) | Vega referansi alindi" % len(jr))

    # bizim katalog
    ours = list(csv.DictReader(open(os.path.join(HERE, "out", "master_catalog.csv"))))
    ora = np.array([float(r["ra_deg"]) for r in ours])
    ode = np.array([float(r["dec_deg"]) for r in ours])

    res = []
    for j in jr:
        d = np.hypot((ora - j["ra"]) * np.cos(np.radians(j["de"])), ode - j["de"]) * 3600.0
        i = int(np.argmin(d))
        if d[i] > MATCH_ARCSEC:
            continue
        o = ours[i]
        # OLCULEN SED -> gercek dM
        f_meas = f_vega * sed_ratio(j["cols"], vega_cols, grid)
        Ii = PCI(grid, f_meas, r_inst, grid)
        Iv = PCI(grid, f_meas, r_v, grid)
        if Ii <= 0 or Iv <= 0:
            continue
        dm_true = -2.5 * np.log10((Ii / Iv) / vega_ratio)
        res.append(dict(hr=j["hr"], spt=j["spt"], bv=(float(o["bv"]) if o["bv"] else np.nan),
                        dm_true=dm_true, dm_ours=float(o["delta_v"]),
                        sptype_ours=o["sptype"], matched=o["matched_sptype"],
                        varflag=o.get("varflag", ""), multflag=o.get("multflag", "")))

    if not res:
        raise RuntimeError("eslesen yildiz yok")

    dt = np.array([r["dm_true"] for r in res])
    do = np.array([r["dm_ours"] for r in res])
    err = do - dt
    print("eslesen yildiz: %d (<%.0f arcsec)" % (len(res), MATCH_ARCSEC))
    print()
    print("=" * 66)
    print("  SONUC — bizim yontem vs OLCULEN 13-renk fotometri")
    print("=" * 66)
    print("  ortalama hata : %+.4f kadir" % err.mean())
    print("  std sapma     :  %.4f kadir" % err.std())
    print("  medyan |hata| :  %.4f" % np.median(np.abs(err)))
    print("  maks |hata|   :  %.4f" % np.abs(err).max())

    # Ayni gercek degere karsi B-V regresyonu (makale Tablo 1 karsiligi)
    m = np.isfinite([r["bv"] for r in res])
    bv = np.array([r["bv"] for r in res])[m]
    dt_m, do_m = dt[m], do[m]
    print()
    print("  --- ayni %d yildizda B-V regresyonu (makale kiyasi) ---" % m.sum())
    print("  %-22s %10s %10s" % ("yontem", "ort", "std"))
    print("  %-22s %+10.4f %10.4f" % ("BIZIM (sentetik)", (do_m - dt_m).mean(), (do_m - dt_m).std()))
    for deg, nm in ((1, "B-V lineer"), (2, "B-V kuadratik"), (4, "B-V 4.derece")):
        c = np.polyfit(bv, dt_m, deg)          # GERCEK degere uydur (en iyi ihtimal)
        e = np.polyval(c, bv) - dt_m
        print("  %-22s %+10.4f %10.4f" % (nm, e.mean(), e.std()))

    # ── Makale ornekleme filtresi: degisken + cift yildizlari ELE ──
    # Lu&Wu: "Variable stars, double stars and stars with incomplete information
    # are deleted. Only 827 stars are used." -> adil kiyas icin ayni filtre.
    clean = [r for r in res
             if str(r.get("varflag", "")).strip() in ("", "--")
             and str(r.get("multflag", "")).strip() in ("", "--")]
    if clean:
        ec = np.array([r["dm_ours"] - r["dm_true"] for r in clean])
        print()
        print("  --- MAKALE FILTRESI (degisken+cift elendi): %d yildiz ---" % len(clean))
        print("  ortalama hata : %+.4f | std : %.4f | medyan|hata| : %.4f | maks : %.4f"
              % (ec.mean(), ec.std(), np.median(np.abs(ec)), np.abs(ec).max()))
        bvc = np.array([r["bv"] for r in clean])
        dtc = np.array([r["dm_true"] for r in clean])
        mm = np.isfinite(bvc)
        for deg, nm in ((1, "B-V lineer"), (4, "B-V 4.derece")):
            c = np.polyfit(bvc[mm], dtc[mm], deg)
            e = np.polyval(c, bvc[mm]) - dtc[mm]
            print("  %-22s std : %.4f" % (nm, e.std()))

    print()
    print("  --- en buyuk 6 sapma ---")
    res.sort(key=lambda r: -abs(r["dm_ours"] - r["dm_true"]))
    print("  %-6s %-12s %-10s %8s %8s %8s" % ("HR", "SpT(J13)", "eslesen", "gercek", "bizim", "hata"))
    for r in res[:6]:
        print("  %-6s %-12s %-10s %8.3f %8.3f %+8.3f"
              % (r["hr"].strip(), r["spt"][:12], r["matched"][:10],
                 r["dm_true"], r["dm_ours"], r["dm_ours"] - r["dm_true"]))
    print("=" * 66)


if __name__ == "__main__":
    main()
