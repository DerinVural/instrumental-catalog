"""Instrumental kadir katalogu ureteci (ana akis).

Girdi : hipparcos.fits (astrometri+V+B-V) + hip_sptype.csv (SpType) + Pickles + tepki egrileri
Cikti : out/master_catalog.csv, out/hip_instrumental.txt, out/report.md
"""
import os
import csv
import sys
import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import response as resp
import spectra as spec
from photometry import Photometry
import radiometry as radio

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "out")
HIP_FITS = r"C:\Users\derin\star-sim-tu\data\catalogs\hipparcos.fits"

from config import MAG_LIMIT, COLOR_MARGIN, VMAG_PREFILTER
MAS_TO_RAD = 4.84813681e-9      # milli-arcsec -> radyan
DEC_CLIP_DEG = 89.5             # cos(dec) bolmesi icin guvenlik siniri
TARGET_EPOCH = 2000.0           # cikti epogu (J2000, mevcut hygdata konvansiyonu)

# [B3] Kaynak epogu BASLIKTAN OKUNMAZ, VERIDEN OLCULUR.
#
# UYARI (kendi hatamiz, 2026-08-31): FITS basligi EPOCH='J2000' diyor, ama bu
# EKINOKS'u (ICRS/J2000 koordinat cercevesi) belirtiyor; konumlarin EPOGUNU
# degil. Hipparcos konumlari J1991.25 epogundadir. Basliga guvenip PM
# tasimasini kaldirmak, yuksek oz-hareketli yildizlarda 60-90 arcsec hataya
# yol acar (Barnard 90.6", Groombridge 1830 61.8").
#
# Bu yuzden epok, BILINEN J2000 konumlariyla kiyaslanarak belirlenir.
# Referanslar: SIMBAD ICRS J2000 (epok 2000.0), yuksek oz-hareketli yildizlar.
EPOCH_REF_J2000 = {
    # HIP: (ad, RA_deg, Dec_deg)
    87937: ("Barnard", 269.4520769, 4.693391),
    57939: ("Groombridge 1830", 178.2448713, 37.7186775),
    24186: ("Kapteyn", 77.919083, -45.018417),
    54035: ("Lalande 21185", 165.834125, 35.969889),
}
EPOCH_CANDIDATES = (1991.25, 2000.0)      # denenecek kaynak epoklari
EPOCH_TOL_ARCSEC = 0.5                    # kabul esigi


def measure_source_epoch(tab, target_epoch=2000.0):
    """Kaynak epogunu VERIDEN olcer: her aday epok icin referans yildizlari
    hedef epoga tasiyip bilinen J2000 konumlariyla kiyaslar; en kucuk artigi
    veren adayi secer. Hicbiri esigi saglamazsa KOSUYU DURDURUR."""
    hips = np.asarray(tab["HIP"])
    best = None
    detail = []
    for cand in EPOCH_CANDIDATES:
        dt = target_epoch - cand
        errs = []
        for hip, (nm, ra0, de0) in EPOCH_REF_J2000.items():
            idx = np.where(hips == hip)[0]
            if not len(idx):
                continue
            i = idx[0]
            ra, de = float(tab["RA"][i]), float(tab["DEC"][i])
            pra, pde = float(tab["PMRA"][i]), float(tab["PMDEC"][i])
            cosd = np.cos(np.radians(de))
            ra_t = ra + dt * (pra / 1000.0 / 3600.0) / max(cosd, 1e-6)
            de_t = de + dt * (pde / 1000.0 / 3600.0)
            errs.append(np.hypot((ra_t - ra0) * np.cos(np.radians(de0)),
                                 de_t - de0) * 3600.0)
        if not errs:
            continue
        med = float(np.median(errs))
        detail.append((cand, med, len(errs)))
        if best is None or med < best[1]:
            best = (cand, med, len(errs))

    if best is None:
        raise SystemExit("[B3] Epok olculemedi: referans yildizlarin hicbiri bulunamadi.")
    for cand, med, n in detail:
        print("     epok adayi J%.2f -> medyan artik %8.3f arcsec (%d yildiz)"
              % (cand, med, n))
    if best[1] > EPOCH_TOL_ARCSEC:
        raise SystemExit(
            "[B3] Epok belirlenemedi: en iyi aday J%.2f bile %.2f arcsec artik "
            "birakiyor (esik %.2f). Konumlar/PM beklenenden farkli — VARSAYIM "
            "YAPILMAYACAK." % (best[0], best[1], EPOCH_TOL_ARCSEC))
    return best[0], best[1]


def load_sptype_map():
    """hip -> dict(sptype, plx, e_plx, varflag, multflag)"""
    path = os.path.join(DATA, "hip_sptype.csv")
    m = {}

    def fnum(s):
        try:
            return float(s)
        except (TypeError, ValueError):
            return float("nan")

    with open(path) as f:
        for row in csv.DictReader(f):
            try:
                hip = int(row["hip"])
            except (ValueError, KeyError):
                continue
            m[hip] = dict(
                sptype=row.get("sptype", "").strip(),
                plx=fnum(row.get("plx", "")),
                e_plx=fnum(row.get("e_plx", "")),
                varflag=row.get("varflag", "").strip(),
                multflag=row.get("multflag", "").strip(),
            )
    return m


def main():
    os.makedirs(OUT, exist_ok=True)

    # ── 1. Tepki egrileri + acikllik ──
    grid = resp.GRID_NM
    r_inst = resp.system_response(grid)
    r_v = resp.v_band_response(grid)
    area_cm2, f_length, d_m = resp.aperture_area_cm2()
    print("[1] tepki egrileri hazir | aciklik D=%.2f mm, A=%.3f cm^2, f=%.2f mm"
          % (d_m * 1e3, area_cm2, f_length * 1e3))

    # ── 2. Pickles + Vega (A0V) ──
    lib = spec.PicklesLibrary()
    vfn, vsp, _ = lib.match(20.0, "V")           # A0V = Vega vekili
    vega_lam, vega_flux = lib.spectrum(vfn)
    print("[2] Pickles: %d tip | Vega vekili: %s (%s)" % (len(lib.entries), vsp, vfn))

    phot = Photometry(grid, r_inst, r_v, vega_lam, vega_flux, area_cm2)
    s0 = phot.s0_electrons_per_s()
    print("[3] S0 (m_inst=0 -> e-/s) = %.4g" % s0)

    # ── 3. Katalog verisi ──
    _hdu = fits.open(HIP_FITS)[1]
    hip_tab = _hdu.data
    print("[3b] epok VERIDEN olculuyor (baslik EPOCH='%s' — EKINOKS, epok DEGIL):"
          % str(_hdu.header.get("EPOCH", "?")).strip())
    src_epoch, resid = measure_source_epoch(hip_tab, TARGET_EPOCH)
    dt_epoch = TARGET_EPOCH - src_epoch
    print("     SECILEN: J%.2f (artik %.3f arcsec) -> PM tasima %+.2f yil"
          % (src_epoch, resid, dt_epoch))
    sptype_map = load_sptype_map()
    print("[4] hipparcos.fits: %d yildiz | SpType kaydi: %d" % (len(hip_tab), len(sptype_map)))

    rows = []
    stats = {"flags": {}, "skipped_v": 0, "skipped_calc": 0, "dec_clip": 0}

    for rec in hip_tab:
        vmag = float(rec["VMAG"])
        if not np.isfinite(vmag):
            stats["skipped_v"] += 1
            continue
        if vmag > VMAG_PREFILTER:          # [B1] turetilmis esik (config.py)
            continue

        hip = int(rec["HIP"])
        bv = float(rec["BV"]) if np.isfinite(rec["BV"]) else None
        meta = sptype_map.get(hip, {})
        sptype = meta.get("sptype", "")

        fn, matched_sp, flag = spec.resolve_spectrum(
            lib, sptype, bv, vmag=vmag,
            plx_mas=meta.get("plx", float("nan")),
            e_plx_mas=meta.get("e_plx", float("nan")))
        lam, flux = lib.spectrum(fn)
        mi = phot.m_inst(lam, flux, vmag)
        if mi is None:
            stats["skipped_calc"] += 1
            continue
        if mi > MAG_LIMIT:
            continue

        # ── astrometri: derece -> radyan, pm mas/yr -> rad/yr ──
        ra_deg = float(rec["RA"])
        dec_deg = float(rec["DEC"])
        pmra_mas = float(rec["PMRA"]) if np.isfinite(rec["PMRA"]) else 0.0   # mu_alpha*
        pmdec_mas = float(rec["PMDEC"]) if np.isfinite(rec["PMDEC"]) else 0.0

        # YI_SIL: ra += dt*pmra  -> cos(dec)'SIZ mu_alpha bekler; Hipparcos cos(dec) ICERIR
        cosd = np.cos(np.radians(min(abs(dec_deg), DEC_CLIP_DEG)))
        if abs(dec_deg) > DEC_CLIP_DEG:
            stats["dec_clip"] += 1
            flag += "+decclip"
        pmra_rad = (pmra_mas / cosd) * MAS_TO_RAD       # rad/yr, cos(dec)'siz
        pmdec_rad = pmdec_mas * MAS_TO_RAD

        # kaynak epok -> J2000 (dt_epoch basliktan turetildi; 0 ise tasima yok)
        ra_rad = np.radians(ra_deg) + dt_epoch * pmra_rad
        dec_rad = np.radians(dec_deg) + dt_epoch * pmdec_rad

        e_per_s = s0 * 10.0 ** (-0.4 * mi)
        stats["flags"][flag] = stats["flags"].get(flag, 0) + 1
        rows.append(dict(
            hip=hip, ra_deg=ra_deg, dec_deg=dec_deg,
            pmra_mas=pmra_mas, pmdec_mas=pmdec_mas,
            vmag=vmag, bv=(bv if bv is not None else ""), sptype=sptype,
            pickles_file=fn, matched_sptype=matched_sp,
            m_inst=mi, delta_v=mi - vmag, e_per_s=e_per_s,
            varflag=meta.get("varflag", ""), multflag=meta.get("multflag", ""),
            ra_rad=ra_rad, dec_rad=dec_rad,
            pmra_rad=pmra_rad, pmdec_rad=pmdec_rad, flags=flag,
        ))

    rows.sort(key=lambda r: r["m_inst"])
    print("[5] hesaplandi: %d yildiz (m_inst <= %.1f)" % (len(rows), MAG_LIMIT))

    # [B1] EKSIKLIK GUARD'i: renk terimi payi yetmediyse katalog EKSIKTIR.
    observed_min = min(r["delta_v"] for r in rows)
    if observed_min < -COLOR_MARGIN:
        raise RuntimeError(
            "COLOR_MARGIN=%.2f YETERSIZ: olculen min delta_v=%.3f. "
            "V<=%.2f on-elemesi, m_inst<=%.1f kosulunu saglayan bazi kirmizi "
            "yildizlari disarida birakmis olabilir -> KATALOG EKSIK. "
            "config.py'de COLOR_MARGIN'i buyutup yeniden calistir."
            % (COLOR_MARGIN, observed_min, VMAG_PREFILTER, MAG_LIMIT))
    print("[5b] renk payi: |min delta_v|=%.3f / COLOR_MARGIN=%.2f -> kalan pay %.3f kadir"
          % (abs(observed_min), COLOR_MARGIN, COLOR_MARGIN - abs(observed_min)))

    # ── 4. master ──
    master = os.path.join(OUT, "master_catalog.csv")
    cols = ["hip", "ra_deg", "dec_deg", "pmra_mas", "pmdec_mas", "vmag", "bv",
            "sptype", "pickles_file", "matched_sptype", "m_inst", "delta_v",
            "e_per_s", "varflag", "multflag", "flags"]
    with open(master, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # ── 5. onboard (YI_SIL drop-in; hygdata.txt ile BIREBIR format) ──
    onboard = os.path.join(OUT, "hip_instrumental.txt")
    with open(onboard, "w", newline="") as f:
        f.write("MAG\tRARAD\tDECRAD\tPMRARAD\tPMDECRAD\n")
        for r in rows:
            # .12g: RA/Dec radyanda ~0.001 mas cozunurluk (.9g 0.25 mas biraktiriyordu)
            f.write("%.4f\t%.12g\t%.12g\t%.9g\t%.9g\n"
                    % (r["m_inst"], r["ra_rad"], r["dec_rad"],
                       r["pmra_rad"], r["pmdec_rad"]))

    # ── 6. rapor ──
    dv = np.array([r["delta_v"] for r in rows])
    bvs = np.array([r["bv"] for r in rows if r["bv"] != ""], float)
    dv_bv = np.array([r["delta_v"] for r in rows if r["bv"] != ""])
    corr = float(np.corrcoef(bvs, dv_bv)[0, 1]) if len(bvs) > 2 else float("nan")

    rep = os.path.join(OUT, "report.md")
    with open(rep, "w") as f:
        f.write("# Instrumental Kadir Katalogu — Rapor\n\n")
        f.write("## Radyometri\n\n")
        f.write("| | |\n|---|---|\n")
        f.write("| Aciklik capi | %.3f mm |\n" % (d_m * 1e3))
        f.write("| Aciklik alani | %.4f cm^2 |\n" % area_cm2)
        f.write("| Odak uzakligi | %.3f mm |\n" % (f_length * 1e3))
        f.write("| **S0 (m_inst=0 -> e-/s)** | **%.6g** |\n" % s0)
        f.write("| Vega vekili | %s (%s) |\n\n" % (vsp, vfn))
        f.write("> **Optik konvansiyon:** Radyometri, YI_SIL `star_tracker_params.h`\n"
                "> degerlerinden turetilir (FOV 14.7 deg -> f=43.662 mm, f/1.585,\n"
                "> A=5.960 cm^2). Gercek Zemax tasarimi farklidir (f=45.237 mm LUT\n"
                "> basligindan, f/1.685) -> aciklik alani 5.661 cm^2. Ucus degerleri\n"
                "> istenirse S0'i **0.947** ile carpin (=%.6g e-/s).\n"
                "> Karar: sim ile tutarlilik esas alindi.\n\n" % (s0 * 0.947))
        radio.write_report_section(f, s0, MAG_LIMIT)

        f.write("## On-eleme ve renk payi [B1]\n\n")
        f.write("| | |\n|---|---|\n")
        f.write("| m_inst limiti | %.1f |\n" % MAG_LIMIT)
        f.write("| COLOR_MARGIN | %.2f |\n" % COLOR_MARGIN)
        f.write("| V on-eleme esigi | %.2f (= limit + pay) |\n" % VMAG_PREFILTER)
        f.write("| olculen delta_v min / max | %.3f / %.3f |\n" % (dv.min(), dv.max()))
        f.write("| **kalan pay** | **%.3f kadir** |\n\n" % (COLOR_MARGIN - abs(dv.min())))

        f.write("## Yildiz sayilari\n\n")
        f.write("- m_inst <= %.1f: **%d yildiz**\n" % (MAG_LIMIT, len(rows)))
        f.write("- V gecersiz (atlandi): %d\n" % stats["skipped_v"])
        f.write("- hesap basarisiz (atlandi): %d\n" % stats["skipped_calc"])
        f.write("- |dec| > %.1f (pmRA clip): %d\n\n" % (DEC_CLIP_DEG, stats["dec_clip"]))
        f.write("## Spektrum eslestirme (flag dagilimi)\n\n| flag | adet |\n|---|---|\n")
        for k, v in sorted(stats["flags"].items(), key=lambda x: -x[1]):
            f.write("| %s | %d |\n" % (k, v))
        f.write("\n## Renk terimi (m_inst - V)\n\n")
        f.write("- ortalama: %.3f kadir\n" % dv.mean())
        f.write("- min / max: %.3f / %.3f\n" % (dv.min(), dv.max()))
        f.write("- **B-V ile korelasyon: %.3f** (negatif beklenir: kirmizi yildiz daha parlak)\n\n" % corr)
        f.write("## Cikti dosyalari\n\n")
        f.write("- `master_catalog.csv` — izlenebilirlik (tum alanlar)\n")
        f.write("- `hip_instrumental.txt` — YI_SIL drop-in (MAG kolonu = m_inst)\n")

    print("[6] yazildi:\n    %s\n    %s\n    %s" % (master, onboard, rep))
    print("    S0=%.6g | renk terimi ort=%.3f | B-V korelasyon=%.3f" % (s0, dv.mean(), corr))


if __name__ == "__main__":
    main()
