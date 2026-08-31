"""Radyometri / tespit butcesi — m_inst'ten SNR, tespit limiti ve attitude hatasina.

Katalog yalnizca m_inst tutar (goreli). Ama sim ve ucus yazilimi "bu pikselde kac
elektron birikecek" sorusunu sorar. S0 bu koprunun tek sabitidir:

    N (e-/s) = S0 * 10^(-0.4 * m_inst)

Bu modul zinciri sonuna kadar goturur:
    m_inst -> e-/s -> (poz) -> e- -> (gurultu) -> SNR -> centroid -> ATTITUDE

Cikti report.md'ye yazilir; boylece S0 degisince (or. AR kaplama duzeltmesi)
tespit limitinin ne kadar kaydigi SESSIZ KALMAZ.
"""
import numpy as np

# ── Dedektor / sahne parametreleri — YI_SIL star_tracker_params.h ile AYNI ──
READ_NOISE_E = 13.0        # e- rms
DARK_E_PX_S = 7.5          # e-/px/s
SKY_E_PX_S = 0.3           # e-/px/s
FULL_WELL_E = 13500.0      # tam kuyu
T_EXP_S = (0.015, 0.050)   # min, varsayilan poz
PLATE_ARCSEC_PX = 25.84    # SIM optigi (FOV 14.7 -> f=43.662mm); bkz. README optik konvansiyon

# ── Fotometri penceresi ──
N_PIX = 9                  # 3x3 merkez penceresi (PSF EE80 yaricapi ~1.35 px)
SIGMA_PSF_PX = 0.70        # ZEMAX LUT on-axis 2. moment (olculdu)
PEAK_FRAC = 0.27           # toplam akinin tepe pikseldeki payi (7x7 kernel, olculdu)

SNR_THRESHOLDS = (6.0, 10.0, 15.0)


def electrons(s0, m_inst, t_exp):
    """Toplam sinyal [e-]."""
    return s0 * 10.0 ** (-0.4 * m_inst) * t_exp


def snr(s0, m_inst, t_exp, n_pix=N_PIX):
    """Sok + okuma + karanlik + gokyuzu gurultusuyle SNR."""
    e = electrons(s0, m_inst, t_exp)
    var = e + n_pix * (READ_NOISE_E ** 2 + (DARK_E_PX_S + SKY_E_PX_S) * t_exp)
    return e / np.sqrt(var)


def centroid_sigma_px(s0, m_inst, t_exp):
    """Merkez bulma hatasi ~ sigma_PSF / SNR."""
    s = snr(s0, m_inst, t_exp)
    return SIGMA_PSF_PX / s if s > 0 else np.inf


def limiting_mag(s0, t_exp, snr_thr):
    """Verilen SNR esigini saglayan en sonuk m_inst (ikili arama)."""
    lo, hi = -5.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if snr(s0, mid, t_exp) > snr_thr:
            lo = mid
        else:
            hi = mid
    return lo


def saturation_mag(s0, t_exp):
    """Tepe pikselin tam kuyuyu doldurdugu m_inst (bundan PARLAK yildizlar doyar)."""
    e_sat = FULL_WELL_E / PEAK_FRAC          # doyma icin gereken TOPLAM sinyal
    return -2.5 * np.log10(e_sat / (s0 * t_exp))


def write_report_section(f, s0, mag_limit):
    """report.md'ye radyometri butcesi bolumunu yazar."""
    t_def = T_EXP_S[1]
    f.write("## Radyometri butcesi (S0 -> tespit)\n\n")
    f.write("Zincir: `m_inst -> e-/s -> poz -> e- -> SNR -> centroid -> attitude`. "
            "Poz %.0f ms, okuma gurultusu %.0f e-, pencere %d px, sigma_PSF %.2f px, "
            "plaka olcegi %.2f as/px.\n\n"
            % (t_def * 1000, READ_NOISE_E, N_PIX, SIGMA_PSF_PX, PLATE_ARCSEC_PX))

    f.write("| m_inst | e-/s | e- (%.0f ms) | SNR | centroid px | tek yildiz (as) |\n"
            % (t_def * 1000))
    f.write("|---|---|---|---|---|---|\n")
    for m in (2, 4, 6, mag_limit, mag_limit + 1):
        e_s = s0 * 10.0 ** (-0.4 * m)
        e = electrons(s0, m, t_def)
        s = snr(s0, m, t_def)
        c = centroid_sigma_px(s0, m, t_def)
        f.write("| %.1f | %.3g | %.0f | %.1f | %.4f | %.3f |\n"
                % (m, e_s, e, s, c, c * PLATE_ARCSEC_PX))

    f.write("\n### Tespit limiti\n\n")
    f.write("| poz | " + " | ".join("SNR>%.0f" % t for t in SNR_THRESHOLDS) + " |\n")
    f.write("|---|" + "---|" * len(SNR_THRESHOLDS) + "\n")
    for t_exp in T_EXP_S:
        vals = " | ".join("%.2f" % limiting_mag(s0, t_exp, thr) for thr in SNR_THRESHOLDS)
        f.write("| %.0f ms | %s |\n" % (t_exp * 1000, vals))

    m_sat = saturation_mag(s0, t_def)
    f.write("\n### Doyma\n\n")
    f.write("- Tam kuyu %.0f e-, tepe piksel payi %.0f%% -> **m_inst < %.2f olan yildizlar "
            "DOYAR** (%.0f ms pozda).\n" % (FULL_WELL_E, PEAK_FRAC * 100, m_sat, t_def * 1000))
    f.write("- Doymus yildizin merkezi kayar; algoritma tarafina 'bu listeye merkez guvenme' "
            "notu gerekir.\n")

    # 15 yildizla attitude
    c6 = centroid_sigma_px(s0, 6.0, t_def) * PLATE_ARCSEC_PX
    f.write("\n### Attitude (kaba kestirim)\n\n")
    f.write("- m_inst=6 yildizlar, 15 yildizli cozum: **%.3f arcsec** "
            "(tek yildiz %.3f / sqrt(15)).\n" % (c6 / np.sqrt(15.0), c6))
    f.write("- NOT: yalnizca foton/okuma gurultusu. Smear, alan-bagimli PSF, kalibrasyon "
            "ve katalog hatalari DAHIL DEGIL.\n\n")
