"""Sistem spektral tepkisi: R_inst(lambda) = QE(lambda) x T_optics(lambda), ve R_V.

Tum egriler ortak 1 nm izgaraya (LAM_MIN..LAM_MAX) interpole edilir.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")

LAM_MIN = 380.0     # nm
LAM_MAX = 1000.0    # nm
LAM_STEP = 1.0

GRID_NM = np.arange(LAM_MIN, LAM_MAX + LAM_STEP, LAM_STEP)


def _load_csv(name):
    """iki kolonlu CSV (lambda_nm, deger) -> (lam, val).

    '#' ile baslayan yorum satirlari ve basliksatiri atlanir (veri dosyalarinda
    proveniyans/uyari notlari bulunabilir — bkz. t_optics_zemax.csv).
    """
    path = os.path.join(DATA, name)
    lam, val = [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split(",")
            if len(parts) < 2:
                continue
            try:
                lam.append(float(parts[0]))
                val.append(float(parts[1]))
            except ValueError:
                continue                      # baslik satiri
    if not lam:
        raise ValueError("%s: sayisal veri bulunamadi" % name)
    return np.asarray(lam, float), np.asarray(val, float)


def _interp_clamped(lam_src, val_src, grid):
    """Izgaraya interpole; olculen aralik disinda UC DEGERE sabitlenir (clamp).

    QE 1000 nm'de %2'ye dustugu icin ust-uc sabitlemenin etkisi ihmal edilebilir;
    T_optics 900 nm ustu sabitlenir (ZEMAX 900 nm'de bitiyor).
    """
    order = np.argsort(lam_src)
    lam_src, val_src = lam_src[order], val_src[order]
    return np.interp(grid, lam_src, val_src, left=val_src[0], right=val_src[-1])


def system_response(grid=GRID_NM):
    """R_inst(lambda) = QE x T_optics  (dolgu faktoru QE egrisinde ZATEN dahil)."""
    lam_q, qe = _load_csv("qe_cmv4000_e5.csv")
    lam_t, t_opt = _load_csv("t_optics_zemax.csv")
    qe_i = _interp_clamped(lam_q, qe, grid)
    t_i = _interp_clamped(lam_t, t_opt, grid)
    # QE olculen aralik disinda (400 nm alti) hizla duser -> 400 nm altini sifirla
    qe_i = np.where(grid < lam_q.min(), 0.0, qe_i)
    return qe_i * t_i


def v_band_response(grid=GRID_NM):
    """Johnson-Cousins V bant tepkisi; olculen aralik disi SIFIR (bant sinirli)."""
    lam_v, rv = _load_csv("johnson_v.csv")
    r = np.interp(grid, lam_v, rv, left=0.0, right=0.0)
    return np.clip(r, 0.0, None)


def aperture_area_cm2(pixel_num=2048, pixel_size_m=5.5e-6, fov_deg=14.7, f_num=1.585):
    """YI_SIL star_tracker_params.h ile AYNI zincir: f_length -> D -> alan."""
    f_length = (pixel_num * pixel_size_m) / (2.0 * np.tan(np.radians(fov_deg) / 2.0))
    d_m = f_length / f_num
    area_m2 = np.pi * (d_m / 2.0) ** 2
    return area_m2 * 1.0e4, f_length, d_m
