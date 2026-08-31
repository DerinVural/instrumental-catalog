"""Pickles spektrum kutuphanesi: SpType ayristirma, eslestirme, fallback zinciri."""
import os
import re
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PICKLES = os.path.join(HERE, "data", "pickles")

# Tayf sinifi -> sicaklik indeksi tabani (sicak -> soguk artan)
_CLASS_BASE = {"O": 0, "B": 10, "A": 20, "F": 30, "G": 40, "K": 50, "M": 60}

# Isitma sinifi normalizasyonu (Ia/Iab/Ib -> I)
_LUM_ORDER = ["V", "IV", "III", "II", "I"]

# B-V -> anakol (V) tayf indeksi. Kaynak: standart anakol renkleri
# (Johnson; Mamajek tablosu ile uyumlu). SpType yoksa fallback icin.
_BV_TO_IDX = [
    (-0.33, 5.0),   # O5
    (-0.30, 10.0),  # B0
    (-0.17, 15.0),  # B5
    (-0.05, 19.0),  # B9
    (0.00, 20.0),   # A0
    (0.15, 25.0),   # A5
    (0.30, 30.0),   # F0
    (0.44, 35.0),   # F5
    (0.58, 40.0),   # G0
    (0.63, 42.0),   # G2
    (0.68, 45.0),   # G5
    (0.81, 50.0),   # K0
    (1.15, 55.0),   # K5
    (1.40, 60.0),   # M0
    (1.49, 62.0),   # M2
    (1.64, 65.0),   # M5
]


def _parse_subnum(digits):
    """Alt-sinif sayisi.

    '5'    -> 5.0
    '57'   -> 6.0    (B5-B7 ARALIGI: ilk basamak < ikinci)
    '10'   -> 10.0   (ARALIK DEGIL: ilk >= ikinci -> tek sayi; Pickles 'M10III')
    '2.5'  -> 2.5
    ''     -> 5.0    (sinif ortasi varsayilan)
    """
    digits = digits.strip()
    if not digits:
        return 5.0
    if "." in digits:
        try:
            return float(digits)
        except ValueError:
            return 5.0
    if len(digits) == 2 and digits.isdigit():
        a, b = int(digits[0]), int(digits[1])
        if a < b:
            return (a + b) / 2.0        # gercek aralik (57, 12, 01, 47, ...)
        return float(digits)            # tek sayi (10) — aralik degil
    try:
        return float(digits)
    except ValueError:
        return 5.0


def parse_sptype(sp):
    """SpType metni -> (idx, lum) veya None.

    idx: surekli tayf indeksi (O0=0 ... M9=69), lum: 'V','IV','III','II','I'
    Ornekler: 'B3Vp'->(13,'V'), 'K5III'->(55,'III'), 'F7:Ib-IIv SB'->(37,'I'),
              'B9'->(19,'V' varsayilan), 'M2.5V'->(62.5,'V')
    """
    if sp is None:
        return None
    s = str(sp).strip()
    if not s or s.lower() in ("--", "nan", "none"):
        return None
    m = re.match(r"^\s*([a-z]{0,2})\s*([OBAFGKM])\s*([\d.]*)(.*)$", s)
    if not m:
        return None
    prefix, letter, digits, rest = m.group(1), m.group(2), m.group(3), m.group(4)
    idx = _CLASS_BASE[letter] + _parse_subnum(digits)

    lum = None
    lm = re.search(r"(IV|III|II|Iab|Ia|Ib|I|V)", rest)
    if lm:
        tok = lm.group(1)
        lum = "I" if tok in ("Ia", "Ib", "Iab", "I") else tok
    if lum is None and prefix:
        # ESKI TARZ onekler (Johnson/Mitchell, HD): g=dev, d=cuce, sd=altcuce, c=superdev
        lum = {"g": "III", "d": "V", "sd": "V", "c": "I"}.get(prefix)
    if lum is None:
        lum = "V"          # isitma sinifi belirtilmemis -> anakol varsay
    return idx, lum


def bv_to_index(bv):
    """B-V -> anakol tayf indeksi (lineer interpolasyon, uclarda sabit)."""
    if bv is None or not np.isfinite(bv):
        return None
    xs = [p[0] for p in _BV_TO_IDX]
    ys = [p[1] for p in _BV_TO_IDX]
    return float(np.interp(bv, xs, ys))


class PicklesLibrary:
    """Pickles UVILIB: tip eslestirme + spektrum yukleme (onbellekli)."""

    def __init__(self, path=PICKLES):
        self.path = path
        idx = fits.open(os.path.join(path, "pickles.fits"))[1].data
        self.entries = []          # (filename, sptype_str, idx, lum)
        for row in idx:
            fn = row["FILENAME"].strip()
            sp = row["SPTYPE"].strip()
            # w*/r* = metalce fakir/zengin varyantlar -> gunes metalligini kullan, atla
            if sp[0] in ("w", "r"):
                continue
            p = parse_sptype(sp)
            if p is None:
                continue
            self.entries.append((fn, sp, p[0], p[1]))
        if not self.entries:
            raise RuntimeError("Pickles indeksi bos/okunamadi")
        self._cache = {}

    def match(self, idx, lum):
        """En yakin Pickles tipi: once AYNI isitma sinifi, yoksa en yakin sinif.

        Donus: (filename, sptype_str, exact_lum_bool)
        """
        same = [e for e in self.entries if e[3] == lum]
        if same:
            best = min(same, key=lambda e: abs(e[2] - idx))
            return best[0], best[1], True
        # isitma sinifi yok -> tum kutuphaneden en yakin (once idx, sonra sinif mesafesi)
        li = _LUM_ORDER.index(lum) if lum in _LUM_ORDER else 0
        best = min(
            self.entries,
            key=lambda e: (abs(e[2] - idx), abs(_LUM_ORDER.index(e[3]) - li)),
        )
        return best[0], best[1], False

    def spectrum(self, filename):
        """(lambda_nm, flux) — goreli akı, artan dalgaboyu."""
        if filename in self._cache:
            return self._cache[filename]
        d = fits.open(os.path.join(self.path, filename + ".fits"))[1].data
        lam_nm = np.asarray(d["WAVELENGTH"], float) / 10.0    # Angstrom -> nm
        flux = np.asarray(d["FLUX"], float)
        order = np.argsort(lam_nm)
        out = (lam_nm[order], flux[order])
        self._cache[filename] = out
        return out


def lum_class_from_absmag(vmag, plx_mas, e_plx_mas=None):
    """Mutlak kadirden isitma sinifi (SpType'ta sinif YOKKEN kullanilir).

    M_V = V + 5 + 5 log10(plx_arcsec).  Paralaks guvenilir degilse None.
    Esikler kaba ama ayirici: dev/cuce farki bu bantta ~1 kadire kadar dM degistirir.
    """
    if plx_mas is None or not np.isfinite(plx_mas) or plx_mas <= 1.0:
        return None                      # >1 kpc veya negatif -> guvenilmez
    if e_plx_mas is not None and np.isfinite(e_plx_mas) and e_plx_mas > 0:
        if plx_mas / e_plx_mas < 5.0:    # <5 sigma -> guvenme
            return None
    if vmag is None or not np.isfinite(vmag):
        return None
    m_v = vmag + 5.0 + 5.0 * np.log10(plx_mas / 1000.0)
    if m_v < 0.0:
        return "I"
    if m_v < 3.5:
        return "III"
    if m_v < 4.5:
        return "IV"
    return "V"


def resolve_spectrum(lib, sptype, bv, vmag=None, plx_mas=None, e_plx_mas=None):
    """Fallback zinciri: SpType(+paralaks) -> B-V -> G2V.

    Donus: (filename, matched_sptype, flag)
      flag: 'sptype' | 'sptype_lum_plx' | 'sptype_lum_approx' | 'bv' | 'default_g2v'
    """
    p = parse_sptype(sptype)
    if p is not None:
        idx, lum = p
        # SpType'ta isitma sinifi YOKSA parse_sptype 'V' varsayar. Parlak katalogda
        # bu gec-M yildizlarda YANLIS (M6-M8 cuceler cok sonuk; bunlar dev).
        # Paralaks varsa mutlak kadirden gercek sinifi tureti.
        explicit = re.search(r"(IV|III|II|Iab|Ia|Ib|I|V)", str(sptype or ""))
        if not explicit:
            lum_plx = lum_class_from_absmag(vmag, plx_mas, e_plx_mas)
            if lum_plx is not None and lum_plx != lum:
                fn, sp, exact = lib.match(idx, lum_plx)
                return fn, sp, "sptype_lum_plx"
        fn, sp, exact = lib.match(idx, lum)
        return fn, sp, ("sptype" if exact else "sptype_lum_approx")

    idx = bv_to_index(bv)
    if idx is not None:
        fn, sp, _ = lib.match(idx, "V")
        return fn, sp, "bv"

    fn, sp, _ = lib.match(42.0, "V")      # G2V
    return fn, sp, "default_g2v"
