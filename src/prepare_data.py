"""ZEMAX transmission dosyasini parse eder, QE ve Johnson V tablolarini yazar.

Cikti: data/t_optics_zemax.csv, data/qe_cmv4000_e5.csv, data/johnson_v.csv
"""
import os
import re
import csv

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
ZEMAX_SRC = r"C:\Users\derin\Desktop\datas\Transmission Data.txt"

# ── CMV4000 Mono E5 QE (datasheet DS000728 v8-01, Figure 7'den digitize) ──
# NOT: bu egri dolgu faktorunu (FF) ICERIR -> fill_factor ayrica carpilmamali.
QE_CMV4000_E5 = [
    (400, 0.33), (450, 0.50), (500, 0.60), (520, 0.64), (550, 0.60),
    (600, 0.58), (650, 0.58), (700, 0.52), (750, 0.44), (800, 0.35),
    (850, 0.24), (900, 0.15), (950, 0.08), (1000, 0.02),
]

# ── Johnson-Cousins V bandi (Bessell 1990, PASP 102:1181) ──
JOHNSON_V = [
    (470, 0.000), (480, 0.030), (490, 0.163), (500, 0.458), (510, 0.780),
    (520, 0.967), (530, 1.000), (540, 0.973), (550, 0.898), (560, 0.792),
    (570, 0.684), (580, 0.574), (590, 0.461), (600, 0.359), (610, 0.270),
    (620, 0.197), (630, 0.135), (640, 0.081), (650, 0.045), (660, 0.025),
    (670, 0.017), (680, 0.013), (690, 0.009), (700, 0.006), (710, 0.003),
    (720, 0.001), (730, 0.000),
]


def parse_zemax(path):
    """ZEMAX Polarization Transmission Data (UTF-16) -> {field_deg: {lambda_um: T}}"""
    with open(path, "rb") as f:
        raw = f.read()
    for enc in ("utf-16", "utf-16-le", "utf-8-sig", "latin-1"):
        try:
            text = raw.decode(enc)
            if "Transmission" in text:
                break
        except UnicodeDecodeError:
            continue
    else:
        raise RuntimeError("ZEMAX dosyasi cozulemedi")

    fields = {}
    cur = None
    for line in text.splitlines():
        s = " ".join(line.split())          # coklu bosluklari sadelestir
        # "Chief Ray ... Surface By Surface" bolumu Field Pos basliklarini TEKRAR eder
        # ama "Transmission at" satiri icermez -> sozlugu sifirlar. Orada dur.
        if s.startswith("Chief Ray Transmission Surface"):
            break
        m = re.match(r"Field Pos\s*:\s*([-\d.]+)\s*\(deg\)", s)
        if m:
            cur = float(m.group(1))
            fields[cur] = {}
            continue
        m = re.match(r"Transmission at\s+([\d.]+)\s*:\s*([\d.]+)", s)
        if m and cur is not None:
            fields[cur][float(m.group(1))] = float(m.group(2))
    return fields


def main():
    os.makedirs(DATA, exist_ok=True)

    # 1) ZEMAX
    fields = parse_zemax(ZEMAX_SRC)
    if not fields:
        raise RuntimeError("ZEMAX'tan hic alan okunamadi")
    onaxis_key = min(fields, key=abs)          # en kucuk |field| = on-axis
    onaxis = fields[onaxis_key]
    out = os.path.join(DATA, "t_optics_zemax.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lambda_nm", "T_optics"])
        for um in sorted(onaxis):
            w.writerow([um * 1000.0, onaxis[um]])
    print("yazildi: %s  (on-axis field=%.3f deg, %d dalgaboyu)"
          % (out, onaxis_key, len(onaxis)))
    print("  tum alan pozisyonlari: %s" % sorted(fields))

    # 2) QE
    out = os.path.join(DATA, "qe_cmv4000_e5.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lambda_nm", "QE"])
        w.writerows(QE_CMV4000_E5)
    print("yazildi: %s  (%d nokta)" % (out, len(QE_CMV4000_E5)))

    # 3) Johnson V
    out = os.path.join(DATA, "johnson_v.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lambda_nm", "R_V"])
        w.writerows(JOHNSON_V)
    print("yazildi: %s  (%d nokta)" % (out, len(JOHNSON_V)))


if __name__ == "__main__":
    main()
