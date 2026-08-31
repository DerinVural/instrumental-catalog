"""VizieR I/239/hip_main -> HIP basina SpType ceker, data/hip_sptype.csv yazar.

Yalnizca Vmag <= VMAG_CUT olan yildizlar (m_inst, V'den en fazla ~1 kadir parlak
olabilir; m_inst<=7.0 hedefi icin V<=8.0 guvenli pay birakir).
"""
import os
import csv
import warnings

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(DATA, "hip_sptype.csv")
VMAG_CUT = 8.0


def main():
    from astroquery.vizier import Vizier

    v = Vizier(columns=["HIP", "Vmag", "B-V", "SpType",
                        "Plx", "e_Plx", "VarFlag", "MultFlag"], row_limit=-1)
    print("VizieR sorgusu: I/239/hip_main, Vmag < %.1f ..." % VMAG_CUT)
    res = v.query_constraints(catalog="I/239/hip_main", Vmag="<%.1f" % VMAG_CUT)
    if not len(res):
        raise RuntimeError("VizieR bos dondu")
    t = res[0]
    print("  gelen satir: %d" % len(t))

    os.makedirs(DATA, exist_ok=True)
    n_sp = 0
    def g(row, key):
        v = row[key]
        try:
            return "" if v is None or (hasattr(v, "mask") and v is None) else str(v).strip()
        except Exception:
            return ""

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hip", "vmag", "bv", "sptype", "plx", "e_plx", "varflag", "multflag"])
        for row in t:
            hip = row["HIP"]
            sp = g(row, "SpType")
            if sp and sp.lower() != "--":
                n_sp += 1
            w.writerow([hip, row["Vmag"], row["B-V"], sp,
                        g(row, "Plx"), g(row, "e_Plx"),
                        g(row, "VarFlag"), g(row, "MultFlag")])
    print("yazildi: %s  (%d satir, %d tanesinde SpType var)" % (OUT, len(t), n_sp))


if __name__ == "__main__":
    main()
