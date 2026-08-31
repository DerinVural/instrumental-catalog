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

    v = Vizier(columns=["HIP", "Vmag", "B-V", "SpType"], row_limit=-1)
    print("VizieR sorgusu: I/239/hip_main, Vmag < %.1f ..." % VMAG_CUT)
    res = v.query_constraints(catalog="I/239/hip_main", Vmag="<%.1f" % VMAG_CUT)
    if not len(res):
        raise RuntimeError("VizieR bos dondu")
    t = res[0]
    print("  gelen satir: %d" % len(t))

    os.makedirs(DATA, exist_ok=True)
    n_sp = 0
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hip", "vmag", "bv", "sptype"])
        for row in t:
            hip = row["HIP"]
            sp = str(row["SpType"]).strip() if row["SpType"] is not None else ""
            if sp and sp.lower() != "--":
                n_sp += 1
            w.writerow([hip, row["Vmag"], row["B-V"], sp])
    print("yazildi: %s  (%d satir, %d tanesinde SpType var)" % (OUT, len(t), n_sp))


if __name__ == "__main__":
    main()
