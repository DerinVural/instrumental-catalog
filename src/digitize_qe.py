"""CMV4000 Mono E5 QE egrisini datasheet PDF'inin VEKTOR verisinden digitize eder.

Goz karariyla okuma yerine, Figure 7'deki siyah polyline'in gercek koordinatlari
eksen etiketleriyle kalibre edilerek cikarilir.
Dogrulama capalari (datasheet METNI): 550 nm ~ %60, 900 nm ~ %8.
"""
import os
import numpy as np
import pymupdf

PDF = r"C:\Users\derin\Desktop\datas\CMV4000_DS000728_8-01.pdf"
PAGE = 19            # 0-indexed -> sayfa 20, Figure 7
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "data", "qe_cmv4000_e5.csv")

# Eksen kalibrasyonu: etiket metin kutularinin merkezi (PDF koordinati)
# y ekseni: QE [%]      x ekseni: dalgaboyu [nm]
Y_CAL = [(355.6, 0.0), (171.5, 70.0)]        # '0' ve '70' etiket merkezleri
X_CAL = [(224.0, 400.0), (441.2, 1000.0)]    # '400' ve '1000' tik merkezleri


def main():
    pg = pymupdf.open(PDF)[PAGE]
    black = None
    for g in pg.get_drawings():
        if g.get("color") != (0.0, 0.0, 0.0):
            continue
        pts = []
        for it in g["items"]:
            if it[0] == "l":
                pts += [it[1], it[2]]
            elif it[0] == "c":
                pts += [it[1], it[4]]
        if len(pts) > 100:                    # egri (izgara/cerceve degil)
            black = pts
            break
    if black is None:
        raise RuntimeError("Mono E5 (siyah) egrisi bulunamadi")

    (yp0, yv0), (yp1, yv1) = Y_CAL
    (xp0, xv0), (xp1, xv1) = X_CAL
    lam = np.array([xv0 + (p.x - xp0) * (xv1 - xv0) / (xp1 - xp0) for p in black])
    qe = np.array([yv0 + (p.y - yp0) * (yv1 - yv0) / (yp1 - yp0) for p in black]) / 100.0

    order = np.argsort(lam)
    lam, qe = lam[order], qe[order]

    # Girisim salinimlarini (fringing) duzlestir: 10 nm izgaraya ortalama
    grid = np.arange(400, 1001, 25.0)
    sm = []
    for g0 in grid:
        sel = (lam >= g0 - 12.5) & (lam < g0 + 12.5)
        sm.append(qe[sel].mean() if sel.any() else np.interp(g0, lam, qe))
    sm = np.clip(np.array(sm), 0.0, 1.0)

    print("digitize edildi: %d ham nokta -> %d izgara noktasi" % (len(lam), len(grid)))
    print("DOGRULAMA (datasheet metin capalari):")
    for target, exp in ((550, 0.60), (900, 0.08)):
        got = float(np.interp(target, grid, sm))
        ok = abs(got - exp) < 0.03
        print("   %4d nm: okunan %.3f | beklenen ~%.2f  -> %s"
              % (target, got, exp, "UYUMLU" if ok else "SAPMA"))

    with open(OUT, "w", newline="") as f:
        f.write("lambda_nm,QE\n")
        for g0, v in zip(grid, sm):
            f.write("%.0f,%.4f\n" % (g0, v))
    print("yazildi: %s" % OUT)
    print("  " + ", ".join("%.0f:%.2f" % (g0, v) for g0, v in zip(grid, sm)))


if __name__ == "__main__":
    main()
