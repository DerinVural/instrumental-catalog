# Instrumental Kadir Katalogu — Rapor

## Radyometri

| | |
|---|---|
| Aciklik capi | 27.547 mm |
| Aciklik alani | 5.9600 cm^2 |
| Odak uzakligi | 43.662 mm |
| **S0 (m_inst=0 -> e-/s)** | **3.55406e+06** |
| Vega vekili | A0V (pickles_9) |

> **Optik konvansiyon:** Radyometri, YI_SIL `star_tracker_params.h` degerlerinden
> turetilir (FOV 14.7 deg -> f=43.662 mm, f/1.585, A=5.960 cm^2). Gercek Zemax
> tasarimi farklidir (f=45.237 mm LUT basligindan, f/1.685) -> A=5.661 cm^2.
> Ucus degerleri istenirse S0'i **0.947** ile carpin (= 3.36584e+06 e-/s).
> Karar: sim ile tutarlilik esas alindi.

## Yildiz sayilari

- m_inst <= 7.0: **19699 yildiz**
- V gecersiz (atlandi): 0
- hesap basarisiz (atlandi): 0
- |dec| > 89.5 (pmRA clip): 1

## Spektrum eslestirme (flag dagilimi)

| flag | adet |
|---|---|
| sptype | 18748 |
| bv | 950 |
| sptype+decclip | 1 |

## Renk terimi (m_inst - V)

- ortalama: -0.293 kadir
- min / max: -3.832 / 0.030
- **B-V ile korelasyon: -0.738** (negatif beklenir: kirmizi yildiz daha parlak)

## Cikti dosyalari

- `master_catalog.csv` — izlenebilirlik (tum alanlar)
- `hip_instrumental.txt` — YI_SIL drop-in (MAG kolonu = m_inst)
