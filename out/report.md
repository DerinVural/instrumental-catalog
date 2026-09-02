# Instrumental Kadir Katalogu — Rapor

## Radyometri

| | |
|---|---|
| Aciklik capi | 27.547 mm |
| Aciklik alani | 5.9600 cm^2 |
| Odak uzakligi | 43.662 mm |
| **S0 (m_inst=0 -> e-/s)** | **6.8298e+06** |
| Vega vekili | A0V (pickles_9) |

> **Optik konvansiyon:** Radyometri, YI_SIL `star_tracker_params.h`
> degerlerinden turetilir (FOV 14.7 deg -> f=43.662 mm, f/1.585,
> A=5.960 cm^2). Gercek Zemax tasarimi farklidir (f=45.237 mm LUT
> basligindan, f/1.685) -> aciklik alani 5.661 cm^2. Ucus degerleri
> istenirse S0'i **0.947** ile carpin (=6.46782e+06 e-/s).
> Karar: sim ile tutarlilik esas alindi.

## Radyometri butcesi (S0 -> tespit)

Zincir: `m_inst -> e-/s -> poz -> e- -> SNR -> centroid -> attitude`. Poz 50 ms, okuma gurultusu 13 e-, pencere 9 px, sigma_PSF 0.70 px, plaka olcegi 25.84 as/px.

| m_inst | e-/s | e- (50 ms) | SNR | centroid px | tek yildiz (as) |
|---|---|---|---|---|---|
| 2.0 | 1.08e+06 | 54122 | 229.4 | 0.0031 | 0.079 |
| 4.0 | 1.72e+05 | 8578 | 85.3 | 0.0082 | 0.212 |
| 6.0 | 2.72e+04 | 1359 | 25.3 | 0.0277 | 0.715 |
| 7.0 | 1.08e+04 | 541 | 11.9 | 0.0588 | 1.519 |
| 8.0 | 4.31e+03 | 215 | 5.2 | 0.1355 | 3.502 |

### Tespit limiti

| poz | SNR>6 | SNR>10 | SNR>15 |
|---|---|---|---|
| 15 ms | 6.52 | 5.91 | 5.40 |
| 50 ms | 7.83 | 7.22 | 6.71 |

### Doyma

- Tam kuyu 13500 e-, tepe piksel payi 27% -> **m_inst < 2.09 olan yildizlar DOYAR** (50 ms pozda).
- Doymus yildizin merkezi kayar; algoritma tarafina 'bu listeye merkez guvenme' notu gerekir.

### Attitude (kaba kestirim)

- m_inst=6 yildizlar, 15 yildizli cozum: **0.184 arcsec** (tek yildiz 0.715 / sqrt(15)).
- NOT: yalnizca foton/okuma gurultusu. Smear, alan-bagimli PSF, kalibrasyon ve katalog hatalari DAHIL DEGIL.

## On-eleme ve renk payi [B1]

| | |
|---|---|
| m_inst limiti | 7.0 |
| COLOR_MARGIN | 4.50 |
| V on-eleme esigi | 11.50 (= limit + pay) |
| olculen delta_v min / max | -3.787 / 0.026 |
| **kalan pay** | **0.713 kadir** |

## Yildiz sayilari

- m_inst <= 7.0: **18993 yildiz**
- V gecersiz (atlandi): 0
- hesap basarisiz (atlandi): 0
- |dec| > 89.5 (pmRA clip): 1

## Spektrum eslestirme (flag dagilimi)

| flag | adet |
|---|---|
| sptype | 15435 |
| sptype_lum_plx | 3436 |
| bv | 121 |
| sptype+decclip | 1 |

## Renk terimi (m_inst - V)

- ortalama: -0.249 kadir
- min / max: -3.787 / 0.026
- **B-V ile korelasyon: -0.702** (negatif beklenir: kirmizi yildiz daha parlak)

## Cikti dosyalari

- `master_catalog.csv` — izlenebilirlik (tum alanlar)
- `hip_instrumental.txt` — YI_SIL ESKI format (5 kolon, MAG = m_inst)
- `katalog_instrumental.csv` — YENI `katalog_tam` formati (upstream f9e49c4). Vmag kolonu = **m_inst**, Hpmag kolonu = V (izlenebilirlik; fark = renk terimi). Kullanim:
  `YI_CATALOG_FILE=data/katalog_instrumental.csv ./YI_SIL --algo=pyramid`
  Hicbir dosyayi ezmez; konumlar J2000'e tasinmis birim vektor olarak gomulu.
