# İnstrümental Kadir Kataloğu — Tasarım Dokümanı

**Tarih:** 2026-08-28
**Repo:** `C:\Users\derin\instrumental-catalog` (yerel, push yok)
**Tüketici:** YI_SIL star tracker SIL simülasyonu (ileride uçuş yazılımı)

---

## 1. Amaç

Hipparcos kataloğundaki **V kadiri** yerine, bu star tracker'ın **kendi optik + dedektör
zincirinin gerçekte gördüğü parlaklığı** (instrümental kadir, `m_inst`) içeren bir yıldız
kataloğu üretmek.

**Problem:** V kadiri standart bir fotometrik banttır (insan gözü/Johnson V). Silisyum
dedektör (CMV4000) yakın-kızılötesine duyarlıdır ve optik geçirgenlik dalgaboyuna göre
değişir. Sonuç: kırmızı bir yıldız, V kadirinin öngördüğünden **daha fazla** elektron
üretir; mavi bir yıldız daha az. Sim/uçuş yazılımı V kadirini kullandığı sürece
öngörülen sinyal (SNR, pozlama, tespit eşiği) sistematik olarak yanlıştır.

**Çözüm:** Her yıldızın gerçek spektrumunu, sistemin spektral tepkisiyle konvolüe ederek
`m_inst` hesaplamak.

---

## 2. Kapsam

### Dahil
- Hipparcos → `m_inst` hesaplayan offline Python aracı
- Gerçek yıldız spektrumları (Pickles kütüphanesi, spektral tipe göre)
- CMV4000 Mono E5 QE(λ) + ZEMAX optik geçirgenlik T(λ) entegrasyonu
- İki çıktı: izlenebilir **master** dosya + YI_SIL'e **drop-in onboard** dosya
- Sim radyometrisiyle tutarlı `S0` (mag 0 → e⁻/s) hesabı

### Hariç (bilinçli)
- **Alan-bağımlı** vignetting/geçirgenlik → katalog on-axis'tir; alan etkisi render/FSW işi
- **Sıcaklık** bağımlılığı → katalog referans sıcaklıktadır; termal kayma render/FSW işi
- YI_SIL C kodunda değişiklik → çıktı mevcut formata birebir uyar
- Yıldızlararası kızarma (reddening) modellemesi → gözlenen B−V zaten içerir

---

## 3. Girdi verileri

| Veri | Kaynak | Durum |
|---|---|---|
| Astrometri + V + B−V | `star-sim-tu/data/catalogs/hipparcos.fits` (117955 yıldız: HIP, RA, DEC, VMAG, PMRA, PMDEC, BV, PLX) | ✅ mevcut |
| Spektral tip (SpType) | VizieR `I/239/hip_main`, HIP üzerinden JOIN | ✅ erişim doğrulandı |
| Yıldız spektrumları | Pickles UVILIB, `https://ssb.stsci.edu/cdbs/grid/pickles/` | ✅ erişim doğrulandı |
| QE(λ) | CMV4000 datasheet DS000728 v8-01, Figure 7, **Mono E5** eğrisi | ✅ digitize edildi |
| T_optics(λ) | `Desktop/datas/Transmission Data.txt` (ZEMAX, StarTracker_F1.68_v2_Final_TUZAY.zmx) | ✅ mevcut |
| Filtre | Yok (mono sensör) → `T_filter = 1` | — |

### Digitize edilmiş QE (CMV4000 Mono E5)

| λ (nm) | 400 | 450 | 500 | 520 | 550 | 600 | 650 | 700 | 750 | 800 | 850 | 900 | 950 | 1000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| QE % | 33 | 50 | 60 | 64 | 60 | 58 | 58 | 52 | 44 | 35 | 24 | 15 | 8 | 2 |

Datasheet metniyle çapraz doğrulama: "QE·FF %60 @ 550 nm with micro lenses" ✓

**DİKKAT:** Bu eğri **dolgu faktörünü (FF) içerir**. Radyometri zincirinde `fill_factor`
ayrıca çarpılmamalı — çift sayım olur.

### ZEMAX geçirgenliği (on-axis, Field Pos 0.0°)

| λ (µm) | 0.400 | 0.480 | 0.550 | 0.6328 | 0.720 | 0.800 | 0.900 |
|---|---|---|---|---|---|---|---|
| T | 0.0908 | 0.2866 | 0.3259 | 0.3417 | 0.3506 | 0.3542 | 0.3581 |

Dosya UTF-16 kodlu; parse edilip `data/t_optics_zemax.csv`'ye normalize edilir.
ZEMAX notu: *"Aperture, Fresnel, coating, vignetting, and internal transmittance effects
are considered."* → **açıklık etkisi dahil olabilir**; mutlak `S0` hesabında açıklık alanıyla
çift sayım riski var (bkz. §8 doğrulama).

---

## 4. Yöntem

### 4.1 Sistem spektral tepkisi
```
R_inst(λ) = QE(λ) × T_optics(λ)        [400–1000 nm, 1 nm'e interpole]
R_V(λ)    = Johnson V bant eğrisi       [referans bant]
```
Her iki tepki de **foton-sayım** modunda kullanılır (CCD böyle çalışır).

### 4.2 İnstrümental kadir (oran yöntemi)

```
                    ⎡  ∫F_yıldız·R_inst·(λ/hc)dλ  /  ∫F_yıldız·R_V·(λ/hc)dλ  ⎤
m_inst = V − 2.5·log⎢ ──────────────────────────────────────────────────────  ⎥
                    ⎣  ∫F_vega·R_inst·(λ/hc)dλ   /  ∫F_vega·R_V·(λ/hc)dλ    ⎦
```

**Neden bu form:** Mutlak akı kalibrasyonu gerekmez — yalnızca **göreli** spektrum şekli
ve iki tepki eğrisi yeter. Pickles spektrumları göreli olduğu için doğrudan kullanılabilir.
Tanım gereği Vega (A0V) için `m_inst = V` çıkar → sıfır noktası kendiliğinden tutarlı.

### 4.3 Mutlak ölçek (S0)
`S0` = `m_inst = 0` olan bir yıldızın ürettiği e⁻/s:
```
S0 = A_açıklık · ∫ F_vega,V=0(λ) · R_inst(λ) · (λ/hc) dλ
```
`F_vega,V=0`: Pickles A0V spektrumu, V bandında `F_λ(5500Å) = 3.63e-9 erg/s/cm²/Å`
olacak şekilde ölçeklenir (Bessell standart sıfır noktası).

Bu `S0`, YI_SIL `star_tracker_params.h`'deki kaba `phi_spectral × band_width × QE_avg`
formülünün yerine geçecek tutarlı değerdir (rapor edilir; C tarafına elle girilir).

### 4.4 Spektrum eşleme ve fallback
```
SpType parse başarılı  →  Pickles spektrumu (tip + ışıtma sınıfı)
SpType yok/bozuk       →  B−V ile en yakın anakol (V) tipi
B−V de yok             →  G2V (güneş benzeri)        [flag]
```
Her yıldız `flags` kolonunda hangi yolu izlediğiyle işaretlenir.

### 4.5 Birim ve konvansiyon dönüşümleri (kritik)

| Alan | Hipparcos | YI_SIL beklentisi | Dönüşüm |
|---|---|---|---|
| RA, Dec | derece | radyan | `× π/180` |
| pmRA | **μ_α\* = μ_α·cos δ** (mas/yr) | **μ_α** (rad/yr) | `÷ cos δ`, sonra `× 4.8481368e-9` |
| pmDec | μ_δ (mas/yr) | rad/yr | `× 4.8481368e-9` |
| Epok | J1991.25 (Hipparcos) | J2000 | PM ile J2000'e taşı |

**pmRA cos δ uyarısı:** YI_SIL `csv_parser.c:49`'da `corrected_ra = ra + Δt·pmra` uygular —
yani `cos δ` çarpanı **beklemez**. Hipparcos `pmRA` ise `cos δ` içerir. Bölünmezse yüksek
deklinasyonda konum hatası oluşur (δ=60°'de 2×). Kutuplarda `cos δ → 0` olduğundan
`|δ| > 89.5°` için bölme sınırlanır (`clip`) ve `flag` konur.

**Epok:** Çıktı **J2000** konumları içerir (mevcut `hygdata.txt` konvansiyonu). YI_SIL kendi
`dt_years` mekanizmasıyla güncel epoğa taşımaya devam eder.

---

## 5. Mimari

```
instrumental-catalog/
├── src/
│   ├── response.py       # QE + T_optics → R_inst(λ); Johnson V → R_V(λ)
│   ├── spectra.py        # Pickles yükleme, SpType→spektrum eşleme, fallback zinciri
│   ├── photometry.py     # integrasyon, Vega ZP, m_inst, e⁻/s, S0
│   └── build_catalog.py  # orkestrasyon: veri çek → hesapla → çıktılar + rapor
├── data/
│   ├── qe_cmv4000_e5.csv      # digitize edilmiş QE (bkz. §3)
│   ├── t_optics_zemax.csv     # ZEMAX txt'den parse (bkz. §3)
│   └── johnson_v.csv          # Johnson-Cousins V bant eğrisi (Bessell 1990,
│                              #   PASP 102:1181, Tablo 2) — repoya gömülü sabit tablo
├── tests/                # birim testler
├── out/                  # üretilen kataloglar + rapor
└── docs/superpowers/specs/
```

**Modül sınırları:**

| Modül | Girdi | Çıktı | Bağımlılık |
|---|---|---|---|
| `response.py` | CSV eğrileri | `R_inst(λ)`, `R_V(λ)` dizileri | numpy |
| `spectra.py` | SpType / B−V | göreli spektrum `F(λ)` + flag | astropy, numpy |
| `photometry.py` | `F(λ)`, `R(λ)`, V | `m_inst`, `e⁻/s`, `S0` | numpy |
| `build_catalog.py` | — | 2 dosya + rapor | yukarıdakiler |

Her modül tek sorumluluk taşır ve bağımsız test edilebilir; `build_catalog.py` yalnızca
akışı yönetir.

---

## 6. Çıktılar

### 6.1 Master (izlenebilirlik) — `out/master_catalog.csv`
```
hip, ra_deg, dec_deg, pmra_mas, pmdec_mas, vmag, bv, sptype,
pickles_file, m_inst, e_per_s, delta_v (=m_inst−V), flags
```
Amaç: doğrulama, hata ayıklama, uçuş konfigürasyon yönetimi, "bu değer neden böyle" sorusu.

### 6.2 Onboard (YI_SIL drop-in) — `out/hip_instrumental.txt`
```
MAG	RARAD	DECRAD	PMRARAD	PMDECRAD        ← mevcut hygdata.txt ile BİREBİR aynı başlık/format
8.73	1.5700e-05	0.019007	-2.52e-08	-9.11e-09
 ↑ m_inst (V DEĞİL)
```
- Tab ayraçlı, radyan, J2000
- **MAG kolonu artık `m_inst`** — YI_SIL'de kod değişikliği gerekmez, yalnızca dosya yolu
- Kadir limiti: `m_inst ≤ 7.0` (YI_SIL kendi `mag_limit`'iyle daha da süzer)

### 6.3 Rapor — `out/report.md`
- Hesaplanan `S0` (e⁻/s, mag 0) → C tarafına girilecek değer
- Yıldız sayıları (toplam, limit sonrası, fallback kullanılanlar)
- `m_inst − V` vs `B−V` istatistiği + grafik
- Doğrulama tablosu (bilinen yıldızlar)

---

## 7. Hata durumları

| Durum | Davranış |
|---|---|
| VizieR erişilemiyor | Hata ver, dur (kısmi katalog üretme) |
| SpType eksik/parse edilemez | B−V fallback, flag |
| B−V de eksik | G2V fallback, flag |
| Pickles dosyası eksik | Hata ver, dur (eksik spektrum sessizce atlanmaz) |
| `\|δ\| > 89.5°` (cos δ ≈ 0) | pmRA bölmesi sınırlanır, flag |
| Negatif/anlamsız V | Yıldız atlanır, rapora yazılır |

İlke: **sessiz veri kaybı yok** — atlanan/fallback uygulanan her yıldız raporlanır.

---

## 8. Doğrulama

| # | Test | Beklenen |
|---|---|---|
| 1 | **Vega yuvarlak-gidiş** | A0V spektrumu → `m_inst = V` (±0.001 kadir) |
| 2 | **Renk terimi yönü** | `m_inst − V` ile `B−V` arasında monoton **negatif** ilişki (kırmızı yıldız daha parlak) |
| 3 | **Bilinen yıldızlar** | Betelgeuse (M2Iab, B−V≈1.85): belirgin negatif Δ. Rigel (B8Ia, B−V≈−0.03): Δ≈0 |
| 4 | **S0 mertebesi** | Yeni `S0`, mevcut kaba `S0` ile aynı mertebede (>10× sapma → açıklık çift-sayımı şüphesi, §3 notu) |
| 5 | **PM konvansiyonu** | Yüksek-δ, yüksek-pm yıldızlarda (ör. Barnard) J2000 konumu bağımsız kaynakla ±50 mas uyum |
| 6 | **Format uyumu** | `hip_instrumental.txt`, YI_SIL `Parse_StarCatalog` ile yüklenir; yıldız sayısı beklenen |
| 7 | **Uçtan uca sim** | YI_SIL yeni katalogla koşar; acquisition + tracking sağlıklı, SNR makul |
| 8 | **Birim testler** | interpolasyon, spektrum ölçekleme (V yuvarlak-gidiş), SpType parse, fallback zinciri, birim dönüşümleri |

---

## 9. Bilinen sınırlamalar

1. **QE eğrisi grafikten digitize** — datasheet sayısal veri vermiyor; ±%2-3 okuma belirsizliği.
   Girişim salınımları (fringing) düzleştirildi.
2. **ZEMAX geçirgenliği 7 dalgaboyunda** — aradaki değerler interpole, 900 nm üstü sabitlenir
   (clamp). 1000 nm'e kadar QE zaten %2'ye düşüyor, etki küçük.
3. **Pickles ışıtma sınıfı kapsamı** — kütüphanede her tip/sınıf kombinasyonu yok; en yakın
   eşleşme kullanılır.
4. **Tek referans sıcaklık** — QE ve optik sıcaklıkla değişir; katalog tek noktada geçerlidir.
5. **Değişken yıldızlar** — Hipparcos ortalama V değeri kullanılır; genlik modellenmez.
6. **Çift yıldızlar** — Hipparcos çözülmüş bileşenleri ayrı satır olarak gelir; birleşik
   fotometri modellenmez.

---

## 10. Başarı ölçütü

- `out/hip_instrumental.txt` YI_SIL'de sorunsuz çalışır (Test 6, 7)
- Vega yuvarlak-gidiş ve renk terimi yönü doğrulanır (Test 1, 2)
- Master dosya her yıldız için hesabın izini taşır (SpType, kullanılan spektrum, flag)
- Rapor, C tarafına girilecek tutarlı `S0` değerini verir
