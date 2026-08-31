# instrumental-catalog

Star tracker için **instrümental kadir** (instrumental magnitude) yıldız kataloğu üreteci.

Hipparcos kataloğundaki standart **V kadiri** yerine, bu star tracker'ın **kendi optik +
dedektör zincirinin gerçekte gördüğü parlaklığı** hesaplar.

> **Neden:** V kadiri Johnson V bandında tanımlıdır. Silisyum dedektör (CMV4000) yakın
> kızılötesine duyarlıdır ve optik geçirgenlik dalgaboyuna göre değişir. Sonuç: kırmızı
> yıldızlar V'nin öngördüğünden **daha fazla** elektron üretir. V kadiri kullanıldığı
> sürece öngörülen sinyal (SNR, pozlama, tespit eşiği) sistematik olarak yanlıştır.

## Yöntem

```
R_inst(λ) = QE(λ) × T_optics(λ)          CMV4000 Mono E5 × ZEMAX (on-axis)

                    ⎡ (∫F·R_inst·λdλ / ∫F·R_V·λdλ)_yıldız ⎤
m_inst = V − 2.5·log⎢ ─────────────────────────────────── ⎥
                    ⎣ (∫F·R_inst·λdλ / ∫F·R_V·λdλ)_Vega   ⎦
```

Oran yöntemi olduğu için **mutlak akı kalibrasyonu gerekmez**; Vega (A0V) tanım gereği
`m_inst = V` verir. Yıldız spektrumları **Pickles UVILIB** kütüphanesinden, her yıldızın
spektral tipine göre seçilir.

## Sonuçlar

| | |
|---|---|
| Yıldız sayısı | 19807 (m_inst ≤ 7.0) |
| S0 (m_inst=0 → e⁻/s) | 3.65×10⁶ |
| Ortalama renk terimi (m_inst − V) | −0.296 kadir |
| corr(B−V, m_inst−V) | −0.742 |

| Yıldız | Tip | V | m_inst | Δ |
|---|---|---|---|---|
| Vega | A0V | 0.03 | 0.030 | 0.000 (referans) |
| Rigel | B8Ia | 0.18 | 0.179 | −0.001 |
| Arcturus | K2III | −0.05 | −0.444 | −0.394 |
| Betelgeuse | M2Ib | 0.45 | −0.588 | −1.038 |

## Kullanım

```bash
python src/prepare_data.py     # ZEMAX/QE/Johnson V tablolarını hazırla
python src/fetch_sptype.py     # VizieR'den SpType çek
python src/build_catalog.py    # katalogları üret
python src/validate.py         # doğrulama testleri
```

Pickles kütüphanesi ilk çalıştırmada `data/pickles/` altına indirilir
(<https://ssb.stsci.edu/cdbs/grid/pickles/dat_uvi/>).

## Çıktılar (`out/`)

| Dosya | İçerik |
|---|---|
| `hip_instrumental.txt` | Onboard katalog — 5 kolon (`MAG RARAD DECRAD PMRARAD PMDECRAD`), **MAG = m_inst** |
| `master_catalog.csv` | İzlenebilirlik — HIP, V, B−V, SpType, kullanılan Pickles spektrumu, m_inst, Δ, e⁻/s, flag |
| `report.md` | S0, sayımlar, flag dağılımı, istatistikler |

Onboard dosya, YI_SIL'in mevcut katalog formatıyla **birebir uyumludur** (drop-in).

## Modüller

| Modül | Sorumluluk |
|---|---|
| `response.py` | QE + T_optics → R_inst(λ); Johnson V → R_V(λ); açıklık alanı |
| `spectra.py` | Pickles kütüphanesi, SpType ayrıştırma/eşleme, fallback zinciri |
| `photometry.py` | Foton-sayım integrali, Vega sıfır noktası, m_inst, S0 |
| `build_catalog.py` | Orkestrasyon → 2 katalog + rapor |
| `validate.py` | İç tutarlılık testleri (7/7 geçiyor) |
| `validate_johnson13.py` | Bağımsız doğruluk — ölçülen 13-renk fotometrisine karşı |
| `cv_compare.py` | Çapraz doğrulanmış B−V karşılaştırması |
| `digitize_qe.py` | QE eğrisini datasheet PDF vektöründen çıkarır |

## Optik konvansiyon

Radyometri (`S0`, açıklık alanı) **YI_SIL `star_tracker_params.h`** değerlerinden türetilir:
FOV 14.7° → f = 43.662 mm, f/1.585, A = 5.960 cm².

Gerçek Zemax tasarımı bundan farklıdır (PSF LUT başlığı: f = 45.237 mm; f/1.685 →
A = 5.661 cm²). Uçuş değerleri istenirse `S0` **0.947** ile çarpılır. Karar: sim ile
tutarlılık esas alındı — kadirler (`m_inst`) bu seçimden **etkilenmez**, yalnızca mutlak
ölçek (`S0`) etkilenir.


## Bağımsız doğrulama (Johnson & Mitchell 13-renk fotometrisi)

`validate_johnson13.py`, VizieR **II/84** (Johnson+ 1975, 1380 parlak yıldız) veri setini
kullanır: her yıldızın **ölçülen** 13-renk fotometrisinden (337–1104 nm) gerçek SED'i
yeniden kurar, sistem tepkisiyle konvolüe eder ve "gerçek" ΔM üretir. Bu, Lu & Wu (2019)
makalesinin doğrulama yöntemidir.

| Örneklem | Bizim (sentetik) | B−V lineer | B−V 4. derece |
|---|---|---|---|
| Tüm eşleşenler (1192) | **0.125** | 0.219 | 0.200 |
| Temiz (530; değişken+çift elendi) | 0.045 | 0.038 | **0.021** |

*(std, kadir. B−V için 5-kat çapraz doğrulama — `cv_compare.py`)*

**Yorum (dürüstlük notu):**

- **Tüm örneklemde** sentetik yöntem net üstün (~1.7×). B−V, düzensiz yıldızlarda
  (Mira değişkenleri, çiftler, S-tipi) çöker — makalenin ana iddiası doğrulandı.
- **Temiz örneklemde** B−V polinomu bizi geçer (0.021 vs 0.045). "Uslu" yıldızlarda B−V
  mükemmel bir göstergedir; bizim artık hatamız **Pickles tip kuantalanmasından** gelir
  (131 tip sonlu; K3.5III yıldıza en yakın K3III verilir).
- **Ama** B−V polinomu kalibrasyon verisi gerektirir — o veriyi üretmenin yolu zaten
  sentetik fotometridir. Yeni bir optik/dedektör için B−V katsayıları **yoktur**;
  sentetik yöntem sıfırdan çalışır (makalenin asıl argümanı).

**İndirilemez hata kaynakları:** Mira değişkenleri (13-renk tek epok, Hipparcos V ortalama),
S/C-tipi yıldızlar (Pickles'ta yok), tip kuantalanması.

## Kapsam sınırı (uçuş)

Katalog **on-axis** ve **referans sıcaklıkta** geçerlidir. Alan-bağımlı vignetting ve
termal kayma katalogda **değil**, render/FSW tarafında uygulanır — bir yıldızın katalog
kadiri, dedektörde nereye düştüğüne bağlı olamaz.

Tasarım ayrıntısı: `docs/superpowers/specs/2026-08-28-instrumental-magnitude-catalog-design.md`

## Veri kaynakları

| Veri | Kaynak |
|---|---|
| Astrometri, V, B−V | Hipparcos (ESA, 1997) |
| Spektral tip | VizieR `I/239/hip_main` |
| Yıldız spektrumları | Pickles (1998) UVILIB, STScI |
| QE(λ) | ams-OSRAM CMV4000 datasheet DS000728 v8-01, Fig. 7 (Mono E5) |
| T_optics(λ) | ZEMAX optik tasarım analizi (on-axis) |
| Johnson V bandı | Bessell (1990), PASP 102:1181 |
