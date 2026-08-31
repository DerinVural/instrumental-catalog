"""Ortak esikler — TEK KAYNAK.

[B1/N2] Onceden iki ayri kesme vardi ve tutarsizdi:
    fetch_sptype.py : VMAG_CUT = 8.0
    build_catalog.py: vmag > 8.5 -> ele
Sonuc: (a) 8.5 < V <= ~10.8 arasindaki KIRMIZI yildizlar m_inst<=7.0'i
sagladiklari halde katalogda yoktu (eksik katalog -> sahnede karsiligi
olmayan tespit -> sahte yildiz), (b) 8.0 < V <= 8.5 arasi SpType alamayip
sessizce B-V dalina dusuyordu (parlakliga bagli sistematik kalite dususu).

Artik ikisi de buradan turetilir.
"""

MAG_LIMIT = 7.0        # katalog m_inst ust siniri (master + onboard)

# Renk terimi payi: m_inst = V + C ve C NEGATIF olabilir (kirmizi yildiz).
# Bir yildiz m_inst <= MAG_LIMIT kosulunu V = MAG_LIMIT + |C| kadar sonukken
# de saglayabilir. Olculen en negatif C ~ -3.83 (gec-M devleri, N8 ile
# fiziksel oldugu dogrulandi) -> emniyet payiyla 4.5.
COLOR_MARGIN = 4.5

# On-eleme esigi: bu degerden sonuk yildizlar hic hesaplanmaz.
VMAG_PREFILTER = MAG_LIMIT + COLOR_MARGIN      # 11.5
