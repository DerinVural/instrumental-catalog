"""Fotometri: foton-sayim integrasyonu, Vega sifir noktasi, m_inst ve S0.

Formul (oran yontemi — mutlak aki kalibrasyonu GEREKMEZ):

              [ N_inst(yildiz) / N_V(yildiz) ]
m_inst = V - 2.5 log10 ------------------------------
              [ N_inst(vega)  / N_V(vega)  ]

N_X = \\int F(lam) R_X(lam) lam dlam   (foton sayim; hc sabiti oranlarda sadelesir)

Tanim geregi Vega (A0V) icin m_inst = V cikar.
"""
import numpy as np

# Fiziksel sabitler (CGS)
HC_ERG_CM = 1.98644586e-16          # h*c [erg*cm]
ANG_TO_CM = 1.0e-8

# Johnson V sifir noktasi: V=0 yildizin V bandindaki akisi (Bessell 1998)
F_LAMBDA_V0 = 3.63e-9               # erg/s/cm^2/Angstrom


def photon_count_integral(lam_nm, flux, response_grid, grid_nm):
    """N = \\int F(lam) R(lam) lam dlam  (goreli foton sayimi).

    lam_nm/flux: kaynak spektrum; response_grid: ortak izgarada tepki egrisi.
    """
    f = np.interp(grid_nm, lam_nm, flux, left=0.0, right=0.0)
    integrand = f * response_grid * grid_nm
    return np.trapezoid(integrand, grid_nm) if hasattr(np, "trapezoid") \
        else np.trapz(integrand, grid_nm)


class Photometry:
    """Sabit tepki egrileriyle m_inst/S0 hesaplayici."""

    def __init__(self, grid_nm, r_inst, r_v, vega_lam_nm, vega_flux, area_cm2):
        self.grid = grid_nm
        self.r_inst = r_inst
        self.r_v = r_v
        self.area_cm2 = area_cm2

        # Vega referans oranlari (m_inst = V normalizasyonu)
        self.n_inst_vega = photon_count_integral(vega_lam_nm, vega_flux, r_inst, grid_nm)
        self.n_v_vega = photon_count_integral(vega_lam_nm, vega_flux, r_v, grid_nm)
        if self.n_inst_vega <= 0 or self.n_v_vega <= 0:
            raise ValueError("Vega integralleri pozitif degil")
        self.vega_ratio = self.n_inst_vega / self.n_v_vega

        self.vega_lam = vega_lam_nm
        self.vega_flux = vega_flux
        self._s0 = None

    def m_inst(self, lam_nm, flux, vmag):
        """Yildizin instrumental kadiri."""
        n_inst = photon_count_integral(lam_nm, flux, self.r_inst, self.grid)
        n_v = photon_count_integral(lam_nm, flux, self.r_v, self.grid)
        if n_inst <= 0 or n_v <= 0:
            return None
        ratio = (n_inst / n_v) / self.vega_ratio
        return float(vmag - 2.5 * np.log10(ratio))

    def s0_electrons_per_s(self):
        """m_inst = 0 olan bir yildizin urettigi e-/s (mutlak olcek).

        Vega spektrumu, V bandinda ortalama aki yogunlugu F_LAMBDA_V0 (V=0) olacak
        sekilde olceklenir; ardindan R_inst ile foton integrali alinir.
        """
        if self._s0 is not None:
            return self._s0

        grid_a = self.grid * 10.0                       # nm -> Angstrom
        f = np.interp(self.grid, self.vega_lam, self.vega_flux, left=0.0, right=0.0)

        # V bandinda foton-agirlikli ortalama aki yogunlugu
        w = self.r_v * self.grid
        trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
        mean_f_v = trapz(f * w, self.grid) / trapz(w, self.grid)
        scale = F_LAMBDA_V0 / mean_f_v                  # V=0 olacak sekilde olcekle

        f_abs = f * scale                               # erg/s/cm^2/Angstrom
        lam_cm = grid_a * ANG_TO_CM
        photons = f_abs * lam_cm / HC_ERG_CM            # foton/s/cm^2/Angstrom
        n = trapz(photons * self.r_inst, grid_a)        # foton/s/cm^2 (QE dahil -> e-)
        self._s0 = float(n * self.area_cm2)
        return self._s0
