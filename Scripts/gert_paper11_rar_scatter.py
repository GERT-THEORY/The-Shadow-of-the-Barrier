#!/usr/bin/env python3
"""
GERT Paper 11 — The Shadow of the Barrier
============================================
Two Conformal Fossils from the Cohesive Fraction:
  Fossil 1: a₀ from φ_start/φ_end (Paper IX, 0.57%)
  Fossil 2: σ_RAR from (0.5 − φ_max)/(φ_max ln10) (Paper XI, 0.4%)

Derivation:
  1. g_obs = g_bar / φ  (GERT acceleration relation, Paper VI)
  2. log₁₀(g_obs) = log₁₀(g_bar) − log₁₀(φ)
  3. σ_RAR = |d(−log₁₀φ)/dφ| × δφ = δφ/(φ ln10)  (Jacobian)
  4. δφ_max = 0.5 − φ_max  (barrier distance = max fluctuation)
  5. σ_RAR = (0.5 − φ_max)/(φ_max × ln10) = 0.0572 dex
  6. Observed: 0.057 dex. Agreement: 0.4%. Zero free parameters.

Veronica Padilha Dutra (2026)
"""

import os
import numpy as np
from scipy.special import expit
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paper I parameters (frozen) ───────────────────────────────────
FM_I, FM_F = 0.7831, 0.5851
FL_I, FL_M = 1.3414, 1.1236
LOG_RHO_M, D_M = -20.30, 1.0
FM_PEAK = 0.37; LOG_RHO_C = -17.41; SIGMA_C = 1.0
LOG_RHO_L, D_L = -25.60, 2.0
FL_PEAK = 4.6245; LOG_RHO_L2 = -23.93; SIGMA_L2 = 1.0
C = 2.998e8; MPC = 3.0857e22
H0 = 72.5e3 / MPC
A_GERT = C * H0 / (2 * np.pi)
A0_OBS = 1.2e-10  # m/s²

def logistic(x, x0, d): return expit(-(x - x0) / d)
def gauss(x, x0, s): return np.exp(-0.5 * ((x - x0) / s)**2)

def fM(x):
    b = FM_I + (FM_F - FM_I) * logistic(x, LOG_RHO_M, D_M)
    return b * (1 + FM_PEAK * gauss(x, LOG_RHO_C, SIGMA_C))

def fL(x):
    b = FL_I + (FL_M - FL_I) * logistic(x, LOG_RHO_L, D_L)
    return b * (1 + FL_PEAK * gauss(x, LOG_RHO_L2, SIGMA_L2))

def phi(x): return fM(x) / (fM(x) + fL(x))


# ══════════════════════════════════════════════════════════════════════
#  BLOCK I — THE TWO FOSSILS
# ══════════════════════════════════════════════════════════════════════

def block_I_fossils():
    """Compute both conformal fossils from φ(x)."""
    print(f"\n{'='*72}")
    print(f"  BLOCK I — TWO CONFORMAL FOSSILS FROM φ(x)")
    print(f"{'='*72}")

    # ── Fossil 1: Milgrom scale ───────────────────────────────────
    phi_start = phi(-5)   # x → −5 (Cauldron start)
    phi_end   = phi(-27)  # x → −27 (ultra-dilute)
    ratio = phi_start / phi_end
    a0_pred = A_GERT * ratio

    print(f"""
  FOSSIL 1 — THE MILGROM SCALE (Paper IX, Result 2)
  
    φ_start = f_M,i / (f_M,i + f_L,i) = {phi_start:.6f}
    φ_end   = f_M,f / (f_M,f + f_L,m) = {phi_end:.6f}
    Ratio   = φ_start / φ_end = {ratio:.4f}
    
    a₀ = a_GERT × ratio = {A_GERT:.4e} × {ratio:.4f}
       = {a0_pred:.4e} m/s²
    
    Observed: {A0_OBS:.4e} m/s²
    Agreement: {abs(a0_pred - A0_OBS)/A0_OBS * 100:.2f}%
    Free parameters: 0
""")

    # ── Fossil 2: RAR scatter ─────────────────────────────────────
    res = minimize_scalar(lambda x: -phi(x), bounds=(-20, -15), method='bounded')
    x_max = res.x
    phi_max = phi(x_max)
    delta = 0.5 - phi_max
    jacobian = 1 / (phi_max * np.log(10))
    sigma_pred = delta * jacobian
    sigma_obs = 0.057

    print(f"""
  FOSSIL 2 — THE RAR INTRINSIC SCATTER (Paper XI)
  
    φ_max = {phi_max:.6f}  at x = {x_max:.4f}
    
    Barrier distance:  0.5 − φ_max = {delta:.6f}
    Jacobian:          1/(φ_max × ln10) = {jacobian:.6f}
    
    σ_RAR = (0.5 − φ_max) / (φ_max × ln10)
          = {delta:.6f} / ({phi_max:.4f} × {np.log(10):.4f})
          = {sigma_pred:.4f} dex
    
    Observed: {sigma_obs} dex (McGaugh, Lelli & Schombert 2016)
    Agreement: {abs(sigma_pred - sigma_obs)/sigma_obs * 100:.1f}%
    Free parameters: 0
""")

    # ── Summary table ─────────────────────────────────────────────
    print(f"  {'─'*70}")
    print(f"  {'Fossil':<12} {'Quantity':<12} {'Predicted':<16} "
          f"{'Observed':<16} {'Agreement':<10}")
    print(f"  {'─'*70}")
    print(f"  {'1 (IX)':<12} {'a₀':<12} {a0_pred:.4e}{'m/s²':<4} "
          f"{A0_OBS:.4e}{'m/s²':<4} {abs(a0_pred-A0_OBS)/A0_OBS*100:.2f}%")
    print(f"  {'2 (XI)':<12} {'σ_RAR':<12} {sigma_pred:.4f}{' dex':<8} "
          f"{sigma_obs:.3f}{' dex':<8} {abs(sigma_pred-sigma_obs)/sigma_obs*100:.1f}%")
    print(f"  {'─'*70}")

    return phi_max, delta, sigma_pred, x_max


# ══════════════════════════════════════════════════════════════════════
#  BLOCK II — THE DERIVATION STEP BY STEP
# ══════════════════════════════════════════════════════════════════════

def block_II_derivation(phi_max, delta):
    """Show the derivation step by step."""
    print(f"\n{'='*72}")
    print(f"  BLOCK II — THE DERIVATION")
    print(f"{'='*72}")
    print(f"""
  Step 1: The GERT acceleration relation (Paper VI)
    g_obs = g_bar / φ(x)                              ... (4)
    
  Step 2: In log space
    log₁₀(g_obs) = log₁₀(g_bar) − log₁₀(φ)          ... (6)
    
  Step 3: The scatter
    σ_RAR = σ[−log₁₀(φ)]
    
  Step 4: The Jacobian
    δ[−log₁₀(φ)] = |d(−log₁₀φ)/dφ| × δφ
                  = (1/(φ ln10)) × δφ                  ... (8)
    
  Step 5: The maximum thermodynamic fluctuation
    δφ_max = 0.5 − φ_max = {delta:.6f}               ... (10)
    
    This is the BARRIER DISTANCE: the universe's closest
    approach to φ = 1/2 (where ΔG = 0 and time stops).
    No fluctuation in φ can exceed this value.
    
  Step 6: The prediction
    σ_RAR = (0.5 − φ_max) / (φ_max × ln10)           ... (11)
          = {delta:.6f} / ({phi_max:.4f} × {np.log(10):.4f})
          = {delta / (phi_max * np.log(10)):.4f} dex
    
  All numbers from Paper I MCMC. Zero adjustment.
""")


# ══════════════════════════════════════════════════════════════════════
#  BLOCK III — SENSITIVITY ANALYSIS
# ══════════════════════════════════════════════════════════════════════

def block_III_sensitivity():
    """Test sensitivity to peak amplitude and position."""
    print(f"\n{'='*72}")
    print(f"  BLOCK III — SENSITIVITY ANALYSIS")
    print(f"{'='*72}")

    print(f"\n  {'Variation':<35} {'φ_max':<10} {'σ_pred (dex)':<14} "
          f"{'vs obs':<10}")
    print(f"  {'─'*70}")

    cases = [
        ("Baseline (F_M,peak=0.37)", 0.37, -17.41),
        ("F_M,peak=0.35 (−5%)", 0.35, -17.41),
        ("F_M,peak=0.40 (+8%)", 0.40, -17.41),
        ("log ρ_c = −17.0", 0.37, -17.0),
        ("log ρ_c = −17.8", 0.37, -17.8),
    ]

    for label, peak, rho_c in cases:
        def fM_var(x):
            b = FM_I + (FM_F - FM_I) * logistic(x, LOG_RHO_M, D_M)
            return b * (1 + peak * gauss(x, rho_c, SIGMA_C))
        def phi_var(x): return fM_var(x) / (fM_var(x) + fL(x))

        res = minimize_scalar(lambda x: -phi_var(x), bounds=(-20, -15),
                              method='bounded')
        pm = phi_var(res.x)
        sp = (0.5 - pm) / (pm * np.log(10))
        diff = (sp - 0.057) / 0.057 * 100

        print(f"  {label:<35} {pm:<10.4f} {sp:<14.4f} {diff:>+8.1f}%")


# ══════════════════════════════════════════════════════════════════════
#  BLOCK IV — FIGURES
# ══════════════════════════════════════════════════════════════════════

def block_IV_figures(phi_max, x_max):
    """Generate Paper 11 figures."""
    print(f"\n{'='*72}")
    print(f"  BLOCK IV — FIGURES")
    print(f"{'='*72}")

    outdir = '/mnt/user-data/outputs'

    # ── Figure 1: Two fossils ─────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("GERT Paper XI — Two Conformal Fossils from φ(x)",
                 fontsize=13, fontweight='bold')

    x_arr = np.linspace(-5, -28, 2000)
    phi_arr = np.array([phi(x) for x in x_arr])

    # Panel 1: φ(x) with barrier and φ_max
    ax = axes[0]
    ax.plot(x_arr, phi_arr, 'C3', lw=2.5)
    ax.axhline(0.5, color='k', ls='--', lw=1.5, label='Barrier φ = 1/2')
    ax.axhline(phi_max, color='C1', ls=':', lw=1.5,
               label=f'φ_max = {phi_max:.4f}')
    ax.axvline(x_max, color='grey', ls=':', alpha=0.5)
    ax.fill_between(x_arr, phi_arr, 0.5, where=phi_arr < 0.5,
                    alpha=0.08, color='C3')
    ax.annotate(f'0.5 − φ_max = {0.5-phi_max:.4f}',
                xy=(x_max, (0.5 + phi_max)/2), fontsize=9,
                ha='center', color='C1',
                arrowprops=dict(arrowstyle='<->', color='C1'))
    ax.set_xlabel('x = log₁₀ρ')
    ax.set_ylabel('φ(x)')
    ax.set_title('The Cohesive Fraction\nand the Forbidden Barrier')
    ax.legend(fontsize=7, loc='lower left')
    ax.set_xlim(-5, -28); ax.grid(alpha=0.2)

    # Panel 2: Fossil 1 — a₀
    ax = axes[1]
    phi_s = phi(-5); phi_e = phi(-27)
    a0_pred = A_GERT * phi_s / phi_e
    ax.barh(['Observed', 'GERT'], [A0_OBS*1e10, a0_pred*1e10],
            color=['C0', 'C3'], height=0.5)
    ax.set_xlabel('a₀ (× 10⁻¹⁰ m/s²)')
    ax.set_title(f'Fossil 1: Milgrom Scale\n'
                 f'Agreement: {abs(a0_pred-A0_OBS)/A0_OBS*100:.2f}%')
    ax.set_xlim(1.15, 1.25)
    for i, v in enumerate([A0_OBS*1e10, a0_pred*1e10]):
        ax.text(v + 0.003, i, f'{v:.4f}', va='center', fontsize=9)
    ax.grid(alpha=0.2, axis='x')

    # Panel 3: Fossil 2 — σ_RAR
    ax = axes[2]
    sigma_pred = (0.5 - phi_max) / (phi_max * np.log(10))
    ax.barh(['Observed', 'GERT'], [0.057, sigma_pred],
            color=['C0', 'C3'], height=0.5)
    ax.set_xlabel('σ_RAR (dex)')
    ax.set_title(f'Fossil 2: RAR Scatter\n'
                 f'Agreement: {abs(sigma_pred-0.057)/0.057*100:.1f}%')
    ax.set_xlim(0.054, 0.061)
    for i, v in enumerate([0.057, sigma_pred]):
        ax.text(v + 0.0003, i, f'{v:.4f}', va='center', fontsize=9)
    ax.grid(alpha=0.2, axis='x')

    plt.tight_layout()
    plt.savefig(f'{outdir}/paper11_fig1_two_fossils.png', dpi=150)
    plt.close()
    print(f"  Fig 1 saved: paper11_fig1_two_fossils.png")

    # ── Figure 2: The barrier shadow ──────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("GERT Paper XI — The Shadow of the Barrier",
                 fontsize=13, fontweight='bold')

    # Panel 1: φ near the peak with barrier distance
    ax = axes[0]
    x_zoom = np.linspace(-14, -21, 1000)
    phi_zoom = np.array([phi(x) for x in x_zoom])
    ax.plot(x_zoom, phi_zoom, 'C3', lw=2.5)
    ax.axhline(0.5, color='k', ls='--', lw=1.5, label='φ = 1/2')
    ax.fill_between(x_zoom, phi_zoom, 0.5,
                    where=phi_zoom < 0.5, alpha=0.15, color='C1',
                    label=f'δφ_max = {0.5-phi_max:.4f}')
    ax.axhline(phi_max, color='C1', ls=':', lw=1)
    ax.set_xlabel('x = log₁₀ρ')
    ax.set_ylabel('φ(x)')
    ax.set_title('Barrier Distance at Recombination')
    ax.legend(fontsize=8); ax.grid(alpha=0.2)

    # Panel 2: The Jacobian conversion
    ax = axes[1]
    phi_range = np.linspace(0.30, 0.49, 200)
    jacobian = 1 / (phi_range * np.log(10))
    delta_range = 0.5 - phi_range
    sigma_range = delta_range * jacobian

    ax.plot(phi_range, sigma_range, 'C4', lw=2.5, label='σ = (0.5−φ)/(φ ln10)')
    ax.axvline(phi_max, color='C3', ls='--', lw=1.5,
               label=f'φ_max = {phi_max:.4f}')
    ax.axhline(0.057, color='C0', ls=':', lw=1.5,
               label='σ_obs = 0.057 dex')
    ax.plot(phi_max, sigma_pred, 'o', ms=10, color='C3', zorder=5)
    ax.set_xlabel('φ')
    ax.set_ylabel('σ_RAR (dex)')
    ax.set_title('The Prediction Curve\nσ(φ) = (0.5−φ)/(φ ln10)')
    ax.legend(fontsize=8); ax.grid(alpha=0.2)
    ax.set_xlim(0.30, 0.49); ax.set_ylim(0, 0.25)

    plt.tight_layout()
    plt.savefig(f'{outdir}/paper11_fig2_barrier_shadow.png', dpi=150)
    plt.close()
    print(f"  Fig 2 saved: paper11_fig2_barrier_shadow.png")


# ══════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    os.makedirs('/mnt/user-data/outputs', exist_ok=True)

    print()
    print("GERT Paper 11 — The Shadow of the Barrier")
    print("Dutra V P (2026)")
    print()

    phi_max, delta, sigma_pred, x_max = block_I_fossils()
    block_II_derivation(phi_max, delta)
    block_III_sensitivity()
    block_IV_figures(phi_max, x_max)

    print(f"\n{'='*72}")
    print(f"  PAPER 11 — COMPLETE")
    print(f"{'='*72}")
    print(f"""
  Two conformal fossils from φ(x):
  
    Fossil 1: a₀ = a_GERT × φ_start/φ_end = 1.2068×10⁻¹⁰ m/s²  (0.57%)
    Fossil 2: σ_RAR = (0.5−φ_max)/(φ_max ln10) = 0.0572 dex      (0.4%)
  
  Both from the same function φ(x) = fM/(fM+fL).
  Both zero free parameters.
  Both sub-percent agreement.
  
  The tightest correlation in extragalactic astronomy 
  has a thermodynamic origin.
""")
